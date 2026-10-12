# Endpoint Extension Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Native execution is recommended for these two closely related tasks.

**Goal:** Add opt-in endpoint extension to road realignment so pathfinding can use context beyond the original endpoints, while returned measurements cover only the trimmed road.

**Architecture:** Keep the original centerline for validation, attributes, and existing fallback behavior. After existing short-road checks and chunk splitting, attempt extraction and pathfinding on an extended copy, then orient and trim the corrected geometry before measurement and spline fitting. If that attempt is unusable, run the original extraction and search once.

**Tech Stack:** R >= 4.1.0; existing sf, lwgeom, lidR, raster and gdistance dependencies; base-R regression assertions.

**Spec:** The bounded design approved in this chat on 2026-10-12, reproduced below, based on [ALSroads issue #46](https://github.com/r-lidar-lab/ALSroads/issues/46). This plan is explicitly requested; there is no separate architectural spec.

## Approved design and global constraints

- Extend both endpoints along their existing directions using `st_extend_line(line, distance, end = "BOTH")`; do not add neighbor selection.
- Add `param$extraction$endpoint_extension = 0`. This is a distance in the projected CRS's linear units, consistent with `section_length` and `road_buffer`. Missing/NULL values also mean zero, preserving older parameter lists.
- Require one finite, non-negative numeric value. Report invalid values before starting extraction. An initial opt-in setting is `3 * param$extraction$section_length` (30 with current defaults).
- Zero extension and `constraint$confidence = 1` preserve the existing flow. Original length checks and chunk splitting remain ahead of extension; apply the same distance to each processed chunk.
- Extend before point-cloud extraction. Pass the extended copy through the existing conductivity generation/cropping, water masking, end-cap construction, and least-cost search.
- Orient the resulting path consistently with the search centerline, then remove the configured distance from each end of that path by arc length. Trimming is approximate relative to the original endpoint positions, as proposed in the issue.
- Measure widths and fit the spline only after trimming. Preserve CRS, the geometry column name, original attributes, and existing class-based geometry retention.
- Retry original extraction and pathfinding if the extended attempt errors, extracts no points, reports zero conductivity, or yields invalid, disconnected, empty, or too-short geometry. Let errors from the original attempt propagate normally.
- Retain the existing LCP `CONDUCTIVITY` estimate; when extension is enabled it describes the full search route, including the extensions. Do not introduce a new scoring algorithm. Document this limitation because it can affect SCORE/CLASS.
- Add no runtime or testing dependencies, public function arguments, changes to snapping, or unrelated refactoring. Leave `st_extend_line()` unchanged for its existing callers.

## Review focus

1. Reversed corrected paths must trim to the same geometry and direction (Task 1).
2. Curved paths must retain interior bends, attributes and CRS after trimming (Task 1).
3. Missing coverage, failed searches and unusable trims must retry the original flow, rather than classify the road from the failed extension (Task 2).
4. Older parameter lists, invalid distances and confidence-one calls must behave as specified (Task 2).
5. Split roads and duplicate endpoint coordinates must avoid corrupting geometry: extend each chunk after splitting; skip the extension attempt when either terminal segment has zero length (Task 2).

## File map

| File | Responsibility |
| --- | --- |
| `R/line_tools.R` | Add one internal orientation-and-trimming helper beside `st_extend_line()` |
| `R/als_road.R` | Validate the distance and integrate the optional extended attempt into `measure_road()` |
| `R/parameters.R` | Default parameter, roxygen documentation, opt-in example |
| `man/alsroads_default_parameters.Rd` | Regenerated parameter documentation |
| `tests/endpoint-extension.R` | One runnable regression script, with geometry and pipeline sections |
| `NEWS.md` | Announce the opt-in behavior and link issue #46 |
| `.Rbuildignore` | Exclude `docs/superpowers` planning artifacts from the R package build |

No edits are needed in `R/lcp_tools.R` or `R/query_tools.R`: their existing centerline arguments carry the extended geometry through the complete search.

### Task 1: Orient and trim an extended path

**Files:** Modify `R/line_tools.R`; create `tests/endpoint-extension.R`.

**Interfaces:**

- Consumes: existing `st_extend_line(line, distance, end = "BOTH")` and `lwgeom::st_linesubstring(x, from, to)`.
- Produces: internal `trim_extended_path(path, reference, distance)`, where `path` and `reference` are single-feature projected sf LINESTRING objects. Return a trimmed sf object with path attributes preserved, or NULL for unusable geometry. `reference` is the extended search centerline, not the original line.
- Test command: `Rscript --vanilla tests/endpoint-extension.R geometry`. Source `R/line_tools.R` into a local environment; this section needs only sf and lwgeom.

- [ ] **Write failing geometry assertions.** Construct reference coordinates `(-30,0), (0,0), (40,0), (40,40), (80,40), (110,40)` and a corrected path `(-30,0), (0,0), (30,20), (40,40), (80,40), (110,40)` in EPSG:32198. Give the path `ID = 46`, `CONDUCTIVITY = 0.75`, and a geometry column named `road`. Define `forward` from this path, `reversed` from reversed coordinates, and `too_short` from a 60-unit straight line, all trimmed by 30. Assert:

```r
expected <- rbind(c(0, 0), c(30, 20), c(40, 40), c(80, 40))
stopifnot(isTRUE(all.equal(sf::st_coordinates(forward)[, 1:2],
                          expected, check.attributes = FALSE, tolerance = 1e-8)))
stopifnot(isTRUE(sf::st_equals_exact(forward, reversed, par = 1e-8, sparse = FALSE)[1, 1]))
stopifnot(identical(sf::st_crs(forward), sf::st_crs(reference)),
          identical(attr(forward, "sf_column"), "road"),
          identical(sf::st_drop_geometry(forward), sf::st_drop_geometry(path)),
          is.null(too_short))
```

Also assert NULL for empty or MULTILINESTRING input, and for path length below twice the extension.

- [ ] **Run the geometry command.** Expected: failure because `trim_extended_path()` does not exist.
- [ ] **Implement the helper in `R/line_tools.R`.** Reject unusable geometry and non-finite/insufficient length. Compare the summed squared distances from path endpoints to the reference endpoints in the forward and reverse pairing; reverse coordinates only when the reverse pairing is closer. Trim using fractions `distance / length` and `1 - distance / length`, updating geometry without rebuilding the attribute table. Require a nonempty LINESTRING result. Keep this helper internal with `@noRd`.
- [ ] **Run the geometry command again.** Expected: exit code 0 and a printed `geometry checks passed` message. Both directions must also have the same start and end coordinates; geometric equality alone does not prove direction.
- [ ] **Commit the helper and regression script.** Suggested message: `feat: add orientation-aware trimming for extended road paths`.

### Task 2: Integrate opt-in extension, fallback, and documentation

**Files:** Modify `R/als_road.R`, `R/parameters.R`, `tests/endpoint-extension.R`, `NEWS.md`, `.Rbuildignore`; regenerate `man/alsroads_default_parameters.Rd`.

**Interfaces:**

- Consumes: Task 1's `trim_extended_path(path, reference, distance)`.
- Produces: unchanged `measure_road()` and `measure_roads()` signatures and result schemas; optional behavior controlled only by `param$extraction$endpoint_extension`.
- Test commands: `Rscript --vanilla tests/endpoint-extension.R pipeline` and `Rscript --vanilla tests/endpoint-extension.R` (both sections). With no arguments, run both sections. In repository runs, source local files; in R CMD check runs, obtain the corresponding functions from the installed ALSroads namespace. Rebind the test copies of both `measure_road()` and `measure_roads()` to the local test environment so recursive chunk calls use the same captured functions. Load bundled fixtures with `system.file("extdata", package = "ALSroads")`.

- [ ] **Write failing pipeline checks.** Use the bundled LAScatalog as a valid input, and deterministic nonempty LAS/path results. Source modified functions into a local environment whose parent is the installed ALSroads namespace; override only `extract_road`, `least_cost_path`, `road_measure`, and `adjust_spline` there to capture inputs and control outcomes. Do not mutate the package namespace. Use `Windex = FALSE` for these calls. Leave real metric aggregation, scoring, and final attribute assembly active. Establish these cases and assertions:

```r
# Records contain result, extracted/search geometries, and measurement/spline inputs.
stopifnot(identical(zero$result, missing_parameter$result),
          length(zero$searched) == 1L,
          length(no_relocation$searched) == 0L)
stopifnot(isTRUE(all.equal(as.numeric(sf::st_length(enabled$extracted[[1]])),
                          as.numeric(sf::st_length(original)) + 60)),
          isTRUE(sf::st_equals_exact(enabled$extracted[[1]],
                                    enabled$searched[[1]], par = 1e-8, sparse = FALSE)[1, 1]))
stopifnot(isTRUE(sf::st_equals_exact(enabled$measured[[1]],
                                    expected_trimmed, par = 1e-8, sparse = FALSE)[1, 1]))
```

Define the named records by calling the actual `measure_road()` orchestration with distance zero, a missing field, distance 30, and confidence one. Use more than four valid slice-metric rows in the enabled case so spline capture proves that no extension points reach smoothing. Assert the same enabled result for reversed LCP output, retained original attributes/geometry-column name/CRS, and zero defaults in both exported parameter lists. Test each fallback trigger separately: extension error, empty extraction, zero conductivity, MULTILINESTRING, and a trimmed length below `4 * section_length`. Each must make a second extraction on the original geometry and produce the same returned result as zero extension. If the original retry errors, assert that its error reaches the caller. For invalid distances (`-1`, NA, Inf, character, numeric(0), and length two), assert the documented validation error and no extraction calls. For duplicate terminal coordinates, assert only the original search is attempted. For a straight road over `road_max_len` and a closed loop, assert each search receives an extended chunk, rather than an extension being included in the initial split.

- [ ] **Run pipeline checks against current behavior.** Expected: failure on the extended-extraction assertion or invalid-distance validation. Existing zero-distance baseline checks must pass.
- [ ] **Add the parameter and integration.** Validate the new value near the existing input checks, using error text `'endpoint_extension' must be a finite, non-negative number`. After splitting, and only when relocation is enabled and distance is positive, construct an extended copy by assigning `st_extend_line(sf::st_geometry(centerline), distance)` to its geometry. Skip this attempt for zero-length terminal segments. In a narrowly scoped `tryCatch`, extract points, search, trim, and require positive finite conductivity plus trimmed length at least `length_min`; return a local list containing `las` and `res` on success. On failure, fall through to the original extraction and original LCP call. Reuse the successful attempt's point cloud and trimmed result in the existing downstream flow. Keep original geometry checks, class handling, measurement, smoothing, and final attribute assembly intact.
- [ ] **Run both regression sections.** Expected: exit code 0, geometry and pipeline success messages, no changes to the zero/missing-distance baseline.
- [ ] **Document the option.** Add `endpoint_extension = 0` to the extraction defaults, explain projected units, symmetric trimming and fallback, and show `param$extraction$endpoint_extension <- 3 * param$extraction$section_length`. State that CONDUCTIVITY still summarizes the full search route. Update NEWS with issue #46; add `^docs/superpowers$` to `.Rbuildignore`. Regenerate Rd with `Rscript --vanilla -e "roxygen2::roxygenise('.', roclets = 'rd')"`; inspect and discard unrelated generated changes.
- [ ] **Verify the package and sample data.** Run the commands below. Compare disabled and enabled `measure_road()` results on the bundled j5gr LAScatalog, centerline, and precomputed conductivity using the existing README/example loading pattern. Confirm valid returned LINESTRING, retained CRS/attributes, and usable metrics; record endpoint displacement and SCORE/CLASS changes. This smoke comparison does not establish a universal accuracy improvement, because the bundle has no endpoint ground-truth assertion.

```text
R CMD INSTALL --no-multiarch .
Rscript --vanilla tests/endpoint-extension.R
R CMD build .
R CMD check --no-manual ALSroads_0.2.0.tar.gz
git diff --check
```

Expected: installation/build succeed, regression assertions pass, no new package-check errors or warnings attributable to this change, and no whitespace errors. Report existing or environment-dependent check findings separately. Rscript was not found on PATH during planning; execution needs a working R installation, the existing package dependencies and native build toolchain. roxygen2 is a documentation tool, not a new package dependency.

- [ ] **Commit the integrated feature, tests and documentation.** Suggested message: `feat: add opt-in endpoint extension with legacy fallback`.

## Completion and handoff

- [ ] Confirm every approved requirement and review-focus item has its implementation and regression check.
- [ ] Review the final diff for unrelated changes and accidental edits to `st_extend_line()` or snapping.
- [ ] Report test outcomes and the observed sample-data comparison, including the full-search conductivity limitation.

Plan self-review completed against the approved chat design. Implementation has not started. Recommended execution: native, because both tasks share one small geometry contract and one regression script. Review and execution-method selection are required before implementation under the writing-plans skill.
