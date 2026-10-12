# Run from the repository: Rscript --vanilla tests/endpoint-extension.R geometry
# R CMD check uses the installed namespace instead of repository source files.
test_env <- new.env(parent = globalenv())
if (file.exists("R/line_tools.R")) {
  for (file in list.files("R", pattern = "[.]R$", full.names = TRUE))
    sys.source(file, envir = test_env)
} else {
  ns <- asNamespace("ALSroads")
  parent.env(test_env) <- ns
  test_env$trim_extended_path <- get("trim_extended_path", ns)
  for (name in c("measure_road", "measure_roads")) {
    test_env[[name]] <- get(name, ns)
    environment(test_env[[name]]) <- test_env
  }
  test_env$alsroads_default_parameters <- get("alsroads_default_parameters", ns)
  test_env$alsroads_default_parameters2 <- get("alsroads_default_parameters2", ns)
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 0L || (length(args) == 1L && args %in% c("geometry", "pipeline")))

geometry_checks <- function() {
# Losing an interior bend, reversing the result, or dropping attributes must fail.
reference <- sf::st_sf(road = sf::st_sfc(sf::st_linestring(rbind(
  c(-30, 0), c(0, 0), c(40, 0), c(40, 40), c(80, 40), c(110, 40)
)), crs = 32198))
xy <- rbind(c(-30, 0), c(0, 0), c(30, 20), c(40, 40), c(80, 40), c(110, 40))
path <- sf::st_sf(ID = 46, CONDUCTIVITY = 0.75,
                  road = sf::st_sfc(sf::st_linestring(xy), crs = 32198))
forward <- test_env$trim_extended_path(path, reference, 30)
reverse_path <- path
sf::st_geometry(reverse_path) <- sf::st_sfc(sf::st_linestring(xy[nrow(xy):1, ]), crs = 32198)
reversed <- test_env$trim_extended_path(reverse_path, reference, 30)
expected <- rbind(c(0, 0), c(30, 20), c(40, 40), c(80, 40))
expected_line <- sf::st_sfc(sf::st_linestring(expected), crs = 32198)
ends <- sf::st_coordinates(forward)[c(1, nrow(sf::st_coordinates(forward))), 1:2]
stopifnot(as.numeric(sf::st_distance(forward, expected_line, which = "Hausdorff")) < 1e-8,
          isTRUE(all.equal(ends, expected[c(1, 4), ],
                          check.attributes = FALSE, tolerance = 1e-8)),
          isTRUE(all.equal(sf::st_coordinates(reversed), sf::st_coordinates(forward))),
          identical(sf::st_crs(forward), sf::st_crs(reference)),
          identical(attr(forward, "sf_column"), "road"),
          identical(sf::st_drop_geometry(forward), sf::st_drop_geometry(path)))

# An extension that consumes the whole path cannot be used for measurement.
for (length in c(59, 60)) {
  short <- path
  sf::st_geometry(short) <- sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(length, 0))), crs = 32198)
  stopifnot(is.null(test_env$trim_extended_path(short, reference, 30)))
}
empty <- path
sf::st_geometry(empty) <- sf::st_sfc(sf::st_linestring(matrix(numeric(), ncol = 2)), crs = 32198)
disconnected <- path
sf::st_geometry(disconnected) <- sf::st_sfc(sf::st_multilinestring(list(xy[1:3, ], xy[4:6, ])), crs = 32198)
stopifnot(is.null(test_env$trim_extended_path(empty, reference, 30)),
          is.null(test_env$trim_extended_path(disconnected, reference, 30)))
cat("geometry checks passed\n")
}

pipeline_checks <- function() {
  options(ALSroads.debug.finding = FALSE, ALSroads.debug.measuring = FALSE,
          ALSroads.debug.metrics = FALSE, ALSroads.debug.verbose = FALSE,
          ALSroads.debug.progress = FALSE)
  dir <- system.file("extdata", package = "ALSroads")
  if (!nzchar(dir)) dir <- "inst/extdata"
  ctg <- lidR::readLAScatalog(dir)
  data <- data.frame(X = c(0, 1), Y = c(0, 1), Z = c(0, 0))
  las <- lidR::LAS(data, lidR::LASheader(data))
  original <- sf::st_sf(ID = 46, road = sf::st_sfc(sf::st_linestring(
    rbind(c(0, 0), c(40, 0), c(40, 40), c(80, 40))), crs = 32198))

  # Capture the real orchestration; replace only expensive extraction/search/measurement.
  run_case <- function(param, road = original, outcome = "ok", reverse = FALSE,
                       retry_error = FALSE) {
    env <- new.env(parent = test_env)
    for (name in c("measure_road", "measure_roads")) {
      env[[name]] <- test_env[[name]]
      environment(env[[name]]) <- env
    }
    observed <- new.env()
    observed$extracted <- observed$searched <- observed$measured <- observed$spline_inputs <- list()
    observed$error <- NULL
    observed$warnings <- character()
    env$extract_road <- function(las, road, param) {
      observed$extracted[[length(observed$extracted) + 1L]] <- road
      if (outcome == "empty" && length(observed$extracted) == 1L) return(las_fixture[0L])
      las_fixture
    }
    las_fixture <- las
    env$least_cost_path <- function(las, centerline, dtm, conductivity, water, param) {
      observed$searched[[length(observed$searched) + 1L]] <- centerline
      if (retry_error && length(observed$searched) == 2L) stop("original retry unavailable")
      path <- sf::st_sf(CONDUCTIVITY = 0.75, geometry = sf::st_geometry(centerline))
      if (length(observed$extracted) == 1L) {
        if (outcome == "error") stop("extension unavailable")
        if (outcome == "zero") path$CONDUCTIVITY <- 0
        if (outcome == "disconnected") {
          xy <- sf::st_coordinates(centerline)[, 1:2]
          sf::st_geometry(path) <- sf::st_sfc(sf::st_multilinestring(
            list(xy[1:2, ], xy[(nrow(xy) - 1):nrow(xy), ])), crs = 32198)
        }
        if (outcome == "short")
          sf::st_geometry(path) <- sf::st_sfc(sf::st_linestring(
            rbind(c(-30, 0), c(50, 0))), crs = 32198)
      }
      if (reverse) {
        xy <- sf::st_coordinates(path)[, 1:2]
        sf::st_geometry(path) <- sf::st_sfc(sf::st_linestring(xy[nrow(xy):1, ]), crs = 32198)
      }
      path
    }
    env$road_measure <- function(las, road, param) {
      observed$measured[[length(observed$measured) + 1L]] <- road
      points <- sf::st_cast(sf::st_line_sample(sf::st_geometry(road),
                                              sample = seq(0, 1, length.out = 6)), "POINT")
      sf::st_sf(road_width = rep(8, 6), drivable_width = rep(8, 6),
                pzabove05 = rep(0, 6), pzabove2 = rep(0, 6),
                number_accotements = rep(2, 6), zroad = rep(100, 6),
                distance_to_start = seq(0, as.numeric(sf::st_length(road)), length.out = 6),
                find_score = rep(1, 6), geometry = points)
    }
    env$adjust_spline <- function(points) {
      observed$spline_inputs[[length(observed$spline_inputs) + 1L]] <- points
      sf::st_geometry(observed$measured[[length(observed$measured)]])
    }
    invisible(capture.output(result <- tryCatch(withCallingHandlers(
      env$measure_road(ctg, road, conductivity = 1, param = param, Windex = FALSE),
      warning = function(w) {
        observed$warnings <- c(observed$warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        observed$error <- conditionMessage(e)
        NULL
      })))
    c(list(result = result), as.list(observed))
  }

  param <- test_env$alsroads_default_parameters
  param$extraction$endpoint_extension <- 0
  zero <- run_case(param)
  stopifnot(is.null(zero$error), zero$result$CLASS == 1,
            identical(sf::st_drop_geometry(zero$result)["ID"], sf::st_drop_geometry(original)["ID"]),
            identical(attr(zero$result, "sf_column"), "road"))
  missing <- param
  missing$extraction$endpoint_extension <- NULL
  missing_parameter <- run_case(missing)
  stopifnot(identical(zero$result, missing_parameter$result), length(zero$searched) == 1L)
  param$extraction$endpoint_extension <- 30
  enabled <- run_case(param)
  stopifnot(is.null(enabled$error),
            as.numeric(sf::st_length(enabled$extracted[[1]])) == 180,
            identical(sf::st_geometry(enabled$extracted[[1]]), sf::st_geometry(enabled$searched[[1]])),
            as.numeric(sf::st_distance(enabled$measured[[1]], original, which = "Hausdorff")) < 1e-8,
            identical(sf::st_crs(enabled$result), sf::st_crs(original)),
            identical(attr(enabled$result, "sf_column"), "road"),
            enabled$result$ID == 46, enabled$result$CLASS == 1,
            as.numeric(sf::st_length(enabled$measured[[1]])) < 121,
            length(enabled$spline_inputs) == 1L,
            max(enabled$spline_inputs[[1]]$distance_to_start) < 121)
  reversed <- run_case(param, reverse = TRUE)
  stopifnot(identical(enabled$result, reversed$result))
  no_relocation_param <- param
  no_relocation_param$constraint$confidence <- 1
  no_relocation <- run_case(no_relocation_param)
  stopifnot(is.null(no_relocation$error), length(no_relocation$searched) == 0L,
            identical(sf::st_geometry(no_relocation$extracted[[1]]), sf::st_geometry(original)))

  for (outcome in c("error", "empty", "zero", "disconnected", "short")) {
    fallback <- run_case(param, outcome = outcome)
    stopifnot(is.null(fallback$error), length(fallback$extracted) == 2L,
              identical(sf::st_geometry(fallback$extracted[[2]]), sf::st_geometry(original)),
              identical(fallback$result, zero$result))
  }
  failed_retry <- run_case(param, outcome = "error", retry_error = TRUE)
  stopifnot(identical(failed_retry$error, "original retry unavailable"))
  for (distance in list(-1, NA_real_, NaN, Inf, "30", numeric(), c(0, 30), TRUE)) {
    invalid <- param
    invalid$extraction$endpoint_extension <- distance
    rejected <- run_case(invalid)
    stopifnot(identical(rejected$error, "'endpoint_extension' must be a finite, non-negative number"),
              length(rejected$extracted) == 0L)
  }

  # Existing angle validation rejects duplicated terminal coordinates before extraction.
  duplicate <- original
  sf::st_geometry(duplicate) <- sf::st_sfc(sf::st_linestring(rbind(
    c(0, 0), c(0, 0), c(40, 0), c(80, 40))), crs = 32198)
  duplicate_zero <- run_case(missing, road = duplicate)
  duplicate_enabled <- run_case(param, road = duplicate)
  stopifnot(!is.null(duplicate_zero$error),
            identical(duplicate_zero$error, duplicate_enabled$error),
            length(duplicate_enabled$extracted) == 0L)
  long <- original
  sf::st_geometry(long) <- sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(240, 0))), crs = 32198)
  chunk_param <- param
  chunk_param$extraction$road_max_len <- 100
  chunked <- run_case(chunk_param, road = long)
  stopifnot(is.null(chunked$error), length(chunked$searched) == 3L,
            all(vapply(chunked$searched, function(x) abs(as.numeric(sf::st_length(x)) - 140) < 1e-8, TRUE)))
  loop <- original
  sf::st_geometry(loop) <- sf::st_sfc(sf::st_linestring(rbind(
    c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0))), crs = 32198)
  looped <- run_case(param, road = loop)
  stopifnot(is.null(looped$error), length(looped$searched) == 2L,
            all(vapply(looped$searched, function(x) abs(as.numeric(sf::st_length(x)) - 260) < 1e-8, TRUE)))
  stopifnot(test_env$alsroads_default_parameters$extraction$endpoint_extension == 0,
            test_env$alsroads_default_parameters2$extraction$endpoint_extension == 0)
  cat("pipeline checks passed\n")
}

if (length(args) == 0L || args == "geometry") geometry_checks()
if (length(args) == 0L || args == "pipeline") pipeline_checks()
