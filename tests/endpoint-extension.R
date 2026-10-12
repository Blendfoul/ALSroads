# Run from the repository: Rscript --vanilla tests/endpoint-extension.R geometry
# R CMD check uses the installed namespace instead of repository source files.
test_env <- new.env(parent = globalenv())
if (file.exists("R/line_tools.R")) {
  sys.source("R/line_tools.R", envir = test_env)
} else {
  ns <- asNamespace("ALSroads")
  test_env$trim_extended_path <- get("trim_extended_path", ns)
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 0L || identical(args, "geometry"))

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
