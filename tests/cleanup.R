# Run from the repository: Rscript --vanilla tests/cleanup.R
test_env <- new.env(parent = globalenv())
if (file.exists("R/als_road.R")) {
  for (file in c("als_road.R", "snap_tools.R", "lcp_tools.R"))
    sys.source(file.path("R", file), envir = test_env)
} else {
  ns <- asNamespace("ALSroads")
  for (name in c("road_class4", "make_bridge", "sobel.matrix"))
    test_env[[name]] <- get(name, ns)
}

# Class-4 fallbacks preserve attributes, geometry, CRS and column ordering.
centerline <- sf::st_sf(ID = 46, road = sf::st_sfc(sf::st_linestring(
  rbind(c(0, 0), c(10, 0))), crs = 32198), label = "original", ROADWIDTH = 9)
result <- test_env$road_class4(centerline)
stopifnot(result$ID == 46, result$label == "original",
          result$ROADWIDTH == 0, result$DRIVABLEWIDTH == 0,
          result$PERCABOVEROAD == 100, result$SHOULDERS == 0,
          is.na(result$SINUOSITY), result$CONDUCTIVITY == 0,
          result$SCORE == 0, result$CLASS == 4,
          identical(sf::st_geometry(result), sf::st_geometry(centerline)),
          identical(sf::st_crs(result), sf::st_crs(centerline)),
          identical(attr(result, "sf_column"), "road"),
          tail(names(result), 1) == "road")

# Every endpoint pairing makes the same continuous bridge.
pairings <- expand.grid(end_1 = c(FALSE, TRUE), end_2 = c(FALSE, TRUE))
for (i in seq_len(nrow(pairings))) {
  end_1 <- pairings$end_1[i]
  end_2 <- pairings$end_2[i]
  first <- rbind(c(0, 0), c(10, 0))
  second <- rbind(c(12, 0), c(20, 0))
  if (!end_1) first <- first[2:1, ]
  if (end_2) second <- second[2:1, ]
  roads <- sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(first),
                                         sf::st_linestring(second), crs = 32198))
  node <- data.frame(pos = 1:2, end = c(end_1, end_2), X = c(10, 12), Y = c(0, 0))
  bridge <- test_env$make_bridge(roads, node)
  stopifnot(isTRUE(all.equal(sf::st_coordinates(bridge$bridge)[, 1:2],
                            rbind(c(0, 0), c(11, 0), c(20, 0)),
                            check.attributes = FALSE)),
            identical(bridge$reversed, c(!end_1, end_2)),
            identical(sf::st_crs(bridge$bridge), sf::st_crs(roads)))
}

# The supported Sobel kernel keeps its edge magnitudes and missing cells.
image <- matrix(0, 5, 5)
image[3, 3] <- 1
image[1, 1] <- NA
edges <- test_env$sobel.matrix(image)
expected <- matrix(c(sqrt(2), 2, sqrt(2), 2, 0, 2, sqrt(2), 2, sqrt(2)), 3, 3)
stopifnot(is.na(edges[1, 1]),
          isTRUE(all.equal(edges[2:4, 2:4], expected, tolerance = 1e-8)))
cat("cleanup checks passed\n")
