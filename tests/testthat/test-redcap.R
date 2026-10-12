# REDCAP (full-order Ward). Needs the compiled Rust backend, so skipped on
# CRAN like the other regionalization tests.

nc_data <- function() {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")
  skip_on_cran()
  sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
}

regions_contiguous <- function(result) {
  nb <- spdep::poly2nb(result, queen = TRUE)
  all(vapply(split(seq_len(nrow(result)), result$.region), function(idx) {
    length(idx) == 1 || suppressWarnings(spdep::n.comp.nb(spdep::subset.nb(nb, seq_len(nrow(result)) %in% idx))$nc) == 1
  }, logical(1)))
}

test_that("redcap returns the requested number of contiguous regions", {
  nc <- nc_data()
  result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6)
  expect_s3_class(result, "sf")
  expect_true(".region" %in% names(result))
  expect_equal(length(unique(result$.region)), 6)
  expect_true(regions_contiguous(result))
  meta <- attr(result, "spopt")
  expect_equal(meta$algorithm, "redcap")
  expect_equal(meta$method, "fullorder-ward")
  expect_gt(meta$objective, 0)
})

test_that("redcap respects a floor constraint", {
  nc <- nc_data()
  result <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 5,
                   floor = "BIR74", floor_value = 50000)
  expect_true(all(tapply(nc$BIR74, result$.region, sum) >= 50000))
  expect_true(regions_contiguous(result))
})

test_that("redcap is deterministic and accepts the rgeoda method name", {
  nc <- nc_data()
  a <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6)
  b <- redcap(nc, attrs = c("SID74", "SID79"), n_regions = 6, method = "fullorder-wardlinkage")
  expect_identical(a$.region, b$.region)
  expect_error(redcap(nc, attrs = "SID74", n_regions = 4, method = "single"), "should be")
})

test_that("redcap stays close to the constrained Ward partition at the same k", {
  # REDCAP full-order Ward cuts the tree built from the same merge sequence
  # that ward_spatial() uses. Each cut re-optimizes where a region splits, so
  # the two usually land within a few percent of each other (REDCAP is
  # typically, not always, the more homogeneous). This guards against a
  # broken tree, not against small differences.
  nc <- nc_data()
  for (k in c(4, 6, 10)) {
    r <- redcap(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = k)
    w <- ward_spatial(nc, attrs = c("SID74", "SID79", "BIR74"), n_regions = k)
    expect_equal(length(unique(r$.region)), k)
    expect_true(regions_contiguous(r))
    expect_lte(attr(r, "spopt")$objective, attr(w, "spopt")$objective * 1.10)
  }
})
