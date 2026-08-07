test_that("kIndex and ijIndex are inverse of each other", {

  nr <- 3
  nc <- 5
  ij <- expand.grid(i = 1:nr, j = 1:nc)

  # Column-stacking ordering
  kCols <- kIndex(i = ij[, 1], j = ij[, 2], nr = nr, nc = nc)
  expect_equal(unname(ijIndex(kCols, nr = nr, nc = nc)),
               unname(as.matrix(ij)))

  # Row-stacking ordering
  kRows <- kIndex(i = ij[, 1], j = ij[, 2], nr = nr, nc = nc, byRows = TRUE)
  expect_equal(unname(ijIndex(kRows, nr = nr, nc = nc, byRows = TRUE)),
               unname(as.matrix(ij)))

})

test_that("kColToRow and kRowToCol permute linear indexes consistently", {

  nr <- 2
  nc <- 5
  k <- 1:(nr * nc)
  expect_equal(kRowToCol(kColToRow(k, nr = nr, nc = nc), nr = nr, nc = nc), k)
  expect_equal(kColToRow(kRowToCol(k, nr = nr, nc = nc), nr = nr, nc = nc), k)

  # kColToRow permutes a column-stacked vector into its row-stacked ordering
  A <- matrix(10 * seq_len(nr * nc), nrow = nr, ncol = nc)
  expect_equal(as.vector(A)[kColToRow(k, nr = nr, nc = nc)], as.vector(t(A)))

})

test_that("repRow and repCol replicate as documented", {

  expect_equal(repRow(1:5, 2), matrix(1:5, nrow = 2, ncol = 5, byrow = TRUE))
  expect_equal(repCol(1:5, 2), matrix(1:5, nrow = 5, ncol = 2))

  A <- rbind(1:5, 5:1)
  expect_equal(dim(repRow(A, 2)), c(nrow(A) * 2, ncol(A)))
  expect_equal(dim(repCol(A, 2)), c(nrow(A), ncol(A) * 2))

})

test_that("matMatch finds the correct rows and columns", {

  A <- rbind(5:6, repRow(1:2, 3), 3:4)
  B <- unique(A)
  ind <- matMatch(x = A, mat = B)
  expect_equal(B[ind, ], A)

  # By columns and via base::match
  A2 <- cbind(5:6, repCol(1:2, 3), 3:4)
  B2 <- t(unique(t(A2)))
  ind2 <- matMatch(x = A2, mat = B2, rows = FALSE)
  expect_equal(B2[, ind2], A2)
  expect_equal(matMatch(x = A, mat = B, useMatch = TRUE), ind)

})
