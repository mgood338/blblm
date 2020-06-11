test_that("blblog() computes logistic regression with bag of little bootstraps without parallelization", {
  bi_iris <- sample(iris[1:100,], 10000, replace = TRUE) # binomial iris

  fit1 <- blblog(Species ~ Sepal.Width+Sepal.Length, data = bi_iris, m = 2, B = 100, parallel = FALSE)

  expect_s3_class(fit1, "blblog")
})


test_that("blblog() computes logistic regression with bag of little bootstraps with parallelization", {
  bi_iris <- sample(iris[1:100,], 10000, replace = TRUE) # binomial iris

  fit2 <- blblog(Species ~ Sepal.Width+Sepal.Length, data = bi_iris, m = 2, B = 100, parallel = TRUE)

  expect_s3_class(fit2, "blblog")
})