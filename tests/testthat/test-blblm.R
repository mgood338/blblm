test_that("blblm() computes linear regression with bag of little bootstraps without parallelization", {

  fit4 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)

  expect_s3_class(fit4, "blblm")
})


test_that("blblm() computes linear regression with bag of little bootstraps with parallelization", {

  fit5 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)

  expect_s3_class(fit5, "blblm")
})