context("test-dof")

test_that("DOF computation", {
  ## Unit test: DoF
  # Model that does not lose a DoF
  testModel1 <- mxModel("noRegModel", type="RAM", manifestVars = "X",
                        mxPath("X", arrows=2, values=1, free=FALSE),
                        mxPath("one", "X", values=.5, free=TRUE, label="Mean"),
                        # True mean very non-zero: DoF not recovered
                        mxData(data.frame(X=rnorm(10, 10, 1)), type="raw"))
  regTest1 <- suppressWarnings(summarizeRegularized(mxRun(regularizeMxModel(testModel1, "Mean"))))


  testModel2 <- mxModel(testModel1, name="regModel",
                        # Alter to make True Mean Zero: one DoF recovered
                        mxData(data.frame(X=rnorm(10, 0, 1)), type="raw"))
  regTest2 <- suppressWarnings(summarizeRegularized(mxRun(regularizeMxModel(testModel2, "Mean"))))

  expect_equal(regTest1$degreesOfFreedom, 9)
  expect_equal(regTest2$degreesOfFreedom, 10)
})
