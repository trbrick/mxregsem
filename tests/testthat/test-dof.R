context("test-dof")

test_that("DOF computation", {
  ## Unit test: DoF
  # Model that does not lose a DoF
  testModel1 <- mxModel("noRegModel", type="RAM", manifestVars = "X",
                        mxPath("X", arrows=2, values=1, free=FALSE),
                        mxPath("one", "X", values=.5, free=TRUE, label="Mean"),
                        # True mean very non-zero: DoF not recovered
                        mxData(data.frame(X=rnorm(10, 10, 1)), type="raw"))
  testModel1 <- mxModel(testModel1, mxRegularizeLASSO("Mean", "MeanLASSO", 10))
  regTest1 <- mxRun(testModel1)


  testModel2 <- mxRename(testModel1, oldname = "noRegModel", newname = "RegModel")
  testModel2 <- mxModel(testModel2,
                        # Alter to make True Mean Zero: one DoF recovered
                        mxData(data.frame(X=rnorm(10, 0, 1)), type="raw"))
  regTest2 <- suppressWarnings(mxRun(testModel2))

  expect_equal(summary(regTest1)$degreesOfFreedom, 9)
  expect_equal(summary(regTest2)$degreesOfFreedom, 10)
})
