box::use(
  testthat[expect_equal, test_that],
  utils[read.csv],
)
box::use(
  app/logic/data_input[example_file],
  app/logic/mock_data[make_mock_data],
)

test_that("the generator reproduces the shipped example", {
  expect_equal(make_mock_data(), read.csv(example_file, check.names = FALSE))
})
