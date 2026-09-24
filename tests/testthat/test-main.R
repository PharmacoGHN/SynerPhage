box::use(
  shiny[testServer],
  testthat[expect_equal, expect_match, test_that],
)
box::use(
  app/main[server],
)

test_that("a typed pair flows through the engine to the headline numbers", {
  testServer(server, {
    session$setInputs(mode = "pair", pair_strain = "Strain 1", pair_name_a = "P1",
                      pair_name_b = "P2", pair_a = "0.2 0.25", pair_b = "0.4, 0.45",
                      pair_ab = "0.9 0.88", scale = "auto", auto_rope = TRUE)
    expect_equal(output$vb_strains, "1")
    expect_equal(output$vb_combos, "1")
    expect_equal(output$vb_calls, "1 / 1")
  })
})

test_that("the example table is read and analysed", {
  testServer(server, {
    session$setInputs(mode = "table", source = "example", layout = "wide",
                      strain_col = "Bacteria", rep_col = "Replica", sep = "",
                      ignore_cols = "Nophage", scale = "auto", auto_rope = TRUE)
    expect_equal(output$vb_strains, "8")
    expect_equal(output$vb_combos, "4")
    expect_match(output$vb_rope, "pts$")
  })
})
