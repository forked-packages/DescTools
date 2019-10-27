context("Tests: ZTest()")

test_that("ZTest() works", {
  # Unit tests
  expect_error(ZTest())
  expect_error(ZTest(1))
  expect_error(ZTest(-5:5))

  my_rez <- ZTest(-5:5, mu = 0, sd_pop = 1)
  expect_is(my_rez, class = "htest")
  expect_length(my_rez, n = 10)
  # Expected names are the same as for t.test:
  # c("statistic", "parameter", "p.value",  "conf.int", "estimate", "null.value", 
  #   "stderr", "alternative", "method", "data.name")
  expect_named(my_rez, names(t.test(-5:5)))
  
  expect_equal(ZTest(-5:5, sd_pop = 1)$p.value,           expected = 1)
  expect_equal(ZTest(-5:5, mu = 100, sd_pop = 1)$p.value, expected = 0)
  
  v_1 <- rep(96, times = 55)
  z_1 <- (96 - 100) / (12/sqrt(55))
  expect_equal(
    ZTest(v_1, mu = 100, sd_pop = 12)$statistic[["z"]], 
    expected = z_1
  )
  
  expect_equal(
    ZTest(v_1, mu = 100, sd_pop = 12, alternative = "less")$p.value, 
    expected = pnorm(z_1)
  )
    
})
