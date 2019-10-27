context("Tests: TTestA()")

test_that("TTestA() works", {
  # Prepare data
  x = with(sleep, extra[group == 1])
  y = with(sleep, extra[group == 2])
  
  mx = mean(x)
  sx = sd(x)
  nx = length(x)
  
  my = mean(y)
  sy = sd(y)
  ny = length(y)

  # Unit tests
  expect_error(TTestA())
  expect_error(TTestA(1, 1))
  expect_error(TTestA(1, 1, 1))
  
  my_rez <- TTestA(1, 1, 10)
  expect_is(my_rez, class = "htest")
  expect_length(my_rez, n = 10)
  expect_named(my_rez, names(t.test(x)))
  
  expect_identical(
    TTestA(mx = mx, sx = sx, nx = nx), 
    t.test(x)
  )
  
  expect_identical(
    TTestA(mx = mx, sx = sx, nx = nx, conf.level = .99, mu = .75,
      alternative = "greater"), 
    t.test(x, conf.level = .99, mu = .75, alternative = "greater")
  )
  
  expect_identical(
    TTestA(mx = mx, sx = sx, nx = nx, my = my, sy = sy, ny = ny),
    t.test(x, y)
  )
    
})
