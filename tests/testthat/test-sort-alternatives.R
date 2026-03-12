test_that("all sort alternatives produce correct output", {
  results <- robscale:::sort_shootout_correctness()
  failures <- results[!results$matches_reference, ]
  expect_equal(nrow(failures), 0L,
    info = paste("Failures:", paste(capture.output(print(failures)), collapse = "\n")))
})

test_that("sort shootout benchmark runs without error", {
  bm <- robscale:::sort_shootout_benchmark(rounds = 10L, batch = 100L)
  expect_true(is.data.frame(bm))
  expect_true(all(c("method", "n", "median_ns") %in% names(bm)))
  # At least 4 methods x 7 sizes = 28 rows
  expect_gte(nrow(bm), 28L)
})
