# test-robscale-nr.R — Removal guards for diagnostic NR exports.
#
# WU-PROD-4: rob_scale_nr_impl and rob_scale_nr_iters were diagnostic exports
# added for the NR A/B test phase (WU-NR-1 through WU-NR-4). NR is now the
# production algorithm. These exports are dead weight and have been removed.
#
# These guards prevent re-addition.

test_that("diagnostic NR exports removed from namespace (guard against re-addition)", {
  expect_false("rob_scale_nr_impl"  %in% ls(envir = asNamespace("robscale")))
  expect_false("rob_scale_nr_iters" %in% ls(envir = asNamespace("robscale")))
})
