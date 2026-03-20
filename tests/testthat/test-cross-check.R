test_that("robscale matches revss baseline", {
  skip_if_not_installed("revss")
  tol <- 1e-4
  # revss v3.0.0 (2026-03-18) changed robScale() and robLoc() defaults:
  #   - robScale() uses bias-corrected "AA" constants (madn) as starting value
  #   - robLoc() uses the new madn-based scale internally
  # Both changes produce different output for all n. Only adm() is unaffected.
  compare_iterative <- packageVersion("revss") < "3.0.0"
  set.seed(42)
  for (n in 3:20) {
    for (rep in 1:100) {
      x <- runif(n, -100, 100)
      expect_equal(robscale::adm(x), revss::adm(x), tolerance = tol)
      if (compare_iterative) {
        expect_equal(robscale::robLoc(x), revss::robLoc(x), tolerance = tol)
        expect_equal(robscale::robScale(x), revss::robScale(x), tolerance = tol)
      }
    }
  }
})

test_that("robLoc reproduces golden reference values", {
  # Golden values computed with robscale v0.2.1, verified against revss v2.0.0.
  # Decoupled from upstream revss — catches regressions in our own algorithm.
  golden <- c(
     82.96120869927108288,
     28.34910377860069275,
     31.39845808036625385,
     41.01295680738985538,
     -7.54143549129366875,
    -30.60884221170773145,
     18.60681201118527639,
    -27.92670497395137019,
    -16.44768253778534017,
     44.61202669309833624,
     29.70421395501008277,
    -17.35552615062394466,
     52.85815256497575376,
    -14.24107006138247655,
    -35.64729924577549980,
     14.65914841827957105,
     39.31696859241033337,
    -17.61513735996397756,
    -17.54360035055339750,
     38.01664944755422226,
     14.53988554388956445,
     -8.25708276819363007,
      7.57430703442530628,
     -0.40615615866561922,
     13.54013021964397723,
      3.76169471177400316,
      0.52287570009941731,
    -15.14688748187069933,
    -20.46209027590011331,
      7.12914617182963894,
    -15.89020205088178095,
     -3.45778719631556619,
     46.05800347891764801,
     29.42221170015346132,
    -17.73478756001277290,
     -2.04266042679924276,
     -8.02465194351456468,
    -18.60825557960116683,
    -31.00009231623163330,
     -3.54494254953130428,
    -13.53466281911471469,
      2.73649562863284235,
    -11.36847055762331848,
      9.09189576030746416,
      3.77623708249237433,
      7.61408334865766712,
      1.17751694394114814,
      3.06425511295394948,
     -0.55113303590639540,
    -22.80421680971541321,
     -5.12935674967489597,
     22.28807344057985418,
      3.53658571710340253,
     -7.85758346205312197,
      1.49903758730641412,
     -0.84182850963680345,
    -12.21142816372030993,
    -32.19388410119852750,
    -22.97860853195781417,
      4.66338270452734083,
    -23.77788324172414747,
     -0.60613787330851632,
     31.70928530303438819,
    -10.25573953706127028,
      0.87882902485949321,
    -19.03436430048506267,
     10.13598728325839282,
     -6.31939983053498455,
     -4.49275911154324703,
    -15.66596525452547439,
     14.43699766434143150,
     28.32193377864704686,
    -30.08891918345009486,
      1.75138190786630021,
     -4.78367079774357595,
    -20.62943640264032297,
     -8.29530766507329353,
     17.40507000133415261,
     -9.25243749809783544,
    -13.05948525293304741,
    -12.50728211401600909,
      0.24700454009225337,
    -24.65607955212726310,
     -2.87494082318324473,
      0.85427544591204629,
     11.57162780055238116,
    -23.85327496830368688,
     28.83522969413644077,
     14.20892044129155707,
     -4.94908029037604535
  )
  set.seed(42)
  idx <- 1L
  for (n in 3:20) {
    for (rep in 1:5) {
      x <- runif(n, -100, 100)
      expect_equal(robscale::robLoc(x), golden[idx], tolerance = 1e-12,
                   label = paste0("robLoc n=", n, " rep=", rep))
      idx <- idx + 1L
    }
    for (rep in 6:100) runif(n, -100, 100)
  }
})

test_that("robScale converges to true fixed point within 2*sqrt(eps)", {
  # Correctness test: when robScale() converges tightly, its result must lie
  # within 2*sqrt(eps) of the true M-scale fixed point s*, computed to near-
  # machine-precision by the same iterative formula (tol = 1e-14, maxit = 10000).
  #
  # Tolerance derivation (exact, via mean-value theorem on g(s) = s * v(s)):
  #   Let r = g'(s*) be the local contraction rate and tol_impl = sqrt(eps).
  #   At the stopping point robScale has |v_prev - 1| <= tol_impl, and
  #   fp_residual(s_impl) = |v(s_impl) - 1| ≈ r * |v_prev - 1| ≈ r * tol_impl.
  #   Using fpr as a certificate of r: r ≤ fpr / tol_impl.
  #   Error bound (first-order): |s_impl - s*| / s* ≤ r/(1-r) * tol_impl.
  #   We skip when fpr > tol_impl/2, which certifies r ≤ 0.5 and thus
  #   error ≤ r/(1-r) * tol_impl ≤ 1.0 * tol_impl = sqrt(eps).
  #   tolerance = 2*sqrt(eps) absorbs O(r^2) corrections from the linearisation.
  #   Reference error (tol = 1e-14) is negligible.
  #
  # Skip conditions:
  #   (a) fpr > tol_impl  — robScale exhausted maxit without converging.
  #   (b) fpr > tol_impl/2 — contraction rate r > 0.5; error bound exceeds tol.
  #   Both are expected for pathologically slow inputs, not correctness failures.
  #   The skip criterion is algorithm-independent (no assumptions about Aitken).
  C     <- 0.37394112142347236
  MAD_C <- 1.482602218505602

  rob_scale_ref <- function(x, tol = 1e-14, maxit = 10000L) {
    center <- median(x)
    s <- MAD_C * median(abs(x - center))
    if (s == 0) return(0)
    for (k in seq_len(maxit)) {
      u <- (x - center) * (0.5 / (C * s))
      v <- sqrt(2 * mean(tanh(u)^2))
      s <- s * v
      if (abs(v - 1) <= tol) break
    }
    s
  }

  # fp_residual: |v(s) - 1|. Zero at s*; certificates r via fpr ≈ r * tol_impl.
  fp_residual <- function(x, s) {
    center <- median(x)
    u <- (x - center) * (0.5 / (C * s))
    abs(sqrt(2 * mean(tanh(u)^2)) - 1)
  }

  tol_impl <- sqrt(.Machine$double.eps)     # robScale's convergence criterion
  tol2     <- 2 * sqrt(.Machine$double.eps) # test tolerance (see derivation above)

  set.seed(42)
  for (n in c(4L, 8L, 16L, 32L, 64L, 100L, 500L)) {
    for (rep in 1:10) {
      x      <- rnorm(n)
      s_impl <- robscale::robScale(x)
      fpr    <- fp_residual(x, s_impl)
      # Skip (a): maxit reached without convergence.
      # Skip (b): r > 0.5; error bound exceeds the asserted tolerance.
      if (fpr > tol_impl / 2) next
      expect_equal(s_impl, rob_scale_ref(x), tolerance = tol2,
                   label = paste0("robScale n=", n, " rep=", rep))
    }
  }
})
