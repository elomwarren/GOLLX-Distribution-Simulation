# ==============================================================================
# STAT 702 : ADVANCED STATISTICAL COMPUTING
# Topic : GOLLX Simulation
# @Elom Kodjoh-Kpakpassou
# ==============================================================================

# Reproducibility
set.seed(9148)

# ==========
# GOLLX
# ==========

################
# PDF of GOLLX #
################
dgollx <- function(x, Lambda, Alpha, Beta) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(all(p > 0 & p < 1), Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)
    phi_x <- (0.5 * w * x^2) + (w * x) + 1
    Mu_x <- phi_x / exp(Lambda * x)
    PDF <- (
        ((Beta * Alpha * Lambda^2) / (Lambda + 1)) *
            (((1 + 0.5 * Lambda * x^2) *
                phi_x^(Beta - 1)) / (1 - Mu_x)^(1 - Beta * Alpha)) *
            exp(Lambda * Beta * x) *
            (
                (1 - Mu_x)^Alpha + Mu_x^Alpha
            )^(-Beta - 1)
    )
    return(PDF)
}

################
# CDF of GOLLX #
################
pgollx <- function(x, Lambda, Alpha, Beta) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(all(p > 0 & p < 1), Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)

    phi_x <- (0.5 * w * x^2) + (w * x) + 1
    Mu_x <- phi_x / exp(Lambda * x)
    CDF <- (
        ((1 - Mu_x)^Alpha + (Mu_x^Alpha))^(-Beta)
    ) * (1 - Mu_x)^(Alpha * Beta)
    return(CDF)
}

############################################
# Quantile function (inverse CDF) of GOLLX
############################################
qgollx <- function(p, Lambda, Alpha, Beta) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(all(p > 0 & p < 1), Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)

    # Routine to compute inverse CDF
    sapply(p, function(u) {
        # Step 1: Analytically invert the GOLL-G layer to get target Mu
        k <- (u^(-1 / Beta) - 1)^(1 / Alpha)
        M_target <- k / (1 + k)

        # Step 2: numerically solve Mu(x) = M_target
        # where Mu(x) = (0.5*w*x^2 + w*x + 1) * exp(-Lambda*x)
        obj <- function(x) {
            (0.5 * w * x^2 + w * x + 1) * exp(-Lambda * x) - M_target
        }

        # Bracket the root: obj(0) = 1 - M_target > 0, need obj(upper) < 0
        upper <- 1
        while (obj(upper) > 0) upper <- upper * 2

        stats::uniroot(
            obj,
            lower = 0,
            upper = upper,
            tol = .Machine$double.eps^0.5
        )$root
    })
}

# Trial
# qgollx(c(0.25, 0.50, 0.75), Lambda = 2, Alpha = 1.5, Beta = 0.8)

###########################
# Generation of GOLLX rvs
###########################
rgollx <- function(n, Lambda, Alpha, Beta) {
    qgollx(runif(n), Lambda, Alpha, Beta)
}

# Trial
# rgollx(n = 100, Lambda = 2, Alpha = 1.5, Beta = 0.8)
