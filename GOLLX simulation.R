# ==============================================================================
# STAT 702 : ADVANCED STATISTICAL COMPUTING
# Topic : GOLLX Simulation
# Elom L. Warren Kodjoh-Kpakpassou
# @ elwkodjoh-kpakpassou@st.ug.edu.gh
# ==============================================================================

# Reproducibility
set.seed(9148)

# ==========
# GOLLX
# ==========

################
# PDF of GOLLX #
################
dgollx <- function(x, Lambda, Alpha = 9, Beta = 0.51) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(x > 0, Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)
    phi_x <- (0.5 * w * Lambda * x^2) + (w * x) + 1
    Mu_x <- phi_x * exp(-Lambda * x)
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
pgollx <- function(x, Lambda, Alpha = 9, Beta = 0.51) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(x > 0, Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)

    phi_x <- (0.5 * w * Lambda * x^2) + (w * x) + 1
    Mu_x <- phi_x * exp(-Lambda * x)
    CDF <- (
        ((1 - Mu_x)^Alpha + (Mu_x^Alpha))^(-Beta)
    ) * (1 - Mu_x)^(Alpha * Beta)
    return(CDF)
}

############################################
# Quantile function (inverse CDF) of GOLLX
############################################
qgollx <- function(p, Lambda, Alpha = 9, Beta = 0.51) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(all(p > 0 & p < 1), Lambda > 0, Alpha > 0, Beta > 0)

    w <- Lambda / (Lambda + 1)

    # Routine to compute inverse CDF
    sapply(p, function(u) {
        # Step 1: Analytically invert the GOLL-G layer to get target Mu
        k <- (u^(-1 / Beta) - 1)^(1 / Alpha)
        M_target <- k / (1 + k)

        # Step 2: numerically solve Mu(x) = M_target
        # where Mu(x) = (0.5*w*Lambda*x^2 + w*x + 1) * exp(-Lambda*x)
        obj <- function(x) {
            (0.5 * w * Lambda * x^2 + w * x + 1) * exp(-Lambda * x) - M_target
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
rgollx <- function(n, Lambda, Alpha = 9, Beta = 0.51) {
    # n > 0
    stopifnot(n > 0)
    qgollx(runif(n), Lambda, Alpha, Beta)
}

# Trial
# rgollx(n = 100, Lambda = 2, Alpha = 1.5, Beta = 0.8)

###########################################################
# GOLLX Negative log-likelihood for MLE optimization in R #
###########################################################
negllgollx <- function(Lambda, Alpha = 9, Beta = 0.51, x) {
    -sum(dgollx(x, Lambda, Alpha, Beta))
}

# ==========
# XGAMMA
# ==========


#################
# PDF of XGAMMA #
#################

dxgamma <- function(x, Lambda) {
    # x > 0, Lambda > 0
    stopifnot(x > 0, Lambda > 0)

    ((Lambda + 1)^(-1)) * exp(-Lambda * x) * (1 + 0.5 * Lambda * x^2) * Lambda^2
}

#################
# CDF of XGAMMA #
#################
pxgamma <- function(x, Lambda) {
    # x > 0, Lambda > 0
    stopifnot(x > 0, Lambda > 0)

    w <- Lambda / (Lambda + 1)
    phi_x <- (0.5 * w * x^2) + (w * x) + 1

    CDF <- 1 - exp(-Lambda * x) * phi_x
    return(CDF)
}

############################################
# Quantile function (inverse CDF) of XGAMMA
############################################
qgamma <- function(p, Lambda) {
    # x > 0, Lambda > 0, Alpha, Beta > 0
    stopifnot(all(p > 0 & p < 1), Lambda > 0)

    w <- Lambda / (Lambda + 1)
    # Routine to compute inverse CDF
    sapply(p, function(u) {
        # Define target for inversion
        target <- u - 1

        # numerically solve Mu(x) = M_target
        # where Mu(x) = (0.5*w*Lambda*x^2 + w*x + 1) * exp(-Lambda*x)
        obj <- function(x) {
            (0.5 * w * Lambda * x^2 + w * x + 1) * exp(-Lambda * x) - target
        }

        # Bracket the root: obj(0) = 1 - M_target > 0, need obj(upper) < 0
        upper <- 1
        while (obj(upper) < 0) upper <- upper * 2

        stats::uniroot(
            obj,
            lower = 0,
            upper = upper,
            tol = .Machine$double.eps^0.5
        )$root
    })
}

###########################
# Generation of XGAMMA rvs
###########################
rxgamma <- function(n, Lambda) {
    # n > 0, Lambda > 0
    stopifnot(n > 0, Lambda > 0)

    t <- Lambda / (Lambda + 1)
    U <- runif(n)
    V <- rexp(n, rate = Lambda)
    W <- rgamma(n, shape = 3, rate = Lambda)
    XGAMMA <- ifelse(U <= t, V, W)
    return(XGAMMA)
}

############################################################
# XGAMMA Negative log-likelihood for MLE optimization in R #
############################################################
negllxgamma <- function(Lambda, x) {
    -sum(dgamma(x, Lambda))
}

# =============================
# SIMULATION HELPER FUNCTIONS
# =============================

###########################
# MLE 1D optimization helper
###########################
FitMLE <- function(x, llfn, Start = 1 / mean(x), upper = 1e6) {
    # x :: vector, of data
    # llfn :: function string, e.g. llfn = "negllxgamma"

    fit <- optim(
        par = Start,
        fn = match.fun(llfn),
        x = x,
        method = "Brent",
        lower = 0,
        upper = upper,
        hessian = FALSE
    )
    if (fit$convergence != 0) {
        warning("optim did not converge: ", fit$message)
    }
    fit$par
}

######################
# BIAS & MSE helper
######################
BiasMSE <- function(estim, param) {
    # estim :: vector, estimates
    # param :: scalar, parameter value
    Bias <- mean(estim - param)
    MSE <- mean((estim - param)^2)
    tab <- round(c(Bias, MSE), 4)
    names(tab) <- c("Bias", "MSE")
    return(tab)
}

##################################
# Batch MLE across sample matrix
##################################
FitMLEBatch <- function(samples, llfns) {
    # samples :: matrix, columns = independent samples/replicates
    # llfns   :: named list or character vector of function-name strings,
    #            e.g. c(gollx = "negllgollx", xgamma = "negllxgamma")
    if (is.null(names(llfns))) {
        names(llfns) <- llfns
    }
    lapply(llfns, function(fn) {
        apply(samples, MARGIN = 2, function(X) FitMLE(X, llfn = fn))
    })
}

#####################################
# Generate R samples of size n
#####################################
GenerateSamples <- function(R, n, param, rvfn) {
    # R     :: number of repetitions/replicates
    # n     :: sample size per replicate
    # param :: scalar, parameter value for the RNG function
    # rvfn  :: function string, e.g. rvfn = "rexp"
    replicate(
        n = R,
        match.fun(rvfn)(n, param)
    )
}

########################################
# MONTE CARLO SIMULATION HELPER
########################################
MonteCarloSim <- function(R, n,
                          param, rvfn = "rgollx",
                          llfns = c(gollx = "negllgollx", xgamma = "negllxgamma")) {
    # R :: Number of Repetitions
    # n :: sample size
    # param :: scalar, value of parameter assessed, 1D parameter
    # rvfn :: function to generate random variates, e.g. rvfn = "rexp"

    # Generate R random samples of size n from distribution of interest
    samples <- GenerateSamples(R, n, param, rvfn)

    # Maximum Likelihood Estimates for generated samples
    param.hat <- FitMLEBatch(samples, llfns)

    # param.hat$gollx    # vector of lambda.hat across columns, GOLLX
    # param.hat$xgamma   # vector of lambda.hat across columns, XGAMMA

    # Compute Bias and MSE
    lapply(
        param.hat,
        function(par) BiasMSE(estim = par, param = param)
    )
}

# ==============================================================================
# SIMULATION
# ==============================================================================
# Simulation parameters
R <- 10000
n <- 100
Lambda <- 0.5

MonteCarloSim(R, n, Lambda)
