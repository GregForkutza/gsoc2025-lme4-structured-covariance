---
title: "Design Document: Structured Covariance Matrices for `lme4`"
format: html
engine: knitr
bibliography: gsoc2025.bib
---

## Task 0: Class Structure

::: {.callout-note}
**BMB**: It may be worth thinking about several (potential) characteristics of a covariance-matrix class in addition to the mathematical structure of the covariance matrix itself:

- how is it parameterized? (e.g. cf. `nlme`'s `pdLogChol` vs `pdNatural`, which are two different parameterizations of general-positive-definite matrices; there are several more options here, see @pinheiroUnconstrained1996 (perhaps these are potential sub-classes of a class defined by structure alone? Check the class structure underneath `nlme`?)
- how are any auxiliary data passed to the model? For example, we could fit a phylogenetic model with `propto` by computing the correlation structure from the tree ourselves (see `ape::vcv.phylo`), or by passing a phylogenetic tree object and asking the package to compute the correlation for us. Similarly, spatial models can take their data by taking a distance matrix, or by taking `x`, `y` coordinates and computing Euclidean distances internally (or by taking an arbitrary set of coordinates and a distance function?). I find both `nlme` and `glmmTMB` a bit clunky in this respect: `nlme` hard-codes a Euclidean distance, while `glmmTMB` has a `numFactor` argument for packaging coordinate information with values. In `glmmTMB`, `ou` and `exp` both represent exponentially decaying correlation functions, with `ou` intended for 1D/time-series and `exp` intended for 2D/spatial problems ...

:::


### Hierarchical Class Structure

-   **StructuredCov**
    -   **IdentityStructuredCov**
        -   `diag` (diagonal, homogeneous): $\theta = \sigma^2$
        -   `homdiag` (diagonal, heterogeneous): $\theta = (\sigma_1^2, \sigma_2^2, \ldots)$
    -   **UnstructuredCov**
        -   `us` (unstructured): $\theta = \mathrm{Cholesky}(\Sigma)$
        -   `rr` (reduced rank): $\theta = \Lambda \Lambda^\top$, $\dim(\Lambda) = q \times r$

::: {.callout-note}
**BMB**: I wouldn't say reduced rank is "unstructured". There is really only one type of "unstructured" covariance matrix, i.e. one whose only constraint is that it is positive semidefinite.
:::
<!-- have to hack the list structure in order to have callouts in the middle ... -->

-
    -   **CorrelationStructuredCov**
        -   `cs` (compound symmetric): $\theta = (\sigma^2, \rho)$
        -   `hetcs` (heterogeneous compound symmetric): $\theta = \{\sigma_1^2, \ldots, \sigma_n^2\},\ \rho$
        -   `propto` (fixed correlation matrix): fixed correlation $+$ variance
        -   `gau` (Gaussian correlation): $\theta = \sigma^2 + \text{shape}$
::: {.callout-note}
**BMB**: I don't understand the justification for this category (maybe you could include some notes explaining the scope of each category?) (1) at least in `glmmTMB`, `gau()` is a spatial correlation structure --- although arguably you could use it (a Gaussian process with a 'squared exponential'/'Gaussian'/'radial basis function' kernel) in 1D, as a time-series model as well. In any case, it's distance-dependent in a way that doesn't naturally seem to fit with the other cases.
:::
        
        
-
    -   **TimeSeriesStructuredCov**
        -   **ARStructuredCov**
            -   `ar1` (autoregressive order 1, homogeneous): $\theta = (\sigma^2, \rho)$
            -   `hetar1` (autoregressive order 1, heterogeneous): $\theta = (\sigma_1^2, \ldots, \sigma_n^2, \rho)$
            -   `ante` (ante-dependence): $\theta = (\sigma^2, \rho_1, \ldots, \rho_k)$
        -   `ARMA` (autoregressive moving average): $\theta = (\text{AR}_1, \ldots, \text{AR}_p,\ \text{MA}_1, \ldots, \text{MA}_q)$
        -   `toep` (Toeplitz structure): $\theta =$ banded entries of $\Sigma$
        -   `mat` (Matérn kernel): $\theta =$ scale $+$ shape
        -   `ou` (Ornstein–Uhlenbeck process): $\theta = (\sigma^2, \phi)$
        
::: {.callout-note}
**BMB**: `toep` can be considered a time-series structure (I think `toep` is `ante`/AR with the maximum number of steps?), but is also often used for discrete spatial designs. `mat` should almost certainly be spatial rather than time-series (although see notes above about `ou` vs. `exp`, and `gau`, and 1- vs 2D Gaussian processes)
:::

-
    -   **SpatialStructuredCov**
        -   `exp` (exponential): $\theta = (\sigma^2, \phi)$
        -   `lin` (linear): $\theta = (\sigma^2, \phi)$
        -   `spher` (spherical): $\theta = (\sigma^2, \phi)$
        -   `ratquad` (rational quadratic): $\theta = (\sigma^2, \alpha, \ell)$
    -   **PhyloStructuredCov**
        -   `brownian` (Brownian motion): needs tree
        -   `pagel` (Pagel's $\lambda$ transformation): tree $+$ $\lambda$
        -   `blomberg` (Blomberg's $K$ statistic): tree $+$ $K$
        -   `grafen` (Grafen’s transformation): tree $+$ $\alpha$

### Subclass Analysis

#### `SpatialStructuredCov`

##### Mathematical Structure of `SpatialStructuredCov`

All spatial covariance structures share a common mathematical form:

$$
\text{Cov}(x_i, x_j) = \sigma^2 \cdot \rho_\theta(d_{ij})
$$

where $d_{ij}$ is the distance between observations $x_i$ and $x_j$, $\sigma^2$ is a marginal variance, and $\rho_\theta$ is a structure-specific kernel function parameterized by one or more values in the vector $\theta$. Each subclass of `SpatialStructuredCov`, such as exponential, linear, spherical, or rational quadratic, corresponds to a different form of $\rho_\theta$.

Despite their differences in functional form, all spatial covariance structures share assumptions. First, they assume stationarity, meaning that correlations between observations depend only on the distance between them. Second, they are typically isotropic, so that the correlation depends only on the magnitude of the distance (e.g., Euclidean norm) and not its direction. Finally, there is a degree of parameter homogeneity: nearly all spatial structures rely on a scale or range parameter $\phi$, often in conjunction with a variance parameter $\sigma^2$. Some structures, such as the rational quadratic kernel, also introduce additional shape parameters like $\alpha$ or $\ell$.

::: {.callout-note}
**BMB**: pretty much everything we're considering will assume stationarity. (There are some 'adaptive' splines in `mgcv` that allow the smoothing to vary as a function of `x`, but I don't think we're going there.) Isotropy is a good point to make, although obviously we *could* construct anisotropic versions of these models if we wanted to (with rotation/scaling parameters for the `x` and `y` coordinates ...)
:::

##### Constructor Logic and Validation

###### S3 Representation of `SpatialStructuredCov`

In an S3 system, we can define an exponential covariance structure as a plain list. The kernel function is stored directly, and a character vector lists the required parameters. We simulate inheritance by giving the object multiple class tags.

```{r s3-spatial}
exp_struct <- structure(
  list(
    kernelFn = function(d, theta) exp(-d / theta["phi"]),  # Defines how correlations decay with distance
    paramNames = c("sigma", "phi")                         # Declares which parameters the kernel expects
  ),
  class = c("exp", "SpatialStructuredCov", "StructuredCov") # Class tags for method dispatch
)
```

S3 has no built-in validity mechanism. Validation must be done manually by a constructor function:

```{r s3-spatial-valid}
make_spatial_cov <- function(type = "exp", theta = list(), kernelFn = NULL, paramNames = NULL, distanceFn = NULL) {
  # Validate required parameters
    test_d <- matrix(c(0, 1, 1, 0), nrow = 2)  # symmetric distance matrix
    test_theta <- list(phi = 1, sigma = 1)

  if (!is.function(kernelFn)) stop("kernelFn must be a valid function")

    test_kernel <- try(kernelFn(test_d, test_theta), silent = TRUE)
  if (inherits(test_kernel, "try-error") || !is.numeric(test_kernel)) {
      stop("kernelFn must accept (d, theta) and return a numeric matrix")
    }

  if (!is.null(distanceFn)) {
    test_mat <- distanceFn(matrix(c(0, 0), nrow = 1))  # minimal test input
    if (!is.matrix(test_mat) || !isSymmetric(test_mat)) stop("distanceFn must return a square symmetric matrix")
  }
  
  if (is.null(paramNames)) stop("paramNames must be provided")
  missing_params <- setdiff(paramNames, names(theta))
  if (length(missing_params) > 0) {
    stop(sprintf("theta must contain parameters: %s", paste(missing_params, collapse = ", ")))
  }
  
  # Parameter-specific constraints
  if (!is.null(theta$phi) && theta$phi <= 0) stop("phi must be > 0")
  if (!is.null(theta$sigma) && theta$sigma < 0) stop("sigma must be >= 0")
  if (!is.null(theta$alpha) && theta$alpha <= 0) stop("alpha must be > 0")
  if (!is.null(theta$ell) && theta$ell <= 0) stop("ell must be > 0")

  structure(
    list(
      theta = theta,
      type = type,
      kernelFn = kernelFn,
      paramNames = paramNames,
      distanceFn = distanceFn
    ),
    class = c(type, "SpatialStructuredCov", "StructuredCov")
  )
}
```

::: {.callout-note}
**BMB**: I might try to avoid repeating by writing functions like this:
```{r assert-function}
strip_parent <- function(x) gsub("^[^$]\\$", "", x)
check_positive <- function(x) {
  nm <- deparse(substitute(x)) |> strip_parent()
  if (!is.null(x) && x <= 0) stop(sprintf("%s must be > 0", nm),
                                  .call = FALSE)
}
```
(and similarly for `check_nonnegative`, etc.) (PS we should probably avoid pipes in production; in principle `lme4` still works with R versions < 4 and we don't want to break back-compatibility more than necessary)
:::

###### S4 Representation of `SpatialStructuredCov`

In an S4 system, we define a virtual class `SpatialStructuredCov` that extends a general `StructuredCov` base class. This class declares formal slots for the kernel function and parameter names, ensuring every instance has the expected components. Specific structures (like exponential, linear or spherical) then subclass this definition.

```{r s4-basics}
library(stats4)

setClass("StructuredCov", slots = list()) 

setClass("SpatialStructuredCov",
         contains = "StructuredCov",
         slots = list(
           kernelFn = "function",         # Defines how correlation decays with distance
           paramNames = "character",      # Declares which parameters the kernel expects
           theta = "list",                # Stores parameter values like sigma, phi, etc.
           distanceFn = "functionOrNULL"  # Optional user-supplied distance function
         ))

# Helper class union to allow NULL as a valid value for distanceFn
setClassUnion("functionOrNULL", c("function", "NULL"))

setClass("ExpCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = function(d, theta) exp(-d / theta[["phi"]]),
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

setClass("LinCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = function(d, theta) pmax(1 - d / theta[["phi"]], 0),
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

setClass("SpherCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = function(d, theta) {
             r <- d / theta[["phi"]]
             out <- ifelse(r < 1, 1 - 1.5 * r + 0.5 * r^3, 0)
             return(out)
           },
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

setClass("RatQuadCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = function(d, theta) {
             alpha <- theta[["alpha"]]
             ell <- theta[["ell"]]
             (1 + (d^2) / (2 * alpha * ell^2))^(-alpha)
           },
           paramNames = c("sigma", "alpha", "ell"),
           theta = list(sigma = 1, alpha = 1, ell = 1),
           distanceFn = NULL
         ))
```

``` r
setValidity("SpatialStructuredCov", function(object) {
  msg <- character()

  # Check paramNames are present in theta
  missing_params <- setdiff(object@paramNames, names(object@theta))
  if (length(missing_params) > 0) {
    msg <- c(msg, paste("Missing parameters in theta:", paste(missing_params, collapse = ", ")))
  }

  # Validate each parameter constraint if present
  if (!is.null(object@theta[["phi"]]) && object@theta[["phi"]] <= 0) {
    msg <- c(msg, "phi must be > 0")
  }

  if (!is.null(object@theta[["sigma"]]) && object@theta[["sigma"]] < 0) {
    msg <- c(msg, "sigma must be >= 0")
  }

  if (!is.null(object@theta[["alpha"]]) && object@theta[["alpha"]] <= 0) {
    msg <- c(msg, "alpha must be > 0")
  }

  if (!is.null(object@theta[["ell"]]) && object@theta[["ell"]] <= 0) {
    msg <- c(msg, "ell must be > 0")
  }

  # Optionally test the kernel function
  test_d <- matrix(c(0, 1, 1, 0), nrow = 2)
  test_theta <- object@theta
  test_kernel <- try(object@kernelFn(test_d, test_theta), silent = TRUE)

  if (inherits(test_kernel, "try-error") || !is.numeric(test_kernel)) {
    msg <- c(msg, "kernelFn must accept (d, theta) and return a numeric matrix")
  }

  # Optional distanceFn test if present
  if (!is.null(object@distanceFn)) {
    test_input <- matrix(c(0, 0), nrow = 1)
    test_mat <- try(object@distanceFn(test_input), silent = TRUE)
    if (inherits(test_mat, "try-error") || !is.matrix(test_mat) || !isSymmetric(test_mat)) {
      msg <- c(msg, "distanceFn must return a square symmetric matrix")
    }
  }

  if (length(msg) == 0) TRUE else msg
})
```

##### Inheritable Methods

Spatial covariance structures differ primarily in the form of their kernel function, but they otherwise share a common computational interface.

###### Shared Functional Requirements

All spatial structures should support: 

- `makeTheta()`: constructs the initial parameter vector $\theta$. 
- `nTheta()`: returns the number of free parameters. 
- `makeLambda(structure, theta, distmat)`: builds the relative covariance factor $\Lambda_\theta$ using the kernel function and a distance matrix. 
- `validate(structure)`: checks that all parameter values are within valid bounds.

These methods can be defined once for the `SpatialStructuredCov` superclass and reused by all subclasses, except when a structure has special behavior (e.g., `ratquad` has more parameters than `exp`).

###### S3 Representation of `SpatialStructuredCov`

In the S3 system, method dispatch is controlled using `UseMethod()` and class vectors. We define generics and then create default methods for `"SpatialStructuredCov"`.

```r
# Generic for makeTheta
makeTheta <- function(structure, ...) {
  UseMethod("makeTheta", structure)
}

# Default method for spatial structures
makeTheta.SpatialStructuredCov <- function(structure) {
  # Create named vector with default values
  theta <- setNames(rep(1, length(structure$paramNames)), structure$paramNames)
  return(theta)
}

# Generic and method for nTheta
nTheta <- function(structure, ...) {
  UseMethod("nTheta", structure)
}

nTheta.SpatialStructuredCov <- function(structure) {
  length(structure$paramNames)
}

# Generic and default makeLambda
makeLambda <- function(structure, theta, distmat) {
  UseMethod("makeLambda", structure)
}

makeLambda.SpatialStructuredCov <- function(structure, theta, distmat) {
  kernelMat <- structure$kernelFn(distmat, theta)
  chol(kernelMat)  # returns Cholesky factor 
}
```

###### S4 Representation of `SpatialStructuredCov`

In the S4 system, we define formal generics and provide superclass-level methods that apply to any subclass of SpatialStructuredCov.

``` r
# Define S4 generics
setGeneric("makeTheta", function(structure) standardGeneric("makeTheta"))
setGeneric("nTheta", function(structure) standardGeneric("nTheta"))
setGeneric("makeLambda", function(structure, theta, distmat) standardGeneric("makeLambda"))

# Default method for all SpatialStructuredCov objects
setMethod("makeTheta", "SpatialStructuredCov", function(structure) {
  theta <- setNames(rep(1, length(structure@paramNames)), structure@paramNames)
  return(theta)
})

setMethod("nTheta", "SpatialStructuredCov", function(structure) {
  length(structure@paramNames)
})

setMethod("makeLambda", "SpatialStructuredCov", function(structure, theta, distmat) {
  kernelMat <- structure@kernelFn(distmat, theta)
  chol(kernelMat)
})
```

##### Validation Testing
You can run each script to test the class structure, parameter validation, and kernel behavior for both the S3 and S4 implementations.

###### S3
```r
#--- Constructor and validation function ---#
make_spatial_cov <- function(type = "exp", theta = list(), kernelFn = NULL, paramNames = NULL, distanceFn = NULL) {
  test_d <- matrix(c(0, 1, 1, 0), nrow = 2)
  test_theta <- switch(
    type,
    exp     = list(phi = 1, sigma = 1),
    ratquad = list(sigma = 1, alpha = 1, ell = 1),
    stop("Unsupported covariance type")
  )
  
  if (!is.function(kernelFn)) stop("kernelFn must be a valid function")
  
  test_kernel <- try(kernelFn(test_d, test_theta), silent = TRUE)
  if (inherits(test_kernel, "try-error") || !is.numeric(test_kernel)) {
    stop("kernelFn must accept (d, theta) and return a numeric matrix")
  }
  
  if (!is.null(distanceFn)) {
    test_mat <- distanceFn(matrix(c(0, 0), nrow = 1))
    if (!is.matrix(test_mat) || !isSymmetric(test_mat)) stop("distanceFn must return a square symmetric matrix")
  }
  
  if (is.null(paramNames)) stop("paramNames must be provided")
  missing_params <- setdiff(paramNames, names(theta))
  if (length(missing_params) > 0) {
    stop(sprintf("theta must contain parameters: %s", paste(missing_params, collapse = ", ")))
  }
  
  if (!is.null(theta$phi) && theta$phi <= 0) stop("phi must be > 0")
  if (!is.null(theta$sigma) && theta$sigma < 0) stop("sigma must be >= 0")
  if (!is.null(theta$alpha) && theta$alpha <= 0) stop("alpha must be > 0")
  if (!is.null(theta$ell) && theta$ell <= 0) stop("ell must be > 0")
  
  structure(
    list(
      theta = theta,
      type = type,
      kernelFn = kernelFn,
      paramNames = paramNames,
      distanceFn = distanceFn
    ),
    class = c("SpatialStructuredCov", "StructuredCov", type)
  )
}

#--- Shared kernel function definitions ---#
ratquad_kernel <- function(d, theta) {
  alpha <- theta[["alpha"]]
  ell   <- theta[["ell"]]
  (1 + (d^2) / (2 * alpha * ell^2))^(-alpha)
}

exp_kernel <- function(d, theta) {
  exp(-d / theta[["phi"]])
}

#--- S3 method definitions ---#

makeTheta <- function(structure, ...) UseMethod("makeTheta")
nTheta    <- function(structure, ...) UseMethod("nTheta")
makeLambda <- function(structure, theta, distmat, ...) UseMethod("makeLambda")
makeCovMatrix <- function(structure, theta, distmat, ...) UseMethod("makeCovMatrix")

makeTheta.SpatialStructuredCov <- function(structure, ...) {
  setNames(rep(1, length(structure$paramNames)), structure$paramNames)
}

nTheta.SpatialStructuredCov <- function(structure, ...) {
  length(structure$paramNames)
}

makeLambda.SpatialStructuredCov <- function(structure, theta, distmat, ...) {
  kernelMat <- structure$kernelFn(distmat, theta)
  if (!isSymmetric(kernelMat)) stop("Kernel matrix must be symmetric")
  chol(kernelMat)
}

makeCovMatrix.SpatialStructuredCov <- function(structure, theta, distmat, ...) {
  sigma2 <- theta[["sigma"]]
  kernelMat <- structure$kernelFn(distmat, theta)
  if (!isSymmetric(kernelMat)) stop("Kernel matrix must be symmetric")
  sigma2 * kernelMat
}

#--- Post-construction validator ---#
validate.SpatialStructuredCov <- function(object) {
  result <- try(make_spatial_cov(
    type = object$type,
    theta = object$theta,
    kernelFn = object$kernelFn,
    paramNames = object$paramNames,
    distanceFn = object$distanceFn
  ), silent = TRUE)
  
  if (inherits(result, "try-error")) {
    return(conditionMessage(attr(result, "condition")))
  } else TRUE
}


#--- Test cases ---#
# Example spatial structures 
exp_struct <- make_spatial_cov(
  type = "exp",
  theta = list(sigma = 1, phi = 1),
  kernelFn = exp_kernel,  
  paramNames = c("sigma", "phi")
)

ratquad_struct <- make_spatial_cov(
  type = "ratquad",
  theta = list(sigma = 1, alpha = 1, ell = 1),
  kernelFn = ratquad_kernel,
  paramNames = c("sigma", "alpha", "ell")
)


# Kernel matrix input
dmat <- matrix(c(0, 1, 1, 0), nrow = 2)

# Print results
print(makeTheta(exp_struct))       # Expected: sigma=1, phi=1
print(nTheta(exp_struct))          # Expected: 2
print(makeLambda(exp_struct, exp_struct$theta, dmat))

print(makeTheta(ratquad_struct))   # Expected: sigma=1, alpha=1, ell=1
print(nTheta(ratquad_struct))      # Expected: 3
print(makeLambda(ratquad_struct, ratquad_struct$theta, dmat))

print(makeCovMatrix(exp_struct, exp_struct$theta, dmat))
print(makeCovMatrix(ratquad_struct, ratquad_struct$theta, dmat))

print(validate.SpatialStructuredCov(exp_struct))     
print(validate.SpatialStructuredCov(ratquad_struct)) 


# TESTING MANUAL VALIDATION FAILURES  

# Missing required param (sigma)
print(try(make_spatial_cov(
  type = "exp",
  theta = list(phi = 1),
  kernelFn = function(d, theta) exp(-d / theta[["phi"]]),
  paramNames = c("sigma", "phi")
), silent = TRUE))

# Negative phi
print(try(make_spatial_cov(
  type = "exp",
  theta = list(sigma = 1, phi = -1),
  kernelFn = function(d, theta) exp(-d / theta[["phi"]]),
  paramNames = c("sigma", "phi")
), silent = TRUE))

# Negative sigma
print(try(make_spatial_cov(
  type = "exp",
  theta = list(sigma = -1, phi = 1),
  kernelFn = function(d, theta) exp(-d / theta[["phi"]]),
  paramNames = c("sigma", "phi")
), silent = TRUE))

# Negative alpha
print(try(make_spatial_cov(
  type = "ratquad",
  theta = list(sigma = 1, alpha = -1, ell = 1),
  kernelFn = function(d, theta) {
    alpha <- theta[["alpha"]]
    ell <- theta[["ell"]]
    (1 + (d^2) / (2 * alpha * ell^2))^(-alpha)
  },
  paramNames = c("sigma", "alpha", "ell")
), silent = TRUE))

# Negative ell
print(try(make_spatial_cov(
  type = "ratquad",
  theta = list(sigma = 1, alpha = 1, ell = -0.1),
  kernelFn = function(d, theta) {
    alpha <- theta[["alpha"]]
    ell <- theta[["ell"]]
    (1 + (d^2) / (2 * alpha * ell^2))^(-alpha)
  },
  paramNames = c("sigma", "alpha", "ell")
), silent = TRUE))

# Invalid kernelFn (not a function)
print(try(make_spatial_cov(
  type = "exp",
  theta = list(sigma = 1, phi = 1),
  kernelFn = "not_a_function",
  paramNames = c("sigma", "phi")
), silent = TRUE))

# Invalid kernelFn (wrong output type)
print(try(make_spatial_cov(
  type = "exp",
  theta = list(sigma = 1, phi = 1),
  kernelFn = function(d, theta) "not a matrix",
  paramNames = c("sigma", "phi")
), silent = TRUE))

# Invalid distanceFn (not symmetric)
print(try(make_spatial_cov(
  type = "exp",
  theta = list(sigma = 1, phi = 1),
  kernelFn = function(d, theta) exp(-d / theta[["phi"]]),
  paramNames = c("sigma", "phi"),
  distanceFn = function(x) matrix(1:3, nrow = 1)  # not square/symmetric
), silent = TRUE))

```
###### S4
```r
library(methods)
library(Matrix)
library(stats4)

# Allow optional distanceFn slot
setClassUnion("functionOrNULL", c("function", "NULL"))

# Base and spatial class
setClass("StructuredCov", slots = list())

setClass("SpatialStructuredCov",
         contains = "StructuredCov",
         slots = list(
           kernelFn = "function",
           paramNames = "character",
           theta = "list",
           distanceFn = "functionOrNULL"
         ))

# Kernel functions
exp_kernel <- function(d, theta) exp(-d / theta[["phi"]])

ratquad_kernel <- function(d, theta) {
  alpha <- theta[["alpha"]]
  ell <- theta[["ell"]]
  (1 + (d^2) / (2 * alpha * ell^2))^(-alpha)
}

lin_kernel <- function(d, theta) {
  pmax(1 - d / theta[["phi"]], 0)
}

spher_kernel <- function(d, theta) {
  r <- d / theta[["phi"]]
  ifelse(r < 1, 1 - 1.5 * r + 0.5 * r^3, 0)
}

# Spatial structure subclasses
setClass("ExpCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = exp_kernel,
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

setClass("RatQuadCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = ratquad_kernel,
           paramNames = c("sigma", "alpha", "ell"),
           theta = list(sigma = 1, alpha = 1, ell = 1),
           distanceFn = NULL
         ))

setClass("LinCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = lin_kernel,
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

setClass("SpherCov",
         contains = "SpatialStructuredCov",
         prototype = list(
           kernelFn = spher_kernel,
           paramNames = c("sigma", "phi"),
           theta = list(sigma = 1, phi = 1),
           distanceFn = NULL
         ))

# Validity checker for all SpatialStructuredCovs
setValidity("SpatialStructuredCov", function(object) {
  msg <- character()
  missing_params <- setdiff(object@paramNames, names(object@theta))
  if (length(missing_params) > 0)
    msg <- c(msg, paste("Missing parameters in theta:", paste(missing_params, collapse = ", ")))
  if (!is.null(object@theta[["phi"]]) && object@theta[["phi"]] <= 0)
    msg <- c(msg, "phi must be > 0")
  if (!is.null(object@theta[["sigma"]]) && object@theta[["sigma"]] < 0)
    msg <- c(msg, "sigma must be >= 0")
  if (!is.null(object@theta[["alpha"]]) && object@theta[["alpha"]] <= 0)
    msg <- c(msg, "alpha must be > 0")
  if (!is.null(object@theta[["ell"]]) && object@theta[["ell"]] <= 0)
    msg <- c(msg, "ell must be > 0")

  test_d <- matrix(c(0, 1, 1, 0), nrow = 2)
  test_theta <- object@theta
  test_kernel <- try(object@kernelFn(test_d, test_theta), silent = TRUE)
  if (inherits(test_kernel, "try-error") || !is.numeric(test_kernel))
    msg <- c(msg, "kernelFn must accept (d, theta) and return a numeric matrix")

  if (!is.null(object@distanceFn)) {
    test_input <- matrix(c(0, 0), nrow = 1)
    test_mat <- try(object@distanceFn(test_input), silent = TRUE)
    if (inherits(test_mat, "try-error") || !is.matrix(test_mat) || !isSymmetric(test_mat))
      msg <- c(msg, "distanceFn must return a square symmetric matrix")
  }

  if (length(msg) == 0) TRUE else msg
})

# S4 generics
setGeneric("makeLambda", function(structure, theta, distmat) standardGeneric("makeLambda"))
setGeneric("makeCovMatrix", function(structure, theta, distmat) standardGeneric("makeCovMatrix"))

# Default methods for SpatialStructuredCov
setMethod("makeLambda", "SpatialStructuredCov", function(structure, theta, distmat) {
  kernelMat <- structure@kernelFn(distmat, theta)
  if (!isSymmetric(kernelMat)) stop("Kernel matrix must be symmetric")
  chol(kernelMat)
})

setMethod("makeCovMatrix", "SpatialStructuredCov", function(structure, theta, distmat) {
  sigma2 <- theta[["sigma"]]
  kernelMat <- structure@kernelFn(distmat, theta)
  if (!isSymmetric(kernelMat)) stop("Kernel matrix must be symmetric")
  sigma2 * kernelMat
})

#--- Run validations and tests ---#
dmat <- matrix(c(0, 1, 1, 0), nrow = 2)

# Valid objects
exp_obj <- new("ExpCov")
rat_obj <- new("RatQuadCov")
lin_obj <- new("LinCov")
spher_obj <- new("SpherCov")

print(validObject(exp_obj))  
print(validObject(rat_obj))  
print(validObject(lin_obj))  
print(validObject(spher_obj))

# Output tests
print(exp_obj@kernelFn(dmat, exp_obj@theta))  
print(rat_obj@kernelFn(dmat, rat_obj@theta))  
print(lin_obj@kernelFn(dmat, lin_obj@theta))  
print(spher_obj@kernelFn(dmat, spher_obj@theta))  

print(makeLambda(exp_obj, exp_obj@theta, dmat))
print(makeCovMatrix(exp_obj, exp_obj@theta, dmat))

print(makeLambda(rat_obj, rat_obj@theta, dmat))
print(makeCovMatrix(rat_obj, rat_obj@theta, dmat))

print(makeLambda(lin_obj, lin_obj@theta, dmat))
print(makeCovMatrix(lin_obj, lin_obj@theta, dmat))

print(makeLambda(spher_obj, spher_obj@theta, dmat))
print(makeCovMatrix(spher_obj, spher_obj@theta, dmat))

# Failure tests
exp_obj2 <- new("ExpCov")
exp_obj2@theta <- list(sigma = 1, phi = -0.5)
print(try(validObject(exp_obj2), silent = TRUE))  # phi < 0

exp_obj3 <- new("ExpCov")
exp_obj3@theta <- list(phi = 1)
print(try(validObject(exp_obj3), silent = TRUE))  # missing sigma

rat_obj2 <- new("RatQuadCov")
rat_obj2@theta <- list(sigma = 1, alpha = -0.5, ell = 1)
print(try(validObject(rat_obj2), silent = TRUE))

lin_bad <- new("LinCov")
lin_bad@theta <- list(sigma = 1, phi = -1)
print(try(validObject(lin_bad), silent = TRUE))

spher_bad <- new("SpherCov")
spher_bad@theta <- list(phi = 1)
print(try(validObject(spher_bad), silent = TRUE))

bad_dist_fn <- function(x) matrix(c(1, 2, 3), nrow = 1)  # Not symmetric
bad_dist <- new("ExpCov")
bad_dist@distanceFn <- bad_dist_fn
print(try(validObject(bad_dist), silent = TRUE))
```

## references
