
# Designing Class Structures for Structured Covariance in lme4: S3 vs S4
This document compares the S3 and S4 object systems in R, focusing on how they support structured covariance matrices in `lme4`. The project aims to implement a variety of structures, such as diagonal, compound symmetric, autoregressive, and user-defined forms, that differ in internal logic but require a common interface. Each must support methods for constructing the relative covariance factor (`makeLambda`), accessing and updating parameters (`makeTheta`, `updateTheta`), and optionally building the full covariance matrix (`makeSigma`). The central design question is how to organize these components cleanly and extensibly using S3 or S4.

## S3 Class System

S3 lets us tag objects with a class and use naming patterns to run the right method. We can create a class for each covariance type and write specific methods for each one.

```r
# Example: Diagonal Covariance (S3)

    # Constructor
    makeDiagCov <- function(sd) {
      structure(list(sd = sd), class = c("DiagCov", "StructuredCov"))
    }

    # Method to construct relative covariance factor
    makeLambda <- function(x, ...) UseMethod("makeLambda")

    makeLambda.DiagCov <- function(x, ...) {
      diag(x$sd)
    }

    # Method to extract parameters
    makeTheta <- function(x, ...) UseMethod("makeTheta")

    makeTheta.DiagCov <- function(x, ...) {
      x$sd
    }
```

Here, the object created by `makeDiagCov` has two classes: `"DiagCov"` and `"StructuredCov"`. This allows fallback dispatch to `makeLambda.StructuredCov()` if a specific method is not found.


```r
# Usage Pattern
cov1 <- makeDiagCov(c(1, 2, 3))
makeLambda(cov1)
makeTheta(cov1)
```

This allows you to call `makeLambda()` on any structure, and it will automatically use the right version for that type.


## S4 Class System

S4 lets us define classes more formally, with typed fields and support for dispatching on multiple arguments.

```r
### Example: Diagonal Covariance (S4)
    # Define a virtual superclass
    setClass("StructuredCovariance", contains = "VIRTUAL")

    # Define a concrete subclass
    setClass(
      "DiagonalCovariance",
      slots = list(sd = "numeric"),
      contains = "StructuredCovariance",
      validity = function(object) {
        if (any(object@sd < 0)) "Standard deviations must be non-negative" else TRUE
      }
    )

    # Constructor
    DiagonalCovariance <- function(sd) {
      new("DiagonalCovariance", sd = sd)
    }

    # Define generic and method
    setGeneric("makeLambda", function(x, ...) standardGeneric("makeLambda"))

    setMethod("makeLambda", "DiagonalCovariance", function(x, ...) {
      diag(x@sd)
    })
```

```r
# Usage Pattern
cov1 <- DiagonalCovariance(c(1, 2, 3))
makeLambda(cov1)
```

S4 ensures that all `DiagonalCovariance` objects have a properly structured `sd` slot, and it can validate input at construction time.


## Comparison of S3 and S4 
The following is me thinking through the pros and cons of S3 versus S4 from a design and development perspective. That said, this is my first time doing serious development with either system, so my view is probably naive when it comes to the actual experience of working with them and the downstream effects that might come up in practice when developing in `lme4`. The point here is to make a case for S4 as the preferred approach, while setting that up as a contrast to explore why S3 might actually be more practical.

The downside is that S3 does not check types or structure, and it only dispatches on the first argument. There’s no guarantee that an object passed to a method has the components the method needs. If something is missing or mislabeled, the function might break in unexpected ways, or worse, return the wrong result without warning. Single-argument dispatch can also become a problem when a method needs to respond to how multiple objects interact, which happens in lme4’s modular design. In those cases, S3 forces workarounds like checking types manually inside each function or repeating the same logic across different methods. This clutters the code and makes the programmer responsible for catching errors. Since lme4 is built from separate modules that pass data through different stages of model setup and optimization, this kind of fragility can spread. If we have to validate everything manually, it becomes harder to understand the system as a whole, and it’s possible that a small change in one part could silently break another. Without formal class definitions, shared methods across covariance structures have to be managed by hand, which makes it easier to introduce mistakes or break things without noticing. As more structures are added, it may become harder to keep things consistent with just S3. If we want others to add their own structures, the lack of formal rules could make mistakes more likely and coordination more difficult.

S4 gives us a way to define classes formally. We can specify what slots an object should have, what types they must be, and how they relate to other classes. This adds some overhead, but it helps catch problems early. If someone passes the wrong kind of object or leaves out a required field, we’ll get an error right away instead of chasing a bug later. S4 also lets us dispatch on multiple arguments, which can be useful when method behavior depends on how different parts of the model interact. Since `lme4` already uses S4 in places, especially through the Matrix package, this approach would fit well into the existing system. It could also help us design an interface for structured covariance, especially if we want to support more structures over time. We can define a virtual superclass for all covariance types and write shared methods that apply to anything that inherits from it. This gives us a way to extend the system while keeping the core logic consistent. The main drawback is that S4 takes more time to write and is harder for new contributors to learn. 


## File and Function Organization

```
R/
  diag_cov.R        # Diagonal covariance structure (S3 or S4)
  ar1_cov.R         # AR(1) covariance structure
  structured_cov.R  # Shared methods and virtual classes
  generics.R        # Definitions for makeLambda, makeTheta, updateTheta, etc.
```

Each file should include:

- Constructor
- Class definition (S3 or S4)
- Methods implementing `makeLambda`, `makeTheta`, etc.
- Optional helper functions (e.g. matrix builders)

# Bibliography

- Bates, Douglas. “Converting Packages to S4.” R News, vol. 3, no. 1, 2003, pp. 20–24. https://cran.r-project.org/doc/Rnews/Rnews_2003-1.pdf.

- R Core Team. “Methods for S3 Generic Functions.” R Documentation, R Development Version. https://stat.ethz.ch/R-manual/R-devel/library/methods/html/Methods_for_S3.html.

- Stack Overflow. “Class in R: S3 vs. S4.” Stack Overflow, 23 Jun. 2011. https://stackoverflow.com/questions/6450803/class-in-r-s3-vs-s4.

- Wickham, Hadley. Advanced R. 2nd ed., CRC Press, 2019. https://adv-r.hadley.nz.
