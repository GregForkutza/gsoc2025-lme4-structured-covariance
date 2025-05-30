# Designing Class Structures for Structured Covariance in lme4: S3 vs S4
This document compares the S3 and S4 object systems in R, focusing on how they support structured covariance matrices in `lme4`. The project aims to implement a variety of structures, such as diagonal, compound symmetric, autoregressive, and user-defined forms, that differ in internal logic but require a common interface. Each must support methods for constructing the relative covariance factor (`makeLambda`), accessing and updating parameters (`makeTheta`, `updateTheta`), and optionally building the full covariance matrix (`makeSigma`). The central design question is how to organize these components using S3 or S4.

## S3 Class System
In S3, an object is considered to belong to a class simply because it has a class attribute.There’s no enforcement that the object has any specific fields (slots), nor that those fields are of the right type. S3 lets us tag objects with a class and use naming patterns to run the right method. We can create a class for each covariance type and write specific methods for each one.

S3 method dispatch looks only at the class of the first argument to decide which method to use. If your function depends on the interaction between two or more objects, S3 can’t help you automatically choose the right method.

```r
foo <- function(x, y) UseMethod("foo")

foo.A <- function(x, y) cat("Called foo for class A\n")
foo.B <- function(x, y) cat("Called foo for class B\n")

objA <- structure(list(), class = "A")
objB <- structure(list(), class = "B")

foo(objA, objB)  # Calls foo.A
foo(objB, objA)  # Calls foo.B
```
If the class is misapplied or a field is missing, S3 will not stop you at construction time. Errors only appear during execution.

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

### Type and Structure Checking

S3 doesn’t check object structure. If you assign the right class name, R assumes the object is valid, even if it’s missing components the method expects.

```r
# S3: no structure enforcement
makeS3Thing <- function(x) structure(list(x = x), class = "S3Thing")

describe <- function(obj) UseMethod("describe")
describe.S3Thing <- function(obj) paste("This is an S3Thing with x =", obj$x)

obj1 <- makeS3Thing(42)
describe(obj1)  # Works

badObj <- list(y = 99)
class(badObj) <- "S3Thing"
describe(badObj)  # Error: object 'x' not found
````

In contrast, S4 enforces structure at creation time.

```r
# S4: formal structure checking
setClass("S4Thing", slots = list(x = "numeric"))
S4Thing <- function(x) new("S4Thing", x = x)

setGeneric("describe", function(obj) standardGeneric("describe"))
setMethod("describe", "S4Thing", function(obj) paste("This is an S4Thing with x =", obj@x))

obj2 <- S4Thing(10)
describe(obj2)  # Works

# Fails at creation time
# badObj <- new("S4Thing", y = 99)
```

S4 guarantees that required components are present and correctly typed, so methods can safely access them.

### Method Dispatch

S3 only dispatches on the first argument, which becomes limiting when method behavior depends on how multiple objects interact. This forces manual type checks and code duplication.

```r
# S3: simulate multi-dispatch with manual type check
combine <- function(a, b) UseMethod("combine")
combine.S3Thing <- function(a, b) {
  if (!inherits(b, "OtherThing")) stop("Expected OtherThing")
  paste("Combining x =", a$x, "with z =", b$z)
}
a <- structure(list(x = 5), class = "S3Thing")
b <- structure(list(z = "hello"), class = "OtherThing")
combine(a, b)  # Works
```

S4 handles this naturally with multi-argument dispatch.

```r
# S4: multi-dispatch example
setClass("OtherThing", slots = list(z = "character"))
setGeneric("combine", function(a, b) standardGeneric("combine"))

setMethod("combine", signature("S4Thing", "OtherThing"), function(a, b) {
  paste("Combining x =", a@x, "with z =", b@z)
})

combine(S4Thing(1), new("OtherThing", z = "hello"))  # Works
```
Here’s a clear and natural way to write that section. It builds directly off your comparison section, explains what a virtual superclass is in plain terms, and then proposes a design pattern that uses it effectively in your project.

### Inheritance/interface consistency
To add.

## Using a Virtual Superclass for Structured Covariance

In S4, a virtual superclass is a class that you never create directly. It just defines a common structure for other classes to inherit from. In the example above, `StructuredCovariance` is a virtual class, and `DiagonalCovariance` is a subclass that inherits from it. The idea is to group all structured covariance types under a single parent class—like diagonal, AR(1), compound symmetric, or user-defined structures.

This lets us write generic functions like `makeLambda()`, `makeTheta()`, and `updateTheta()` that apply to any subclass. Each structure can override the methods it needs, but everything shares the same interface. For example, `DiagonalCovariance` has its own version of `makeLambda()`, but it still inherits from `StructuredCovariance`, so it fits into the same system.

Once that structure is in place, downstream code can work with any object that inherits from `StructuredCovariance`, without needing to know the specific type. This makes it easier to add new covariance types later, or let contributors define their own, as long as they follow the expected interface.

`lme4` is built in a modular way. Different parts of the model, like the design matrices, covariance factors, and parameter vectors are created and updated in separate stages. By using a virtual superclass for structured covariance, we can make sure each structure supports the same basic methods. This means the rest of the code can call things like `makeLambda()` without needing to know what type of structure it’s working with. That keeps the pieces separate and makes it easier to add new structures later without changing how the rest of the system works.



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
