
# Work Log and Schedule

## Bonding Period

### May
In May I accomplished:

- **May 15:** 5 hours – research and writing  
- **May 23:** 3 hours – research and writing  
- **May 27:** 2 hours – research and writing  
- **May 30:** 6 hours – research and writing  

### June
- **June 1:**  
    - Created this work log to document weekly goals and daily tasks throughout the GSoC working period. This will track specific research, design, coding, and discussions, and help align with the proposed schedule.


## Work Period

### Week 1: June 3 – June 8 (And most likely Week 2)

 Finalize the S4 (or possibly S3—though I’m assuming S4 for now to be able to write things out concretely) class hierarchy and method interface for structured covariance matrices in lme4, and begin building the formula parsing and model construction infrastructure that links user-specified structured random effects to the internal model-fitting process.

The first objective this week is to build out a basic class structure using S4. This involves setting up a general parent class called `StructuredCov`, and then adding more specific types to cover different kinds of covariance structures. These include a diagonal type (`DiagonalStructuredCov`) with examples like `IndepCov` and `HeteroDiagCov`, a time series type (`TimeSeriesStructuredCov`) for things like `AR1Cov` or `ToeplitzCov`, a block-diagonal type (`BlockStructuredCov`), spatial types like `ExpCov` and `RatQuadCov` under `SpatialStructuredCov`, and a category for correlation-based structures (`CorrelationStructuredCov`) that would include things like `CSCov` and `UnstructuredCorrCov`. For each of these, each class should store its parameters (in a slot called `theta`), list the parameter names (in something like `paramNames`), and include any other data needed to build the structure—for example, a distance function for spatial structures or information about groupings for block-diagonal ones. I’ll also make sure there are basic checks to validate the input, and that each type can respond to core methods like `makeLambda()`, `makeCovMatrix()`, and `validObject()`. My goal here is to make sure each structure has the right shape and tools to plug into the modeling pipeline.

The second task is to implement a user-facing formula interface that allows users to specify structured random effects using functions like `ar1()`, `cs()`, and `diag()` within a formula. 

Next, the formula parsing infrastructure needs to be able to detect and extract these structured terms. A new function, `findbars_x()`, will be written as a wrapper around `reformulas::findbars()`. This function should preserve the annotations on calls like `ar1(Time | Subject)` and carry along any attached metadata (such as `covtype` and `args`). I may add a tag like `"StructuredTerm"` to help identify these terms later, though I’m still working out the right way to think about how tagging fits into the bigger picture. The function will also include basic checks to make sure each structured term is syntactically valid, uses an acceptable grouping factor, and includes any necessary arguments. In parallel, I will use `splitForm()` to separate the formula into fixed effects, standard random effects, and the structured terms, which will make downstream processing more manageable.

After that, I’ll begin modifying the part of the code that builds the model components from the formula. This involves extending or wrapping the function `mkReTrms()`. In this modified version, I’ll loop over the structured terms, extract their metadata (like the `covtype` and any user-supplied arguments), and then use that information to create the right kind of structured covariance object. At this stage, the idea is to use a function that acts like a constructor: it takes the structure name and arguments, and uses them to create an object of the appropriate S4 class using `new()`. These structured objects will be stored in a new field, probably `covstruc`, in the same object that normally holds all the random effect structures. The grouping design matrices will still be built here, but the actual structured covariance won’t be applied until later, when I call something like `makeLambda()` inside the deviance function.

Finally, I plan(well see...) to wrap up the week with basic prototyping and testing. I’ll generate simple simulated datasets to check that structured terms in the formula are being recognized correctly, that their metadata is preserved and passed through, and that the right kind of structure object is created. I’ll also check that `makeLambda()` and `makeCovMatrix()` are being called where expected and that the model optimization works for at least one structure type, such as a diagonal or AR(1) structure. 

### Week 2: June 9 – June 15 
This week, my main goal is to start building the S4 class structure for covariance matrices, focusing on separating how they're defined from how their parameters are handled for optimization, and then begin connecting this to how formulas are processed. So far this week (Monday and Tuesday), I've successfully refined the S4 class hierarchy. This includes designing parallel virtual classes for parameterization types like LogChol and LogScaleBoundedCor, which will allow concrete classes for specific structures (such as Diagonal, Unstructured, Compound Symmetry, and AR(1)) to inherit from both their mathematical type and their parameter handling method.  While I don't plan to implement multiple different parameterizations for each structure right now (except possibly for unstructured covariance), this architectural choice  will easily allow for future development of alternative parameterizations without rebuilding the core hierarchy.

I've also established all the main S4 functions every covariance object needs to implement, including methods for getting/setting parameters, computing the covariance matrix, finding its inverse, and generating start values. Looking ahead, our next steps for the remainder of the week will involve starting to write the actual S4 class definitions for Diagonal, Unstructured, Compound Symmetry, and AR(1), focusing on setting up their internal storage and the functions that convert parameters for optimization. In parallel, I'll continue planning for the user-facing formula input using functions like ar1(), cs(), and diag(), and how to extract these special terms, including reviewing examples to map formula tags to the correct S4 class names.

### Week 3: June 16 - 22
Complete implementation and testing for `ar()` and `cs()` concrete covariance classes. Then, begin developing the core logic for formula parsing to return the correct S4 covariance object with its associated methods based on formula input.

**Note:** I will not be working on Tuesday, June 17th and Wednesday, June 18th, as I’ll be preparing for a job interview on Wednesday and have tutoring commitments in the evenings. I plan to make up this time over the weekend and early next week, as this is my final week of tutoring.

### Daily Log
| Date    | Hours | Task Summary                                                                                                                                                                   |
| ------  | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| June 1  | 2 hrs | Created and structured the work log; reviewed `glmmTMB` formula syntax documentation                                                                                           |
| June 2  | 4 hrs | Zoom meeting with Mikael on S3 vs S4 design. Drafted architectural notes for formula parsing and user interface                                                                |
| June 3  | 8 hrs | Implemented class structure for `SpatialStructuredCov` as part of S3 vs S4 design document. Began defining remaining classes, including diagonal and block-diagonal structures |
| June 4  | 2 hrs | Read S3 vs S4 chapter in *Advanced R* to deepen understanding of class system differences                                                                                      |
| June 6  | 6 hrs | Researched covariance structure design patterns and rewrote class hierarchy document in response to mentor feedback                                                            |
| June 8  | 2 hrs | Reviewed flexLambda branch and to identify key methods for the S4 covariance class abstract interface.                                                                         |
| June 9  | 2 hrs | Began refining the S4 class hierarchy for covariance structures, focusing on the abstract interface and initial parameterization virtual classes. |
| June 10 | 3 hrs | Continued refining the S4 class hierarchy, setting up concrete classes for AR, Unstructured, Diagonal, and CS, incorporating parallel virtual classes for parameterization. |
| June 11 | 6 hrs | Committed "feat: Implemented foundational S4 covariance structure classes and DiagonalCov." Committed core S4 covariance structure definitions and DiagonalCov implementation to covariance_structures.R (feat: Implemented foundational S4 covariance structure classes and DiagonalCov.). Started designing and prototyping tests for test-covariance_structures.R, planning to push with initial tests today. |
| June 12 | 6 hrs | Finalized and pushed the complete testthat unit test suite for the DiagonalCov class, including edge cases. Changed project strategy to build all classes first and confirmed this with Ben. Analyzed Bens technical feedback on reformulas integration which will be tomorrows task. | 
| June 13 | 8 hrs | Committed: "feat: Implemented Unstructured covariance class definition and methods." Also created a separate script containing helper functions for the `get_parameters` and `set_parameters` methods. Spent a significant portion of the day deepening my understanding of the numerical linear algebra involved in the internals of `lme4`, which meant studying the `Matrix` package in greater detail and researching optimal strategies for matrix computations. Tested the methods informally and plan to use a rough script to begin drafting a `testthat` test suite tomorrow. 1 hour weekly meeting with mentors.  |
| June 14 | 4 hrs | Committed: "refactor(cov): Simplify S4 class structure based on mentor feedback." Refactored the covariance S4 classes to align with best practices and feedback from the weekly meeting. Removed redundant slots to establish `@parameters` as the only source. All computation methods now calculate matrices and values on the fly. Replaced the custom `validate_parameters()` generic with standard `setValidity()`. Implemented multiple dispatch on `set_parameters()` for `UnstructuredCov` using the signature argument for `setMethod`. |
| June 16 | 6 hrs | Debugged UnstructuredCov and DiagonalCov test failures, identified and resolved `diag()` behavior bug for `d=1` log-determinant calculation; refactored helper functions, methods, and test scripts to ensure all `UnstructuredCov` and `DiagonalCov` tests pass. |
| June 19 | 3 hrs | Commited: "wip: feat(cov): Implement `CompoundSymmetryCov` class and methods." and "feat(cov): Add helper function for compound symmetry params". Began working on methods for CS and added a helper function for accessing the sd and correlation parameter from an instantiation. |
| June 20 | 6 hrs | Commited: "feat(cov): Implement `CompoundSymmetryCov` class and methods" and "test: Add tests for `CompoundSymmetryCov`". Finished writing the methods for the `CS` class. Wrote accompanying test suite for `CS`. Debugged test suite and now all tests pass. |


