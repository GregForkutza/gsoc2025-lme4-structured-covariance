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

### Week 1: June 3 – June 7 (And most likely Week 2)

 Finalize the S4 (or possibly S3—though I’m assuming S4 for now to be able to write things out concretely) class hierarchy and method interface for structured covariance matrices in lme4, and begin building the formula parsing and model construction infrastructure that links user-specified structured random effects to the internal model-fitting process.

The first objective this week is to build out a basic class structure using S4. This involves setting up a general parent class called `StructuredCov`, and then adding more specific types to cover different kinds of covariance structures. These include a diagonal type (`DiagonalStructuredCov`) with examples like `IndepCov` and `HeteroDiagCov`, a time series type (`TimeSeriesStructuredCov`) for things like `AR1Cov` or `ToeplitzCov`, a block-diagonal type (`BlockStructuredCov`), spatial types like `ExpCov` and `RatQuadCov` under `SpatialStructuredCov`, and a category for correlation-based structures (`CorrelationStructuredCov`) that would include things like `CSCov` and `UnstructuredCorrCov`. For each of these, each class should store its parameters (in a slot called `theta`), list the parameter names (in something like `paramNames`), and include any other data needed to build the structure—for example, a distance function for spatial structures or information about groupings for block-diagonal ones. I’ll also make sure there are basic checks to validate the input, and that each type can respond to core methods like `makeLambda()`, `makeCovMatrix()`, and `validObject()`. My goal here is to make sure each structure has the right shape and tools to plug into the modeling pipeline.


The second task is to implement a user-facing formula interface that allows users to specify structured random effects using functions like `ar1()`, `cs()`, and `diag()` within a formula. 

Next, the formula parsing infrastructure needs to be able to detect and extract these structured terms. A new function, `findbars_x()`, will be written as a wrapper around `reformulas::findbars()`. This function should preserve the annotations on calls like `ar1(Time | Subject)` and carry along any attached metadata (such as `covtype` and `args`). I may add a tag like `"StructuredTerm"` to help identify these terms later, though I’m still working out the right way to think about how tagging fits into the bigger picture. The function will also include basic checks to make sure each structured term is syntactically valid, uses an acceptable grouping factor, and includes any necessary arguments. In parallel, I will use `splitForm()` to separate the formula into fixed effects, standard random effects, and the structured terms, which will make downstream processing more manageable.

After that, I’ll begin modifying the part of the code that builds the model components from the formula. This involves extending or wrapping the function `mkReTrms()`. In this modified version, I’ll loop over the structured terms, extract their metadata (like the `covtype` and any user-supplied arguments), and then use that information to create the right kind of structured covariance object. At this stage, the idea is to use a function that acts like a constructor: it takes the structure name and arguments, and uses them to create an object of the appropriate S4 class using `new()`. These structured objects will be stored in a new field, probably `covstruc`, in the same object that normally holds all the random effect structures. The grouping design matrices will still be built here, but the actual structured covariance won’t be applied until later, when I call something like `makeLambda()` inside the deviance function.

Finally, I plan(well see...) to wrap up the week with basic prototyping and testing. I’ll generate simple simulated datasets to check that structured terms in the formula are being recognized correctly, that their metadata is preserved and passed through, and that the right kind of structure object is created. I’ll also check that `makeLambda()` and `makeCovMatrix()` are being called where expected and that the model optimization works for at least one structure type, such as a diagonal or AR(1) structure. 





### Daily Log

| Date   | Hours | Task Summary                                                                                                                                                                   |
| ------ | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| June 1 | 2 hrs | Created and structured the work log; reviewed `glmmTMB` formula syntax documentation                                                                                           |
| June 2 | 4 hrs | Zoom meeting with Mikael on S3 vs S4 design. Drafted architectural notes for formula parsing and user interface                                                                |
| June 3 | 8 hrs | Implemented class structure for `SpatialStructuredCov` as part of S3 vs S4 design document. Began defining remaining classes, including diagonal and block-diagonal structures |

