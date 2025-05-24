## resources/references

* JSS paper
* flexLambda branch
* ex_flexLambda.html presentation from BIRS meeting
* https://github.com/lme4/lme4/blob/flexLambda/R/reGenerators.R
* https://bbolker.github.io/mixedmodels-misc/notes/varmats.html
* https://github.com/glmmTMB/glmmTMB/blob/604f820c67ed2894625b9c6d1b31cf1ff87aea0a/glmmTMB/R/VarCorr.R

## questions

* backward compatibility is critical
* some covariance structures (e.g diagonal, Toeplitz) can easily be implemented by just changing `Lambda`, `Lind`, `theta`. Others (such as AR1) need a function hook that's applied to `theta` before substituting into `Lambda` ... 
* appropriate level of object-orientation: S3, S4???
* leveraging reformulas for formula processing
    * change in formula formatting (more like glmmTMB)
    * can glmmTMB/lme4 share machinery in reformulas?
    * the machinery for constructing Lambdat/Lind, which is inside `mkReTrms`, now lives inside the `reformulas` package. This might need to be refactored. `reformulas` is used by both `lme4` and `glmmTMB`, as well as a handful of downstream packages. I doubt any of the downstream packages (or any un-packaged code in the wild) other than `lme4` actually uses the Lambdat/Lind components (they probably all use only the Z construction machinery). We need to proceed carefully here, but refactoring will almost certainly improve things. (We may also want to move more of the machinery tagging different covariance types from glmmTMB into reformulas ...)
* can we re-use/abstract/refactor any of the (relatively new) machinery from `glmmTMB` for formatting VarCorr objects?
* do any other R packages have this type of functionality? What machinery do they use?
* what's the best way of organizing het/hom pairs? (What's there in `reGenerators.R` is ugly - way too much repeated code)
