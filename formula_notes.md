Some thoughts on what you will need to adjust/modify for the covariance structures you have in mind:

* diagonal: the 'set up Lambdat with integer indices and stick the values of theta in' approach will still work.
* ar1, cs: you'll also have to implement the 'hook' function that transforms theta parameters to elements of the Cholesky factor ...

Here's the stuff that comes out of MkReTrms:

* Zt, Ztlist, flist, cnms, Gp, nl: unchanged  (this is generally true, as
these don't have to do specifically with the covariance structure or how
it is parameterized, but with the mapping from model formula to latent
variables and the mapping from latent variables to linear predictor ...)

* Lambdat: needs to be diagonal (duh)

* Lind: pattern matching Lambdat (if `n` is the number of latent variables
per block and `nb` is the number of blocks, then this will be something
like `rep(seq(n), nb)` [for each random effect term - these are then
concatenated)
   (see mkBlist in reformulas)

theta: needs to change in length (should probably be rep(1, ntheta)
lower: vector of all zeros

The big architectural question is how much we should try to change
things in reformulas::mkReTrms ...   The tricky thing is that
reformulas::mkReTrms is used by a variety of packages, although I
*think* that most of them are only using it to get Z-components. In any
case, we want to maintain back-compatibility with the current version of
lme4, so I think that should work ...

Reverse imports: buildmer, clustTMB, glmmTMB, gmvjoint, lme4, mbest,
nimbleMacros, predictmeans, robustlmm, tramME, ubms, unmarked
Reverse suggests:       performance

## updates

Suggested strategy: copy the machinery needed for constructing Lambdat, Lind out of `reformulas::mkReTrms`, make a separate `mkLambdat` (or whatever) function that works how you need it to work. You can call `mkReTrms` first so that you have cnms, Gp, nl to work with (indexing/dimension info)
