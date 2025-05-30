##
library(reformulas)
## from misc/ex_flexLambda.Rmd in flexLambda brancha
## form <- y ~ time + ar1d(~(time|id))
form <- y ~ time + ar1(time|id) + rr(time|id, d = 2) + (1|id)

fr <- expand.grid(time = 1:10, id = factor(1:10))
findbars_x(form)
tt <- mkReTrms(findbars(form), fr)

splitForm(form, specials = c(names(glmmTMB:::.valid_covstruct), "s"))
