setClass("VirtualCorrelationModel", contains = "VIRTUAL",
         slots = c(parameters = "numeric"))

#' Represents a Compound Symmetry (CS) correlation structure.
setClass("CSCorrelationModel", contains = "VirtualCorrelationModel")

#' Represents a first-order Autoregressive (AR1) correlation structure.
setClass("AR1CorrelationModel", contains = "VirtualCorrelationModel")


# Composer class()
#' A flexible covariance structure composed of a variance and a correlation model.
setClass("FlexibleCovariance",
         slots = c(
             dimension = "integer",
             variance_model = "VirtualVarianceModel",
             correlation_model = "VirtualCorrelationModel"
         )
)

# Methods for Correlation Models

setGeneric("compute_correlation_matrix", function(object, d) standardGeneric("compute_correlation_matrix"))

#' Method for Compound Symmetry Correlation.
setMethod("compute_correlation_matrix", "CSCorrelationModel", function(object, d) {
    rho <- tanh(object@parameters[1])
    R <- Matrix(rho, nrow = d, ncol = d)
    diag(R) <- 1
    return(R)
})

#' Method for AR(1) Correlation.
setMethod("compute_correlation_matrix", "AR1CorrelationModel", function(object, d) {
    rho <- tanh(object@parameters[1])
    time_diffs <- abs(outer(1:d, 1:d, "-"))
    R <- rho^time_diffs
    return(R)
})

# Method for the Composer class()
setGeneric("compute_covariance_matrix", function(object) standardGeneric("compute_covariance_matrix"))

setMethod("compute_covariance_matrix", "FlexibleCovariance", function(object) {
    R <- compute_correlation_matrix(object@correlation_model, d = object@dimension)
    if (is(object@variance_model, "HomogeneousVarianceModel")) {
        sigma_sq <- exp(object@variance_model@parameters[1])
        return(sigma_sq * R)
    } else {
        st_devs <- exp(0.5 * object@variance_model@parameters)
        D <- Diagonal(object@dimension, x = st_devs)
        return(D %*% R %*% D)
    }
})

# Example

homo_cs <- new("FlexibleCovariance",
               dimension = 3L,
               variance_model = new("HomogeneousVarianceModel", parameters = log(4)),
               correlation_model = new("CSCorrelationModel", parameters = atanh(0.5))
)

hetero_cs <- new("FlexibleCovariance",
                 dimension = 3L,
                 variance_model = new("HeterogeneousVarianceModel", parameters = log(c(4, 9, 16))),
                 correlation_model = new("CSCorrelationModel", parameters = atanh(0.5))
)
