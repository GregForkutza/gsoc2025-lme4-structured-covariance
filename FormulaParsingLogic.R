
# Goal: Translate a formula with covariance terms like `y ~ x + ar1(1 | group)` 
# into specific data strucutres: `Lambda`, `theta`, and `Lind`. 

# 1. user call lmer(y ~ x + ar1(1 | g), data = df)
# 2. lmer Entry Point: The modified lmer function immediately uses 
# reformulas to parse the ar1(1|g) term and creates a FlexibleCovariance object 
# with an AR1CorrelationModel.
# 3. Call to mkReTrms: lmer eventually calls mkReTrms, passing it the parsed term
# and the newly created S4 object.
# 4. mkReTrms Dispatch: Inside mkReTrms, an if statement sees that the S4 object
# is not the standard UnstructuredCovariance.
# 5. S4 Method Calls: It calls get_lambda(s4_object) and get_lind(s4_object) to generate 
# the specialized AR1 correlation structure and index vector.
# 6. Assembly: mkReTrms takes these custom components and assembles them into the final ReTrms object. 
# 7. The rest of the lmer function proceeds as normal, passing this ReTrms object to the 
# C++ optimization engine, which now has the correct recipe to fit an AR1 model.




# Step 1: Define the user facing syntax
cs <- function(term) term
ar1 <- function(term) term
us <- function(term) term     
diag <- function(term) term   
sqexp <- function(term) term  
matern <- function(term) term 
# ... and so on 

# Step 2: Implement the Formula Parser 
# Deconstruct the formula and using a lookup table build the 
# correct S4 object. 

# a) Look up table: maps the string tag from the formula to the corresponding S4 class name. 

cov_struct_table <- list(
    # tag    # Correlation S4 Class        # Variance S4 Class (default)
    cs       = list(cor_class = "CSCorrelationModel"),
    ar1      = list(cor_class = "AR1CorrelationModel"),
    sqexp    = list(cor_class = "SquaredExpKernelModel"),
    matern   = list(cor_class = "MaternKernelModel"),
    diag     = list(cor_class = "IdentityCorrelationModel"), 
    us       = list(is_unstructured = TRUE) # Flag for the default case. 
    # Add new structures here, e.g., toeplitz = list(...)
)

# b) Parser Function: this is called by the lmer wrapper . Uses reformulas to find the random terms 
# and the lookup table to build S4 object. 

#' Create Covariance Object from a Parsed Formula Term
#' @param term A single random effect term parsed by `reformulas`.
#' @param data The model data frame.
#' @return An S4 object inheriting from `VirtualCovariance`.
create_covariance_object_from_term <- function(term, data) {
    
    structure_name <- term$type 
    
    spec <- cov_struct_table[[structure_name]]
    if (is.null(spec)) {
        stop("Unknown covariance structure provided: ", structure_name)
    }

    random_effects_matrix <- model.matrix(term$term_form, data)
    dimension <- ncol(random_effects_matrix)

     if (!is.null(spec$is_unstructured) && spec$is_unstructured) {
        return(new("UnstructuredCovariance", dimension = dimension))
    } else {
        # All other cases use the FlexibleCovariance composer
        variance_model <- new("HeterogeneousVarianceModel", dimension = dimension)
        correlation_model <- new(spec$cor_class) # Dynamically create class from string

        return(new("FlexibleCovariance",
            dimension = dimension,
            variance_model = variance_model,
            correlation_model = correlation_model
        ))
    }
}

# Step 3: Add methods to generate lme4 components 

# Generic to get the relative covariance factor (Cholesky of the correlation matrix)
setGeneric("get_lambda", function(object) standardGeneric("get_lambda"))

# Generic to get the index mapping vector from theta to Lambda
setGeneric("get_lind", function(object) standardGeneric("get_lind"))

# We already have a method for getting theta. 


# Step 4: Modify the lmer function entry point.
# Insert logic to parse the formula 

# Inside the existing lmer function...
lmer <- function(formula, data, ...) {    
     
    # 1. Use reformulas to find and parse all random effects terms.
    parsed_formula <- reformulas::parse_formula(formula, ...)
    
    # 2. For each parsed random term, create the corresponding S4 object.
      s4_object_list <- lapply(
        parsed_formula$random,
        create_covariance_object_from_term,
        data = data
    )
    
}

# Step 5: New Logic inside mkReTerms 
# It receives the parsed random term and its corresponding S4 object.
# Make mkReTrms a dispatcher that either runs the original lme4 code
# for unstructured models or delegates creation of Lambda and Lind 
# to the S4 objects. 
# 

# We need to modify `mkReTrms to also accept a s4_object_list, which is generated
# at the beginning of modifed lmer function.
# Conceptual Sketch of logic inside `mkReTerms`. 

for (i in seq_along(random_effects)) {
    
    # random_effects is the relevant component taken from the list that is
    # returned by the reformulas::parse_formula() function.
    term <- random_effects[[i]]
    s4_object <- s4_object_list[[i]]

    # Build the model matrix Zt using `term`. 
        
    if (is(s4_object, "UnstructuredCovariance")) {
        
        # If it's the default case, use the existing lme4 code
        # to generate the standard unstructured Lambda and Lind.
        # call original lme4 logic
        
    } else {
        
        # If it's a FlexibleCovariance object,
        # generate the components by calling the S4 methods.
        lambda <- get_lambda(s4_object)
        lind <- get_lind(s4_object)
        theta <- get_start_values(s4_object)
        
        # now we can use these new components to build this part of the ReTrms list
        # ReTrms is a list that hold all the essential matrices and vectors related to 
        # the random effects part of a mixed model. It contains everything the lme4 
        # optimization engine needs to know about the random effects structure. 
    }
}
