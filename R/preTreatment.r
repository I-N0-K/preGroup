# utils::globalVariables("%dopar%")

#' PreTreatment Function for Building Ensemble Models
#'
#' This function facilitates the creation of ensemble models by grouping data based on 
#' treatment indicators and other factors. It allows for various configurations and customizations
#' to tailor the model fitting process to specific needs, including handling of different 
#' statistical family distributions.
#'
#' @param formula an object of class 'formula' describing the model to be fitted.
#' @param data a data frame containing the variables specified in the formula.
#' @param treatment_indicator a character string specifying the column used for treatment division.
#' @param alpha.mvs a dummay vector of length two, indicating the alpha values used
#' within  and between group penalization. A value of 1 yields lasso, a value of 0 ridge, a value
#' between 0 and 1 elastic net. The default uses ridge for feature selection and lasso
#' for group selection. 
#' @param family a description of the error distribution to be used in the model. It can be one 
#'        of 'gaussian', 'binomial', etc.
#' @param ad.alpha optional alpha parameter for adjustment.
#' @param ad.penalty the penalty type to be used in adjustment, default is 'lambda.min'.
#' @param use.grad logical, indicating whether to use gradients.
#' @param weights an optional vector of 'weights' to be used in the fitting process.
#' @param type the type of model components to include: 'rules', 'linear', or 'both'.
#' @param sampfrac the fraction of the data to be used in each bootstrap sample.
#' @param maxdepth maximum depth of the trees in the model.
#' @param learnrate learning rate for boosting.
#' @param mtry number of variables randomly sampled as candidates at each split.
#' @param ntrees number of trees to grow.
#' @param confirmatory optional vector of confirmatory conditions.
#' @param singleconditions logical, to include only single conditions in rules.
#' @param winsfrac the fraction of extreme values to winsorize.
#' @param normalize logical, whether to normalize predictor variables.
#' @param standardize logical, whether to standardize predictor variables.
#' @param ordinal logical, to treat ordered factors as numeric.
#' @param nfolds number of folds for cross-validation.
#' @param tree.control list of control parameters for tree construction.
#' @param tree.unbiased logical, indicating whether to use unbiased trees (conditional inference trees; ctrees) or biased trees (recursive partitionaing; rpart).
#' @param removecomplements logical, whether to remove complement rules.
#' @param removeduplicates logical, whether to remove duplicate rules.
#' @param verbose logical, indicating whether to print detailed output during fitting.
#' @param par.init logical, indicating whether to parallelize the initial fit.
#' @param par.final logical, indicating whether to parallelize the final model fitting.
#' @param sparse logical, to use sparse model matrices.
#' @param ... additional arguments affecting the model fitting.
#'
#' @return A list containing the fitted model object, the original call, and a classification 
#' of the model rules into types such as 'linear', 'prognostic', and 'prescriptive'.
#' @importFrom pre pre
#' @importFrom mvs MVS
# ' @method preTreatment
#' @examples
#' set.seed(123)  # For reproducibility
#' # Number of rows
#' n <- 200
#' 
#' # Generate 5 binary variables
#' binary_vars <- as.data.frame(replicate(5, sample(0:1, n, replace = TRUE)))
#' names(binary_vars) <- paste0("X", 1:5)
#' 
#' # Generate 10 categorical variables with 3 levels each
#' categorical_vars <- as.data.frame(replicate(10, sample(letters[1:3], n, replace = TRUE)))
#' names(categorical_vars) <- paste0("X", 6:15)
#' 
#' # Generate 10 continuous variables
#' continuous_vars <- as.data.frame(replicate(10, rnorm(n)))
#' names(continuous_vars) <- paste0("X", 16:25)

#' # Combine all predictor variables
#' predictors <- cbind(binary_vars, categorical_vars, continuous_vars)

# # Create a non-linear relationship for the response variable y
# # For simplicity, we'll use a combination of sine, cosine, and polynomial terms
#' y <- with(predictors, {
#'   2 * sin(X16) + 3 * cos(X17) + X1 * X18^2 - X2 * X19 + X3 * X20^3 +
#'     as.numeric(X1 == 1 & X9 == "b") + as.numeric(X1 == 1 & X16 < 0) +
#'     as.numeric(X6 == "a") * X21 + as.numeric(X7 == "b") * X22 +
#'     as.numeric(X8 == "c") * X23 +
#'     25 * as.numeric(X8 == "c") + 2 * as.numeric(X9 == "b") + X16 + X18 + X19 + rnorm(n)  # Adding some noise
#' })

#' #Combine predictors and response into a single data frame
#' simdata1 <- data.frame(y = y, predictors)
#' result <- preTreatment(y ~ ., treatment_indicator = "X1", data = simdata1, alpha.mvs = c(0,1))
#' result
#' 
#' @export

preTreatment <- function(formula, data, treatment_indicator, alpha.mvs = c(0,1), family = "gaussian", ad.alpha = NA, 
                ad.penalty = "lambda.min",
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500,
                confirmatory = NULL, singleconditions = FALSE,
                winsfrac = .025, normalize = TRUE, standardize = FALSE,
                ordinal = TRUE, nfolds = 10L, tree.control, tree.unbiased = TRUE, 
                removecomplements = TRUE, removeduplicates = TRUE, 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, 
                sparse = FALSE, ...) {
    
    if(is.null(treatment_indicator)) {
        stop("The treatment_indicator cannot be null to fit the preTreatment model.")
    }

    if(!is.factor(data[ ,treatment_indicator])) {
        message("The treatment_indicator is not a factor. It will be converted to a factor.")
        data[ ,treatment_indicator] <- as.factor(data[ ,treatment_indicator])
    }

    # Convert chacracter variables to factor to fix issues like:
    # Error in inum.default(x, nmax = nmax, ...) : 
    # cannot handle objects of class ‘character’
    data <- data.frame(lapply(data, function(x) if (is.character(x)) as.factor(x) else x))

    if(length(levels(data[ ,treatment_indicator])) != 2) {
        stop("The treatment_indicator must have only 2 levels.")
    }

    # Extract variable names from the formula
    formula_vars <- all.vars(formula)

    if(family %in% c("binomial", "Binomial")) {

      data[ , formula_vars[1]] <- as.factor(data[ , formula_vars[1]])

      if(length(levels(data[ , formula_vars[1]])) != 2) {
        stop("The response variable must have only 2 levels for binomial family.")
      }

      message("The family is binomial but the response is not binary factor. It will be converted to a binary factor.")
    }

    # Check if treatment indicator is in the formula
    if (!(treatment_indicator %in% formula_vars)) {
        message("Treatment indicator is not included in the formula. It will be automatically added to the formula to fit the preTreatment model.")
        # Modify the formula to include the treatment indicator
        formula <- as.formula(paste(deparse(formula), "+", treatment_indicator))
    }


    # fit the pre funtion, fit.final = FALSE
    pre_fit <- pre(formula = formula, data = data, fit.final = FALSE, family = family, ad.alpha = ad.alpha, 
               ad.penalty = ad.penalty, use.grad = use.grad, weights = weights, type = type, sampfrac = sampfrac, 
               maxdepth = maxdepth, learnrate = learnrate, mtry = mtry, ntrees = ntrees,
               confirmatory = confirmatory, singleconditions = singleconditions,
               winsfrac = winsfrac, normalize = normalize, standardize = standardize,
               ordinal = ordinal, nfolds = nfolds, tree.control = tree.control, tree.unbiased = tree.unbiased, 
               removecomplements = removecomplements, removeduplicates = removeduplicates, 
               verbose = verbose, par.init = par.init, par.final = par.final, 
               sparse = sparse,  ...)
    
    # MVS only support numeric group ids
    group_id <- make_group_id(pre_fit, treatment_indicator)

    # MVS does not run when the base level learners have less than 2 columns
    if(min(table(group_id)) < 2) {

            # stop("The mininum number of linear terms / prognostic rules / prescriptive terms is less than 2. The preTreatment model cannot be fitted. The pre model will be fitted.")
            stop("The mininum number of linear terms / prognostic rules / prescriptive terms is less than 2. The preTreatment model cannot be fitted.")
            # return(pre(formula = formula, data = data, family = family, ad.alpha = ad.alpha, 
            #         ad.penalty = ad.penalty, use.grad = use.grad, weights = weights, type = type, sampfrac = sampfrac, 
            #         maxdepth = maxdepth, learnrate = learnrate, mtry = mtry, ntrees = ntrees,
            #         confirmatory = confirmatory, singleconditions = singleconditions,
            #         winsfrac = winsfrac, normalize = normalize, standardize = standardize,
            #         ordinal = ordinal, nfolds = nfolds, tree.control = tree.control, tree.unbiased = tree.unbiased, 
            #         removecomplements = removecomplements, removeduplicates = removeduplicates, 
            #         verbose = verbose, par.init = par.init, par.final = par.final, 
            #         sparse = sparse,  ...))

    } else {
        y_var <- if(family == "binomial") {
          as.factor(pre_fit$data[, 1])
      } else {
          pre_fit$data[, 1]
      }

      mvs_fit <- MVS(x = as.matrix(pre_fit$modmat), y = y_var, alpha = alpha.mvs, nnc = c(0,1),
                view = as.matrix(group_id), family = family, cvloss = "mse")
      
      # now convert the group id to meaningful rule types
      rule_type <- c("Intercept", group_id)

      rule_type[rule_type == 1] <- "linear"
      rule_type[rule_type == 2] <- "prognostic"
      rule_type[rule_type == 3] <- "prescriptive"
      
      result <- list(mvs_fit = mvs_fit, pre_fit = pre_fit, call = match.call(), rule_type = rule_type, treatment_indicator = treatment_indicator)

      class(result) <- "preTreatment"
      return(result)
    }
}

# rules is a vector comes from pre_model$rules$description
split_treatment_indicator_rule <- function(rules, treatment_indicator) {

  if(is.null(treatment_indicator)) {
    stop("The treatment_indicator must be specified.")
  }
    # if(treatment_indicator == "linear") {
    #     stop("The label 'linear' is reserved for linear terms. Please change the treatment_indicator to another name.")
    # }
  treatment_indicator <- as.character(treatment_indicator)

  # Use word boundaries in the regex pattern to ensure the treatment indicator is matched as a whole word
  pattern <- paste0("\\b", treatment_indicator, "\\b")
  # Create a boolean vector indicating whether each rule contains the treatment indicator
  has_treatment_indicator <- grepl(pattern = pattern, x = rules)

  rule_index <- rep(NA, length(rules))
#   prognostic
  rule_index[!has_treatment_indicator] <- 2	
#   prescriptive
  rule_index[has_treatment_indicator] <- 3

  return(rule_index)
}

make_group_id <- function(pre_model, treatment_indicator) {

    if(!inherits(pre_model, "pre")) {
        stop("The model must belong to class 'pre'. ")
    }

    column_names <- colnames(pre_model$modmat)
    rule_presence <- sapply(column_names, function(x) grepl("rule", x))
    # Finally, sum the TRUE values to get the count of column names containing "rule"
    linear_count <- sum(!rule_presence)

    group_ids <- c(rep(1, linear_count), split_treatment_indicator_rule(pre_model$rules$description, treatment_indicator))

}

# Helper function to extract alphas from alpha.mvs
extract_numeric_vector <- function(str) {
  # Remove 'c(' and ')'
  str <- gsub("^c\\(|\\)$", "", str)
  # Convert to numeric
  numeric_vector <- as.numeric(str[2:3])
  return(numeric_vector)
}

#' Summarize a preTreatment Model
#'
#' @description
#' Provides a summary of the preTreatment model, including the tree generating algorithm, variable selection algorithm, and the number of nonzero coefficients for linear, prognostic, and prescriptive predictors.
#'
#' @param object An object of class 'preTreatment'.
#' @return A summary list of the preTreatment model. Can extract the the tree generating algorithm, variable selection algorithm. and number of nonzero coefficients
#' @export
summary.preTreatment <- function(object) {

    preTreatment_model <- object

    if (!inherits(preTreatment_model, "preTreatment")) {
        stop("The model must belong to class 'preTreatment'.")
    }

    # Extract the call and formals
    model_call <- as.list(match.call(preTreatment, preTreatment_model$call))
    model_formals <- formals(preTreatment)

    # Extract parameters from the call or use default values from formals
    tree_unbiased <- if (!is.null(model_call$tree.unbiased)) {
        eval(model_call$tree.unbiased)
    } else {
        model_formals$tree.unbiased
    }

    use_grad <- if (!is.null(model_call$use.grad)) {
        eval(model_call$use.grad)
    } else {
        model_formals$use.grad
    }

    alpha_mvs <- if (!is.null(model_call$alpha.mvs)) {
        eval(model_call$alpha.mvs)
    } else {
        model_formals$alpha.mvs
    }

    alpha_mvs <- extract_numeric_vector(as.character(alpha_mvs))
    # Extract coefficients data frame
    coefficients_df <- coef(preTreatment_model)

    # Calculate the number of nonzero coefficients for each type
    if ("linear" %in% coefficients_df$rule_type) {
        linear_nonzero <- sum(coefficients_df$coefficient != 0 & coefficients_df$rule_type == "linear")
    } else {
        linear_nonzero <- 0
    }

    if ("prognostic" %in% coefficients_df$rule_type) {
        prognostic_nonzero <- sum(coefficients_df$coefficient != 0 & coefficients_df$rule_type == "prognostic")
    } else {
        prognostic_nonzero <- 0
    }

    if ("prescriptive" %in% coefficients_df$rule_type) {
        prescriptive_nonzero <- sum(coefficients_df$coefficient != 0 & coefficients_df$rule_type == "prescriptive")
    } else {
        prescriptive_nonzero <- 0
    }

    # Determine the tree generating algorithm and type
    tree_algorithm <- if (!is.null(use_grad) && !use_grad) {
        "GLM Trees"
    } else if (!is.null(tree_unbiased) && tree_unbiased) {
        "Conditional Inference Trees"
    } else {
        "Recursive Partitioning Trees"
    }

        # Check the first element of alpha_mvs
    level1_algorithm <- if (alpha_mvs[1] == 0) {
        "Ridge for individual feature selection."
    } else if (alpha_mvs[1] == 1) {
        "Lasso for individual feature selection."
    } else {
        sprintf("Elastic net penalty for individual feature selection. Alpha: %f.", alpha_mvs[1])
    }

    # Check the second element of alpha_mvs
    level2_algorithm <- if (alpha_mvs[2] == 1) {
        "Lasso for grouped features selection."
    } else if (alpha_mvs[2] == 0) {
        "Ridge for grouped features selection."
    } else {
        sprintf("Elastic net penalty for grouped features selection. Alpha: %f.", alpha_mvs[2])
    }

    # Print the summary
    cat("Summary of preTreatment Model:\n")
    cat("1. Tree Generating Algorithm:", tree_algorithm, "\n")
    cat("2. Variable Selection Algorithm:", "\n")
    cat("   - Level 1:", level1_algorithm, "\n")
    cat("   - Level 2:", level2_algorithm, "\n")
    cat("3. Number of Nonzero Coefficients:\n")
    cat("   - Linear Predictors:", linear_nonzero, "\n")
    cat("   - Prognostic Predictors:", prognostic_nonzero, "\n")
    cat("   - Prescriptive Predictors:", prescriptive_nonzero, "\n")

        # Create a summary list
    summary_list <- list(
        tree_generating_algorithm = tree_algorithm,
        level_1_algorithm = level1_algorithm,
        level_2_algorithm = level2_algorithm,
        n_linear_nonzero = linear_nonzero,
        n_prognostic_nonzero = prognostic_nonzero,
        n_prescriptive_nonzero = prescriptive_nonzero
    )

    # Return the summary list
    return(summary_list)
}

# print.preTreatment <- function(preTreatment_model) {

# }

#' Extract Grouped Coefficients from preTreatment Model
#'
#' This function calculates the aggregated coefficients for a `preTreatment` model,
#' combining intercepts and coefficients from different hierarchical levels.
#'
#' @param object A model object, which must be of class `preTreatment`.
#' @param ... Additional arguments passed to or from other methods.
#' @return A data frame of combined coefficients and descriptions, sorted with respect to origin.
#' @details The function first verifies the class of the model, then extracts and processes
#' coefficients from multiple levels of the model to compute the grand intercept and other
#' subgroup coefficients. It finally merges these results with the rule descriptions from
#' the `pre_fit` component of the model, handles missing descriptions, and sorts the resulting data frame.
#' @method coef preTreatment
#' @importFrom mvs coef.MVS
#' @export
#' @examples
#' set.seed(123)  # For reproducibility
#' # Number of rows
#' n <- 200
#' 
#' # Generate 5 binary variables
#' binary_vars <- as.data.frame(replicate(5, sample(0:1, n, replace = TRUE)))
#' names(binary_vars) <- paste0("X", 1:5)
#' 
#' # Generate 10 categorical variables with 3 levels each
#' categorical_vars <- as.data.frame(replicate(10, sample(letters[1:3], n, replace = TRUE)))
#' names(categorical_vars) <- paste0("X", 6:15)
#' 
#' # Generate 10 continuous variables
#' continuous_vars <- as.data.frame(replicate(10, rnorm(n)))
#' names(continuous_vars) <- paste0("X", 16:25)

#' # Combine all predictor variables
#' predictors <- cbind(binary_vars, categorical_vars, continuous_vars)

# # Create a non-linear relationship for the response variable y
# # For simplicity, we'll use a combination of sine, cosine, and polynomial terms
#' y <- with(predictors, {
#'   2 * sin(X16) + 3 * cos(X17) + X1 * X18^2 - X2 * X19 + X3 * X20^3 +
#'     as.numeric(X1 == 1 & X9 == "b") + as.numeric(X1 == 1 & X16 < 0) +
#'     as.numeric(X6 == "a") * X21 + as.numeric(X7 == "b") * X22 +
#'     as.numeric(X8 == "c") * X23 +
#'     25 * as.numeric(X8 == "c") + 2 * as.numeric(X9 == "b") + X16 + X18 + X19 + rnorm(n)  # Adding some noise
#' })

#' #Combine predictors and response into a single data frame
#' simdata1 <- data.frame(y = y, predictors)
#' 
#' result <- preTreatment(y ~ ., treatment_indicator = "X1", data = simdata1, alpha.mvs = c(0,1))
#' 
#' coef(result)
#' 


coef.preTreatment <- function(object, ...) {

    # Ensure the object is of the correct class
    if (!inherits(object, "preTreatment")) {
        stop("The model must belong to class 'preTreatment'.")
    }

    preTreatment_model <- object

    mvs_coef_list <- coef(preTreatment_model$mvs_fit)

    level_2_intercept <- mvs_coef_list$"Level 2"[[1]][1]

    level_2_coef <- mvs_coef_list$"Level 2"[[1]][-1]

    grand_intercept <- level_2_intercept

    coefs <- c()

    for(i in seq_along(level_2_coef)) {
        # print("level 2 coef")
        # print(level_2_coef[i])

        # print("level 1 coef")
        # print(mvs_coef_list$"Level 1"[[i]])

        # print("level 2 coef * level 1 coef")
        # print(level_2_coef[i] * mvs_coef_list$"Level 1"[[i]])

        subgroup_intercept <- as.matrix(level_2_coef[i] * mvs_coef_list$"Level 1"[[i]])[1]
        grand_intercept  <- grand_intercept + subgroup_intercept

        x_names <- dimnames(mvs_coef_list$"Level 1"[[i]])[[1]][-1]
        x_coefs <- as.matrix(level_2_coef[i] * mvs_coef_list$"Level 1"[[i]])[-1]
        
        # print(length(x_names))
        # print(length(x_coefs))
        subgroup_coef <- cbind(rule = x_names, coefficient = x_coefs)

        coefs <- as.data.frame(rbind(coefs, subgroup_coef))

    }
    grand_intercept <- cbind(rule = "(Intercept)", coefficient = as.numeric(grand_intercept))

    coefs <- rbind(grand_intercept, coefs)

    coefs <- cbind(coefs, rule_type = preTreatment_model$rule_type)

    # Merge the coefs
    coefs <- merge(coefs, preTreatment_model$pre_fit$rules, by = "rule", all.x = TRUE)

    # Add a helper column to flag rows that originally have NA descriptions
    coefs$na_originally <- is.na(coefs$description)

    # Fill missing descriptions with the rule name for intercept and linear terms
    coefs$description[is.na(coefs$description)] <- coefs$rule[is.na(coefs$description)]

    # Sort the dataframe to have NA originally at the top
    coefs <- coefs[order(coefs$na_originally, decreasing = TRUE), ]

    # Remove the helper column
    coefs$na_originally <- NULL

    coefs$coefficient <- as.numeric(coefs$coefficient)

    return(coefs)
}

#' Predicts responses based on a `preTreatment` model and new data.
#'
#' @param object A `preTreatment` model object from which predictions will be made.
#' @param newdata New data on which predictions are to be made.
#' @param type The type of prediction to be made: 'response' or 'HTE'. Type 'response' returns the predicted 
#' response values, while 'HTE' returns the heterogeneous treatment effects (for causal modelling).
#' For gaussian outcome, the predicted result is on the original scale of the response variable. 
#' For binomial outcome, the predicted result is the p(Y = 1|Xs) for type "response" and the 
#' odds ratio of the treatment condition against the control condition for type "HTE":
#' p1 = p(Y = 1|Treatment, Xs)
#' p2 = p(Y = 1|Control, Xs)
#' Odds ratio = (p1 / (1-p1)) / (p2 / (1 - p2)).
#' @param ... Additional arguments passed to or from other methods.
#' @return A vector of predicted responses.
#' @details The function checks if the model object is a `preTreatment` and then uses the
#' `get_new_X` function to transform new data into the model matrix format. Predictions
#' are then made using the multivariate structure of the `preTreatment` model.
#' @method predict preTreatment
#' @importFrom mvs predict.MVS
#' @export
#' @examples
#' # Assuming `model` is a preTreatment model and `new_data` is available
#' \dontrun{
#' mod <- preTreatment(Species ~ ., data = iris, treatment_indicator = "Petal.Width")
#' new_data <- iris[1:5, ]
#' predictions <- predict(mod, new_data)
#' }

predict.preTreatment <- function(object, newdata, type = "response", ...) {
    # Validate model class
    if (!inherits(object, "preTreatment")) {
        stop("The model must belong to class 'preTreatment'.")
    }

    # Define and check permissible types
    valid_types <- c("response", "HTE", "hte", "Hte")
    if (!type %in% valid_types) {
        stop("type must be one of 'response', 'HTE'.")
    }

    preTreatment_model <- object

    new_x <- get_new_X(preTreatment_model, newdata)

    # Process prediction based on type
    if (type  %in% c("response", "Response", "RESPONSE")) {
        return(predict(preTreatment_model$mvs_fit, new_x, predtype = "response", cvlambda = "lambda.min"))

    } else if (type %in% c("HTE", "hte", "Hte")) {
        # Avoid repeating predictions by computing them once
        return(calculate_hte(preTreatment_model, newdata))
    }
}

# Helper function for HTE prediction
calculate_hte <- function(preTreatment_model, newdata) {
    raw_treatment_indicator <- preTreatment_model$treatment_indicator
    # print(raw_treatment_indicator)
    # print(preTreatment_model$pre_fit$data[, raw_treatment_indicator])
    levels <- tryCatch(sort(unique(as.numeric(levels(as.factor(preTreatment_model$pre_fit$data[, raw_treatment_indicator])))), decreasing = TRUE),
                       error = function(e) {levels(as.factor(preTreatment_model$pre_fit$data[, raw_treatment_indicator]))})

    # Check if enough levels exist
    if (length(levels) != 2) {
        stop("Only 2 different treatment levels should be used to compute HTE. The current model has ", length(levels), " levels.")
    }

    message("HTE will be calculated for treatment ", levels[1], " against treatment ", levels[2], ".\n")

    newdata_g1 <- newdata
    newdata_g1[, raw_treatment_indicator] <- factor(as.numeric(levels[1]), levels = levels)
    # print(newdata_g1[, raw_treatment_indicator])
    # print(colnames(newdata_g1))

    newdata_g2 <- newdata
    newdata_g2[, raw_treatment_indicator] <- factor(as.numeric(levels[2]), levels = levels)
    # print(newdata_g2[, raw_treatment_indicator])
    # print(headings(newdata_g2))
    # print("starting HTE calculation")

    # Calculate HTE
    if (preTreatment_model$call$family %in% c("binomial", "Binomial")) {

      # Obtain probabilities for each group in binomial family
      p1 <- predict(preTreatment_model, newdata = newdata_g1, type = "response")
      p2 <- predict(preTreatment_model, newdata = newdata_g2, type = "response")

      # Use odds ratio: treatment vs control as HTE for binomial family
      hte <- (p1/(1 - p1)) / (p2/(1 - p2))

    } else {
      hte <- predict(preTreatment_model, newdata = newdata_g1, type = "response") - 
             predict(preTreatment_model, newdata = newdata_g2, type = "response")
    }

    return(hte)
}

calculate_hte_pre <- function(pre_model, newdata, treatment_indicator) {

    raw_treatment_indicator <- treatment_indicator
    newdata[ ,raw_treatment_indicator] <- as.factor(newdata[ ,raw_treatment_indicator])
    levels <- sort(unique(as.numeric(levels(newdata[ ,raw_treatment_indicator]))), decreasing = TRUE)

    # Check if enough levels exist
    if (length(levels) < 2) {
        stop("Not enough treatment levels to compute HTE.")
    }

    message("HTE will be calculated for treatment ", levels[1], " against treatment ", levels[2], ".\n")

    # Set up data for two groups
    newdata_g1 <- newdata
    newdata_g1$raw_treatment_indicator <-  factor(as.numeric(levels[1]), levels = levels)
    # print(newdata_g1$raw_treatment_indicator)
    newdata_g2 <- newdata
    newdata_g2$raw_treatment_indicator <-  factor(as.numeric(levels[2]), levels = levels)
    # print(newdata_g2$raw_treatment_indicator)
    # within(newdata, {
    #     raw_treatment_indicator <- factor(levels[1], levels = levels)
    # })
    # newdata_g2 <- within(newdata, {
    #     raw_treatment_indicator <- factor(levels[2], levels = levels)
    # })
    
    # Calculate HTE
    hte <- predict(pre_model, newdata = newdata_g1) - 
           predict(pre_model, newdata = newdata_g2)
    return(hte)
}

#' Generate Model Matrix from New Data for preTreatment Model
#'
#' This function processes new data to fit the structure of a `preTreatment` model,
#' facilitating the prediction process.
#'
#' @param preTreatment_model A `preTreatment` model object, specifically its `pre_fit` component.
#' @param new_data New data to be transformed into model matrix format.
#' @return A list containing the transformed new data as a model matrix.
#' @details The function retrieves settings and parameters from the `pre_fit` component of
#' the `preTreatment` model and applies these to the new data to generate a suitable model matrix.
#' This involves applying transformation rules and handling different types of predictors as
#' specified in the model.
#' @export
# ' @method get_new_X preTreatment
#' @examples
#' # Assuming `model` is a preTreatment model and `new_data` is ready
#' \dontrun{
#' mod_matrix <- get_new_X(model, new_data)
#' }

get_new_X <- function(preTreatment_model, new_data) {

    pre_model <- preTreatment_model$pre_fit

    treatment_indicator <- preTreatment_model$treatment_indicator

    if(!is.factor(new_data[ ,treatment_indicator])) {
        message("The treatment_indicator is not a factor. It will be converted to a factor.")
        new_data[ ,treatment_indicator] <- as.factor(new_data[ ,treatment_indicator])
    }

    # get_modmat returns a list instead of the matrix
    new_modlist <- get_modmat(
      wins_points = pre_model$wins_points, 
      x_scales = pre_model$x_scales, 
      formula = pre_model$formula, 
      data = new_data, 
      rules = if (pre_model$type == "linear" || is.null(pre_model$rules)) {NULL} else {
        structure(pre_model$rules$description, names = pre_model$rules$rule)}, 
      type = pre_model$type, 
      winsfrac = if (is.null(pre_model$winsfrac)) {
        formals(pre)$winsfrac
      } else {
        pre_model$winsfrac
      },
      x_names = pre_model$x_names, 
      normalize = pre_model$normalize,
      y_names = NULL,
      confirmatory = tryCatch(
          eval(pre_model$call$confirmatory),
          error = function(e) NULL
      )
    )
    
    return(new_modlist$x)
}


#' @importFrom utils head
#' @importFrom stats sd quantile model.frame quantile terms model.matrix predict coef
#' @importFrom Formula Formula
#' @importFrom survival is.Surv
#' @importFrom MatrixModels model.Matrix
#' @importFrom pre pre
#' @importFrom graphics barplot par text
get_modmat <- function(
  # Pass these if you already have an object
  wins_points = NULL, x_scales = NULL,
  # These should be passed in all calls
  formula, data, rules, type, x_names, winsfrac, normalize, 
  # Response variable names is optional:
  y_names = NULL, sparse = FALSE, rulevars = NULL,
  confirmatory = NULL) {
  
  ## TODO: is next line necessary? Only for converting ordered categorical vars to linear terms?
  ## Is only used twice, in get_rulemat function below.
  ## Perhaps should be conditional on ordinal argument of pre()?
  data_org <- data # needed to evaluate rules later
  
  # convert ordered categorical predictor variables to linear terms:
  data[ , sapply(data, is.ordered)] <- 
    as.numeric(data[ , sapply(data, is.ordered)])
  
  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    conf_vars <- confirmatory[confirmatory %in% x_names]
    conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_rules) > 0L) {
      ## Evaluate rules (add to model matrix later):
      eval_rules <- function(data, rule_descr) {
        1L * with(data, eval(parse(text = rule_descr)))
      }
      conf_rules <- mapply(FUN = eval_rules, rule_descr = conf_rules, 
                           MoreArgs = list(data = data))
    }
  }
  
  #####
  # Perform winsorizing and normalizing
  ## TODO: Allow for not supplying all variables, but only variables with non-zero importances:
  if (type != "rules" && any(sapply(data[ , x_names], is.numeric))) {
    #####
    # if type is not rules, linear terms should be prepared:
    
    # Winsorize numeric variables (section 5 of F&P(2008)):
    if (winsfrac > 0) {
      
      if (miss_wins_points <- is.null(wins_points)) {
        wins_points <- data.frame(varname = x_names, value = NA, lb = NA, ub = NA, 
                                  stringsAsFactors = FALSE)
      }
      
      j <- 0L
      tol <- sqrt(.Machine$double.eps)
      for (i in x_names) {
        j <- j + 1L
        if (is.numeric(data[[i]])) {
          if (miss_wins_points) {
            lim <- quantile(data[ , i], probs = c(winsfrac, 1 - winsfrac))
            wins_points$value[j] <- paste(lim[1L], "<=", i, "<=", lim[2L])
            lb <- lim[1L]
            ub <- lim[2L]
            if (ub - lb < tol) {
              ## If lower and upper bound are equal, do not winsorize and issue warning:
              warning("Variable ", x_names[j], " will be winsorized employing winsfrac = 0, to prevent reducing the variance of its linear term to 0.", immediate. = TRUE)
              wins_points$lb[j] <- min(data[ , i])
              wins_points$ub[j] <- max(data[ , i])
            } else {
              wins_points$lb[j] <- lb
              wins_points$ub[j] <- ub
              data[ , i][data[ , i] < lb] <- lb
              data[ , i][data[ , i] > ub] <- ub
            }
          } else {
            data[ , i][data[ , i] < wins_points$lb[j]] <- wins_points$lb[j]
            data[ , i][data[ , i] > wins_points$ub[j]] <- wins_points$ub[j]
          }
        }
      }
    }
    
    # normalize numeric variables:
    if (normalize) { 
      # Normalize linear terms (section 5 of F&P08), if there are any:
      needs_scaling <- x_names[sapply(data[ , x_names, drop = FALSE], is.numeric)]
      if (length(needs_scaling) > 0) {
        if (is.null(x_scales)) {
          x_scales <- apply(
            data[ , needs_scaling, drop = FALSE], 2L, sd, na.rm = TRUE) / 0.4
        }
        ## check if variables have zero variance (if so, do not scale):
        tol <- sqrt(.Machine$double.eps)
        almost_zero_var_inds <- which(x_scales < tol)
        if (length(almost_zero_var_inds) > 0) {
          # print warning and set all those x_scales to 1
          warning("A total of ", length(almost_zero_var_inds), " variable(s) (e.g.,", paste0(head(needs_scaling[almost_zero_var_inds]), collapse = ", "), ") have sd < ", tol, " and will not be normalized.")  
          # omit from needs_scaling:
          x_scales[almost_zero_var_inds] <- 1
        }
        data[ , needs_scaling] <- scale(
          data[ , needs_scaling, drop = FALSE], center = FALSE, scale = x_scales)
      }
    }
  }

  ## Combine rules and variables:
  if (length(rules) == 0L) rules <- NULL
  if (type == "linear" || is.null(rules)) {
    if (sparse) {
      mf <- model.frame(Formula(formula), data)
      x <- model.Matrix(terms(mf), mf, sparse = TRUE)
      rm(mf)
      
    } else {
      x <- model.matrix(Formula(formula), data = data)
    }
  } else if (type %in% c("both", "rules") && !is.null(rules)) {
    if (sparse) {
      x <- if (is.null(rulevars)) { 
        .get_rules_mat_sparse(data_org, rules) 
        } else {
          rulevars
        }
      if (type == "both") {
        mf <- model.frame(Formula(formula), data)
        x <- cbind(model.Matrix(terms(mf), mf, sparse = TRUE), x)
        rm(mf)
      }
    } else { 
      x <- if (is.null(rulevars)) {
        .get_rules_mat_dense(data_org, rules)
        } else { 
          rulevars
        }
      if (type == "both") {
        x <- cbind(model.matrix(Formula(formula), data = data), x)
      }
    }
    
  } else { 
    stop("not implemented with type ", sQuote(type), " and is.null(rules) is ", 
         sQuote(is.null(rules)))
  }
  
  #####
  # Remove intercept
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  ## Add confirmatory terms:
  if (!is.null(confirmatory)) {
    ## Get confirmatory rules and variables:
    #conf_vars <- confirmatory[confirmatory %in% x_names]
    #conf_rules <- confirmatory[!(confirmatory %in% conf_vars)]
    if (length(conf_vars) > 0L) {
      ## If type == "rules", add variables to model matrix:
      if (type == "rules") {
        x <- cbind(x, as.matrix(data[conf_vars], rownames.force = FALSE))
        ## TODO: Then check if it should be winsorized or normalized, or issue warning if not
      }
    }
    if (is.matrix(conf_rules)) {
      #if (!sparse) {
        x <- cbind(x, conf_rules)
        ## I believe both sparse and non-sparse matrices can use cbind:
        ## Matrix::cBind returns:
        ## .Defunct(msg = "'cBind' is defunct; 'base::cbind' handles S4 objects since R 3.2.0")
      #} else {
        ## TODO: implement something for sparse matrices.
      #}
    }
  }
  
  if (is.null(y_names)) {
    y <- NULL
  } else {
    y <- data[ , y_names]
    if (is.Surv(y) || length(y_names) > 1L) {
      y <- as.matrix(y)
    }
  }
  
  if (!exists("wins_points", inherits = FALSE)) { wins_points <- NULL }
  if (!exists("x_scales", inherits = FALSE)) { x_scales <- NULL }
  
  list(x = x, y = y, x_scales = x_scales, wins_points = wins_points)
}

.get_rules_mat_dense <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr, data)
  colnames(x) <- names(rules)
  x
}

.get_rules_mat_sparse <- function(data, rules){
  if(length(rules) == 0)
    return(NULL)
  
  # See https://stackoverflow.com/a/8844057/5861244. 
  #
  # if all rules where binary then we could use the `lsparseMatrix-classes`. 
  # However, this will not work when we call `glmnet` as it requires a double
  # matrix`
  expr <- paste0("cbind_sparse_vec(", paste0(
    'as(as.numeric(', rules, '), "sparseVector")', collapse = ", "), ")")
  x <- eval(parse(text = expr), data)
  colnames(x) <- names(rules)
  x
}

##' @export
importance <- function(x, ...)  UseMethod("importance")

#' Calculate Variable Importances for a preTreatment Model
#'
#' Computes variable importances based on a preTreatment model, offering options for both global and local importance, 
#' standardization of importances, and various additional parameters.
#'
#' @param preTreatment_model An object of class 'preTreatment'.
#' @param standardize Logical; if TRUE, importance scores are standardized.
#' @param global Logical; if TRUE, global importances are calculated, otherwise local importances are provided.
#' @param penalty.par.val A character string specifying the penalty parameter value to be used in calculation.
#' @param gamma The gamma parameter for regularization paths in models; if not NULL, it modifies the penalties.
#' @param quantprobs A numeric vector specifying the quantiles to be used in local importance calculation.
#' @param round Numeric; the number of decimal places to round the importances to (if NA, no rounding is performed).
#' @param plot Logical; if TRUE, a plot of the importances is generated.
#' @param ylab The y-axis label in the plot.
#' @param main The main title for the plot.
#' @param abbreviate An integer specifying how much to abbreviate variable names in the plot.
#' @param diag.xlab Logical; if TRUE, labels are diagonal on the x-axis of the plot.
#' @param diag.xlab.hor Horizontal adjustment for diagonal labels.
#' @param diag.xlab.vert Vertical adjustment for diagonal labels.
#' @param cex.axis The character expansion size for axis names.
#' @param legend A character string or logical; position of the legend in the plot or FALSE to exclude it.
#' @param ... Additional arguments passed to underlying functions.
#' @return A list containing the variable importances as data frames.
#' @method importance preTreatment
#' @export
#' @aliases importance
#' @examples
#' result <- importance(preTreatment_model)

importance.preTreatment <- function(preTreatment_model, standardize = FALSE, global = TRUE,
                           penalty.par.val = "lambda.1se", gamma = NULL,
                           quantprobs = c(.75, 1),
                           round = NA, plot = TRUE, ylab = "Importance",
                           main = "Variable importances", abbreviate = 10L,
                           diag.xlab = TRUE, diag.xlab.hor = 0, diag.xlab.vert = 2,
                           cex.axis = 1, legend = "topright", ...)
{

  if (!inherits(preTreatment_model, what = "preTreatment")) {
    stop("Specified object is not of class 'preTreatment'.")
  }
  
  x <- preTreatment_model$pre_fit

  if (!global) {
    if (x$family %in% c("mgaussian", "multinomial")) {
      warning("Local importances cannot be calculated for multivariate and multinomial outcomes. Global importances will be returned.")
      global <- TRUE 
    }
  }
  
  if (standardize && x$family %in% c("multinomial", "binomial", "cox")) {
    warning("Standardized importances cannot be calculated for binary, multinomial or survival responses. Unstandardized importances will be returned.")
    standardize <- FALSE
  }
  
  ## Step 1: Calculate the importances of the base learners:
  
  ## get coefficients:
  if (is.null(gamma)) {
    coefs <- coef(preTreatment_model, penalty.par.val = penalty.par.val)
  } else {
    coefs <- coef(preTreatment_model, penalty.par.val = penalty.par.val, gamma = gamma)    
  }
  if (x$family %in% c("mgaussian", "multinomial")) {
    # Returns the y variable names
    coef_inds <- names(coefs)[!names(coefs) %in% c("rule", "description")]
  }

  ## continue only when there are nonzero terms besides intercept:
  if ((x$family %in% c("gaussian", "binomial", "poisson") && 
      sum(coefs$coefficient != 0) > 1L ) || 
      (x$family %in% c("mgaussian", "multinomial") && 
       sum(rowSums(coefs[,coef_inds]) != 0) > 1L) ||
      (x$family == "cox" && sum(coefs$coefficient != 0) > 0)) { 
    ## give factors a description:
    if (any(is.na(coefs$description))) {
      coefs$description[is.na(coefs$description)] <-
        paste0(as.character(coefs$rule)[is.na(coefs$description)], " ")
    }
    coefs <- coefs[order(coefs$rule), ]
    
    ## Get SDs for every baselearner:
    if (global) {
      ## Get SDs (x$x_scales should be used to get correct SDs for linear terms)
      if (x$family == "cox") {
        sds <- apply(x$modmat, 2, sd, na.rm = TRUE)  
      } else {
        sds <- c(0, apply(x$modmat, 2, sd, na.rm = TRUE))          
      }
      if (standardize) {
        if (x$family == "mgaussian") {
          sd_y <- sapply(x$data[ , x$y_names], sd)
        } else if (x$family %in% c("gaussian", "poisson")) { 
          sd_y <- sd(as.numeric(x$data[ , x$y_names]))
        }
      }
    } else {
      preds <- predict.pre(x, newdata = x$data, type = "response",
                           penalty.par.val = penalty.par.val, ...)
      local_modmat <- x$modmat[preds >= quantile(preds, probs = quantprobs[1]) &
                                      preds <= quantile(preds, probs = quantprobs[2]),]
      if (nrow(local_modmat) < 2) {stop("Selected subregion contains less than 2 observations, importances cannot be calculated")}
      ## x$x_scales should be used to get correct SDs for linear terms:
      if (x$family == "cox") {
        ## cox prop hazard model has no intercept, so should be omitted
        sds <- apply(local_modmat, 2, sd, na.rm = TRUE)
      } else {
        sds <- c(0, apply(local_modmat, 2, sd, na.rm = TRUE))
      }

      if (standardize) {
        sd_y <- sd(x$data[preds >= quantile(preds, probs = quantprobs[1]) & 
                                 preds <= quantile(preds, probs = quantprobs[2]),
                               x$y_names])
      }
    }
    
    ## Check if there are any " ` " marks in sd names, if so remove:
    if (any(grepl("`", names(sds), fixed = TRUE))) {
      names(sds) <- gsub("`", "", names(sds), fixed = TRUE)
    }
    
    if(x$normalize) {
      sds[names(x$x_scales)] <- sds[names(x$x_scales)] * x$x_scales
    }
    
    if (x$family != "cox") {
      names(sds)[1] <- "(Intercept)"
    }
    
    sds <- sds[order(names(sds))]
  
    if (any(names(sds) != coefs$rule)) {
      warning("There seems to be a problem with the ordering or size of the coefficient and sd vectors. Importances cannot be calculated.")
    }
    
    ## baselearner importance is given by abs(coef*SD) (F&P section 6):
    if (x$family %in% c("multinomial", "mgaussian")) {
      baseimps <- data.frame(coefs, sd = sds)
      baseimps[,gsub("coefficient", "importance", coef_inds)] <- abs(sapply(baseimps[,coef_inds], function(x) x*sds))
    } else {
      baseimps <- data.frame(coefs, sd = sds, imp = abs(coefs$coefficient)*sds)
    }
    
    if (standardize) {
      if (x$family == "mgaussian") {
        for (i in gsub("coefficient", "importance", coef_inds)) {
          baseimps[,i] <- baseimps[,i] / sd_y[gsub("importance.", "", i)]
        }
      } else if (x$family %in% c("gaussian", "poisson")) {
        baseimps$imp <- baseimps$imp / sd_y
      }
    }

    ## Remove nonzero terms:
    if (x$family %in% c("mgaussian", "multinomial")) {
      baseimps <- baseimps[rowSums(baseimps[,coef_inds]) != 0, ]   
    } else {
      baseimps <- baseimps[baseimps$coefficient != 0,]
    }
    
    ## Omit intercept:
    baseimps <- baseimps[baseimps$description != "1",]
    baseimps <- baseimps[baseimps$description != "(Intercept)",]
    
    ## Calculate the number of conditions in each rule:
    baseimps$nterms <- NA
    for(i in 1:nrow(baseimps)) {
      ## If there is " & " in description, there are at least 2 conditions/variables 
      ## in the base learner:
      if (grepl(" & ", baseimps$description[i])) {
        baseimps$nterms[i] <- length(gregexpr("&", baseimps$description)[[i]]) + 1L
      } else {
        baseimps$nterms[i] <- 1L # if not, the number of terms equals 1
      }
    }
    
    ## if no winsorizing is performed, descriptions look different then with winsorizing
    ##    so add a temporary space AFTER description
    if (!is.null(x$call$winsfrac)) {
      if (x$call$winsfrac == 0) {
        linear_term_ids <- which(baseimps$rule == baseimps$description)
        for (i in linear_term_ids) {
          baseimps$description[i] <- paste0(baseimps$description[i], " ")
        }
      }
    }
    
    
    ## Step 2: Calculate variable importances:
    if (x$family %in% c("mgaussian", "multinomial")) {
      varimps <- data.frame(varname = x$x_names, stringsAsFactors = FALSE)
      varimps[,gsub("coefficient", "importance", coef_inds)] <- 0
    } else {
      varimps <- data.frame(varname = x$x_names, imp = 0,
                            stringsAsFactors = FALSE)
    }
    
    for(i in 1:nrow(varimps)) {
      ## Get imps from rules and linear functions
      for(j in 1:nrow(baseimps)) {
        
        ## if the variable name appears in the description (of rule or linear term):
        ##   (Note: EXACT matches are needed, so 1) there should be a space before 
        ##     and after the variable name in the rule and thus 2) there should be 
        ##     a space added before the description of the rule)
        if (grepl(paste0(" ", varimps$varname[i], " "), paste0(" ", baseimps$description[j]))) {
          ## Count the number of times it appears in the rule
          n_occ <- length(gregexpr(paste0(" ", varimps$varname[i], " "),
                                   paste0(" ", baseimps$description[j]), fixed = TRUE)[[1]])
          ## Add to the importance of the variable
          if (x$family %in% c("mgaussian", "multinomial")) {
            varimps[i, gsub("coefficient", "importance", coef_inds)] <- 
              varimps[i, gsub("coefficient", "importance", coef_inds)] + 
            (n_occ * baseimps[j, gsub("coefficient", "importance", coef_inds)] / baseimps$nterms[j])
          } else {
            varimps$imp[i] <- varimps$imp[i] + (n_occ * baseimps$imp[j] / baseimps$nterms[j])
          }
        }
      }
  
      ## Get imps for factors
      if (is.factor(x$data[ , varimps$varname[i]])) { # check if variable is a factor and add importance
          # !is.ordered(x$data[ , varimps$varname[i]])) {
        ## Sum baseimps$imp for which baseimps$rule has varimps$varname[i] as _part_ of its name
        if (x$family %in% c("mgaussian", "multinomial")) {
          varimps[i, gsub("coefficient", "importance", coef_inds)] <-
            varimps[i, gsub("coefficient", "importance", coef_inds)] +
            colSums(baseimps[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE), 
                             gsub("coefficient", "importance", coef_inds)])
        } else { # otherwise, simply add importance(s)
          varimps$imp[i] <- varimps$imp[i] + 
            sum(baseimps$imp[grepl(varimps$varname[i], baseimps$rule, fixed = TRUE)])
        }
      } 
    }
    
    
    ## Step 3: Return (and plot) importances:
    
    if (x$family %in% c("mgaussian", "multinomial")) {
      varimps <- varimps[rowSums(varimps[ , gsub("coefficient", "importance", coef_inds)]) != 0, ]
      ord <- order(rowSums(varimps[ , gsub("coefficient", "importance", coef_inds)]), 
                   decreasing = TRUE, method = "radix")
      varimps <- varimps[ord, ]
    } else {
      baseimps <- baseimps[order(baseimps$imp, decreasing = TRUE, method = "radix"), ]
      varimps <- varimps[order(varimps$imp, decreasing = TRUE, method = "radix"), ]
      # varimps <- varimps[varimps$imp != 0, ]
    }
    
    if (!is.na(round)) {
      baseimps[,sapply(baseimps, is.numeric)] <- round(baseimps[,sapply(baseimps, is.numeric)], digits = round)
      varimps[,sapply(varimps, is.numeric)] <- round(varimps[,sapply(varimps, is.numeric)], digits = round)
    }
    
    if (x$family %in% c("mgaussian","multinomial")) {
      keep <- c("rule", "description", gsub("coefficient", "importance", coef_inds), "rule_type", 
                coef_inds, "sd")
    } else {
      keep <- c("rule", "description", "imp", "coefficient", "sd", "rule_type")
    }

    baseimps <- data.frame(baseimps[, keep], stringsAsFactors = FALSE)

    row.names(baseimps) <- row.names(varimps) <- NULL
    ## Remove added space AFTER description if winsorizing was performed
    if (!is.null(x$call$winsfrac)) {
      if (x$call$winsfrac == 0L) {
        for (i in linear_term_ids) {
          baseimps$description[i] <- substring(baseimps$description[i], first = 2L)
        }
      }
    }
    
    imp <- list(varimps = varimps, baseimps = baseimps)

    # Upgrade the imp object to specify linear, prognostic, and prescriptive importances
    # Only support univariate outcome for now.
    imp <- specify_imps(imp)

    # Create plotting matrix, variable needs to be in columns
    imp_plot <- t(as.matrix(imp$varimps[, c("imp_linear", "imp_prognostic", "imp_prescriptive")]))
    rownames(imp_plot) <- gsub("^imp_", "", rownames(imp_plot))
    colnames(imp_plot) <- imp$varimps$varname

    # Sort importances by total importance
    imp_plot <- imp_plot[, order(colSums(imp_plot), decreasing = TRUE)]

    # Save this matrix for further plotting on the user side
    imp$imp_plot <- imp_plot

    if (plot & nrow(varimps) > 0) {
      if (is.character(legend)) {
        args.legend <- list(x = legend)
        legend.text <- TRUE
      } else {
        legend.text <- NULL
        args.legend <- NULL
      }

      if (x$family %in% c("mgaussian", "multinomial")) {
        # plot_varimps <- t(varimps[ , gsub("coefficient", "importance" , coef_inds)])
        plot_varimps_list <- create_importance_matrix_multivariate(imp)

        plot_varimps <- NA
        
        for(i in seq_along(plot_varimps_list)) {

          plot_varimps  <- plot_varimps_list[[i]]

          if (diag.xlab) {
            xlab.pos <- barplot(plot_varimps, beside = FALSE, ylab = ylab, 
                                names.arg = rep("", times = ncol(plot_varimps)), 
                                main = paste(main, names(plot_varimps_list)[i]), cex.axis = cex.axis, 
                                legend.text = legend.text, 
                                args.legend = args.legend, ...)

            # xlab.pos <- xlab.pos[nrow(xlab.pos),]
            ## add specified number of trailing spaces to variable names:
            plotnames <- varimps$varname
            if (is.numeric(abbreviate) && abbreviate > 0) {
              plotnames <- abbreviate(plotnames, minlength = abbreviate)
            }
            if (diag.xlab.vert > 0) {
              for (i in 1:diag.xlab.vert) {
                plotnames <- paste0(plotnames, " ")
              }
            }
            text(xlab.pos + diag.xlab.hor, par("usr")[3], srt = 45, adj = 1, xpd = TRUE, 
                labels = plotnames, cex = cex.axis)
          } else {
            barplot(plot_varimps, beside = FALSE, main = paste(main, names(plot_varimps_list)[i]), ylab = ylab, 
                    legend.text = legend.text, args.legend = args.legend,
                    cex.axis = cex.axis, ...)

          }
        }
      } else {
        if (diag.xlab) {

          xlab.pos <- barplot(imp_plot, xlab = "", ylab = ylab, 
                              main = main, cex.axis = cex.axis, legend.text = legend.text, 
                                args.legend = args.legend, ...)
          ## add specified number of trailing spaces to variable names:
          plotnames <- varimps$varname
          if (is.numeric(abbreviate) && abbreviate > 0) {
            plotnames <- abbreviate(plotnames, minlength = abbreviate)
          }
          if (diag.xlab.vert > 0) {
            for (i in 1:diag.xlab.vert) {
              plotnames <- paste0(plotnames, " ")
            }
          }
          # text(xlab.pos + diag.xlab.hor, par("usr")[3], srt = 45, adj = 1, xpd = TRUE, 
          #      labels = plotnames, cex = cex.axis)
        } else {
          plotnames <- varimps$varname
          if (is.numeric(abbreviate) && abbreviate > 0) {
            plotnames <- abbreviate(plotnames, minlength = abbreviate)
          }
          barplot(height = imp_plot, names.arg = plotnames, ylab = ylab,
                main = main, cex.axis = cex.axis, legend.text = legend.text, 
                                args.legend = args.legend, ...)
        }
      }
    }

    return(invisible(imp))
   
  } else {
    warning("No non-zero terms in the ensemble. All importances are zero.")
    return(invisible(NULL))
  }
}

#' Specify linear, prognostic, and prescriptive importances
#'
#' Generates a matrix of importances for each variable based on the importances calculated in `importance.preTreatment` 
#' for models with a single response variable.
#'
#' @param imp_in A list containing the importance data generated by `importance.preTreatment`.
#' @return A list: Element varimps is a matrix with importance type (linear, prognostic, prescriptive) as columns and variable names as rows. Element baseimps is a data frame containing rule importances and variable importances in each rule.
#' @details This function only supports univariate outcome. Multivariate version needs to be further implemented.
#' 
specify_imps <- function(imp_in) {
  
  baseimps_in <- imp_in$baseimps
  # Remove rows where rule_type is "Intercept"
  baseimps_in <- baseimps_in[baseimps_in$rule_type != "(Intercept)", ]

  varimps_in <- imp_in$varimps
  
  # # Function to extract variables and their frequencies from a description
  # extract_variables_with_frequency <- function(description) {
  #   # Split the description by '&' to handle subrules
  #   subrules <- strsplit(description, "&")[[1]]
    
  #   # Initialize an empty vector to store variable names
  #   variables <- c()
    
  #   # Define a regular expression pattern to match variable names
  #   pattern <- "\\s*([a-zA-Z0-9_\\.]+)\\s*(<|>|<=|>=|%in%)"
    
  #   # Loop through each subrule and extract variable names
  #   for (subrule in subrules) {
  #     matches <- regmatches(subrule, gregexpr(pattern, subrule, perl = TRUE))
  #     if (length(matches[[1]]) > 0) {
  #       # Extract the variable names from the matches
  #       vars <- gsub(pattern, "\\1", matches[[1]])
  #       variables <- c(variables, vars)
  #     }
  #   }
    
  #   # Count the frequency of each variable
  #   variable_frequency <- table(variables)
    
  #   return(variable_frequency)
  # }
  
    # Function to extract variables and their frequencies from a description
  extract_variables_with_frequency <- function(description, rule_type) {
    if (rule_type %in% c("prognostic", "prescriptive")) {
      # Split the description by '&' to handle subrules
      subrules <- strsplit(description, "&")[[1]]
      
      # Initialize an empty vector to store variable names
      variables <- c()
      
      # Define a regular expression pattern to match variable names
      pattern <- "\\s*([a-zA-Z0-9_\\.]+)\\s*(<|>|<=|>=|%in%)"
      
      # Loop through each subrule and extract variable names
      for (subrule in subrules) {
        matches <- regmatches(subrule, gregexpr(pattern, subrule, perl = TRUE))
        if (length(matches[[1]]) > 0) {
          # Extract the variable names from the matches
          vars <- gsub(pattern, "\\1", matches[[1]])
          variables <- c(variables, vars)
        }
      }
      
      # Count the frequency of each variable
      variable_frequency <- table(variables)
      
      return(variable_frequency)
    } else if (rule_type == "linear") {

      # for linear varaibles, the pattern contains 2 parts: first part the varaible name; length not clear. 
      # Second part the category; length not clear. 
      # to make a safe match, the fucntion should match all possible linear variables 
      # to the description from left to right, exact match. 
      # the variable that has the most length of exact match is the target variable 
      # that the descrtiption matches.

      # Extract the variable name and category from the description for linear variables.
      # Match all possible linear variables to the description from left to right.
      # Identify the variable that has the longest exact match with the description.

      # For linear rules, the description is the variable name
      variable <- description
      
      # Initialize an empty vector to store possible matches
      possible_matches <- c()
      
      # Loop through all variables to find the best match
      for (var in all_variables) {
        if (startsWith(variable, var)) {
          possible_matches <- c(possible_matches, var)
        }
      }
      
     # Find the variable with the longest exact match
      if (length(possible_matches) > 0) {
        exact_matches <- sapply(possible_matches, function(v) {
          match_length <- nchar(v)
          if (substr(variable, 1, match_length) == v) {
            return(match_length)
          } else {
            return(0)
          }
        })
        best_match <- possible_matches[which.max(exact_matches)]
        variable_frequency <- table(best_match)
      } else {
        variable_frequency <- table(variable)
      }
      
      return(variable_frequency)
    }
  }

  all_variables <- varimps_in$varname
  
  # Initialize baseimps_out with baseimps_in and add columns for each variable with zeros
  baseimps_out <- baseimps_in
  for (var in all_variables) {
    baseimps_out[[var]] <- 0
  }
  
  # Populate the columns with the frequency of each variable in the description
  for (i in 1:nrow(baseimps_out)) {
    variable_frequency <- extract_variables_with_frequency(baseimps_out$description[i], baseimps_out$rule_type[i])
    for (var in names(variable_frequency)) {
        baseimps_out[i, var] <- variable_frequency[[var]]
    }
  }
  
  # Calculate the length of the description for each row
  baseimps_out$DescriptionLength <- rowSums(baseimps_out[, all_variables])
  
  # Calculate variable importance for each variable
  for (var in all_variables) {
    imp_col <- paste0(var, "_imp")
    baseimps_out[[imp_col]] <- (baseimps_out$imp / baseimps_out$DescriptionLength) * baseimps_out[[var]]
  }
  
  # Remove the columns corresponding to the variables and DescriptionLength
  baseimps_out[, c(all_variables, "DescriptionLength")] <- NULL
  # Replace NaN values with zeros
  baseimps_out[is.na(baseimps_out)] <- 0
  
  varimps_out <- varimps_in
  
  # Rename the linear imps
  varimps_out$imp <- NULL
  
  # Copy the total imp value
  varimps_out$imp_total <- varimps_in$imp
  
  #
   varimps_out$imp_linear <- sapply(varimps_out$varname, function(v) {
    if (any(baseimps_out$rule_type == "linear")) {
      sum(baseimps_out[baseimps_out$rule_type == "linear", paste0(v, "_imp")], na.rm = TRUE)
    } else { 0 }
  })

  varimps_out$imp_prognostic <- sapply(varimps_out$varname, function(v) {
    if (any(baseimps_out$rule_type == "prognostic")) {
      sum(baseimps_out[baseimps_out$rule_type == "prognostic", paste0(v, "_imp")], na.rm = TRUE)
    } else { 0 }
  })
  
  varimps_out$imp_prescriptive <- sapply(varimps_out$varname, function(v) {
    if (any(baseimps_out$rule_type == "prescriptive")) {
      sum(baseimps_out[baseimps_out$rule_type == "prescriptive", paste0(v, "_imp")], na.rm = TRUE)
    } else { 0 }
  })
  
  varimps_out$imp_total <- varimps_out$imp_linear + varimps_out$imp_prognostic + varimps_out$imp_prescriptive
  
  varimps_out <- varimps_out[, c("varname", "imp_total", "imp_linear", "imp_prognostic", "imp_prescriptive")]
  
  return(list(varimps = varimps_out, baseimps = baseimps_out))
}
