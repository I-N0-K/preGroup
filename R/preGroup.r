require(rpart)
require(partykit)
require(glmnet)
require(pre)
require(mvs)
require(Formula)
require(survival)
# utils::globalVariables("%dopar%")

#' Pre-Group Function for Building Ensemble Models
#'
#' This function facilitates the creation of ensemble models by grouping data based on 
#' treatment indicators and other factors. It allows for various configurations and customizations
#' to tailor the model fitting process to specific needs, including handling of different 
#' statistical family distributions.
#'
#' @param formula an object of class 'formula' describing the model to be fitted.
#' @param data a data frame containing the variables specified in the formula.
#' @param treatment_indicator a character string specifying the column used for treatment division.
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
#' @param tree.unbiased logical, indicating whether to use unbiased trees.
#' @param removecomplements logical, whether to remove complement rules.
#' @param removeduplicates logical, whether to remove duplicate rules.
#' @param verbose logical, indicating whether to print detailed output during fitting.
#' @param par.init logical, indicating whether to parallelize the initial fit.
#' @param par.final logical, indicating whether to parallelize the final model fitting.
#' @param sparse logical, to use sparse model matrices.
#' @param ... additional arguments affecting the model fitting.
#'
#' @return A list containing the fitted model object, the original call, and a classification 
#'         of the model rules into types such as 'linear', 'prognostic', and 'prescriptive'.
#'
#' @examples
#' \dontrun{
#'  data(iris)
#'  result <- preGroup(Species ~ ., data = iris, treatment_indicator = "Petal.Width")
#'  print(result)
#' }
#'
#' @export

preGroup <- function(formula, data, treatment_indicator, family = "gaussian", ad.alpha = NA, 
                ad.penalty = "lambda.min",
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500,
                confirmatory = NULL, singleconditions = FALSE,
                winsfrac = .025, normalize = TRUE, standardize = FALSE,
                ordinal = TRUE, nfolds = 10L, tree.control, tree.unbiased = TRUE, 
                removecomplements = TRUE, removeduplicates = TRUE, 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, 
                sparse = FALSE, ...) {
    
    if(!is.factor(data[ ,treatment_indicator])) {
        message("The treatment_indicator is not a factor. It will be converted to a factor.")
        data[ ,treatment_indicator] <- as.factor(data[ ,treatment_indicator])
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

            warning("The mininum number of linear terms / prognostic rules / prescriptive terms is less than 2. The preGroup model cannot be fitted. The pre model will be fitted.")

            return(pre(formula = formula, data = data, family = family, ad.alpha = ad.alpha, 
                    ad.penalty = ad.penalty, use.grad = use.grad, weights = weights, type = type, sampfrac = sampfrac, 
                    maxdepth = maxdepth, learnrate = learnrate, mtry = mtry, ntrees = ntrees,
                    confirmatory = confirmatory, singleconditions = singleconditions,
                    winsfrac = winsfrac, normalize = normalize, standardize = standardize,
                    ordinal = ordinal, nfolds = nfolds, tree.control = tree.control, tree.unbiased = tree.unbiased, 
                    removecomplements = removecomplements, removeduplicates = removeduplicates, 
                    verbose = verbose, par.init = par.init, par.final = par.final, 
                    sparse = sparse,  ...))

    } else {
        y_var <- if(family == "binomial") {
          as.factor(pre_fit$data[, 1])
      } else {
          pre_fit$data[, 1]
      }

      mvs_fit <- MVS(x = as.matrix(pre_fit$modmat), y = y_var, alpha = c(1,0), nnc = c(0,1),
                view = as.matrix(group_id), family = family, cvloss = "mse")
      
      # now convert the group id to meaningful rule types
      rule_type <- c("Intercept", group_id)

      rule_type[rule_type == 1] <- "linear"
      rule_type[rule_type == 2] <- "prognostic"
      rule_type[rule_type == 3] <- "prescriptive"
      
      result <- list(mvs_fit = mvs_fit, pre_fit = pre_fit, call = match.call(), rule_type = rule_type, treatment_indicator = treatment_indicator)

      class(result) <- "preGroup"
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

    if(!class(pre_model) == "pre") {
        stop("The model must belong to class 'pre'. ")
    }

    column_names <- colnames(pre_model$modmat)
    rule_presence <- sapply(column_names, function(x) grepl("rule", x))
    # Finally, sum the TRUE values to get the count of column names containing "rule"
    linear_count <- sum(!rule_presence)

    group_ids <- c(rep(1, linear_count), split_treatment_indicator_rule(pre_model$rules$description, treatment_indicator))

}

#' Extract Grouped Coefficients from preGroup Model
#'
#' This function calculates the aggregated coefficients for a `preGroup` model,
#' combining intercepts and coefficients from different hierarchical levels.
#'
#' @param preGroup_model A model object, which must be of class `preGroup`.
#' @param ... Additional arguments passed to or from other methods.
#' @return A data frame of combined coefficients and descriptions, sorted with respect to origin.
#' @details The function first verifies the class of the model, then extracts and processes
#' coefficients from multiple levels of the model to compute the grand intercept and other
#' subgroup coefficients. It finally merges these results with the rule descriptions from
#' the `pre_fit` component of the model, handles missing descriptions, and sorts the resulting data frame.
#' @export
#' @examples
#' # Assuming `model` is a preGroup model
#' result <- coef.preGroup(model)
#' print(result)

coef.preGroup <- function(preGroup_model, ...) {

    # Ensure the object is of the correct class
    if (!inherits(preGroup_model, "preGroup")) {
        stop("The model must belong to class 'preGroup'.")
    }

    mvs_coef_list <- coef(preGroup_model$mvs_fit)

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

    coefs <- cbind(coefs, rule_type = preGroup_model$rule_type)

    # Merge the coefs
    coefs <- merge(coefs, preGroup_model$pre_fit$rules, by = "rule", all.x = TRUE)

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


#' Predict Responses Using preGroup Model
#'
#' Predicts responses based on a `preGroup` model and new data.
#'
#' @param preGroup_model A `preGroup` model object from which predictions will be made.
#' @param newdata New data on which predictions are to be made.
#' #' @param type The type of prediction to be made: 'response' or 'HTE'. Type 'response' returns the predicted 
#' response values, while 'HTE' returns the heterogeneous treatment effects (for causal modelling).
#' @param ... Additional arguments passed to or from other methods.
#' @return A vector of predicted responses.
#' @details The function checks if the model object is a `preGroup` and then uses the
#' `get_new_X` function to transform new data into the model matrix format. Predictions
#' are then made using the multivariate structure of the `preGroup` model.
#' @export
#' @examples
#' # Assuming `model` is a preGroup model and `new_data` is available
#' predictions <- predict.preGroup(model, new_data)


predict.preGroup <- function(preGroup_model, newdata, type = "response", ...) {
    # Validate model class
    if (!inherits(preGroup_model, "preGroup")) {
        stop("The model must belong to class 'preGroup'.")
    }

    # Define and check permissible types
    valid_types <- c("response", "HTE", "hte", "Hte")
    if (!type %in% valid_types) {
        stop("type must be one of 'response', 'HTE'.")
    }

    # Process prediction based on type
    if (type == "response") {
        new_x <- get_new_X(preGroup_model, newdata)
        return(predict(preGroup_model$mvs_fit, new_x, predtype = "response", cvlambda = "lambda.min"))
    } else {
        # Avoid repeating predictions by computing them once
        return(calculate_hte(preGroup_model, newdata))
    }
}

# Helper function for HTE prediction
calculate_hte <- function(preGroup_model, newdata) {
    raw_treatment_indicator <- preGroup_model$treatment_indicator
    # print(raw_treatment_indicator)
    # print(preGroup_model$pre_fit$data[, raw_treatment_indicator])
    levels <- tryCatch(sort(unique(as.numeric(levels(as.factor(preGroup_model$pre_fit$data[, raw_treatment_indicator])))), decreasing = TRUE),
                       error = function(e) {levels(as.factor(preGroup_model$pre_fit$data[, raw_treatment_indicator]))})

    # Check if enough levels exist
    if (length(levels) != 2) {
        stop("Only 2 different treatment levels should be used to compute HTE. The current model has ", length(levels), " levels.")
    }

    message("HTE will be calculated for treatment ", levels[1], " against treatment ", levels[2], ".\n")

    # Set up data for two groups
    # newdata_g1 <- within(newdata, {
    #     raw_treatment_indicator <- factor(levels[1], levels = levels)
    # })
    # newdata_g2 <- within(newdata, {
    #     raw_treatment_indicator <- factor(levels[2], levels = levels)
    # })
    newdata_g1 <- newdata
    newdata_g1[, raw_treatment_indicator] <- factor(as.numeric(levels[1]), levels = levels)
    # print(newdata_g1[, raw_treatment_indicator])
    # print(colnames(newdata_g1))

    newdata_g2 <- newdata
    newdata_g2[, raw_treatment_indicator] <- factor(as.numeric(levels[2]), levels = levels)
    # print(newdata_g2[, raw_treatment_indicator])
    # print(headings(newdata_g2))
    print("starting Hte calculation")

    # Calculate HTE
    hte <- predict(preGroup_model, newdata = newdata_g1, type = "response") - 
           predict(preGroup_model, newdata = newdata_g2, type = "response")
    return(hte)
}

# calculate_hte_pre <- function(pre_model, newdata, treatment_indicator) {

#     raw_treatment_indicator <- treatment_indicator
#     levels <- sort(unique(as.numeric(levels(pre_model$data[,raw_treatment_indicator]))), decreasing = TRUE)

#     # Check if enough levels exist
#     if (length(levels) < 2) {
#         stop("Not enough treatment levels to compute HTE.")
#     }

#     message("HTE will be calculated for treatment ", levels[1], " against treatment ", levels[2], ".\n")

#     # Set up data for two groups
#     newdata_g1 <- newdata
#     newdata_g1$raw_treatment_indicator <-  factor(levels[1], levels = levels)
#     print(newdata_g1$raw_treatment_indicator)
#     newdata_g2 <- newdata
#     newdata_g2$raw_treatment_indicator <-  factor(levels[2], levels = levels)
#     print(newdata_g2$raw_treatment_indicator)
#     # within(newdata, {
#     #     raw_treatment_indicator <- factor(levels[1], levels = levels)
#     # })
#     # newdata_g2 <- within(newdata, {
#     #     raw_treatment_indicator <- factor(levels[2], levels = levels)
#     # })
    
#     # Calculate HTE
#     hte <- predict(pre_model, newdata = newdata_g1, type = "response") - 
#            predict(pre_model, newdata = newdata_g2, type = "response")
#     return(hte)
# }

#' Generate Model Matrix from New Data for preGroup Model
#'
#' This function processes new data to fit the structure of a `preGroup` model,
#' facilitating the prediction process.
#'
#' @param preGroup_model A `preGroup` model object, specifically its `pre_fit` component.
#' @param new_data New data to be transformed into model matrix format.
#' @return A list containing the transformed new data as a model matrix.
#' @details The function retrieves settings and parameters from the `pre_fit` component of
#' the `preGroup` model and applies these to the new data to generate a suitable model matrix.
#' This involves applying transformation rules and handling different types of predictors as
#' specified in the model.
#' @export
#' @examples
#' # Assuming `model` is a preGroup model and `new_data` is ready
#' mod_matrix <- get_new_X(model, new_data)

get_new_X <- function(preGroup_model, new_data) {

    pre_model <- preGroup_model$pre_fit

    treatment_indicator <- preGroup_model$treatment_indicator
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

#' Calculate Variable Importances for a preGroup Model
#'
#' Computes variable importances based on a preGroup model, offering options for both global and local importance, 
#' standardization of importances, and various additional parameters.
#'
#' @param preGroup_model An object of class 'preGroup'.
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
#' @export
#' @examples
#' # Assuming `model` is a preGroup model
#' importance_data <- importance.preGroup(model)

importance.preGroup <- function(preGroup_model, standardize = FALSE, global = TRUE,
                           penalty.par.val = "lambda.1se", gamma = NULL,
                           quantprobs = c(.75, 1),
                           round = NA, plot = TRUE, ylab = "Importance",
                           main = "Variable importances", abbreviate = 10L,
                           diag.xlab = TRUE, diag.xlab.hor = 0, diag.xlab.vert = 2,
                           cex.axis = 1, legend = "topright", ...)
{

  if (!inherits(preGroup_model, what = "preGroup")) {
    stop("Specified object is not of class 'preGroup'.")
  }
  
  x <- preGroup_model$pre_fit

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
    coefs <- coef(preGroup_model, penalty.par.val = penalty.par.val)
  } else {
    coefs <- coef(preGroup_model, penalty.par.val = penalty.par.val, gamma = gamma)    
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
      varimps <- varimps[varimps$imp != 0, ]
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

          imp_plot <- create_importance_matrix(imp)

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
          barplot(height = varimps$imp, names.arg = plotnames, ylab = ylab,
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

#' Create Importance Matrix for Single Response Models
#'
#' Generates a matrix of importances for each variable based on the importances calculated in `importance.preGroup` 
#' for models with a single response variable.
#'
#' @param imp A list containing the importance data generated by `importance.preGroup`.
#' @return A matrix with variable names as columns and types of importance as rows (linear, prognostic, prescriptive).
#' @details This function organizes importances into a matrix format for easier visualization and analysis, 
#'          distinguishing between different types of importances.
#' @export
#' @examples
#' # Assuming `imp` is the list returned from `importance.preGroup`
#' importance_matrix <- create_importance_matrix(imp)

create_importance_matrix <- function(imp) {
  # Extract varnames and the importance values from varimps
  varnames <- imp$varimps$varname
  linear_imps <- imp$varimps$imp
  
  # Initialize the matrix with NA values
  importance_matrix <- matrix(NA, nrow = 3, ncol = length(varnames))
  rownames(importance_matrix) <- c("linear", "prognostic", "prescriptive")
  colnames(importance_matrix) <- varnames
  
  # Calculate the values for each column
  for (var in varnames) {
    # First value: the importance of the varname in varimps
    importance_matrix["linear", var] <- linear_imps[imp$varimps$varname == var]
    
    # Second value: sum of baseimps$imp for rules containing the varname and "prognostic"
    prognostic_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$description) & grepl("prognostic", imp$baseimps$rule_type)])
    importance_matrix["prognostic", var] <- prognostic_sum
    
    # Third value: sum of baseimps$imp for rules containing the varname and "prescriptive"
    prescriptive_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$description) & grepl("prescriptive", imp$baseimps$rule_type)])
    importance_matrix["prescriptive", var] <- prescriptive_sum
  }

  overall_importance <- colSums(importance_matrix)
  sorted_importance_matrix <- importance_matrix[, order(overall_importance, decreasing = TRUE)]

  return(sorted_importance_matrix)
}

#' Create Importance Matrices for Multivariate Response Models
#'
#' Generates a list of matrices of importances for each variable and each response in a multivariate model 
#' based on the importances calculated in `importance.preGroup`.
#'
#' @param imp A list containing the importance data generated by `importance.preGroup`.
#' @return A list of matrices, each corresponding to a different response variable in the model.
#' @details Each matrix in the list organizes importances for a specific response, formatted similarly to `create_importance_matrix`.
#' @export
#' @examples
#' # Assuming `imp` is the list returned from `importance.preGroup`
#' importance_matrices <- create_importance_matrix_multivariate(imp)

create_importance_matrix_multivariate <- function(imp) {
  # Get the response names by extracting column names that start with "importance."
  response_names <- grep("^importance\\.", names(imp$varimps), value = TRUE)
  response_names <- sub("^importance\\.", "", response_names)  # Remove "importance." prefix

  # Initialize a list to store the matrices for each response
  importance_matrices <- list()

  for (response in response_names) {
    # Extract varnames and the importance values for the current response
    varnames <- imp$varimps$varname
    importance_col <- paste0("importance.", response)
    linear_imps <- imp$varimps[[importance_col]]

    # Initialize the matrix for the current response
    importance_matrix <- matrix(NA, nrow = 3, ncol = length(varnames))
    rownames(importance_matrix) <- c("linear", "prognostic", "prescriptive")
    colnames(importance_matrix) <- varnames

    # Calculate the values for each column
    for (var in varnames) {
      # First value: the importance of the varname in varimps for the current response
      importance_matrix["linear", var] <- linear_imps[imp$varimps$varname == var]

      # Second value: sum of baseimps$imp for rules containing the varname and "prognostic" for the current response
      prognostic_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$description) & grepl("prognostic", imp$baseimps$rule_type)])
      importance_matrix["prognostic", var] <- prognostic_sum

      # Third value: sum of baseimps$imp for rules containing the varname and "prescriptive" for the current response
      prescriptive_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$description) & grepl("prescriptive", imp$baseimps$rule_type)])
      importance_matrix["prescriptive", var] <- prescriptive_sum
    }

    # Store the matrix in the list with the response name as the key
    importance_matrices[[response]] <- importance_matrix
  }
  return(importance_matrices)
}