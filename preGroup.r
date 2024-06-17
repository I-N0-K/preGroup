# require(rpart)
# require(partykit)
# require(glmnet)
require(pre)
require(mvs)
# utils::globalVariables("%dopar%")

# Problems with calling pre directly (tree.control is not defined). 
# the new problem is that the checking procedure is duplicated twice.

preGroup <- function(formula, data,  treatment_indicator, family = "gaussian", ad.alpha = NA, 
                ad.penalty = "lambda.min",
                use.grad = TRUE, weights, type = "both", sampfrac = .5, 
                maxdepth = 3L, learnrate = .01, mtry = Inf, ntrees = 500,
                confirmatory = NULL, singleconditions = FALSE,
                winsfrac = .025, normalize = TRUE, standardize = FALSE,
                ordinal = TRUE, nfolds = 10L, tree.control, tree.unbiased = TRUE, 
                removecomplements = TRUE, removeduplicates = TRUE, 
                verbose = FALSE, par.init = FALSE, par.final = FALSE, 
                sparse = FALSE, ...) {
                    
    ## Save call:
  cl <- match.call()

    # Check if the data is of type data.frame(
    if (!is.data.frame(data)) {
        stop("Argument 'data' must be of type data.frame")
    }

  ## Check if proper formula argument is specified 
  if (!(inherits(formula, "formula"))) {
    stop("Argument formula should specify and object of class 'formula'.\n")
  }
  ## Check if dot and functions are simultaneously used in formula
  form <- as.character(formula[3])
  for (i in names(data)) {
    form <- gsub(pattern = i, replacement = "", x = form)
  }
  if (any(grepl(".", form, fixed = TRUE))) {
    if (any(grepl("(", form, fixed = TRUE))) {
      if (any(grepl(")", form, fixed = TRUE))) {
        warning("Argument formula contains both one or more functions of predictor variables, as well as a dot ('.'), which should be avoided. Model fitting may fail, and/or both the original variable(s) and their functions may be included as predictor variables.\n", immediate. = TRUE)  
      }
    }
  }
  if (any(grepl("-", as.character(formula), fixed = TRUE))) {
    warning("Argument formula contains a minus sign. Note that the minus sign should not be used to omit the intercept or variables from the ensemble. To omit the intercept from the final ensemble, specify intercept = FALSE\n", immediate. = TRUE)
  }

  ## Check if proper data argument is specified:
  if (!is.data.frame(data)) {
    stop("Argument data should specify a data frame.\n")
  }

  ## Check and set up family argument: 
  if (length(family) > 1L) {
    warning("Argument family has length > 1, only first element will be used.\n")
    family <- family[1L]
  }
  if (is.function(family)) { 
    family <- family()
    if (inherits(family, "family")) {
      link <- family$link
      family <- family$family
      if (family == "gaussian" && link != "identity") {
        warning("The link function specified is currently not supported; identity link will be employed.\n", immediate. = TRUE)
      } else if (family == "binomial" && link != "logit") {
        warning("The link function specified is currently not supported; logit link will be employed.\n", immediate. = TRUE)
      } else if (family == "poisson" && link != "log") {
        warning("The link function specified is currently not supported; log link will be employed.\n", immediate. = TRUE)
      }
    } else {
      stop("Argument family should specify a family object.\n")
    }
  }
  if (is.character(family)) {
    if (!(family %in% c("gaussian", "binomial", "poisson", "mgaussian", "multinomial", "cox"))) {
      stop("Argument family should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
    }
  } else {
    stop("Argument family should be equal to 'gaussian', 'binomial', 'poisson', 'multinomial', 'mgaussian', 'cox', or a corresponding family object.\n")
  }
  
  ## Check if proper weights argument is specified, if specified:
  if (missing(weights)) {
    weights <- rep(1L, times = nrow(data))
  } else {
    if (!is.numeric(weights)) {
      warning("Argument weights should, but does not specify a vector of numeric values.\n", immediate. = TRUE)
    }
    if (length(weights) != nrow(data)) {
      warning("Length of argument weights should be, but is not equal to nrow(data).\n", immediate. = TRUE)
    }
    if (!all(weights > 0)) {
      warning("All weights should be positive but some are not.\n", immediate. = TRUE)
    }
  }
  
  ## Check if proper type argument is specified:
  if (!(length(type) == 1L && type %in% c("rules", "both", "linear"))) {
    stop("Argument type should be 'rules', 'linear' or 'both'.\n")
  }
  
  ## Check if proper sampfrac argument is specified:
  if (!is.function(sampfrac)) {
    if (!(length(sampfrac) == 1 && is.numeric(sampfrac) && sampfrac > 0 && 
          sampfrac <= 1)) {
      stop("Argument sampfrac should be a single numeric value > 0 and <= 1, or a sampling function.\n")
    }
  }

  ## Check if proper maxdepth argument is specified:
  if (is.function(maxdepth)) {
    maxdepth <- maxdepth(ntrees = ntrees)
  } else if (!is.numeric(maxdepth)) {
    stop("Argument maxdepth should be either a numeric vector of length 1 or ntrees, or a random number generating function.\n")
  } else if (!(length(maxdepth) %in% c(1L, ntrees))) {
    warning("Argument maxdepth should be either a numeric vector of length 1 or ntrees, only first value of maxdepth will be used.\n")
    maxdepth <- maxdepth[1L]
  } 
  if (!all(maxdepth > 0)) {
    stop("All values of maxdepth should be > 0.\n")
  } 
  if (!all(maxdepth == suppressWarnings(as.integer(maxdepth)) | is.infinite(maxdepth))) {
    stop("Argument maxdepth should consist of  of integer values or Inf), or a random number generating function.\n")
  }
  
  ## Check if proper learnrate argument is specified:
  if (!(length(learnrate) == 1L && is.numeric(learnrate) && 
        (learnrate >= 0 || learnrate <= 1))) {
    stop("Argument learnrate shoud be a single numeric value >= 0 and <= 1.\n")
  }
  
  ## Check if proper mtry argument is specified:
  if (!(length(mtry) == 1L && mtry > 0 && 
        (mtry == suppressWarnings(as.integer(mtry)) || is.infinite(mtry)))) {
    stop("Argument mtry should be a single integer value, or Inf.\n")
  }
  
  ## Check if proper ntrees argument is specified:
  if (!(length(ntrees) == 1L && ntrees == as.integer(ntrees) && ntrees > 0)) {
    stop("Argument ntrees should be a single positive integer.\n")
  }
  
  ## Check if proper winsfrac argument is specified:
  if (!(length(winsfrac == 1L) && is.numeric(winsfrac) && winsfrac >= 0 && 
        winsfrac < 1)) {
    stop("Argument winsfrac should be a numeric value >= 0 and < 1.\n")
  }

  ## Check if proper nfolds argument is specified:
  if (!(length(nfolds) == 1L && nfolds > 0 && nfolds == as.integer(nfolds))) {
    stop("Argument nfolds should be a positive integer.\n")
  }
  
  ## Check if proper confirmatory argument is specified:
  if (!is.null(confirmatory)) {
    if (!is.character(confirmatory)) {
      stop("Argument confirmatory should specify a character vector.\n")
    }
  }
  
  ## Check if logical arguments of length 1 are properly specified:
  is_logical_and_length_one <- function(x) {is.logical(x) && length(x) == 1L}
  for (i in c(use.grad, removeduplicates, removecomplements, normalize, 
              standardize, ordinal, verbose, tree.unbiased, par.init, 
              par.final)) {
    if (!is_logical_and_length_one(i)) {
      stop("Argument ", i, "should be TRUE or FALSE.\n")
    }
  }
  
  if (par.final || par.init) {
    if(!requireNamespace("foreach")) {
      warning("Parallel computation requires package foreach. Arguments par.init and par.final are set to FALSE.\n")   
      par.init <- par.final <- FALSE
    }
  }

  ## Check if proper tree.control argument is specified:
  if (missing(tree.control)) {
    if (tree.unbiased && use.grad) {
      tree.control <- ctree_control(maxdepth = maxdepth[1], mtry = mtry)
    } else if (tree.unbiased && !use.grad) {
      tree.control <- mob_control(maxdepth = maxdepth[1] + 1, mtry = mtry)
    } else if (!tree.unbiased) {
      if (any(maxdepth > 29)) {
        maxdepth[maxdepth > 29] <- 29L
        warning("If tree.unbiased = FALSE, max(maxdepth) is 29.\n")
      }
      tree.control <- rpart.control(maxdepth = maxdepth[1L])
      if (!is.infinite(mtry)) {
        warning("Value specified for mtry will be ignored if tree.unbiased = FALSE.\n")
      }
    }
  } else {
    if (!is.list(tree.control)) {
      stop("Argument tree.control should specify a list of control parameters.\n")
    }
    if (use.grad && tree.unbiased) {
      if (!setequal(names(ctree_control()), names(tree.control))) {
        args_not_in_ctree_control <- names(tree.control)[which(
          !names(tree.control) %in% names(ctree_control()))]
        warning("Argument tree.control contains named elements: ", 
                paste0(args_not_in_ctree_control, collapse = ", "), 
                ", which should probably not be there, see ?ctree_control.\n")
      }
    } else if (!use.grad && tree.unbiased) { 
      if (!setequal(names(mob_control()), names(tree.control))) {
        args_not_in_glmtree_control <- names(tree.control)[which(
          !names(tree.control) %in% names(mob_control()))]
        warning("Argument tree.control has named elements ", 
                 paste0(args_not_in_glmtree_control, collapse = ", "), 
                ", which should probably not be there, see ?mob_control.\n")
      }
    } else if (!tree.unbiased) {
      if (!setequal(names(rpart.control()), names(tree.control))) {
        warning("Argument tree.control should be a list containing names elements ", paste(names(rpart.control()), collapse = ', '), "\n")
      }
    }
    if (use.grad) { ## if ctree or rpart are employed:
      tree.control$maxdepth <- maxdepth[1L]      
    } else if (tree.unbiased) { ## if glmtree is employed:
      tree.control$maxdepth <- maxdepth[1L] + 1L
    }
    if (tree.unbiased) {
      tree.control$mtry <- mtry
    } else if (mtry != Inf) {
      warning("Argument tree.unbiased was set to FALSE, so rpart is employed for tree induction, and value specified for mtry will be ignored.\n")
      mtry <- Inf
    }
  }
  
  if (!tree.unbiased && !use.grad && learnrate > 0) {
    stop("Employing the rpart algorithm with a learnrate > 0 without gradient boosting is not supported.\n")
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
    mvs_fit <- MVS(x = as.matrix(pre_fit$modmat), y = pre_fit$data[, 1], alpha = c(1,0), nnc = c(0,1),
               view = as.matrix(group_id), family = family, cvloss = "mse")
    
    # now convert the group id to meaningful rule types
    rule_type <- c("Intercept", group_id)
    rule_type[rule_type == 1] <- "linear"
    rule_type[rule_type == 2] <- "prognostic"
    rule_type[rule_type == 3] <- "prescriptive"
    
    result <- list(mvs_fit = mvs_fit, pre_fit = pre_fit, call = match.call(), rule_type = rule_type)

  
  class(result) <- "preGroup"
  return(result)
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
    rule_count <- sum(!rule_presence)

    group_ids <- c(rep(1, rule_count), split_treatment_indicator_rule(pre_model$rules$description, treatment_indicator))

}

coef.preGroup <- function(preGroup_model, ...) {

    if(!class(preGroup_model) == "preGroup") {
        stop("The model must belong to class 'preGroup'. ")
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
        subgroup_coef <- cbind(rule = x_names, coefficients = x_coefs)

        coefs <- rbind(coefs, subgroup_coef)

    }
    grand_intercept <- cbind(rule = "(Intercept)", coefficients = as.numeric(grand_intercept))

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

    return(coefs)
}



importance.preGroup <- function(x, standardize = FALSE, global = TRUE,
                           penalty.par.val = "lambda.1se", gamma = NULL,
                           quantprobs = c(.75, 1),
                           round = NA, plot = TRUE, ylab = "Importance",
                           main = "Variable importances", abbreviate = 10L,
                           diag.xlab = TRUE, diag.xlab.hor = 0, diag.xlab.vert = 2,
                           cex.axis = 1, legend = "topright", ...)
{

  if (!inherits(x, what = "preGroup")) {
    stop("Specified object is not of class 'preGroup'.")
  }
  
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
    coefs <- coef(x, penalty.par.val = penalty.par.val)
  } else {
    coefs <- coef(x, penalty.par.val = penalty.par.val, gamma = gamma)    
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
      keep <- c("rule", "description", gsub("coefficient", "importance", coef_inds), 
                coef_inds, "sd")
    } else {
      keep <- c("rule", "description", "imp", "coefficient", "sd")
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
    prognostic_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$rule) & grepl("prognostic", imp$baseimps$rule)])
    importance_matrix["prognostic", var] <- prognostic_sum
    
    # Third value: sum of baseimps$imp for rules containing the varname and "prescriptive"
    prescriptive_sum <- sum(imp$baseimps$imp[grepl(var, imp$baseimps$rule) & grepl("prescriptive", imp$baseimps$rule)])
    importance_matrix["prescriptive", var] <- prescriptive_sum
  }

  return(importance_matrix)
}

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
      prognostic_sum <- sum(imp$baseimps[[paste0("importance.", response)]][grepl(var, imp$baseimps$rule) & grepl("prognostic", imp$baseimps$rule)])
      importance_matrix["prognostic", var] <- prognostic_sum

      # Third value: sum of baseimps$imp for rules containing the varname and "prescriptive" for the current response
      prescriptive_sum <- sum(imp$baseimps[[paste0("importance.", response)]][grepl(var, imp$baseimps$rule) & grepl("prescriptive", imp$baseimps$rule)])
      importance_matrix["prescriptive", var] <- prescriptive_sum
    }

    # Store the matrix in the list with the response name as the key
    importance_matrices[[response]] <- importance_matrix
  }

  return(importance_matrices)
}