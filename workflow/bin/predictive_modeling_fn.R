############################################################
## Function to extract logistic regression results
############################################################

extract_model_results <- function(
  model,
  model_name
) {

  coef_table <- summary(
    model
  )$coefficients

  coef_table <- coef_table[
    rownames(coef_table) != "(Intercept)",
    ,
    drop = FALSE
  ]

  results <- data.frame(

    model = model_name,

    predictor = rownames(
      coef_table
    ),

    beta = coef_table[
      ,
      "Estimate"
    ],

    SE = coef_table[
      ,
      "Std. Error"
    ],

    OR = exp(
      coef_table[
        ,
        "Estimate"
      ]
    ),

    CI_lower = exp(
      coef_table[
        ,
        "Estimate"
      ] -
        1.96 *
        coef_table[
          ,
          "Std. Error"
        ]
    ),

    CI_upper = exp(
      coef_table[
        ,
        "Estimate"
      ] +
        1.96 *
        coef_table[
          ,
          "Std. Error"
        ]
    ),

    p_value = coef_table[
      ,
      "Pr(>|z|)"
    ],

    row.names = NULL
  )

  return(
    results
  )
}


############################################################
## Model-based ROC analysis
############################################################


############################################################
## Helper function: compare nested models
############################################################

compare_nested_roc <- function(
    base_model,
    extended_model,
    model_name
) {

  ## Make sure both models were fitted on the same observations
  if (!identical(
    rownames(model.frame(base_model)),
    rownames(model.frame(extended_model))
  )) {
    stop(
      paste(
        "Base and extended models use different observations:",
        model_name
      )
    )
  }

  ## Observed response
  y <- model.response(
    model.frame(base_model)
  )

  ## Predicted probabilities
  pred_base <- predict(
    base_model,
    type = "response"
  )

  pred_extended <- predict(
    extended_model,
    type = "response"
  )

  ## ROC curves
  ##
  ## Because these are probabilities of response,
  ## larger probabilities should correspond to response = 1.
  roc_base <- pROC::roc(
    response = y,
    predictor = pred_base,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  roc_extended <- pROC::roc(
    response = y,
    predictor = pred_extended,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  ## AUC
  auc_base <- as.numeric(
    pROC::auc(roc_base)
  )

  auc_extended <- as.numeric(
    pROC::auc(roc_extended)
  )

  ## DeLong comparison of paired ROC curves
  delong_test <- pROC::roc.test(
    roc_base,
    roc_extended,
    paired = TRUE,
    method = "delong"
  )

  ## Likelihood-ratio test
  lrt <- anova(
    base_model,
    extended_model,
    test = "LRT"
  )

  lrt_p <- lrt$`Pr(>Chi)`[2]

  ## Return summary
  data.frame(
    comparison = model_name,
    AUC_base = auc_base,
    AUC_extended = auc_extended,
    delta_AUC = auc_extended - auc_base,
    DeLong_p = delong_test$p.value,
    LRT_p = lrt_p,
    stringsAsFactors = FALSE
  )
}
