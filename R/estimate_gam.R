#' Estimate smooth changes in seropositivity and infection risk over time
#'
#' @description Fits a generalised additive model with monotone increasing SCOP-splines to binary 
#' serological data and returns predictions about seropositivity by age and annual risk of historical infection
#'
#' @param data_in A `<data.frame>` containing individual-level serological data. The required columns
#' are "age" (in years) and "outcome" (0 or 1), representing a seropositive/seronegative result 
#' according to the tested biomarker.
#' 
#' @return A `<data.frame>` with the probability of seropositivity by age and 95 CI, as well as annual proportion of the 
#' population infected by year (where "age" represents the number of years into the past).
#'
#' @export
#'

estimate_gam <- function(data_in) {
  
  # input checking
  checkmate::assert_data_frame(
    data_in,
    min.rows = 1, min.cols = 2
  )
  # check that input `<data.frame>` has columns age and output
  checkmate::assert_names(
    colnames(data_in),
    must.include = c("age", "outcome")
  )
  
  # check for any NAs among data
  checkmate::assert_data_frame(
    data_in[, c("age", "outcome")],
    any.missing = FALSE
  )

  # fit monotonoic gam to data
  b_model <- scam(outcome~s(age,bs="mpi"),family=binomial(link="logit"),data=data_in)
  
  # generate prediction across observed range, rounded for fitting
  age_range_smooth <- round(min(data_in$age)):round(max(data_in$age))
  
  new_dat <- data.frame(age=age_range_smooth)
  
  pred_model_smooth <- predict(b_model,new_dat,type="response",se=TRUE)
  
  # calculate 95% confidence intervals
  ci_lower <- pred_model_smooth$fit - 1.96 * pred_model_smooth$se.fit
  ci_upper <- pred_model_smooth$fit + 1.96 * pred_model_smooth$se.fit
  
  # combine predictions and confidence intervals
  fit_output <- data.frame(
    age = age_range_smooth,
    prob_positive_mid = pred_model_smooth$fit,
    prob_positive_lower = ci_lower,
    prob_positive_upper = ci_upper
  )
  
  # estimate annual infection risk
  vec1 <- head(fit_output$prob_positive_mid,-1) # drop last entry
  vec2 <- tail(fit_output$prob_positive_mid,-1) # drop first entry
  infection_risk <- attack_est(vec1,vec2)
  infection_risk <- c(NA,infection_risk) # pad left hand side to match ages
  
  fit_output$annual_risk <- infection_risk
  
  # return the severity estimate
  fit_output
}



# Calculate infection risk per year
attack_est <- function(pos_1, # positive at age i
                       pos_2 # positive at age i+1
){
  # Calculate remaining susceptibility at age i
  susceptible_1 <- (1-pos_1)
  
  # Calculate proportion of susceptible group infected
  inf_1 <- (pos_2-pos_1)/susceptible_1
  
  return(inf_1)
}

