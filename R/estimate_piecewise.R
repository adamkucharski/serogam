#' Estimate step-wise changes in seropositivity and infection risk over time
#'
#' @description Fits a generalised additive model with step-wise temporal risk (i.e. multiple outbreaks) to binary 
#' serological data and returns predictions about seropositivity by age and annual risk of historical infection
#'
#' @param data_in A `<data.frame>` containing individual-level serological data. The required columns
#' are "age" (in years) and "outcome" (0 or 1), representing a seropositive/seronegative result 
#' according to the tested biomarker.
#' 
#' @param outbreak_years A `<vector>` containing the years in which historical outbreaks occurred.
#' 
#' @param survey_year A `<number>` giving the year in which the ages were reported.
#' 
#' @return A `<data.frame>` with the probability of seropositivity by age and 95 CI, as well as annual proportion of the 
#' population infected by year (where "age" represents the number of years into the past).
#'
#' @export
#'

estimate_piecewise <- function(data_in,
                         outbreak_years,
                         survey_year) {
  
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

  checkmate::assert_vector(
    outbreak_years
  )
  
  checkmate::assert_number(
    survey_year,
    lower = 1800, finite = TRUE
  )

  # convert inputs into age breaks
  knot_ages <- sort(survey_year-outbreak_years) # replace with your actual ages for knots

  # convert age into a factor with levels indicating the interval between knots
  df <- data_in
  df$age_group <- cut(df$age, breaks = c(-Inf, knot_ages, Inf), include.lowest = TRUE, labels = FALSE)
  
  
  # piecewise fit
  b_model <- gam(outcome~factor(age_group),family=binomial(link="logit"),data=df)
  
  # generate prediction
  age_range_smooth <- round(min(data_in$age)):round(max(data_in$age))
  new_dat_group <- cut(age_range_smooth, breaks = c(-Inf, knot_ages, Inf), include.lowest = TRUE, labels = FALSE)
  new_dat_group <- data.frame(age_group=new_dat_group)
  
  pred_model <- predict(b_model,new_dat_group,type="response",se=TRUE)

  # calculate 95% confidence intervals
  ci_lower <- pred_model$fit - 1.96 * pred_model$se.fit
  ci_upper <- pred_model$fit + 1.96 * pred_model$se.fit
  
  # combine predictions and confidence intervals
  fit_output <- data.frame(
    age = age_range_smooth,
    prob_positive_mid = pred_model$fit,
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

