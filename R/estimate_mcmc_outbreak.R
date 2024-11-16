#' Estimate sizes of known outbreaks from seropositivity
#'
#' @description Uses information on outbreak time to fit a model with age specific variation in risk
#'
#' @param data_in A `<data.frame>` containing individual-level serological data. The required columns
#' are "age" (in years at time of survey), "survey_year" (year) and "outcome" (0 or 1), representing a seropositive/seronegative result 
#' according to the tested biomarker.
#' 
#' @param outbreak_years A `<vector>` containing the years in which historical outbreaks occurred.
#' 
#' @return A `<data.frame>` with the probability of seropositivity by age and 95 CI, as well as annual proportion of the 
#' population infected by year (where "age" represents the number of years into the past).
#'
#'
#' @export
#'

estimate_mcmc_outbreak <- function(data_in,
                                   outbreak_years,
                                   age_band
                                   ) {
  
  # DEBUG: age_band = c(0,20);  outbreak_years = c(1972,1997)
  
  # input checking
  checkmate::assert_data_frame(
    data_in,
    min.rows = 1, min.cols = 3
  )
  # check that input `<data.frame>` has columns age and output
  checkmate::assert_names(
    colnames(data_in),
    must.include = c("age", "outcome","survey_year")
  )

  checkmate::assert_vector(
    outbreak_years
  )
  
  # Simulation model
  run_model <- function(data_in,
                        attack_outbreaks,
                        age_risk,
                        age_band = c(0,20)){
    
    # Get age range in data
    age_upper <- ceiling(max(data_in$age))
    age_range <- 0:age_upper
    n_age <- length(age_range)
    
    # Allocate ages to risk groups
    age_group <- findInterval(age_range,age_band)

    # Set up storage of simulation results for serosurveys
    survey_sim_range <- (min(outbreak_years):max(data_in$survey_year))
    n_survey_yr <- length(survey_sim_range)
    
    store_seroprevalence <- data.frame(matrix(0, nrow = n_survey_yr, 
                                                  ncol = n_age),
                                       survey_year = survey_sim_range)
                                       
    # Run annual model
    for(ii in 1:n_survey_yr){
      
      # Get survey year
      survey_yr_ii <- survey_sim_range[ii]
      
      # Check for outbreak in survey year
      match_outbreak <- match(survey_yr_ii,outbreak_years)
      
      if(!is.na(match_outbreak)){ # Add outbreak if relevant
        # Get year of outbreak
        outbreak_yr_ii <- outbreak_years[match_outbreak]
        
        # Calculate proportion susceptible
        seropositive <- store_seroprevalence[ii, 1:n_age]
        susceptible <- 1 - seropositive

        # Calculate age-specific attack rate by group
        age_attack <- age_risk*attack_outbreaks[match_outbreak]
        
        new_infection_group <- age_attack[age_group]*susceptible
        
        # Add to seroprevalence
        store_seroprevalence[ii, 1:n_age] <- store_seroprevalence[ii, 1:n_age] + new_infection_group
        
        # Move demography along 1 year
        if(ii < n_survey_yr){ii
          store_seroprevalence[ii+1, 1:n_age] <- c(0, store_seroprevalence[ii, 1:(n_age-1)])
        }
      
      }else{ # Else move along 1 year
        if(ii < n_survey_yr){
          store_seroprevalence[ii+1, 1:n_age] <- c(0, store_seroprevalence[ii, 1:(n_age-1)])
        }
      }
      
    }
    
    store_seroprevalence
  
  }
  
  
  
  # plot(as.numeric(store_seroprevalence[44,1:n_age]))
  
  # Define likelihood
  logit_prior <- function(x){ifelse(x > 0 & x < 1, 1, 0)}
  pos_prior <- function(x){ifelse(x > 0, 1, 0)}
  
  survey_likelihood <- function(param,
                                data){
    
    # Extract parameters
    false_prob <- 1e-5 # Add term for specificity
    
    # DEBUG: Need to add flexibility here
    attack_outbreaks <- param[1:2]
    attack_outbreaks <- pmax(0,pmin(1,attack_outbreaks)) # Constrain 0<=x<=1 (NOTE: APPROX)
    
    age_risk <- c(1,1) #param[3]) DEBUG
    age_risk <- pmax(0,age_risk) # Constrain 0<=x<=1 (NOTE: APPROX)
    
    #print(age_risk)
    #print(attack_outbreaks)
    
    # Get age range in data
    age_upper <- ceiling(max(data_in$age))
    age_range <- 0:age_upper
    n_age <- length(age_range)
    
    # Run model
    model_out <- run_model(data_in,
                           attack_outbreaks,
                           age_risk,
                           age_band = c(0,20))
    
    # DEBUG: Need to add flexibility here on age_band
    
    
    # Get serosurvey years
    serosurvey_years <- unique(data_in$survey_year)
    
    log_likelihood_store <- 0
    
    # Loop over survey years and add up likelihood
    for(ii in 1:length(sero_years)){
      survey_ii <- serosurvey_years[ii]
      
      # Get simulation and data for that survey year
      get_sim_results <- model_out |> dplyr::filter(survey_year==survey_ii)
      
      get_survey_data_pos <- data_in |> dplyr::filter(survey_year==survey_ii & outcome==1)
      get_survey_data_neg <- data_in |> dplyr::filter(survey_year==survey_ii & outcome==0)
      
      print(sum(is.na(get_sim_results)))
      print(sum(get_sim_results==0))
      
      # Get attack rates by age, calculate Bernouilli probability for
      # positives and negatives and add to overall likelihood
      match_data_ages_pos <- match(round(get_survey_data_pos$age),age_range)
      log_lik_pos_ii <- sum(log(get_sim_results[match_data_ages_pos] + false_prob))
      
      match_data_ages_neg <- match(round(get_survey_data_neg$age),age_range)
      log_lik_neg_ii <- sum(log(1-get_sim_results[match_data_ages_neg]+false_prob))
      
      # NOTE: ADD option for false positives here?
      
      # Add prior
      # log_prior <- sum(log(logit_prior(attack_outbreaks))) + 
      #              sum(log(pos_prior(age_risk))) 
      
      # Calculate likelihood
      log_likelihood_store <- log_likelihood_store + log_lik_pos_ii + log_lik_neg_ii 
                              #log_prior
      
      #print(log_likelihood_store)
      
      # if(log_likelihood_store == "NaN" | is.na(log_likelihood_store)){
      #   log_likelihood_store = -Inf
      # }
    }
    
    log_likelihood_store
    
  }
  
  survey_likelihood(c(0.3,0.3,1.5), data_in)
  
  
  initial_param <- c(attack_outbreaks = c(0.1,0.1), # Initial guess for outbreaks
                     age_risk = c(1.5)) # Initial guess for age risk
  
  n_mcmc <- 50
  
  
  result_mcmcpack <- MCMCpack::MCMCmetrop1R(
    fun = survey_likelihood, 
    theta.init = initial_param,
    mcmc = n_mcmc, # Number of MCMC iterations
    burnin = 0, # Burn-in period
    verbose = FALSE, # Turn off verbose output
    data = data_in
  )
  
  
  ## OUTPUTS
  # Calculate effective sample size (i.e. measure of MCMC mixing)
  ess_mcmcpack <- effectiveSize(result_mcmcpack)
  
  # Plot posterior estimates
  plot(result_mcmcpack)
  
  # Define helper function to calculate median and 95% credible interval from data.frame of MCMC samples
  get_param <- function(x){
    apply(x,2,function(y){val = signif(quantile(y,c(0.5,0.025,0.975)),3);
    val_text <- paste0(val[1]," (95%: CrI: ",val[2],"-",val[3],")")})
  }
  
  # Get posterior median and 95% CrI
  posterior_estimates <- get_param(result_mcmcpack)
  
  # Compile table
  results_table <- data.frame(
    Package = "MCMCpack",
    Posterior_R = posterior_estimates[1],
    Posterior_k = posterior_estimates[2],
    ESS_R = ess_mcmcpack[1],
    ESS_k = ess_mcmcpack[2]
  )
  
  

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

