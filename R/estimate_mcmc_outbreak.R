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
#' @param age_band A `<vector>` containing the cutoffs for age-specific attack rates
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
  
  # Load Lazy MCMC
  devtools::install_github("jameshay218/lazymcmc")
  # NEED TO MOVE THIS TO A DEPENDENCY
  
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

  # Check outbreak years is a vector
  checkmate::assert_vector(
    outbreak_years
  )

  # Load and format data
  data_fiji <- read_csv("https://raw.githubusercontent.com/adamkucharski/fiji-denv3-2014/master/data/serology_inputs.csv")
  data_in <- data.frame(age = data_fiji$AGE_2015+2.5, outcome = data_fiji$DENV2P, survey_year = 2015) # Use midpoint of age groups
  

  # Run MCMC functions --------------------------------------------------
  
  # Define parameters
  ## Update the proposal step size every X iterations (opt_freq) for the first X iterations 
  ## (adaptive_period) to aim for an acceptance rate of 0.44 (popt). After the adaptive period, 
  ## run for a further X (iterations) steps. Save every nth rows, where n is "thin" (ie. 1 here).
  ## Write to disk every X iterations (post thinning). Note that the adaptive period is also saved
  ## to disk
  mcmcPars <- c("iterations"=1000,"popt"=0.44,"opt_freq"=20,
                "thin"=1,"adaptive_period"=200,"save_block"=100)
  
  ## The MCMC code uses the parameter table. Here, we should specify some random starting
  ## points in the "values" column.
  startTab <- data.frame(values=c(0.5,0.5,1),
                         names=c("attack_outbreaks1","attack_outbreaks2","age_risk"),
                         fixed=0,
                         lower_bound=c(0,0,0),
                         upper_bound=c(1,1,2),
                         steps=rep(0.1,3),
                         stringsAsFactors=FALSE)
  
  my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
    
    # Define liklihood function:
    f <- function(pars){
      survey_likelihood(param = pars, data)
    }
    
    return(f)
  }
  
  # Run MCMC fitting
  output <- run_MCMC(parTab=startTab, data=data_in, mcmcPars=mcmcPars, filename="test", 
                     CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=NULL, 
                     PRIOR_FUNC = my_prior, OPT_TUNING=0.2)
  
  # Plot data and binom CI --------------------------------------------------
  
  chain <- read.csv(output$file)

  plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
  
  # Simulate 20 random draws from the posterior (run multiple to get uncertainty)
  age_upper <- ceiling(max(data_in$age)); age_range <- 0:age_upper
  
  
  plot(age_range,-100+0*age_range,ylim=c(0,1),ylab="seropositive",xlab="age")
  
  for(ii in 1:20){
    pick_rand <- sample(1:max(chain$sampno),1)
    param <- chain[,startTab$names][pick_rand,]
    output_model <- run_model(data_in, attack_outbreaks = as.numeric(param), age_risk,age_band = c(0,20))
    sim_postitivity <- output_model |> filter(survey_year == 2015)
    
    lines(age_range,sim_postitivity[1:length(age_range)],col=rgb(0,0,1,0.2))
  }
  
  # Calculate total tallies and outcome = 1 tallies
  tallies <- data_in %>%
    group_by(age) %>%
    summarise(
      total_tallies = n(),
      outcome_1_tallies = sum(outcome == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Convert to vectors
  age_vector <- tallies$age
  total_tallies_vector <- tallies$total_tallies
  outcome_1_vector <- tallies$outcome_1_tallies
  
  plot_CI(age_vector, xx=outcome_1_vector, nn=total_tallies_vector)
  
  
  

}


# Simulation model --------------------------------------------------------

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


# Likelihood function -----------------------------------------------------

survey_likelihood <- function(param, # vector of parameters
                              data # input data
                              ){
  
  # DEBUG: param <- c(0.5,0.5)
  
  # Extract parameters
  false_prob <- 1e-5 # Add small error for specificity (to avoid invalid likelihood)
  
  # DEBUG: Need to add flexibility here
  attack_outbreaks <- param[1:2]
  attack_outbreaks <- pmax(0,pmin(1,attack_outbreaks)) # Constrain 0<=x<=1 (NOTE: APPROX)
  
  age_risk <- c(1, param[3]) # relative age risk
  #age_risk <- pmax(0,age_risk) # Constrain 0<=x<=1 (NOTE: APPROX)
  
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
  for(ii in 1:length(serosurvey_years)){
    survey_ii <- serosurvey_years[ii]
    
    # Get simulation and data for that survey year
    get_sim_results <- model_out |> dplyr::filter(survey_year==survey_ii)
    
    get_survey_data_pos <- data_in |> dplyr::filter(survey_year==survey_ii & outcome==1)
    get_survey_data_neg <- data_in |> dplyr::filter(survey_year==survey_ii & outcome==0)
    
    # print(sum(is.na(get_sim_results)))
    # print(sum(get_sim_results==0))
    
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


# Calculate infection risk per year ---------------------------------------

attack_est <- function(pos_1, # positive at age i
                       pos_2 # positive at age i+1
){
  # Calculate remaining susceptibility at age i
  susceptible_1 <- (1-pos_1)
  
  # Calculate proportion of susceptible group infected
  inf_1 <- (pos_2-pos_1)/susceptible_1
  
  return(inf_1)
}


# Calculate CI ------------------------------------------------------------


plot_CI <- function(ages,xx,nn,colA="black") {
  
  # Remove NA entries in denominators if needed
  ages_p <- ages[!is.na(nn) & nn>0]
  xx_p <- xx[!is.na(nn) & nn>0] %>% round()
  nn_p <- nn[!is.na(nn) & nn>0] %>% round()
  
  for(ii in 1:length(nn_p)){
    test_r <- binom.test(xx_p[ii],nn_p[ii])
    CI1 <- as.numeric(test_r$conf.int)[1]
    CI2 <- as.numeric(test_r$conf.int)[2]
    points(ages_p[ii],xx_p[ii]/nn_p[ii],col=colA,pch=19);
    lines(c(ages_p[ii],ages_p[ii]),c(CI1,CI2),col=colA)
  }
  
}

