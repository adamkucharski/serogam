# serogam
Non-parametric estimation of dynamics from serological data

*Note: this package is under development, so some functionality will not yet be stable*

### Quick start

This package uses the `scam` package for fitting smooth increasing seroprevalence curves (i.e. representing ongoing transmission, as well as the `mgcv` package for fitting piecewise seroprevalence curves (i.e. representing epidemics with periods of zero circulation in between). First, load this and other dependencies:

```
if(!require("devtools","readr","scam")) install.packages(c("devtools","readr","scam"))
library(devtools)
library(readr)
library(scam)
```

Next, install and load the development version of `serogam`:

```
install_github("adamkucharski/serogam")
library(serogam)
```

Next, load age-stratified dengue serosurvey data from [Kucharski et al, eLife 2018](https://elifesciences.org/articles/34848) and fit a model with smoothly varying transmission over time to the dengue-2 positivity data:

```
# Load and format data
data_fiji <- read.csv("https://raw.githubusercontent.com/adamkucharski/fiji-denv3-2014/master/data/serology_inputs.csv")
data_in <- data.frame(age = data_fiji$AGE_2015+2.5, # Use midpoint of age groups
                      outcome = data_fiji$DENV2P,
                      survey_year=rep(2015,nrow(data_fiji)) # Add survey year
                      )
                      
# Fit model and plot
# This model uses monotone increasing SCOP-splines to constrain seropositivity to increase with age
est_dynamics <- estimate_gam(data_in)

plot_estimate(est_dynamics,
	            survey_year=2015, # year in which ages reported
              y_range=c(0,10))

```

![plot_estimate1](https://github.com/adamkucharski/serogam/assets/8329046/da962b8e-7d50-45ec-8678-adaadc6d956f)

This approach allows us to estimate smooth trends in transmission over time, and suggests there was a historical peak in transmission around 1970 and 1996-97.

However, in Pacific Islands dengue is typically epidemic-prone, with large outbreaks followed by periods of limited circulation. As it happens, large dengue-2 outbreaks were reported in Fiji in [1971-3 and 1997-98](https://elifesciences.org/articles/34848). We can therefore use an adapted model that accounts for this epidemic dynamic. Rather than a smooth increase in seropositivity over time, we can specify the epidemic years and use a GAM with factors by age to constraint transmission to be zero between epidemic years, then estimate the size of the epidemics.


```
# Fit model and plot
# This model uses a piecewise risk function to constrain transmission to epidemic years
est_dynamics <- estimate_piecewise(data_in,
				   outbreak_years = c(1972,1997)
				   )

plot_estimate(est_dynamics,
              outbreak_years = c(1972,1997),
              survey_year=2015, # year in which ages reported
              y_range=c(0,50)
              )
```

![plot_estimate2](https://github.com/adamkucharski/serogam/assets/8329046/7eb71910-1cda-4c61-b42e-040abbf6e011)

A [post-epidemic serological survey in 1974](https://www.cambridge.org/core/journals/epidemiology-and-infection/article/mosquitoborne-infections-in-fiji-v-the-197173-dengue-epidemic/FA8AA20BB439CAB922974B498AE0F21C) concluded that around 26% of the Suva population were infected, compared to 16% in the above model (although there will be uncertainty around this estimate that is not yet implemented in this version of the analysis).
