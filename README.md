# serogam
Non-parameteric estimation of dynamics from serological data

### Quick start

First install and load the development version of `serogam`:

```
if(!require("devtools","readr")) install.packages(c("devtools","readr"))
library(devtools)

install_github("adamkucharski/serogam")
library(serogam)
```

Next, load age-stratified dengue serosurvey data from [Kucharski et al, eLife 2018](https://elifesciences.org/articles/34848) and fit a model with smoothly varying transmission over time to dengue-3 positivity prior to the large 2013/14 epidemic:

```
# Load and format data
data_fiji <- read_csv("https://raw.githubusercontent.com/adamkucharski/fiji-denv3-2014/master/data/serology_inputs.csv")
data_in <- data.frame(age = data_fiji$AGE_2015+2.5, outcome = data_fiji$DENV2P) # Use midpoint of age groups

# Fit model and plot
# This model uses monotone increasing SCOP-splines to constrain seropositivity to increase with age
est_dynamics <- estimate_gam(data_in)

plot_estimate(est_dynamics,
			  survey_year=2015, # year in which ages reported
			  y_range=c(0,5))

```

![plot_estimate1](https://github.com/adamkucharski/serogam/assets/8329046/da962b8e-7d50-45ec-8678-adaadc6d956f)

This approach allows us to estimate smooth trends in transmission over time, and suggests there was a historical peak in transmission around 1970 and 1996-97.

However, in Pacific Islands dengue is typically epidemic-prone, with large outbreaks followed by periods of limited circulation. As it happens, large dengue-2 outbreaks were reported in Fiji in [1971-3 and 1997-98](https://elifesciences.org/articles/34848). We can therefore use an adapted model that accounts for this epidemic dynamic. Rather than a smooth increase in seropositivity over time, we can specify the epidemic years and use a GAM with factors by age to constraint transmission to be zero between epidemic years, then estimate the size of the epidemics.


```
# Fit model and plot
# This model uses a piecewise risk function to constrain transmission to epidemic years
est_dynamics <- estimate_piecewise(data_in,
								outbreak_years = c(1972,1997),
								survey_year = 2015 
								)

plot_estimate(est_dynamics,
			  outbreak_years = c(1972,1997),
			  survey_year=2015, # year in which ages reported
			  y_range=c(0,50)
			  )

```

![plot_estimate2](https://github.com/adamkucharski/serogam/assets/8329046/7eb71910-1cda-4c61-b42e-040abbf6e011)
