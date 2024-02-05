#' Plot seropositivity and infection risk over time
#'
#' @description Plots an output object from `estimate_gam` or `estimate_piecewise`
#'
#' @param data_in A `<data.frame>` containing the output of `estimate_gam` or `estimate_piecewise`
#' 
#' @param outbreak_years An optional `<vector>` containing the years in which historical outbreaks occurred.
#' 
#' @param survey_year An optional `<number>` giving the year in which the ages were reported.
#' 
#' @param y_range An optional `<vector>` giving alternative y axis range to plot infection risk on.
#' 
#' @param write_file An optional `<string>` giving the file path to write the plots to.
#' 
#' @return A plot with the probability of seropositivity by age and 95 CI, as well as annual proportion of the 
#' population infected by year.
#'
#' @export
#'

plot_estimate <- function(fit_in,
                          outbreak_years=NULL,
                          survey_year=NULL,
                          write_file=NULL,
                          y_range=c(0,100)) {
  
  # input checking
  checkmate::assert_data_frame(
    fit_in,
    min.rows = 1, min.cols = 5
  )
  # check that input `<data.frame>` has columns age and output
  checkmate::assert_names(
    colnames(fit_in),
    must.include = c("age", "prob_positive_mid", "prob_positive_lower","prob_positive_upper","annual_risk")
  )
  
  
  # set up plot
  par(mfrow=c(2,1),mgp=c(2,0.7,0),mar = c(3.5,3.5,1,1))
  
  # plot seropositivity
  plot(fit_in$age,fit_in$prob_positive_mid,col="white",xlab="age",ylab="% seropositive",ylim=c(0,100),main="seropositivity by age")

  lines(fit_in$age,100*fit_in$prob_positive_mid,lwd=2,col="dark blue")
  polygon(c(fit_in$age,rev(fit_in$age)),100*c(fit_in$prob_positive_lower,rev(fit_in$prob_positive_upper)),
          col=rgb(0,0,1,0.1),border=NA)
  
  if(!is.null(survey_year) & !is.null(outbreak_years)){
    abline(v = (survey_year - outbreak_years), lty = 2, col = "dark gray")
  }
  
  # plot historical infection risk
  if(!is.null(survey_year) ){
    xx_val <- survey_year-fit_in$age +1
    xx_lab <- "year"
  }else{
    xx_val <- fit_in$age
    xx_lab <- "age"
  }
  
  plot(xx_val,fit_in$annual_risk,col="white",xlab=xx_lab,ylab="% infected",ylim=y_range,main="infection risk over time")
  lines(xx_val,100*fit_in$annual_risk,lwd=2,col="dark blue")

  if(!is.null(survey_year)){
    abline(v = outbreak_years, lty = 2, col = "dark gray")
  }

  if(!is.null(write_file)){
    dev.copy(png,paste0(write_file,"plot_estimate.png"),units="cm",width=10,height=15,res=150)
    dev.off()
  }
  
}

