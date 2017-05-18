f <- function(data, i,formula){
  d2 <- data[i,]
  model<-glm(recidiven~cote,family="binomial",data=d2)
  return(model$coefficients[2])
}

bootcorr <- boot(data=cryo, statistic=f, R=500,formula=recidiven~cote)
bootcorr
boot.ci(bootcorr, type = "bca")


# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}
# bootstrapping with 1000 replications
results <- boot(data=mtcars, statistic=rsq,
                R=1000, formula=mpg~wt+disp)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")

modell<-function(x){
  with(cryo,{
    model1<-glm(recidiven~x,family="binomial")
    #p<-round(summary(model1)$coefficient[2,4],3)
    s<-summary(model1)
    return(s)
    return(p)
  })
}


esaispasl<-apply(cryo[c("sex","cote")],2,modell)
