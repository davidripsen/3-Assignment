sde2_multi_room <- function(data){
  data$yT1 <- data$yTi1
  data$yT2 <- data$yTi2
  data$yT3 <- data$yTi3
  data$yT4 <- data$yTi4
  data$Ph <- (data$Ph1 + data$Ph2) / 2
  bs1 = data$bs1; bs2 = data$bs2; bs3 = data$bs3; bs4 = data$bs4; bs5 = data$bs5;
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dT1 ~  1/Ci*(1/Ria*(Ta-T1) + 1/Rim*(T1m-T1) + Ph + (a1*bs1 + a2*bs2 + a3*bs3 + a4*bs4 + a5*bs5)*Gv)*dt + exp(p11a)*dw11)
  model$addSystem(dT1m ~  1/Cm*(1/Rim*(T1-T1m) + 1/Riw*(T2-T1))*dt + exp(p22a)*dw12)
  
  model$addSystem(dT2 ~  1/Ci*(1/Rim*(T2m-T2) + Ph)*dt + exp(p11b)*dw21)
  model$addSystem(dT2m ~  1/Cm*(1/Rim*(T2-T2m) + 1/Riw*(T1-T2) + 1/Riw*(T3-T2))*dt + exp(p22b)*dw22)
  
  model$addSystem(dT3 ~  1/Ci*(1/Rim*(T3m-T3) + Ph)*dt + exp(p11c)*dw31)
  model$addSystem(dT3m ~  1/Cm*(1/Rim*(T3-T3m) + 1/Riw*(T2-T3) + 1/Riw*(T4-T3))*dt + exp(p22c)*dw32)
  
  model$addSystem(dT4 ~  1/Ci*(1/Ria*(Ta-T4) + 1/Rim*(T4m-T4)  + Ph + (a1*bs1 + a2*bs2 + a3*bs3 + a4*bs4 + a5*bs5)*Gv)*dt + exp(p11d)*dw41)
  model$addSystem(dT4m ~  1/Cm*(1/Rim*(T4-T4m) + 1/Riw*(T3-T4))*dt + exp(p22d)*dw42)
  
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph, bs1, bs2, bs3, bs4, bs5)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yT1 ~ T1)
  model$addObs(yT2 ~ T2)
  model$addObs(yT3 ~ T3)
  model$addObs(yT4 ~ T4)
  # Set the variance of the measurement error
  model$setVariance(yT1 ~ exp(e11))
  model$setVariance(yT2 ~ exp(e21))
  model$setVariance(yT3 ~ exp(e31))
  model$setVariance(yT4 ~ exp(e41))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(T1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T1m = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T2m = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T3m = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T4 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T4m = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Riw = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(p11a = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22a = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p11b = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22b = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p11c = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22c = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p11d = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22d = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = 1, lb = -50, ub = 10))
  model$setParameter(e21 = c(init = 1, lb = -50, ub = 10))
  model$setParameter(e31 = c(init = 1, lb = -50, ub = 10))
  model$setParameter(e41 = c(init = 1, lb = -50, ub = 10))
  model$setParameter(a1 = c(init = 4, lb = -500, ub = 1000))
  model$setParameter(a2 = c(init = 4, lb = -500, ub = 1000))
  model$setParameter(a3 = c(init = 4, lb = -500, ub = 1000))
  model$setParameter(a4 = c(init = 4, lb = -500, ub = 1000))
  model$setParameter(a5 = c(init = 4, lb = -500, ub = 1000))
  ##----------------------------------------------------------------    
  # Optimization criteria
  model$options$maxNumberOfEval = 1000
  model$options$eps = 1E-5
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}