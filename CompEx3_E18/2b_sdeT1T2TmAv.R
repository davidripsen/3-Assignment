2b_sdeT1T2TmAv <- function(data, yT1,yT2,Ph){
  data$yT1 <- yT1
  data$Ph <- Ph
  data$T2 <- yT2
  bs1 = data$bs1; bs2 = data$bs2; bs3 = data$bs3; bs4 = data$bs4; bs5 = data$bs5; 
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in T1: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dT1 ~  1/Ci*(1/Ria*(Ta-T1) + 1/Rim*(Tm-T1) + Ph +
                                 (a1*bs1 + a2*bs2 + a3*bs3 + a4*bs4 + a5*bs5)*Gv)*dt
                  + exp(p11)*dw1)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(T1-Tm) + 1/R_21*(T2-T1))*dt + exp(p22)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,Gv,Ph,bs1,bs2,bs3,bs4,bs4)
  # Set the observation equation: T1 is the state, yT1 is the measured output
  model$addObs(yT1 ~ T1)
  # Set the variance of the measurement error
  model$setVariance(yT1 ~ exp(e11))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(T1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(R_21 = c(init = 20, lb = 1E-4, ub = 1E5)) # New
  #model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5))
  model$setParameter(p11 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(p22 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e11 = c(init = -1, lb = -50, ub = 10))
  model$setParameter(a1 = c(init = 0.2, lb = 0, ub = 1))
  model$setParameter(a2 = c(init = 0.2, lb = 0, ub = 1))
  model$setParameter(a3 = c(init = 0.2, lb = 0, ub = 1))
  model$setParameter(a4 = c(init = 0.2, lb = 0, ub = 1))
  model$setParameter(a5 = c(init = 0.2, lb = 0, ub = 1))
  ##----------------------------------------------------------------    
  
  # Run the parameter optimization
  
  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}