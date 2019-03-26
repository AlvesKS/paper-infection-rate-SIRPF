


SIR_filter = function(model, Nparti, measures, time,  guess_r, sd_meas, sd_par, sd_model, dt= 0.5){
  realSmeas <- measures

  # Initial guess for the estimation. Here we set as the first value of the synthetic measures vector
  N = length(realSmeas)
  y <- numeric(N)
  r <- numeric(N)
  dt = dt
  Xinit <- realSmeas[1]
  S <- Xinit
  r[1] = guess_r
  time = time
  # particle number
  
  Nparti <- Nparti
  Nparti1 <- (1 / Nparti)
  
  # Variables' and model's error
  
  
  Smeas <- mean(realSmeas)
  stdmeas <- 0.25 * Smeas*(1-Smeas)
  # Synthetic data error
  stdmodel <- 0.005 * S
  # Model error
  stdmodeldif <- sd_par 
  # "parameter"" error
  
  
  # Creating the weights' vector
  
  wparti <- numeric(Nparti)
  
  
  # Creating the others vectors for the SIR-PF algorithm
  
  sinti <- numeric(N)
  sinti[1] <- r[1] 
  Snew <- numeric(N)
  Snew[1] <- S
  xestsir <- numeric(N)
  xestsir[1] <- S
  loop1 <- c(1:(length(time) - 1))
  loop2 <- c(1:Nparti)
  
  xpartires <- numeric(Nparti)
  sintires <- numeric(Nparti)
  wpartires <- numeric(Nparti)
  xpartinew <- numeric(Nparti)
  sintinew <- numeric(Nparti)
  stdsint <- numeric(N)
  stdsir <- numeric(N)
  
  aa <- -1
  bb <- 1
  for (k in loop1) { # loop to every time point
    
    for (l in loop2) { # loop to create the particles
      
      # Rondomic particle for "intensity" :
      
      Sold <- Snew[k] + rnorm(1)*stdmodel
      
      # Rondomic particle for "parameter":
      
      rr <- aa + (bb - aa) * runif(1)
      sintold <- sinti[k] + rr * stdmodeldif 
      
      # Solving the direct problem for every particle
      InitCond <- c(Sold)
      steps <- seq(time[k], time[k + 1], by = dt)
      parms <- list(r = sintold)
        
      
         # Logistic
         if(model == 1){
           ode_logi <- ode(InitCond, steps, logi_fun, parms)  
           y[k + 1] = last(ode_logi[,2])
           }     
         # gompertz
         if(model == 2) {
           ode_gompi =ode(InitCond, steps, gompi_fun, parms)  
           y[k + 1] = last(ode_gompi[,2]) 
         }
       

       xpartinew[l] <- y[k + 1]
      
       sintinew[l] <- sintold
      
      # Calculating the weigths
      
      wparti[l] <- exp(-((xpartinew[l] - realSmeas[k+1]) / stdmeas)^2)
    }
    
    # Setting the weigths between 0 and 1
    
    wtotal <- sum(wparti)
    
    wpartin <- numeric(Nparti)
    wpartin <- wparti / wtotal
    
    
    # Resampling
    
    cresa <- numeric(Nparti)
    uresa <- numeric(Nparti)
    cresa[1] <- wpartin[1]
    
    for (i in 2:Nparti) {
      cresa[i] <- cresa[i - 1] + wpartin[i]
    }
    iresa <- 1
    uresa[1] <- runif(1) * Nparti1
    
    for (j in 1:Nparti) {
      uresa[j] <- uresa[1] + Nparti1 * (j - 1)
      
      while (uresa[j] > cresa[iresa]) {
        iresa <- iresa + 1
      }
      
      xpartires[j] <- xpartinew[iresa]
      sintires[j] <- sintinew[iresa]
      wpartires[j] <- Nparti1
    }
    
    Snew[k + 1] <- mean(xpartires)
    sinti[k + 1] <- mean(sintires)
    
    stdsint[k + 1] <- sd(sintires)
    xestsir[k + 1] <- mean(xpartires)
    
    stdsir[k + 1] <- sd(xpartires)
    xpartiold <- xpartires
    
    # Error atualization
    
    stdmodeldif <- sd_par 
    stdmodel <- sd_model * Snew[k + 1]
    stdmeas <- sd_meas * realSmeas[k + 1]*(1-realSmeas[k + 1])
  }
  
  lbdsiro <- sinti - 2.576 *  stdsint
  ubdsiro <- sinti + 2.576 * stdsint
  lbdsir <- xestsir - 2.576 * stdsir
  ubdsir <- xestsir + 2.576 * stdsir
  
  final <- data.frame(time, realSmeas, xestsir, lbdsir, ubdsir, sinti, lbdsiro, ubdsiro)
  
  return(final)
  }
