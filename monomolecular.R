mono_fun <- function(t, y, par) {
  
  y <- y[1]
  r <- par$r
  dy = y*r
  return(list(c(dy)))
}

logistic = function(N=10,dt=1, y0=0.01, r, sd=1, inf = 1){
  
  time <- seq(0,N, by=dt)
  w = numeric(length(time))
  y <- numeric(length(time))
  y[1] = y0
  aa <- -1
  bb <- 1
  for (k in 1:(length(time) - 1)) {
    #contant
    if(inf == 1){r[k+1] = r[k]}
    
    #Increassing
    if(inf == 2){r[k+1] = r[k]+0.002}
    
    #Decreasing
    if(inf == 3){r[k+1] = r[k]-0.002}
    
    #Sinusoidal
    if(inf == 4){r[k+1] = r[k]+ ((pi * cos((pi * (k - 1)) / 30)) / 360)}
    
    #Random
    rr <- aa + (bb - aa) * runif(1)
    if(inf == 5){r[k+1] = r[k] + 0.05 * rr}
    
    InitCond <- c(y[k])
    steps = seq(time[k],time[k+1], by=dt)
    parms <- list(r=r[k])
    ode_mono <- ode(InitCond, steps, logi_fun, parms)  
    y[k + 1] = last(ode_mono[,2])
    #  y[k + 1] <- dt * (r[k] * y[k] * (1 - y[k])) + y[k]
  }
  
  
  for(i in 1:length(time)){
    
    w[i] = rnorm(1,y[i],sd = sd*y[i]*(1-y[i]))
    
    if(w[i] > 1){
      w[i]=1
    }
    if(w [i] < 0){
      w[i]=0
    }
  }
  return(data.frame(time, Intensity = y, Randon_intensity = w, inf_rate = r))
}