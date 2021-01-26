########################################################################################
### RoystonParmarBayes function
### contact: harlancampbell@gmail.com

# Fresh start 
rm(list = ls(all = TRUE))

#####################
## Determine missing packages and load them:
required_packages <- c("flexsurv", "rjags")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed, type="source")                                           
invisible(suppressWarnings(lapply(required_packages, require, character.only = TRUE)))
rm("not_installed", "required_packages")

#####################
## RoystonParmarBayes function

RoystonParmarBayes <- function(times, 
                       event, 
                       treat_id, 
                       knots_placements = c(0,0.5,1), 
                       NITER = 4000,
                       NTHIN = 10,
                       makeplot = TRUE){
                       	
  treat_id <- as.numeric(as.factor(treat_id))
  n_knots <- length(knots_placements)
  
  # Now we compute a basis for a natural cubic spline, 
  # using the parameterisation described by Royston and Parmar (2002).
  knots <- rbind(apply(cbind(knots_placements), 1, 
                         function(x) {quantile(log(times), x)} ))
    
    knots.rep <- knots[rep(seq_len(nrow(knots)), each = length(times)), ] 
    b <- basis(knots = knots.rep, log(times)) 
    db <- dbasis(knots = knots.rep, log(times)) 

 trt_mat <- model.matrix(lm(rep(1,length(treat_id))~as.factor(treat_id)))[,-1]
 the_trt <- matrix(c(trt_mat), , 1)

  ### START OF JAGS CODE 
  # This code is based on Bayesian one-step IPD network meta-analysis of 
  # time-to-event data using Royston-Parmar models: Supplementary material
  # by S. C. Freeman, J. R. Carpenter

model_splines <- "model {
  	
  v_eta[1:n, 1] <- gamma[1, 1] + 
  					(gamma[2:(n_knots), 1])%*%t(u[1:n, 1:(n_knots-1)]) + 
  						(beta[1, 1])%*%t(trt[1:n,]) +  						 
  							(((beta[1, 2])%*%t(trt[1:n, ]))*lnt[1:n])	
  	
  d.sp[1:n, 1] <- (gamma[2:(n_knots), 1])%*%t(du[1:n,1:(n_knots-1)]) +   
						((beta[1, 2])%*%t(trt[1:n, ]))


  # Likelihood  Freeman et al. (2017) equantion (4)  	
  for(i in 1:n) {
  	
  lnL[i] <- - (
 		 # for an observed event:
			  d[i] * ( log( max(d.sp[i, 1], 0.0000001) ) +  
		  		v_eta[i, 1] - exp(v_eta[i, 1]) ) - 
		 # for a censored observation:
				  (1 - d[i]) * exp(v_eta[i, 1]) )
  }

  # Use zeros trick to maximise likelihood
  C <- 100
  zeros[1] ~ dpois(sum(lnL[1:n]) + C)

  # Prior Distributions
  for(p in 1:(n_knots)) {gamma[p, 1] ~ dnorm(0, 0.0001)}
  beta[1, 1] ~ dnorm(0, 0.0001)
  beta[1, 2] ~ dnorm(0, 0.0001)  
}"
  ###  END OF JAGS CODE ### 
   cat(model_splines, file = "model_splines.txt")
  

  jags.model <- jags.model(file = "model_splines.txt", 
                           data = list(n = length(times),
                                       n_knots = n_knots,
                                       lnt = log(times),
                                       d = event,
                                       trt = the_trt,
                                       u = b[,-1],
                                       du = db[,-c(1)],
                                       zeros = 0
                           ),
                           n.chains = 3,
                           n.adapt = 3000)
                           
params <- c("gamma", "beta")
samps <- coda.samples(jags.model, params, n.iter = NITER,  thin = NTHIN)
summary_samps <- summary(samps)$statistics     

# for ploting fitted surival curves:
if(makeplot){	
	gamma_est <- summary_samps[rownames(summary_samps)[grep("gamma", rownames(summary_samps))],"Mean"]
	beta_est <- summary_samps[rownames(summary_samps)[grep("beta", rownames(summary_samps))],"Mean"]

	plot(survfit(Surv(times, event) ~ paste(treat_id, sep="_")), main= "Survival")

	for(trt in 0:1){
		times_seq <- seq(0,max(times),0.1); 
		treat_id_seq <- rep(trt, length(times_seq))
		knots.rep <- knots[rep(seq_len(nrow(knots)), each = length(times_seq)), ] 
		b_seq <- basis(knots = knots.rep, log(times_seq)) 

		v_eta_est <- (gamma_est[1:(n_knots)])%*%t(b_seq[,]) + 
			(beta_est[1])%*%t(treat_id_seq)  +   
				((beta_est[2]%*%t(treat_id_seq) )*log(times_seq))
	
		lines(exp(-exp(c(v_eta_est))) ~ times_seq, col = treat_id_seq+1, lwd = 5)
		}	
	}

return(summary_samps)}


#####################
## Example:

n <- 100
treat_id <- c(rep(0, n/2), rep(1, n/2))
times <- rweibull(n , 9 + -7*treat_id, 5 + 1*treat_id)
event <- rep(1, n)

RPB <- RoystonParmarBayes(times, 
                       event, 
                       treat_id, 
                       knots_placements = c(0,0.5,1), 
                       NITER = 2000,
                       NTHIN = 5)

# we can compare this result to frequentist fit:                    
splw <- flexsurvspline(Surv(times, event) ~ as.factor(treat_id) + 
				gamma1(as.factor(treat_id)), k=1, scale="hazard")				
lines(splw, col="blue")                           


cbind(coefficients(splw) [c(4,5,1,2,3)], RPB[,"Mean"])
