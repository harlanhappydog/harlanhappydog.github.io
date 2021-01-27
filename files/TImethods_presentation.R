library(rriskDistributions)
library(rjags)


# 50 confirmed cases out of 3330 tests
c(round((binom.test(50, 3330))$conf.int,3))

# adjusted for spec/sens:
(c(0.011, 0.020) + 0.995 - 1)/(0.844 + 0.995 - 1)



########

set.seed(1234)
library(rjags)
model <- "model {

ir   ~ dunif(0,1);
sens ~ dunif(0,1);
spec ~ dunif(0,1);

y_sample ~ dbin(positivity , n_sample);
y_sens ~ dbin(sens , n_sens);
y_spec ~ dbin(spec , n_spec);

positivity = ir * sens + (1 - ir) * (1 - spec);
}
"

cat(model, file = "model1.txt")


jags.model <- jags.model(file = "model1.txt", 
	data = list(
    y_sample = 50,
    n_sample = 3330,
    n_spec   = 401,
    y_spec   = 399,
    n_sens   = 122,
    y_sens   = 103
    ))

params <- c("ir")
samps <- coda.samples(jags.model, params, n.iter = 500000,  thin = 10)
round(summary(samps)$quantile, 3)


#######
IR_CI <- summary(samps)$quantile[c(1,5)]
ab_param <- get.beta.par(p=c(0.025,0.975),q=IR_CI, plot=FALSE)
c(round(ab_param)[1],round(ab_param)[1]+round(ab_param)[2]+1)


###########

modelIFR <- "model {

cc ~ dbin(ir, tests );
cases ~ dbin(ir, pop);
deaths ~ dbin(ifr, cases);

ifr ~ dunif(0,1);
ir ~ dunif(0,1)
}
"
cat(modelIFR, file = "modelIFR.txt")



###########

# k = 1: Gangelt, Germany

ab_param <- get.beta.par(p=c(0.025,0.975), q=c(0.1231, 0.2440), plot=FALSE)
c(round(ab_param)[1],round(ab_param)[1]+round(ab_param)[2]+1)

jags.model <- jags.model(file = "modelIFR.txt", 
	data = list(
    cc = 27,
    tests = 153,
    pop = 12597,
    deaths = 7
    ))

params <- c("ifr")
samps <- coda.samples(jags.model, params, n.iter = 5000000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)
100*round(summary(samps)$quantile[c(1,5)], 4)


###########

# k = 2: Geneva, Switzeland

jags.model <- jags.model(file = "modelIFR.txt", 
	data = list(
    cc = 48,
    tests = 440,
    pop = 506765,
    deaths = 243
    ))

params <- c("ifr")
samps <- coda.samples(jags.model, params, n.iter = 5000000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)
100*round(summary(samps)$quantile[c(1,5)], 4)


###########

# k = 3: Luxembourg

jags.model <- jags.model(file = "modelIFR.txt", 
	data = list(
    cc = 23,
    tests = 1213,
    pop = 615729,
    deaths = 92
    ))

params <- c("ifr")
samps <- coda.samples(jags.model, params, n.iter = 5000000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)
100*round(summary(samps)$quantile[c(1,5)], 4)


###########

# k = 4: Split-Dalmatia, Croatia

jags.model <- jags.model(file = "modelIFR.txt", 
	data = list(
    cc = 12,
    tests = 937,
    pop = 447723,
    deaths = 29
    ))

params <- c("ifr")
samps <- coda.samples(jags.model, params, n.iter = 5000000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)
100*round(summary(samps)$quantile[c(1,5)], 4)


###########

# k = 5: Zurich, Switzerland

jags.model <- jags.model(file = "modelIFR.txt", 
	data = list(
    cc = 13,
    tests = 1166,
    pop = 1520968,
    deaths = 127
    ))

params <- c("ifr")
samps <- coda.samples(jags.model, params, n.iter = 5000000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)
100*round(summary(samps)$quantile[c(1,5)], 4)


###########
	cloglog <- function(x){log(-log(1-x))}
	icloglog <- function(x){1 - exp(-exp(x))}
	
	
metaIFR <- "model {

# Priors:

	icloglog_theta ~ dunif(0, 1); 
	icloglog_beta  ~ dunif(0, 1);
	theta <- log(-log(1-icloglog_theta));
	beta  <- log(-log(1-icloglog_beta));

	inv.var_sig   <- (1/sd_sig)^2 ;
	inv.var_tau   <- (1/sd_tau)^2 ;
	sd_sig     ~ dnorm(0, 1/1) T(0,);
	sd_tau     ~ dnorm(0, 1/0.01) T(0,);

# Likelihood:
	
for(k in 1:K){
	cc[k] ~ dbin(ir[k], tests[k]);
	cases[k] ~ dbin(ir[k], pop[k]);
	deaths[k] ~ dbin(ifr[k], cases[k]);

	cloglog(ifr[k]) <- cloglog_ifr[k];
	cloglog(ir[k]) <- cloglog_ir[k];

	cloglog_ifr[k] ~ dnorm(theta, inv.var_tau);
	cloglog_ir[k] ~ dnorm(beta, inv.var_sig);
  }
}
"
cat(metaIFR, file = "metaIFR.txt")


jags.model <- jags.model(file = "metaIFR.txt", 
	data = list(
    K = 5,
    tests = c(153, 440, 1213, 937, 1166),
    cc = c(27, 48, 23, 12, 13),    
    pop = c(12597, 506765, 615729, 447723, 1520968),
    deaths = c(7, 243, 92, 29, 127)
    ))

params <- c("icloglog_theta")
samps <- coda.samples(jags.model, params, n.iter = 500000,  thin = 100)
100*round(summary(samps)$statistics[c("Mean")], 4)	
100*round((summary(samps)$quantile[c(1,5)]), 4)

	
	
	
	
	
	
	
	
	
	

