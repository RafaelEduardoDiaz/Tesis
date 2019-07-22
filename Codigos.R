##################
#### Packages ####
##################
library(Bayeshmmcts)
library(bridgesampling)
library(rstan)
library(bayesplot)
library(coda)
library(ziphsmm)

#---------- Data homicidios
data("homicides")

#########################################
##### Poisson - Hidden Markov Model #####
#########################################

homicidios <- homicides
colnames(homicidios) <- c("Año","Homicidios","Población","Tasa")
Homicidios <- ts(data = round(homicidios$Tasa), start = 1960)

# modelo clasico
mod2 <- pois.HMM.mle(o = Homicidios, m = 2,
                     lambda0 = c(30, 63), 
                     A0 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE), 
                     stationary = TRUE)

# Algoritmo viterbi (decodificación Global)
viterbi <- pois.HMM.viterbi(o = Homicidios, mod = mod2)

# Verificación de supuestos
residuales <- pois.HMM.pseudo_residuals(o = Homicidios, mod = mod2)
pois.HMM.plot.residuals(residuales)

# Predicción de los estados
`año` <- homicidios$`Año`
Estad_pred <- data.frame(`Año`=`año`[59]+1:16,
                         pois.HMM.state_prediction(h=16,o=Homicidios,mod=mod2))
colnames(Estad_pred) <- c("Año","Estado 1","Estado 2")
Estad_pred$Estado <- apply(Estad_pred, 1, which.max)

# Distribución de prónostico
delta <- pois.HMM.stadist(mod2)
h <- 16
xf <- 5:75
`año` <- homicidios[,1]
forecasts <- pois.HMM.forecast(xf, h, Homicidios, mod2)
par(mfrow = c(4, 4), las = 1)
for(i in 1:h){
  fc <- forecasts[, i]
  plot(xf, fc, type = "h", main = paste("Dist. pronós.", `año`[59] + i),
       xlim = c(5, max(xf+2)), ylim = c(0 ,0.1), cex.main = 0.85,
       xlab = "conteo", ylab = "probabilidad", lwd = 1)
  rect(par("usr")[1], par("usr")[3],par("usr")[2],
       par("usr")[4],col=gray(.9,.9),border='white')
  grid(lty=1, col='white')
  lines(xf, fc, type = "h", lwd = 1)
  dstat <- numeric(length(xf))
  for(j in 1:mod2$m) dstat <- dstat + delta[j] * dpois(xf, mod2$lambda[j])
  lines(xf, dstat, col = "chartreuse2", lwd = 2)
}

#------------------------------------------------------------------------------
# Ajuste Bayesiano PHMM -------------------------------------------------------
#------------------------------------------------------------------------------

# Modelo bayesiano de 2 estados
PHMM_2states <- bayes.PHMM(y = Homicidios, 
                           m = 2, chains = 3, iter = 2000, 
                           control = list(adapt_delta = 0.99))
print(PHMM_2states, digits = 3)

# Modelo bayesiano de 3 estados
PHMM_3states <- bayes.PHMM(y = Homicidios, 
                           m = 3, chains = 3, iter = 2000, 
                           control = list(adapt_delta = 0.99))


# estimates of the log marginal likelihoods
bridge_H0 <- bridge_sampler(samples = PHMM_2states)
bridge_H1 <- bridge_sampler(samples = PHMM_3states)

error_measures(bridge_H0)$percentage
error_measures(bridge_H1)$percentage

#The Bayes factor in favor of H0 over H1 can then be obtained as follows:
bridge_H0$logml
bridge_H1$logml

# Factor de bayes
bf(bridge_H0, bridge_H1)

# Posterior
posterior <- as.array(PHMM_2states)
lp_cp <- log_posterior(PHMM_2states)
np_cp <- nuts_params(PHMM_2states)
rstan::traceplot(PHMM_2states)

# Gráfica de intevalos de credibilidad
color_scheme_set("red")
mcmc_intervals(posterior, prob_outer = 0.95,
               pars = c("A[1,1]", "A[1,2]", "A[2,1]", "A[2,2]")) 

# Histogragmas univariados y gráfico de dispersión bivariado
color_scheme_set("mix-brightblue-gray")
mcmc_pairs(posterior, 
           np = np_cp, 
           pars = c("A[1,1]","A[1,2]","A[2,1]","A[2,2]",
                    "lambda[1]","lambda[2]"), off_diag_args = list(size = 0.75))

# Prueba de Heibelberg y test medio ancho 
PHMM_mcmc <- as.mcmc(as.matrix(PHMM_2states))
Test_HyW <- heidel.diag(PHMM_mcmc)

# Intervalos de Confianza y de credibilidad
intervalos_cred <- mcmc_intervals_data(posterior, prob_outer = 0.95, point_est = "mean")
intervalos_conf <- pois.HMM.confint(mod = mod2, n = 59, B = 250)

#######################################################
##### Zero Inflated Poisson - Hidden Markov Model #####
#######################################################
rm(list = ls())

#---------- Data GIF #
data("wildfires")
incendios <- wildfires
colnames(incendios) <- c("Fecha", "GIF")
GIF <- ts(data = incendios$GIF, start = c(2002,1), frequency = 12)

#### ZIP HMM clásico de 2 estados
ZIPHMM_2states <- hmmfit(y=incendios$GIF,
                         M = 2, prior_init = c(0.6, 0.4),
                         tpm_init = matrix(c(0.9,0.1,
                                             0.5,0.5),2,2,byrow=TRUE),
                         emit_init=c(7, 45), zero_init = c(0.4,0),
                         method="Nelder-Mead",hessian=TRUE,control=list(maxit=1000,trace=1))

#### ZIP HMM clásico de 4 estados
ZIPHMM_4states <- hmmfit(y=incendios$GIF,
                         M = 4, prior_init = c(0.5,0.2,0.2,0.1),
                         tpm_init = matrix(c(0.80,0.15,0.04,0.01,
                                             0.50,0.30,0.15,0.05,
                                             0.15,0.35,0.45,0.05,
                                             0.15,0.35,0.25,0.25),4,4,byrow=TRUE),
                         emit_init=c(3, 15, 43, 100), zero_init = c(0.45,0,0,0),
                         method="Nelder-Mead",hessian=TRUE,control=list(maxit=1000,trace=1))

# Algoritmo Viterbi para el ZIP HMM
ZIP.viterbi <- hmmviterbi(y = incendios$GIF, 
                          ntimes = length(incendios$GIF), 
                          M = 4, prior_init = ZIPHMM_4states$prior, 
                          tpm_init = ZIPHMM_4states$tpm, 
                          emit_init = ZIPHMM_4states$emit_parm, 
                          zero_init = ZIPHMM_4states$zeroprop)

# Gráfica decodificación algoritmo Viterbi
par(mfrow = c(1,1))
par(mar=c(2,2,1,.5)+.5, mgp=c(1.6,.6,0))
### Plot 1
plot(GIF,xlab="Año",type='o', col=4, ylab = "Número")
rect(par("usr")[1], par("usr")[3],par("usr")[2],
     par("usr")[4],col=gray(.9,.9),border='white');
grid(lty=1, col='white')
lines(GIF, type='o', col=4)
abline(h = ZIPHMM_4states$emit_parm, col = "orange", lty = 2)
points(x = time(GIF), 
       y = ifelse(ZIP.viterbi == 1, ZIPHMM_4states$emit_parm[1],
                  ifelse(ZIP.viterbi == 2, ZIPHMM_4states$emit_parm[2],
                         ifelse(ZIP.viterbi == 3,ZIPHMM_4states$emit_parm[3],
                                ZIPHMM_4states$emit_parm[4]))), 
       pch = 21, bg = "orange", col = "white")

#------------------------------------------------------------------------------
# Ajuste Bayesiano ZIP HMM ----------------------------------------------------
#------------------------------------------------------------------------------

# Bayes ZIP HMM 2 y 4 Estaddos
Bayes_ZIPHMM1_2S <- bayes.ZIPHMM1(y=GIF,m=2,chains=4,iter=2000)
Bayes_ZIPHMM1_4S <- bayes.ZIPHMM1(y=GIF,m=4,chains=4,iter=2000)

# Factor de bayes, extrayendo la log verosimilitud marginal, 
# utilizando muestreador por puente
set.seed(1)
bridge_H0 <- bridge_sampler(samples = Bayes_ZIPHMM1_2S)
bridge_H2 <- bridge_sampler(samples = Bayes_ZIPHMM1_4S)
bf(bridge_H0, bridge_H2) # Evidencia extrema para H2

# Gráfica de las cadenas
posteriorZ <- as.array(Bayes_ZIPHMM1_4S)
lp_cpZ <- log_posterior(Bayes_ZIPHMM1_4S)
np_cpZ <- nuts_params(Bayes_ZIPHMM1_4S)
rstan::traceplot(Bayes_ZIPHMM1_4S, pars = c("lambda[1]","lambda[2]","lambda[3]","lambda[4]","A[1,1]","A[1,2]","A[1,3]","A[1,4]","A[2,1]","A[2,2]","A[2,3]","A[2,4]","A[3,1]","A[3,2]","A[3,3]","A[3,4]","A[4,1]","A[4,2]","A[4,3]","A[4,4]","theta"), ncol = 4)

# Prueba de convergencia Heidelberg y Welch
ZIPHMM_mcmc <- as.mcmc(as.matrix(Bayes_ZIPHMM1_4S))
Test_HyW <- heidel.diag(ZIPHMM_mcmc)

# Intervalos de Credibilidad y de Confianza
intervalos_credZ <- as.data.frame(mcmc_intervals_data(posteriorZ, prob_outer = 0.95, point_est = "mean"))
intervalos_credZ$Ancho <- intervalos_credZ$hh - intervalos_credZ$ll
intervalos_confZIP <- ZIP.HMM.confint(mod=ZIPHMM_4states,n=length(GIF),B=100)
