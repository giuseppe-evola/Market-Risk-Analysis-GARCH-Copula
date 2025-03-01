rm(list=ls())
library(moments)
library(MASS)
library(lmtest)
library(sandwich)
library(rugarch)
library(dplyr)
library(evmix)       
library(ggplot2)  
library(QRM) 

load("Market_risk.RData")

#-----------
# Compute the log-returns of these stocks, to be used as risk factors.
S_t <- prices[2:nrow(prices), 1:4]
S_t_1 <- prices[1:(nrow(prices)-1), 1:4]

X <- log(S_t / S_t_1)
X <- na.omit(X)

date <- as.Date(rownames(X))

#-----------
#  Conduct a preliminary descriptive analysis, highlighting the features of these
#data, including the possible heteroscedastic behaviour and (non-)normality
#of asset returns

# Descriptive analysis
summary(X)
skewness(X)  # slight positive asimmetry for any stock's return
kurtosis(X)  #Very strong kurtosis for any stock's return

# Visualization of returns in time
rets <- c("WAB", "BEN", "ES", "ABBV")
par(mfrow = c(2, 2))

for (ret in rets) {                                # easy to see the volatility clustering phenomenon
  plot(x = date, y = X[[ret]], 
       type = 'l',
       xlab = 'Date',
       ylab = 'Returns',
       main = paste('Daily log-returns of', ret),
       cex.main = .8,
       cex.lab = .7)
}

par(mfrow = c(1, 1))

# Non - Normality (Theroretical vs Empirical distribution - QQ Plot - Kolmogorov/Smirnov)
results <- data.frame(Stock = character(), Mu = numeric(), Sigma = numeric())

for (ret in rets) {                             # try to fit the distribution to normal
  fit.norm <- fitdistr(X[[ret]], "normal")
  
  mu.hat <- as.numeric(fit.norm$estimate[1])   
  sigma.hat <- as.numeric(fit.norm$estimate[2])
  results <- rbind(results, data.frame(Stock = ret, Mu = mu.hat, Sigma = sigma.hat))
}

par(mfrow = c(2, 2))     # grafical rapresentation of real frequencies vs theoretical frequencies (normal)
for (ret in rets) {
  hist(X[[ret]], probability = TRUE, ylim = c(0, 40), 
       nclass = 40, main = paste('Histogram of', ret, 'daily returns'),
       xlab = 'log-Returns', cex.main = .8, cex.lab = .7)
  
  x <- seq(from = min(X[[ret]]), to = max(X[[ret]]), length.out = 1000)
  lines(x = x, y = dnorm(x, mean = mu.hat, sd = sigma.hat), col = 'red')
}

# QQ plots      
standardized_X <- data.frame(                               # NO NORMALITY        
  WAB = (X$WAB - mean(X$WAB)) / sd(X$WAB),
  BEN = (X$BEN - mean(X$BEN)) / sd(X$BEN),
  ES = (X$ES - mean(X$ES)) / sd(X$ES),
  ABBV = (X$ABBV - mean(X$ABBV)) / sd(X$ABBV)
)
par(mfrow = c(2, 2))
for (i in 1:ncol(standardized_X)) {
  qqnorm(standardized_X[[i]], main = paste("Q-Q plot di", colnames(standardized_X)[i], "standardized"))
  qqline(standardized_X[[i]], col = "red")  
}   # From this graph we ca easily see that the distribution is very different in the tails 


# Kolmogorov Smirnov test
ks_test_WAB <- ks.test(X$WAB, "pnorm", mean = mean(X$WAB), sd = sd(X$WAB))
ks_test_BEN <- ks.test(X$BEN, "pnorm", mean = mean(X$BEN), sd = sd(X$BEN))
ks_test_ES <- ks.test(X$ES, "pnorm", mean = mean(X$ES), sd = sd(X$ES))
ks_test_ABBV <- ks.test(X$ABBV, "pnorm", mean = mean(X$ABBV), sd = sd(X$ABBV))

ks_table_pvalues <- data.frame(ks_p_WAB = ks_test_WAB$p.value,   # p-value <0.05 -> reject the null: NOT NORMAL DISTRIBUTIONS
                               ks_p_BEN = ks_test_BEN$p.value,
                               ks_p_ES = ks_test_ES$p.value,
                               ks_p_ABBV = ks_test_ABBV$p.value)


# Heteroskedasticity:  (eventualmente prova a dfare test arch)
# 1. we can see it from the visualization of returns above

par(mfrow = c(2, 2))

for (ret in rets) {                                # easy to see the volatility clustering phenomenon
  plot(x = date, y = abs(X[[ret]]),                # abs value of the returns
       type = 'l',
       xlab = 'Date',
       ylab = 'Returns',
       main = paste('Daily log-returns (Absolute value)', ret),
       cex.main = .8,
       cex.lab = .7)
}
par(mfrow = c(1, 1))

# Autocorrelation (squared). To enhance the volatilty clustering effect

par(mfrow = c(2, 2))
acf(X$WAB^2, main="ACF^2 of WAB")
acf(X$BEN^2, main="ACF^2 of BEN")
acf(X$ES^2, main="ACF^2 of ES")
acf(X$ABBV^2, main="ACF^2 of ABBV")


#-------------------------------------
#For each stock, model the heteroscedasticity of the log-returns through a
#GARCH model. Consider at least 3 possible specifications of the GARCH
#model and select the best one. Justify your choice

spec_garch11_normal <- ugarchspec(                                       # Garch(1,1), Innovations are normal
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "norm")

spec_garch11_t <- ugarchspec(                                           # Garch(1,1), Innovations are t-student
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "std")

spec_garch12_t <- ugarchspec(                                            # Garch(1,2), Innovations are t student
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "std")

spec_garch12_normal <- ugarchspec(                                      
  variance.model = list(model = "sGARCH", garchOrder = c(1, 2)),        # Garch(1,2) Innovations are normal
  mean.model = list(armaOrder = c(0, 0), include.mean = F),  
  distribution.model = "norm")

spec_garch22_normal <- ugarchspec(                                      
  variance.model = list(model = "sGARCH", garchOrder = c(2, 2)),        # Garch(2,2) Innovations are normal
  mean.model = list(armaOrder = c(0, 0), include.mean = F),  
  distribution.model = "norm")

spec_garch22_t <- ugarchspec(                                            # Garch(2,2), Innovations are t student
  variance.model = list(model = "sGARCH", garchOrder = c(2, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  distribution.model = "std")


# WAB
fit.garch11_normal_WAB <- ugarchfit(spec = spec_garch11_normal, data = X$WAB)
fit.garch11_t_WAB <- ugarchfit(spec = spec_garch11_t, data = X$WAB)
fit.garch12_normal_WAB <- ugarchfit(spec = spec_garch12_normal, data = X$WAB)
fit.garch12_t_WAB <- ugarchfit(spec = spec_garch12_t , data = X$WAB)
fit.garch22_normal_WAB <- ugarchfit(spec = spec_garch22_normal, data = X$WAB)
fit.garch22_t_WAB <- ugarchfit(spec = spec_garch22_t, data = X$WAB)

BIC_WAB <- data.frame(garch_11_normal = infocriteria(fit.garch11_normal_WAB)[2],
                      garch_11_t = infocriteria(fit.garch11_t_WAB)[2],
                      garch_12_normal = infocriteria(fit.garch12_normal_WAB)[2],
                      garch_12_t = infocriteria(fit.garch12_t_WAB)[2],
                      garch_22_normal =  infocriteria(fit.garch22_normal_WAB)[2],
                      garch_22_t = infocriteria(fit.garch22_t_WAB)[2]
)


#BEN
fit.garch11_normal_BEN <- ugarchfit(spec = spec_garch11_normal, data = X$BEN)
fit.garch11_t_BEN <- ugarchfit(spec = spec_garch11_t, data = X$BEN)
fit.garch12_normal_BEN <- ugarchfit(spec = spec_garch12_normal, data = X$BEN)
fit.garch12_t_BEN <- ugarchfit(spec = spec_garch12_t , data = X$BEN)
fit.garch22_normal_BEN <- ugarchfit(spec = spec_garch22_normal, data = X$BEN)
fit.garch22_t_BEN <- ugarchfit(spec = spec_garch22_t, data = X$BEN)

BIC_BEN <- data.frame(garch_11_normal = infocriteria(fit.garch11_normal_BEN)[2],
                      garch_11_t = infocriteria(fit.garch11_t_BEN)[2],
                      garch_12_normal = infocriteria(fit.garch12_normal_BEN)[2],
                      garch_12_t = infocriteria(fit.garch12_t_BEN)[2],
                      garch_22_normal =  infocriteria(fit.garch22_normal_BEN)[2],
                      garch_22_t = infocriteria(fit.garch22_t_BEN)[2])

# ES
fit.garch11_normal_ES <- ugarchfit(spec = spec_garch11_normal, data = X$ES)
fit.garch11_t_ES <- ugarchfit(spec = spec_garch11_t, data = X$ES)
fit.garch12_normal_ES <- ugarchfit(spec = spec_garch12_normal, data = X$ES)
fit.garch12_t_ES <- ugarchfit(spec = spec_garch12_t , data = X$ES)
fit.garch22_normal_ES <- ugarchfit(spec = spec_garch22_normal, data = X$ES)
fit.garch22_t_ES <- ugarchfit(spec = spec_garch22_t, data = X$ES)

BIC_ES <- data.frame(garch_11_normal = infocriteria(fit.garch11_normal_ES)[2],
                     garch_11_t = infocriteria(fit.garch11_t_ES)[2],
                     garch_12_normal = infocriteria(fit.garch12_normal_ES)[2],
                     garch_12_t = infocriteria(fit.garch12_t_ES)[2],
                     garch_22_normal =  infocriteria(fit.garch22_normal_ES)[2],
                     garch_22_t = infocriteria(fit.garch22_t_ES)[2])

#ABBV
fit.garch11_normal_ABBV <- ugarchfit(spec = spec_garch11_normal, data = X$ABBV)
fit.garch11_t_ABBV <- ugarchfit(spec = spec_garch11_t, data = X$ABBV)
fit.garch12_normal_ABBV <- ugarchfit(spec = spec_garch12_normal, data = X$ABBV)
fit.garch12_t_ABBV <- ugarchfit(spec = spec_garch12_t , data = X$ABBV)
fit.garch22_normal_ABBV <- ugarchfit(spec = spec_garch22_normal, data = X$ABBV)
fit.garch22_t_ABBV <- ugarchfit(spec = spec_garch22_t, data = X$ABBV)

BIC_ABBV <- data.frame(garch_11_normal = infocriteria(fit.garch11_normal_ABBV)[2],
                       garch_11_t = infocriteria(fit.garch11_t_ABBV)[2],
                       garch_12_normal = infocriteria(fit.garch12_normal_ABBV)[2],
                       garch_12_t = infocriteria(fit.garch12_t_ABBV)[2],
                       garch_22_normal =  infocriteria(fit.garch22_normal_ABBV)[2],
                       garch_22_t = infocriteria(fit.garch22_t_ABBV)[2])

BIC_table <- rbind(WAB = BIC_WAB, BEB = BIC_BEN, ES = BIC_ES, ABBV = BIC_ABBV)

# Only Based on the BIC we have chosen: 
# WAB: garch(1,1) - innovation t student distributed
# BEN: garch(1,1) - innovation t student distributed
# ES: garch(1,1) - innovation t student distributed
# ABBV: garch(1,1) - innovation t student distributed

#-------------------------------------------------------
# Use three different methods to compute the (in-sample) daily Value-at
#Risk (VaR) and Expected Shortfall (ES) over the whole sample period
#(i.e. two measures for each day and each stock in the sample). Use a
#95% confidence level for VaR and a 97.5% level for ES.

# Analytical Approach (normal distribution)
VaR_95_AA <- data.frame(matrix(ncol = length(rets), nrow = 1))
VaR_97.5_AA <- data.frame(matrix(ncol = length(rets), nrow = 1))
ES_97.5_AA <- data.frame(matrix(ncol = length(rets), nrow = 1))
colnames(VaR_95_AA) <- rets
colnames(ES_97.5_AA) <- rets
colnames(VaR_97.5_AA) <- rets

for (ret in 1:length(rets)) {
  fit.norm <- fitdistr(X[[rets[ret]]], "normal")
  
  mu.hat <- as.numeric(fit.norm$estimate[1])   # estimated mean
  sigma.hat <- as.numeric(fit.norm$estimate[2]) # estimated sigma
  
  VaR_95_AA[1, ret] <- -(mu.hat + sigma.hat * qnorm(.05))*100  # VaR 95
  VaR_97.5_AA[1,ret] <- -(mu.hat + sigma.hat * qnorm(.025))*100 # VaR 97.5
  
  #ES at 97.5% 
  z_alpha <- qnorm(.025)   
  phi_z_alpha <- dnorm(z_alpha)  
  ES_97.5_AA[1, ret] <- (mu.hat + (sigma.hat * phi_z_alpha) / (0.025))*100
}

# Historical Approach
VaR_95_HS <- data.frame(matrix(ncol = length(rets), nrow = 1))
VaR_97.5_HS <- data.frame(matrix(ncol = length(rets), nrow = 1))
ES_97.5_HS <- data.frame(matrix(ncol = length(rets), nrow = 1))
colnames(VaR_95_HS) <- rets
colnames(ES_97.5_HS) <- rets
colnames(VaR_97.5_HS) <- rets

for(ret in 1:length(rets)){
  VaR_95_HS[ret] <- -quantile(X[[rets[ret]]], 0.05)*100
  VaR_97.5_HS[ret] <- -quantile(X[[rets[ret]]], 0.025)*100
  perc_2.5 <- quantile(X[[rets[ret]]], 0.025)
  ES_97.5_HS[ret] <- -mean(X[[rets[ret]]][X[[rets[ret]]] <= perc_2.5])*100
}

# FHS approach (garch(1,1), t student innovations)

# WAB
vol_t_WAB <- sigma(fit.garch11_t_WAB)  # conditional volatility
res_WAB <- residuals(fit.garch11_t_WAB, standardize = T) # standard residuals
vol_t_WAB <- coredata(vol_t_WAB) # trasforming and using only the vector of data
res_WAB <- coredata(res_WAB)

n_days <- length(vol_t_WAB)

simulated_returns_WAB <- matrix(NA, nrow = n_days, ncol = n_days) # each column is a day of simualated returns
VaR95_WAB_FHS <- numeric(n_days)  # VaR simulated
ES97.5_WAB_FHS <- numeric(n_days)   # ES simulated
VaR97.5_WAB_FHS <- numeric(n_days) # I am going to use it later on

for (i in 1:n_days) {
  simulated_returns_WAB[, i] <- vol_t_WAB[i] * res_WAB  # product between the volatility at each day times the residuals for all the serie

  VaR95_WAB_FHS[i] <- -quantile(simulated_returns_WAB[, i], 0.05) * 100
  VaR97.5_WAB_FHS[i] <- -quantile(simulated_returns_WAB[, i], 0.025) * 100
  ES97.5_WAB_FHS[i] <- -mean(simulated_returns_WAB[, i][simulated_returns_WAB[, i] <= quantile(simulated_returns_WAB[, i], 0.025)]) * 100
}

# BEN

vol_t_BEN <- sigma(fit.garch11_t_BEN)  # conditional volatility
res_BEN <- residuals(fit.garch11_t_BEN, standardize = T) # standard residuals
vol_t_BEN <- coredata(vol_t_BEN) # trasforming and using only the vector of data
res_BEN <- coredata(res_BEN)

simulated_returns_BEN <- matrix(NA, nrow = n_days, ncol = n_days) # each column is a day of simualated returns
VaR95_BEN_FHS <- numeric(n_days)  # VaR simulated
ES97.5_BEN_FHS <- numeric(n_days)   # ES simulated
VaR97.5_BEN_FHS <- numeric(n_days) #  I am going to use it later on

for (i in 1:n_days) {
  simulated_returns_BEN[, i] <- vol_t_BEN[i] * res_BEN  # poduct between the volatility at each day times the residuals for all the serie
  
  VaR95_BEN_FHS[i] <- -quantile(simulated_returns_BEN[, i], 0.05) * 100
  VaR97.5_BEN_FHS[i] <- -quantile(simulated_returns_BEN[, i], 0.025) * 100
  ES97.5_BEN_FHS[i] <- -mean(simulated_returns_BEN[, i][simulated_returns_BEN[, i] <= quantile(simulated_returns_BEN[, i], 0.025)]) * 100
}

# ES
vol_t_ES <- sigma(fit.garch11_t_ES)  # conditional volatility
res_ES <- residuals(fit.garch11_t_ES, standardize = T) # standard residuals
vol_t_ES <- coredata(vol_t_ES) # trasforming and using only the vector of data
res_ES <- coredata(res_ES)

simulated_returns_ES <- matrix(NA, nrow = n_days, ncol = n_days) # each column is a day of simualated returns
VaR95_ES_FHS <- numeric(n_days)  # VaR simulated
ES97.5_ES_FHS <- numeric(n_days)   # ES simulated
VaR97.5_ES_FHS <- numeric(n_days)  # VaR simulated

for (i in 1:n_days) {
  simulated_returns_ES[, i] <- vol_t_ES[i] * res_ES  # product between the volatility at each day times the residuals for all the serie
  
  VaR95_ES_FHS[i] <- -quantile(simulated_returns_ES[, i], 0.05) * 100
  VaR97.5_ES_FHS[i] <- -quantile(simulated_returns_ES[, i], 0.025) * 100
  ES97.5_ES_FHS[i] <- -mean(simulated_returns_ES[, i][simulated_returns_ES[, i] <= quantile(simulated_returns_ES[, i], 0.025)]) * 100
}

# ABBV
vol_t_ABBV <- sigma(fit.garch11_t_ABBV)  # conditional volatility
res_ABBV <- residuals(fit.garch11_t_ABBV, standardize = T) # standard residuals
vol_t_ABBV <- coredata(vol_t_ABBV) # trasforming and using only the vector of data
res_ABBV <- coredata(res_ABBV)

simulated_returns_ABBV <- matrix(NA, nrow = n_days, ncol = n_days) # each column is a day of simualated returns
VaR95_ABBV_FHS <- numeric(n_days)  # VaR simulated
ES97.5_ABBV_FHS <- numeric(n_days)   # ES simulated
VaR97.5_ABBV_FHS <- numeric(n_days)

for (i in 1:n_days) {
  simulated_returns_ABBV[, i] <- vol_t_ABBV[i] * res_ABBV  # product between the volatility at each day times the residuals for all the serie
  
  VaR95_ABBV_FHS[i] <- -quantile(simulated_returns_ABBV[, i], 0.05) * 100
  VaR97.5_ABBV_FHS[i] <- -quantile(simulated_returns_ABBV[, i], 0.025) * 100
  ES97.5_ABBV_FHS[i] <- -mean(simulated_returns_ABBV[, i][simulated_returns_ABBV[, i] <= quantile(simulated_returns_ABBV[, i], 0.025)]) * 100
}

FHS_risk_measures <- data.frame(date = date, 
                                cbind(VaR95_WAB_FHS, VaR95_BEN_FHS, VaR95_ES_FHS, VaR95_ABBV_FHS,
                                      VaR97.5_WAB_FHS, VaR97.5_BEN_FHS, VaR97.5_ES_FHS, VaR97.5_ABBV_FHS,  
                                      ES97.5_WAB_FHS, ES97.5_BEN_FHS, ES97.5_ES_FHS, ES97.5_ABBV_FHS))

# WAB
par(mfrow = c(2, 1))

plot(x = date, y = X$WAB, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of WAB vs VaR95 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .20))            
abline(h = as.numeric(-VaR_95_AA$WAB/100), col = "blue", lwd = 2)    # VaR analitical
abline(h = as.numeric(-VaR_95_HS$WAB/100), col = "red", lwd = 2)     # VaR HS
lines(date, -VaR95_WAB_FHS/100, col = "green", lwd = 2)                  # VaR FHS

plot(x = date, y = X$WAB, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of WAB vs ES97.5 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .20))            
abline(h = as.numeric(-ES_97.5_AA$WAB/100), col = "blue", lwd = 2)    # ES analitical
abline(h = as.numeric(-ES_97.5_HS$WAB/100), col = "red", lwd = 2)     # ES HS
lines(date, -ES97.5_WAB_FHS/100, col = "green", lwd = 2)                   # ES FHS

# BEN
plot(x = date, y = X$BEN, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of BEN vs VaR95 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.15, .15))            
abline(h = as.numeric(-VaR_95_AA$BEN/100), col = "blue", lwd = 2)    # VaR analitical
abline(h = as.numeric(-VaR_95_HS$BEN/100), col = "red", lwd = 2)     # VaR HS
lines(date, -VaR95_BEN_FHS/100, col = "green", lwd = 2)                  # VaR FHS

plot(x = date, y = X$BEN, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of BEN vs ES97.5 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .15))            
abline(h = as.numeric(-ES_97.5_AA$BEN/100), col = "blue", lwd = 2)    # ES analitical
abline(h = as.numeric(-ES_97.5_HS$BEN/100), col = "red", lwd = 2)     # ES HS
lines(date, -ES97.5_BEN_FHS/100, col = "green", lwd = 2)                   # ES FHS

#ES
plot(x = date, y = X$ES, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of ES vs VaR95 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.15, .15))            
abline(h = as.numeric(-VaR_95_AA$ES/100), col = "blue", lwd = 2)    # VaR analitical
abline(h = as.numeric(-VaR_95_HS$ES/100), col = "red", lwd = 2)     # VaR HS
lines(date, -VaR95_ES_FHS/100, col = "green", lwd = 2)                  # VaR FHS

plot(x = date, y = X$ES, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of ES vs ES97.5 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .15))            
abline(h = as.numeric(-ES_97.5_AA$ES/100), col = "blue", lwd = 2)    # ES analitical
abline(h = as.numeric(-ES_97.5_HS$ES/100), col = "red", lwd = 2)     # ES HS
lines(date, -ES97.5_ES_FHS/100, col = "green", lwd = 2)                   # ES FHS

# ABBV
plot(x = date, y = X$ABBV, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of ABBV vs VaR95 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .20))            
abline(h = as.numeric(-VaR_95_AA$ABBV/100), col = "blue", lwd = 2)    # VaR analitical
abline(h = as.numeric(-VaR_95_HS$ABBV/100), col = "red", lwd = 2)     # VaR HS
lines(date, -VaR95_ABBV_FHS/100, col = "green", lwd = 2)                  # VaR FHS

plot(x = date, y = X$ABBV, 
     type = 'p',          
     xlab = 'Date',
     ylab = 'Returns',
     main = 'Daily log-returns of ABBV vs ES97.5 (different approaches)',
     cex.main = 0.8,
     cex.lab = 0.7,
     pch = 20,
     ylim = c(-.20, .20))            
abline(h = as.numeric(-ES_97.5_AA$ABBV/100), col = "blue", lwd = 2)    # ES analitical
abline(h = as.numeric(-ES_97.5_HS$ABBV/100), col = "red", lwd = 2)     # ES HS
lines(date, -ES97.5_ABBV_FHS/100, col = "green", lwd = 2)                   # ES FHS

par(mfrow = c(1, 1))

# ---------
# Conduct a backtesting analysis of the different methods.

# VaR 

# WAB
VaRTest(alpha = 0.05, actual = X$WAB, VaR = rep(-VaR_95_AA$WAB/100, length(X$WAB))) # UC: H0 , CC: H0 - analitical
VaRTest(alpha = 0.05, actual = X$WAB, VaR = rep(-VaR_95_HS$WAB/100, length(X$WAB)))  # UC: H0 , CC: H1 - historical 
VaRTest(alpha = 0.05, actual = X$WAB, VaR = -VaR95_WAB_FHS/100)  # UC: H0 , CC: H0 - FHS

# BEN
VaRTest(alpha = 0.05, actual = X$BEN, VaR = rep(-VaR_95_AA$BEN/100, length(X$BEN))) # UC: H0 , CC: H1
VaRTest(alpha = 0.05, actual = X$BEN, VaR = rep(-VaR_95_HS$BEN/100, length(X$BEN))) # UC: H0 , CC: H0
VaRTest(alpha = 0.05, actual = X$BEN, VaR = -VaR95_BEN_FHS/100) # UC: H0 , CC: H0

# ES
VaRTest(alpha = 0.05, actual = X$ES, VaR = rep(-VaR_95_AA$ES/100, length(X$ES))) # UC: H1 , CC: H1
VaRTest(alpha = 0.05, actual = X$ES, VaR = rep(-VaR_95_HS$ES/100, length(X$ES))) # UC: H0 , CC: H1
VaRTest(alpha = 0.05, actual = X$ES, VaR = -VaR95_ES_FHS/100) # UC: H0 , CC: H0

#ABBV
VaRTest(alpha = 0.05, actual = X$ABBV, VaR = rep(-VaR_95_AA$ABBV/100, length(X$ABBV))) # UC: H1 , CC: H1
VaRTest(alpha = 0.05, actual = X$ABBV, VaR = rep(-VaR_95_HS$ABBV/100, length(X$ABBV))) # UC: H0 , CC: H0
VaRTest(alpha = 0.05, actual = X$ABBV, VaR = -VaR95_ABBV_FHS/100) # UC: H0 , CC: H0

# Exp Shortfall

# WAB
ESTest(alpha = .025, actual= X$WAB, VaR = rep(-VaR_97.5_AA$WAB/100, length(X$WAB)), ES = rep(-ES_97.5_AA$WAB/100, length(X$WAB)), boot = T) #H1
ESTest(alpha = .025, actual= X$WAB, VaR = rep(-VaR_97.5_HS$WAB/100, length(X$WAB)), ES = rep(-ES_97.5_HS$WAB/100, length(X$WAB)), boot = T) #H0
ESTest(alpha = .025, actual= X$WAB, VaR = -VaR97.5_WAB_FHS/100, ES = -ES97.5_WAB_FHS/100, boot = T) #H0

# BEN
ESTest(alpha = .025, actual= X$BEN, VaR = rep(-VaR_97.5_AA$BEN/100, length(X$BEN)), ES = rep(-ES_97.5_AA$BEN/100, length(X$BEN)),boot = T) #H1
ESTest(alpha = .025, actual= X$BEN, VaR = rep(-VaR_97.5_HS$BEN/100, length(X$BEN)), ES = rep(-ES_97.5_HS$BEN/100, length(X$BEN)), boot = T) # H0
ESTest(alpha = .025, actual= X$BEN, VaR = -VaR97.5_BEN_FHS/100, ES = -ES97.5_BEN_FHS/100, boot = T) # H0

# ES
ESTest(alpha = .025, actual= X$ES, VaR = rep(-VaR_97.5_AA$ES/100, length(X$ES)), ES = rep(-ES_97.5_AA$ES/100, length(X$ES)), boot = T) #H1
ESTest(alpha = .025, actual= X$ES, VaR = rep(-VaR_97.5_HS$ES/100, length(X$ES)), ES = rep(-ES_97.5_HS$ES/100, length(X$ES)), boot = T) #H0
ESTest(alpha = .025, actual= X$ES, VaR = -VaR97.5_ES_FHS/100, ES = -ES97.5_ES_FHS/100, boot = T) #H0

# ABBV
ESTest(alpha = .025, actual= X$ABBV, VaR = rep(-VaR_97.5_AA$ABBV/100, length(X$ABBV)), ES = rep(-ES_97.5_AA$ABBV/100, length(X$ABBV)), boot = T) #H1
ESTest(alpha = .025, actual= X$ABBV, VaR = rep(-VaR_97.5_HS$ABBV/100, length(X$ABBV)), ES = rep(-ES_97.5_HS$ABBV/100, length(X$ABBV)), boot = T) # H0
ESTest(alpha = .025, actual= X$ABBV, VaR = -VaR97.5_ABBV_FHS/100, ES = -ES97.5_ABBV_FHS/100, boot = T) # H0

#----------------------------------
# Use a probability integral transformation of the marginal returns, based
# on your best estimation of the marginal distribution

fit.WAB <- fitdistr(x = X$WAB, densfun = "t")
fit.BEN <- fitdistr(x = X$BEN, densfun = "t") 
fit.ES <- fitdistr(x = X$ES, densfun = "t")
fit.ABBV <- fitdistr(x = X$ABBV, densfun = "t")

U_WAB <- pt((X$WAB - fit.WAB$estimate["m"]) / fit.WAB$estimate["s"], df = fit.WAB$estimate["df"])
U_BEN <- pt((X$BEN - fit.BEN$estimate["m"]) / fit.BEN$estimate["s"], df = fit.BEN$estimate["df"])
U_ES <- pt((X$ES - fit.ES$estimate["m"]) / fit.ES$estimate["s"], df = fit.ES$estimate["df"])
U_ABBV <- pt((X$ABBV - fit.ABBV$estimate["m"]) / fit.ABBV$estimate["s"], df = fit.ABBV$estimate["df"])

U_table <- coredata(cbind(U_WAB, U_BEN, U_ES, U_ABBV))

#----------
#Estimate at least two different copula models from the 4-dimensional
#vectors of pseudo-observations. 

fit.normalCopula <- fitCopula(normalCopula(dim = 4),  data = U_table)
fit.tCopula <- fitCopula(tCopula(dim = 4), data = U_table) 
fit.claytonCopula <- fitCopula(claytonCopula(dim = 4), data = U_table)

copula.N <- fit.normalCopula@copula
copula.t <- fit.tCopula@copula
copula.c <- fit.claytonCopula@copula

coef(fit.normalCopula)
coef(fit.tCopula)
coef(fit.claytonCopula)

lambda(copula.N)
lambda(copula.t)
lambda(copula.c)

BIC(fit.normalCopula)
BIC(fit.tCopula)      #easy to see that BIC tCopula is the lowest. We choose this one
BIC(fit.claytonCopula)

pairs(U_table)

#--------------
# Compute the daily portfolio return over the complete sample period. Plot
#this quantity as a function of time.

omega <- rep(1, ncol(X))
X_matrix <- as.matrix(X)

V <- X_matrix %*% omega  # Portfolio returns are the weighted sum of individual asset returns

# Plot the portfolio returns over time
par(mfrow = c(1, 1))

plot(x = date, y = V,
     type = 'l',          # Line plot
     xlab = 'Date',       # Label for the x-axis
     ylab = 'Portfolio returns',  # Label for the y-axis
     main = 'Daily log-returns of the Portfolio',  # Title of the plot
     cex.main = 0.8,      # Size of the title text
     cex.lab = 0.7)       # Size of the axis labels

#----------
#  You now want to derive risk measures for your portfolio but here using an
#in-sample/out-of-sample approach:

#Consider the data before January 2022 to be an initial in-sample period.
#The out-of-sample period is January 2022 to November 2024. From
#now on, assume that you are interested in modelling the out-of-sample
#distribution of the linearized loss operator (LLO) of this portfolio. That
#is, you want to infer the distribution of  L∆
#t+1 = −∆Vt+1 = −b′Xt+1
#using only past risk factors X1,...,Xt. Assume that estimation for t +
#1 is always performed using the previous 3 years of data (i.e. 750
#observations and the first in-sample period is January 2019 - January 2022).

out_sample_start <- which(date == as.Date("2022-01-03"))
out_sample_end <- length(V)

#---------
# Compute the daily VaR of the LLO for this portfolio, at a 99% confidence
#level. Use two methods: on the one hand, a rolling window FHS
#approach based on a GARCH model of your choice for filtering; on the
#other hand, a rolling window Monte Carlo simulation approach combining
#GARCH and an appropriate copula model. Report the proportion
#of violations suffered by each model.

# Rolling window FHS (VaR 99%)

# I decided to fit a garch(1, 1) with t-student innovations because generally it is the best one in terms of BIC

# Montecarlo rolling window approach (Copula fitting) (VaR 99)
# For the Montecarlo approach with copula, I decided to use T-copula for each rolling window because it fit very well to
# the data in its entirety and because generally it the one which is able to capture the tail dependencies

window_size <- 756
VaR99_FHS_RW <- numeric(out_sample_end - out_sample_start + 1)
VaR99_MC <- numeric(out_sample_end - out_sample_start + 1)

for (t in out_sample_start:out_sample_end){
  
  window_data_WAB <- X$WAB[(t-window_size):(t-1)] # in-sample data (rolling windows)
  window_data_BEN <- X$BEN[(t-window_size):(t-1)]
  window_data_ES <- X$ES[(t-window_size):(t-1)]
  window_data_ABBV <- X$ABBV[(t-window_size):(t-1)] 
  
  fit.garch_WAB <- ugarchfit(spec = spec_garch11_t, data = window_data_WAB)
  fit.garch_BEN <- ugarchfit(spec = spec_garch11_t, data = window_data_BEN)
  fit.garch_ES <- ugarchfit(spec = spec_garch11_t, data = window_data_ES)
  fit.garch_ABBV <- ugarchfit(spec = spec_garch11_t, data = window_data_ABBV)
  
  # Rolling window FHS (VaR 99%)
  simulated_rets_WAB <-residuals(fit.garch_WAB, standardize = TRUE)*rep(tail(sigma(fit.garch_WAB), 1),length(window_data_WAB))
  simulated_rets_BEN <-residuals(fit.garch_BEN, standardize = TRUE)*rep(tail(sigma(fit.garch_BEN), 1),length(window_data_BEN))
  simulated_rets_ES <-residuals(fit.garch_ES, standardize = TRUE)*rep(tail(sigma(fit.garch_ES), 1),length(window_data_ES))
  simulated_rets_ABBV <-residuals(fit.garch_ABBV, standardize = TRUE)*rep(tail(sigma(fit.garch_ABBV), 1),length(window_data_ABBV))
  
  simulated_rets <- cbind(simulated_rets_WAB, simulated_rets_BEN, simulated_rets_ES, simulated_rets_ABBV)
  
  V_FHS <- simulated_rets %*% omega
  VaR99_FHS_RW[t - out_sample_start + 1] <- -quantile(V_FHS, .01)
  
  # Montecarlo rolling window approach (Copula fitting) (VaR 99)
  
  U_WAB1 <- pt(q = residuals(fit.garch_WAB, standardize = T), df = fit.garch_WAB@fit$coef["shape"])
  U_BEN1 <- pt(q = residuals(fit.garch_BEN, standardize = T), df = fit.garch_BEN@fit$coef["shape"])
  U_ES1 <- pt(q = residuals(fit.garch_ES, standardize = T), df = fit.garch_ES@fit$coef["shape"])
  U_ABBV1 <- pt(q = residuals(fit.garch_ABBV, standardize = T), df = fit.garch_ABBV@fit$coef["shape"])
  
  U <- coredata(cbind(U_WAB1, U_BEN1, U_ES1, U_ABBV1))
  
  fit.CopulaT <- fitCopula(tCopula(dim = 4), data = U)
  
  B <- 10000
  rU <- rCopula(n = B, copula = fit.CopulaT@copula)
  
  #WAB
  vol.today_WAB <- tail(sigma(fit.garch_WAB), 1)
  MC_WAB <- qt(p = rU[,1], df = fit.garch_WAB@fit$coef["shape"])*rep(vol.today_WAB, B)
  # BEN
  vol.today_BEN <- tail(sigma(fit.garch_BEN), 1)
  MC_BEN <- qt(p = rU[,2], df = fit.garch_BEN@fit$coef["shape"])*rep(vol.today_BEN, B)
  # ES
  vol.today_ES <- tail(sigma(fit.garch_ES), 1)
  MC_ES <- qt(p = rU[,3], df = fit.garch_ES@fit$coef["shape"])*rep(vol.today_ES, B)
  # ABBV
  vol.today_ABBV <- tail(sigma(fit.garch_ABBV), 1)
  MC_ABBV <- qt(p = rU[,4], df = fit.garch_ABBV@fit$coef["shape"])*rep(vol.today_ABBV, B)
  
  MC_sim_rets <- cbind(MC_WAB, MC_BEN, MC_ES, MC_ABBV) # Montecarlo simulated returns
  
  V_MC <- MC_sim_rets%*%omega     # Portfolio returns 
  VaR99_MC[t - out_sample_start + 1] <- -quantile(V_MC, 0.01) # VaR99
}

# Violations
violations_FHS <- sum(V[out_sample_start:out_sample_end] < -VaR99_FHS_RW) / length(V[out_sample_start:out_sample_end])
violations_MC <- sum(V[out_sample_start:out_sample_end] < -VaR99_MC) / length(V[out_sample_start:out_sample_end])

#------
#Plot the daily VaRs of the portfolio, at a ten-day horizon, using the √h
#rule, where h is the VaR horizon.

scaling_factor <- sqrt(10)
VaR99_FHS_RW_10day <- VaR99_FHS_RW * scaling_factor
VaR99_MC_10day <- VaR99_MC * scaling_factor

out_sample_dates <- date[out_sample_start:out_sample_end]

par(mfrow = c(1, 1))

plot(out_sample_dates, VaR99_FHS_RW_10day, type = "l", col = "blue", 
     ylab = "VaR (10-day)", xlab = "Date", main = "10-Day Horizon VaR")
lines(out_sample_dates, VaR99_MC_10day, col = "red")


#------------
#Starting at day 61 of your forecasts, compute the trading book regulatory
# capital for both methods, assuming k = 3.5 and c = 0 (see formula
#in the course). 
k <- 3.5
scaling_factor <- sqrt(10)
VaR99_FHS_RW_10day <- VaR99_FHS_RW * scaling_factor
VaR99_MC_10day <- VaR99_MC * scaling_factor

# ------ FHS ------
n_FHS <- length(VaR99_FHS_RW_10day)

VaR_moving_avg_FHS <- rep(NA, n_FHS)
for (t in 61:n_FHS) { 
  VaR_moving_avg_FHS[t] <- mean(VaR99_FHS_RW_10day[(t-60):(t-1)])
}

RC_FHS <- rep(NA, n_FHS)
for (t in 61:n_FHS) {  
  RC_FHS[t] <- max(VaR99_FHS_RW_10day[t], k * VaR_moving_avg_FHS[t])
}
RC_FHS <- RC_FHS[61:n_FHS]  

# ------ MC ------
n_MC <- length(VaR99_MC_10day)

VaR_moving_avg_MC <- rep(NA, n_MC)
for (t in 61:n_MC) {  
  VaR_moving_avg_MC[t] <- mean(VaR99_MC_10day[(t-60):(t-1)])
}

RC_MC <- rep(NA, n_MC)
for (t in 61:n_MC) {  
  RC_MC[t] <- max(VaR99_MC_10day[t], k * VaR_moving_avg_MC[t])
}
RC_MC <- RC_MC[61:n_MC]  


ylim_range <- range(c(RC_FHS, RC_MC), na.rm = TRUE)

plot(RC_FHS, type = "l", col = "blue", lwd = 2, 
     main = "Trading Book Regulatory Capital (FHS vs MC)", 
     xlab = "Days", ylab = "RC", ylim = ylim_range)
lines(RC_MC, col = "red", lwd = 2)




