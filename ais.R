##############################################################
library(tidyverse)
library(sn)
library(nnls)
library(latex2exp)
library(e1071)
library(survival)
library(aftgee)
library(robustbase)
library(MASS)  # M-estimate, rlm function
#install.packages("FMsmsnReg")
#install.packages("olsrr")
library(olsrr) # cook's distance plot
library(latex2exp)

options(scipen = 100)
############################################################
#1. AIS Data 
library(FMsmsnReg) #install.packages("FMsmsnReg")
##########################################################

data(ais)
colnames(ais)
str(ais)

#male = ais %>% filter(Sex == 0) #%>% select(BMI, LBM)

attach(ais)

Ferr = scale(Ferr)

par(mar=c(5,5,5,5))
plot(BMI ~ LBM, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))
plot(BMI ~ Bfat, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))
plot(BMI ~ Hc, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))
plot(BMI ~ RCC, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))
plot(BMI ~ WCC, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))
plot(BMI ~ Hg, data = ais, cex.axis = 1.5, cex.lab = 2, col = factor(Sex))

hist(BMI, cex.axis = 1.5, cex.lab = 2, freq = F, main = "")

lse = lm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr)
lts = ltsReg(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, alpha = 0.8)
M = rlm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg +  Ferr , method = "M")
NSNSM= AFT_NSNSM(x = as.matrix(cbind(BMI, Sex, LBM, Bfat, Hc,  Hg, Ferr)), y = as.matrix(RCC))

plot(lse, which = 2)
plot(lse, which = 4)

lse$coefficients
lts$coefficients
M$coefficients
NSNSM$beta

par(mar=c(5,5,5,5))
plot(BMI ~ LBM, data = ais, cex.axis = 1.5, cex.lab = 2)
abline(lse$coefficients, lty = 1)
abline(lts$coefficients, lty = 2)
abline(M$coefficients, lty = 3)
abline(NSNSM$beta, col = "red", lty = 5)

y = as.matrix(RCC)
x = as.matrix(cbind(BMI, Sex, LBM, Bfat, Hc, Hg, Ferr))

eps_lse = y - cbind(1,x) %*% lse$coefficients
eps_lts = y - cbind(1,x) %*% lts$coefficients
eps_M = y - cbind(1,x) %*% M$coefficients
eps_NSNSM = y - cbind(1,x) %*% NSNSM$beta

mean(abs(eps_lse)) 
mean(abs(eps_lts)) 
mean(abs(eps_M)) 
mean(abs(eps_NSNSM)) 

par(mar=c(5,5,5,5))
hist(eps_NSNSM, cex.axis = 1.5, cex.lab = 2.5, xlab = "Residuals", main = "", freq = F, breaks = 5)
density_ssnsm(model = NSNSM, true_epsilon = eps_NSNSM, ylim = c(0, 2.7), breaks = 15)

detach(ais)

############################################################
#2. leave-one out 

load("C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/AIS data/prediction.Rdata")

##########################################################
#cross validation
data(ais)
head(ais)
data = ais %>% dplyr::select(RCC, BMI, Sex, LBM, Bfat, Hc, Hg, Ferr)
k = nrow(data)

set.seed(19873)
data$id <- 1:k #sample(1:k, nrow(data), replace = TRUE)
list <- 1:k
head(data)

training_set = list()
test_set = list()
train_x = list()
train_y = list()
test_x = list()
test_y = list()

for(set in 1 : k){

	training_set[[set]] <- as.data.frame(subset(data, id %in% list[-set]))		# remove rows with id i from dataframe to create training set
	test_set[[set]] <- as.data.frame(subset(data, id %in% c(set)))				# select rows with id i to create test set

	train_x[[set]] = as.matrix(training_set[[set]][,2:8])
	train_y[[set]] = (training_set[[set]][,1])

	test_x[[set]] = as.matrix(test_set[[set]][,2:8])
	test_y[[set]] = (test_set[[set]][,1])

}

#nsnsm
nsnsm_temp = list()
eps_nsnsm = NULL

for(set in 1 : k){
	nsnsm_temp[[set]] = AFT_NSNSM(x = train_x[[set]], y = as.matrix(train_y[[set]]))
	eps_nsnsm = append(eps_nsnsm, test_y[[set]] - cbind(1,test_x[[set]]) %*% nsnsm_temp[[set]]$beta)
	print(set)
}

mean(abs(eps_nsnsm))

#LSE
lse_temp = list()
eps_lse = NULL

for(set in 1 : k){
	lse_temp[[set]] = lm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, data = training_set[[set]])
	eps_lse = append(eps_lse, test_y[[set]] - cbind(1,test_x[[set]]) %*% lse_temp[[set]]$coefficients)
}

mean(abs(eps_lse))

#LTS
lts_temp = list()
eps_lts = NULL

for(set in 1 : k){
	lts_temp[[set]] = ltsReg(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, alpha = 0.8, data = training_set[[set]])
	eps_lts = append(eps_lts, test_y[[set]] - cbind(1,test_x[[set]]) %*% lts_temp[[set]]$coefficients)
}

mean(abs(eps_lts))

#M
m_temp = list()
eps_m = NULL

for(set in 1 : k){
	m_temp[[set]] = rlm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, method = "M", data = training_set[[set]])
	eps_m = append(eps_m, test_y[[set]] - cbind(1,test_x[[set]]) %*% m_temp[[set]]$coefficients)
}

mean(abs(eps_m))

mean(abs(eps_lse))
mean(abs(eps_lts))
mean(abs(eps_m))
mean(abs(eps_nsnsm))


save(nsnsm_temp, lse_temp, lts_temp, m_temp, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/AIS data/prediction.Rdata")

##################################################
### 3. Bootstrap
###################################################

B = 500 # Bootstrap size

ais_boots = list()

set.seed(7821)
boots_num = lapply(1:B, FUN = function(i){ 
						sample(1:nrow(ais), 
						size = nrow(ais), replace = TRUE)}
					)

ais_boots = lapply(1:B, FUN = function (i){
						ais[boots_num[[i]],]
					})


y2 = list()
for(i in 1 : B){y2[[i]] = as.matrix(ais_boots[[i]]$RCC)}


x2 = list()
for(i in 1 : B){
	attach(ais_boots[[i]])
	x2[[i]] = cbind(BMI, Sex, LBM, Bfat, Hc, Hg, Ferr)

	detach(ais_boots[[i]])
}


#OLS
OLS_boots = lapply(1:B, FUN = function(i){
					lm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, data = ais_boots[[i]])

							}
					)

OLS_boots_coeff = sapply(1:B, FUN = function(i){
							OLS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((OLS_boots_coeff - apply(OLS_boots_coeff,1,mean))^2,1,sum)/ (B-1))


#LTS
LTS_boots = lapply(1:B, FUN = function(i){
					ltsReg(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, alpha = 0.8, data = ais_boots[[i]])

							}
					)

LTS_boots_coeff = sapply(1:B, FUN = function(i){
							LTS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((LTS_boots_coeff - apply(LTS_boots_coeff,1,mean))^2,1,sum)/ (B-1))

#Huber
Huber_boots = lapply(1:B, FUN = function(i){
					rlm(RCC ~ BMI + Sex + LBM + Bfat + Hc +  Hg + Ferr, method = "M", data = ais_boots[[i]])

							}
					)

Huber_boots_coeff = sapply(1:B, FUN = function(i){
							Huber_boots[[i]]$coefficient
							}
					) 

sqrt(apply((Huber_boots_coeff - apply(Huber_boots_coeff,1,mean))^2,1,sum)/ (B-1))


#SSNSM
SSNSM_boots = lapply(1:B, FUN = function(i){
		AFT_NSNSM(x = x2[[i]], y = y2[[i]], simulation_index = i)
							}
					)

SSNSM_boots_coeff = sapply(1:B, FUN = function(i){
							SSNSM_boots[[i]]$beta
							}
					) 

sqrt(apply((SSNSM_boots_coeff - apply(SSNSM_boots_coeff,1,mean))^2,1,sum)/ (B-1))


save(OLS_boots,LTS_boots,Huber_boots,SSNSM_boots, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/AIS data/Boots.Rdata")





