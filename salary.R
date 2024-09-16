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

options(scipen = 100)
#############################################################################################################################
#2. Salary Data 
################################################################
load("C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/ex_salary.Rdata")

str(salary)
colnames(salary) = c("age", "pay", "bonus", "hours", "number_of_worker", "career", "sex")
pairs(salary[,c("pay", "bonus", "hours", "number_of_worker")])

#########################################
### Original #########
#########################################

salary3 = salary 
salary3$age
salary3$pay = scale(salary3$pay)
salary3$bonus = scale(salary3[,"bonus"])  
salary3$hours = scale(salary3[,"hours"])
salary3$number_of_worker = scale(salary3[,"number_of_worker"])  #/100

salary3$age = as.character(salary3$age)
salary3$age[salary3$age == "-19"] = "young"
salary3$age[salary3$age == "20-24"] = "young"
salary3$age[salary3$age == "25-29"] = "young"
salary3$age[salary3$age == "30-34"] = "young"
salary3$age[salary3$age == "35-39"] = "young"
salary3$age[salary3$age == "40-44"] = "old"
salary3$age[salary3$age == "45-49"] = "old"
salary3$age[salary3$age == "50-54"] = "old"
salary3$age[salary3$age == "55-59"] = "old"
salary3$age[salary3$age == "60-"] = "old"

salary[31,]
str(salary3)
table(salary3$career)

salary3 = transform(salary3, 
			 sex_man = ifelse(sex == "남", 1, 0), 
			career_10 = ifelse(career == "10년이상", 1, 0),
			age_yo = ifelse(age == "old", 1, 0))

table(salary3$career_10)
table(salary3$age_yo)

salary_data = salary3 %>% dplyr::select(pay, bonus, hours, number_of_worker, sex_man, career_10, age_yo)
str(salary_data)
head(salary_data)

save(salary_data, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/salary_data.Rdata")

attach(salary3)
lse3 = lm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, data = salary3)
lts3 = ltsReg(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, alpha = 0.7, data = salary3)
M3 = rlm(pay ~ bonus + hours + number_of_worker + sex_man + career_10+ age_yo, data = salary3, method = "M")



y3 = as.matrix(pay)

x3 = as.matrix(cbind(bonus, hours, number_of_worker, sex_man,
					 career_10, age_yo))

NSNSM3 = AFT_NSNSM(x = x3, y = y3)


eps_lse3 = y3 -  cbind(1,x3) %*% lse3$coefficients
eps_lts3 = y3 -  cbind(1,x3) %*% lts3$coefficients
eps_M3 = y3 -  cbind(1,x3) %*% M3$coefficients
eps_NSNSM3 = y3 -  cbind(1,x3) %*% NSNSM3$beta

detach(salary3)

plot(lse3, which = 2)
plot(lse3, which = 4)

lse3$coefficients
lts3$coefficients
M3$coefficients
NSNSM3$beta

#MAE
mean(abs(eps_lse3)) 
mean(abs(eps_lts3)) 
mean(abs(eps_M3)) 
mean(abs(eps_NSNSM3)) 

par(mar=c(5,5,5,5))
hist(eps_NSNSM3, cex.axis = 1.5, cex.lab = 2.5, xlab = TeX(r"($/epsilon$)"), main = "", freq = F, breaks = 10)
density_ssnsm(model = NSNSM3, true_epsilon = eps_NSNSM3, ylim = c(0, 1.3), breaks = 10)

############
#########################################
### Outlier is deleted #########
#########################################

salary4 = salary[-31,] 
salary4$age
salary4$pay = scale(salary4$pay)
salary4$bonus = scale(salary4[,"bonus"]) 
salary4$hours = scale(salary4[,"hours"])
salary4$number_of_worker = scale(salary4[,"number_of_worker"])  

salary4 = transform(salary4, 
			 sex_man = ifelse(sex == "남", 1, 0), 
			career_10 = ifelse(career == "10년이상", 1, 0),
			age_yo = ifelse(age == "old", 1, 0))

attach(salary4)
lse4 = lm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, data = salary4)
lts4 = ltsReg(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, alpha = 0.7, data = salary4)
M4 = rlm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, data = salary4, method = "M")

y4 = as.matrix(pay)
x4 = as.matrix(cbind(bonus, hours, number_of_worker, sex_man,
					 career_10, age_yo))

NSNSM4 = AFT_NSNSM(x = x4, y = y4)

eps_lse4 = y4 -  cbind(1,x4) %*% lse4$coefficients
eps_lts4 = y4 -  cbind(1,x4) %*% lts4$coefficients
eps_M4 = y4 -  cbind(1,x4) %*% M4$coefficients
eps_NSNSM4 = y4 -  cbind(1,x4) %*% NSNSM4$beta

detach(salary4)

plot(lse4, which = 2)
plot(lse4, which = 4)

sum(abs(lse3$coefficients - lse4$coefficients))
sum(abs(lts3$coefficients - lts4$coefficients))
sum(abs(M3$coefficients - M4$coefficients))
sum(abs(NSNSM3$beta - NSNSM4$beta))

lse4$coefficients
lts4$coefficients
M4$coefficients
NSNSM4$beta

#MAE
mean(abs(eps_lse4)) 
mean(abs(eps_lts4)) 
mean(abs(eps_M4)) 
mean(abs(eps_NSNSM4)) 

par(mar=c(5,5,5,5))
hist(eps_NSNSM4, cex.axis = 1.5, cex.lab = 2.5, xlab = TeX(r"($/epsilon$)"), main = "", freq = F, breaks = 10)
density_ssnsm(model = NSNSM4, true_epsilon = eps_NSNSM4, ylim = c(0, 1.5), breaks = 10)



############################################################
#2. leave-one out 

load("C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/prediction.Rdata")

##########################################################
#cross validation
data = salary3 
k = nrow(data)

set.seed(19873)
data$id <- 1:k #sample(1:k, nrow(data), replace = TRUE)
list <- 1:k
head(data)

str(data)

training_set = list()
test_set = list()
train_x = list()
train_y = list()
test_x = list()
test_y = list()

for(set in 1 : k){

	training_set[[set]] <- as.data.frame(subset(data, id %in% list[-set]))		# remove rows with id i from dataframe to create training set
#	training_set[[set]]$sex = ifelse(training_set[[set]]$sex == "남", 1, 2)
	test_set[[set]] <- as.data.frame(subset(data, id %in% c(set)))				# select rows with id i to create test set
#	test_set[[set]]$sex = ifelse(test_set[[set]]$sex == "남", 1, 2)

	train_x[[set]] = as.matrix(training_set[[set]][,c(3:5,8:10)])
	train_y[[set]] = (training_set[[set]][,2])

	test_x[[set]] = as.matrix(test_set[[set]][,c(3:5,8:10)])
	test_y[[set]] = (test_set[[set]][,2])

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
	lse_temp[[set]] = lm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo
, data = training_set[[set]])

	eps_lse = append(eps_lse, test_y[[set]] - cbind(1,test_x[[set]]) %*% lse_temp[[set]]$coefficients)
}

mean(abs(eps_lse))

#LTS
lts_temp = list()
eps_lts = NULL

for(set in 1 : k){
	lts_temp[[set]] = ltsReg(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, alpha = 0.7, data = training_set[[set]])
	eps_lts = append(eps_lts, test_y[[set]] - cbind(1,test_x[[set]]) %*% lts_temp[[set]]$coefficients)
}

mean(abs(eps_lts))

#M
m_temp = list()
eps_m = NULL

for(set in 1 : k){
	m_temp[[set]] = rlm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo,  method = "M", data = training_set[[set]])
	eps_m = append(eps_m, test_y[[set]] - cbind(1,test_x[[set]]) %*% m_temp[[set]]$coefficients)
}

mean(abs(eps_m))

mean(abs(eps_lse))
mean(abs(eps_lts))
mean(abs(eps_m))
mean(abs(eps_nsnsm))

save(nsnsm_temp, lse_temp, lts_temp, m_temp, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/prediction.Rdata")

##########################################################
#cross validation2
data = salary4 
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
#	training_set[[set]]$sex = ifelse(training_set[[set]]$sex == "남", 1, 2)
	test_set[[set]] <- as.data.frame(subset(data, id %in% c(set)))				# select rows with id i to create test set
#	test_set[[set]]$sex = ifelse(test_set[[set]]$sex == "남", 1, 2)

	train_x[[set]] = as.matrix(training_set[[set]][,c(2,4,5,8)])
	train_y[[set]] = (training_set[[set]][,3])

	test_x[[set]] = as.matrix(test_set[[set]][,c(2,4,5,8)])
	test_y[[set]] = (test_set[[set]][,3])

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
	lse_temp[[set]] = lm(bonus ~ pay + hours + number_of_worker + sex_man
, data = training_set[[set]])

	eps_lse = append(eps_lse, test_y[[set]] - cbind(1,test_x[[set]]) %*% lse_temp[[set]]$coefficients)
}

mean(abs(eps_lse))

#LTS
lts_temp = list()
eps_lts = NULL

for(set in 1 : k){
	lts_temp[[set]] = ltsReg(bonus ~ pay + hours + number_of_worker + sex_man, alpha = 0.8, data = training_set[[set]])
	eps_lts = append(eps_lts, test_y[[set]] - cbind(1,test_x[[set]]) %*% lts_temp[[set]]$coefficients)
}

mean(abs(eps_lts))

#M
m_temp = list()
eps_m = NULL

for(set in 1 : k){
	m_temp[[set]] = rlm(bonus ~ pay + hours + number_of_worker + sex_man,  method = "M", data = training_set[[set]])
	eps_m = append(eps_m, test_y[[set]] - cbind(1,test_x[[set]]) %*% m_temp[[set]]$coefficients)
}

mean(abs(eps_m))

mean(abs(eps_lse))
mean(abs(eps_lts))
mean(abs(eps_m))
mean(abs(eps_nsnsm))

save(nsnsm_temp, lse_temp, lts_temp, m_temp, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/prediction.Rdata")
 
##################################################
### Bootstrap
###################################################

B = 500 # Bootstrap size

salary3_boots = list()

set.seed(7821)
boots_num = lapply(1:B, FUN = function(i){ 
						sample(1:nrow(salary3), 
						size = nrow(salary3), replace = TRUE)}
					)

salary3_boots = lapply(1:B, FUN = function (i){
						salary3[boots_num[[i]],]
					})


y2 = list()
for(i in 1 : B){y2[[i]] = as.matrix(salary3_boots[[i]]$pay)}


x2 = list()
for(i in 1 : B){
	attach(salary3_boots[[i]])
	x2[[i]] = cbind(bonus, hours, number_of_worker, sex_man, career_10, age_yo)

	detach(salary3_boots[[i]])
}


#OLS
OLS_boots = lapply(1:B, FUN = function(i){
					lm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, data = salary3_boots[[i]])

							}
					)

OLS_boots_coeff = sapply(1:B, FUN = function(i){
							OLS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((OLS_boots_coeff - apply(OLS_boots_coeff,1,mean))^2,1,sum)/ (B-1))


#LTS
LTS_boots = lapply(1:B, FUN = function(i){
					ltsReg(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, alpha = 0.7, data = salary3_boots[[i]])

							}
					)

LTS_boots_coeff = sapply(1:B, FUN = function(i){
							LTS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((LTS_boots_coeff - apply(LTS_boots_coeff,1,mean))^2,1,sum)/ (B-1))

#Huber
Huber_boots = lapply(1:B, FUN = function(i){
					rlm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, method = "M", data = salary3_boots[[i]])

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

save(OLS_boots,LTS_boots,Huber_boots,SSNSM_boots, file = "C:/Users/HP/Desktop/3. AFT using NSNSM/Real data/Robust liner regression/salary/Boots.Rdata")




##################################################
### Bootstrap (deleted)
###################################################

B = 500 # Bootstrap size

salary4_boots = list()

set.seed(7821)
boots_num = lapply(1:B, FUN = function(i){ 
						sample(1:nrow(salary4), 
						size = nrow(salary4), replace = TRUE)}
					)

salary4_boots = lapply(1:B, FUN = function (i){
						salary4[boots_num[[i]],]
					})


y2 = list()
for(i in 1 : B){y2[[i]] = as.matrix(salary4_boots[[i]]$pay)}


x2 = list()
for(i in 1 : B){
	attach(salary4_boots[[i]])
	x2[[i]] = cbind(bonus, hours, number_of_worker, sex_man, career_10, age_yo)

	detach(salary4_boots[[i]])
}


#OLS
OLS_boots = lapply(1:B, FUN = function(i){
					lm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, data = salary4_boots[[i]])

							}
					)

OLS_boots_coeff = sapply(1:B, FUN = function(i){
							OLS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((OLS_boots_coeff - apply(OLS_boots_coeff,1,mean))^2,1,sum)/ (B-1))


#LTS
LTS_boots = lapply(1:B, FUN = function(i){
					ltsReg(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, alpha = 0.7, data = salary4_boots[[i]])

							}
					)

LTS_boots_coeff = sapply(1:B, FUN = function(i){
							LTS_boots[[i]]$coefficient
							}
					) 

sqrt(apply((LTS_boots_coeff - apply(LTS_boots_coeff,1,mean))^2,1,sum)/ (B-1))

#Huber
Huber_boots = lapply(1:B, FUN = function(i){
					rlm(pay ~ bonus + hours + number_of_worker + sex_man + career_10 + age_yo, method = "M", data = salary4_boots[[i]])

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



######################################################################################