#Packages loading

#install.packages("FactoMineR")
#install.packages("factoextra")  
#install.packages("ggplot2")
#install.packages("dplyr") 
#install.packages("lme4")
#install.packages("Matrix")
#install.packages("car")
#install.packages("coin")

library(FactoMineR)
library(dplyr)
library (ggplot2)
library(factoextra)
library(tidyverse)
library(tidyr)
library(lme4)
library(Matrix)
library(car)
library(coin)
library(pwr)

#Datas importation

ddCT_2 <- read_excel("JRL résultats/ddCT_2.xlsx")

#Datas framing 
row_to_suppress <- c(1, 12, 16, 17) #these are the numbers of row with 2 or more NA values
ddCT_2_ <- ddCT_2[-row_to_suppress, ]

#We are creating the variable "group" 
group_of_pots <- data.frame (group = rep(c("1","2","3","4","5","6","7","8"), times = c(3,4,3,3,3,4,4,4)))
ddCT_2_B <- cbind(ddCT_2_, group_of_pots)

#We are creating the variable "level of diversity"
level_of_diversity <- data.frame(levels = rep(c("1", "2", "3", "4"), times = c(3, 10, 11, 4)))
ddCT_M_2 <- cbind(ddCT_2_B ,level_of_diversity)

#We are creating a variable "plate" which specify from which plate the sample compes from
plate <- data.frame(qPCR_plate = rep(c("1", "2"), times = c(13, 15)))
ddCT_M_3 <- cbind(ddCT_M_2, plate)

#We are creating a column to specify each species is present in each situation 
ddCT_M_3$Peas <- rep(c("0","1","0","1"), times = c(3,4,9,12))
ddCT_M_3$Maize <- rep(c("0","1","0","1"), times = c(10,10,4,4))
ddCT_M_3$Arabidopsis <- rep(c("0","1","0","1","0","1"), times = c(7,3,3,3,4,8))

#ACP
acp_data <- PCA(ddCT_1_)
fviz_screeplot(acp_data)

#HOMOSCEDASTICITY TEST - Bartlett test
bartlett_result <- bartlett.test(D14, data = ddCT_M_3)

# Affichage des résultats
print(bartlett_result)

#MIXED MODEL WITH LMER ----> DID NOT USED

#Mixed model (pour le gène D14)
ddCT_model_D14 <- lmer(ddCT_M_3$D14 ~ 1 + levels + group + (1|group) + (1|qPCR_plate) + levels:group + (1 | ddCT_M_3$D14), data = ddCT_M_3)
vif_result_D14 <- vif(ddCT_model)

#Mixed model for TOR
ddCT_model_TOR <- lmer(ddCT_M_3$TOR ~ 1 + levels + group + (1|group) + (1|qPCR_plate) + levels:group + (1 | ddCT_M_3$TOR), data = ddCT_M_3)
vif_result_TOR <- vif(ddCT_model_TOR)

#Mixed model for WRKY
ddCT_model_WRKY <- lmer(ddCT_M_3$WRKY ~ 1 + levels + group + (1|group) + (1|qPCR_plate) + levels:group + (1 | ddCT_M_3$WRKY), data = ddCT_M_3)
vif_result_WRKY <- vif(ddCT_model_WRKY)


#Mixed model with sommer
install.packages("sommer")
library("sommer")

mmod<-mmer(fixed = D14~qPCR_plate+levels,
           rcov=~vsr(dsr(levels),units),
           data=ddCT_M_3)

summary(mmod)

anova(mmod)

#KRUSKAL WALLIS TEST 

#On d14 gene
mod<-lm(D14~qPCR_plate+levels,data=ddCT_M_3)
anova(mod)
plot(mod) 

kruskal.test(D14~levels,data=ddCT_M_3)
kruskal.test(D14~group,data=ddCT_M_3)
kruskal.test(D14~Peas, data=ddCT_M_3)
kruskal.test(D14~Arabidopsis, data=ddCT_M_3)
kruskal.test(D14~Maize, data=ddCT_M_3)
kruskal.test(D14~qPCR_plate, data=ddCT_M_3)

#On TOR gene
kruskal.test(TOR~levels,data=ddCT_M_3)
kruskal.test(TOR~group,data=ddCT_M_3)
kruskal.test(TOR~Peas, data=ddCT_M_3)
kruskal.test(TOR~Arabidopsis, data=ddCT_M_3)
kruskal.test(TOR~Maize, data=ddCT_M_3)
kruskal.test(TOR~qPCR_plate, data=ddCT_M_3)

#On WRKY gene
kruskal.test(WRKY~levels,data=ddCT_M_3)
kruskal.test(WRKY~group,data=ddCT_M_3)
kruskal.test(WRKY~Peas, data=ddCT_M_3)
kruskal.test(WRKY~Arabidopsis, data=ddCT_M_3)
kruskal.test(WRKY~Maize, data=ddCT_M_3)
kruskal.test(WRKY~qPCR_plate, data=ddCT_M_3)
