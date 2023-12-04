
##Données qPCR

# Charger données CT de la plaque avec le package readxl
library(readxl)

setwd("C:/Users/Elève/Desktop/D4")
dataset <- read_excel("dataCT.xlsx")

#Séparer les gènes visés et le gènes de référence

gene_reference <- dataset[, c(5,6,7)]
genes <- dataset[,c(2,3,4,8,9,10,11,12,13)]

#regrouper les conditions -> Créer 8 groupes

library(dplyr)

dataset$Groupe <- cut(seq(nrow(dataset)),
                      breaks = 8, labels = 
                        c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))
               

##Calcul des ddCT
library(readxl)

setwd("C:/Users/Elève/Desktop/D4")
ddCT <- read_excel("ddCT.xlsx")

ddCT$Groupe <- cut(seq(nrow(ddCT)),
                      breaks = 8, labels = 
                        c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))



