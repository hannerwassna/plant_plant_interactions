
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

dataset$Condition <- cut(seq(nrow(dataset)),
                         breaks = 8, labels = 
                           c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))
               
colnames(dataset) <- c("Sample","TOR1","TOR2","TOR3","ABC-like1","ABC-like2","ABC-like3","WRKY1","WRKY2","WRKY3","D141","D142","D143","Conditions")

#Moyenne des triplicats

dataset$TOR <- rowMeans(dataset[, c("TOR1", "TOR2", "TOR3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$`ABC-like`<- rowMeans(dataset[, c("ABC-like1", "ABC-like2", "ABC-like3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$WRKY<- rowMeans(dataset[, c("WRKY1", "WRKY2", "WRKY3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$D14<- rowMeans(dataset[, c("D141", "D142", "D143")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

#Package qPCRtools

###Refaire des tables analysable par le package
library(dplyr)
library(qPCRtools)




##Calcul des ddCT
library(readxl)

setwd("C:/Users/Elève/Desktop/D4")
ddCT <- read_excel("ddCT.xlsx")

ddCT$Groupe <- cut(seq(nrow(ddCT)),
                      breaks = 8, labels = 
                        c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))



