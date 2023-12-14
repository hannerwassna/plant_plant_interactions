
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

ddCT$Group <- cut(seq(nrow(ddCT)),
                 breaks = 8, labels = 
                   c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))

library(ggplot2)

ggplot(ddCT, aes(x = Group, y = TOR, fill=Group)) +
  geom_boxplot() +
  labs(title = "gène TOR",
       x = "Groupe",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = Group, y = TOR, fill=Group)) +
  geom_violin() +
  labs(title = "gène TOR",
       x = "Groupe",
       y = "Expression relative (ddCT)")


ggplot(ddCT, aes(x = Group, y = WRKY, fill=Group)) +
  geom_boxplot() +
  labs(title = "gene WRKY",
       x = "Groupe",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = Group, y = WRKY, fill=Group)) +
  geom_violin() +
  labs(title = "gene WRKY",
       x = "Groupe",
       y = "Expression relative (ddCT)")


ggplot(ddCT, aes(x = Group, y = D14,fill=Group)) +
  geom_violin() +
  labs(title = "gene D14",
       x = "Groupe",
       y = "Expression relative (ddCT)")

mod<-lm(D14~Group,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(TOR~Group,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(WRKY~Group,data=ddCT)
summary(mod)
anova(mod)

ddCT<-ddCT %>% 
  mutate(diversity=1)

ddCT[ddCT$Group%in%c("A/M","P/M","P/A"),"diversity"]<-2
ddCT[ddCT$Group%in%c("P/A/M"),"diversity"]<-3
ddCT[ddCT$Group%in%c("Venezio"),"diversity"]<-0

ddCT$diversity<-as.factor(ddCT$diversity)


mod<-lm(D14~diversity,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(TOR~diversity,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(WRKY~diversity,data=ddCT)
summary(mod)
anova(mod)

ggplot(ddCT, aes(x = diversity, y = TOR, fill= diversity)) +
  geom_boxplot() +
  labs(title = "TOR gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = diversity, y = TOR, fill= diversity)) +
  geom_violin() +
  labs(title = "TOR gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = diversity, y = WRKY, fill= diversity)) +
  geom_boxplot() +
  labs(title = "WRKY  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = diversity, y = WRKY, fill= diversity)) +
  geom_violin() +
  labs(title = "WRKY  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = diversity, y = D14, fill= diversity)) +
  geom_violin() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

ggplot(ddCT, aes(x = diversity, y = D14, fill= diversity)) +
  geom_boxplot() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

#Qui est significative dans D14
library(stats)
library(multcomp)

tukey_results <- glht(mod, linfct = mcp(diversity = "Tukey"), alternative = "two.sided")

summary(tukey_results)
# Afficher les résultats
summary(tukey_results)
test <- TukeyHSD(mod)

##Est ce que la nomalité des résidus est respectée ?
shapiro.test(mod$residuals)

##Est ce que l'homoscdecaticité est respectée ?

library(lmtest)
mod <- lm(D14 ~ diversity, data = ddCT)
bp_result <- bptest(mod)
print(bp_result)

##Test de robusté

library(robustbase)
mod_robuste <- lmrob(D14~diversity, data = ddCT)
summary(mod_robuste)

plot(residuals(mod_robuste))

mod_non_robuste <- lm(D14 ~ diversity, data = ddCT)

             # Comparaison des coefficients et des erreurs standard
summary(mod_non_robuste)
summary(mod_robuste)

plot(residuals(mod_non_robuste))
plot(residuals(mod_robuste))

##Modèle linéaire généralisé
mod <- glm(D14 ~ diversity, data = ddCT, family = poisson)

summary(mod)
##Avec le level 1

level1 <- ddCT[ddCT$diversity %in% c("0", "1"), ]

ggplot(level1, aes(x = diversity, y = D14, fill= diversity)) +
  geom_boxplot() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

mod<-lm(D14~diversity,data=level1)
summary(mod)
anova(mod)

##Avec le level 2

level2 <- ddCT[ddCT$diversity %in% c("0", "2"), ]

ggplot(level2, aes(x = diversity, y = D14, fill= diversity)) +
  geom_boxplot() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

mod<-lm(D14~diversity,data=level2)
summary(mod)
anova(mod)

##Avec le level 3

level3 <- ddCT[ddCT$diversity %in% c("0", "3"), ]

ggplot(level3, aes(x = diversity, y = D14, fill= diversity)) +
  geom_boxplot() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

mod<-lm(D14~diversity,data=level3)
summary(mod)
anova(mod)


##Mesure de la référence
ggplot(dataset, aes(x=Conditions, y= `ABC-like`, fill=Conditions)) +
  geom_boxplot()

mod<-lm(`ABC-like`~Conditions,data=dataset)
summary(mod)
anova(mod)

#Sur le plan d'une espèces

#Peas
data_peas <-ddCT %>% 
  mutate(Species="")

data_peas[data_peas$Group%in%c("Peas","P/A","P/M","P/A/M"),"Species"]<-"Peas" 
data_peas[data_peas$Group%in%c("Venezio"),"Species"]<-"Venezio"

data_Peas <- data_peas[data_peas$Species %in% c("Venezio", "Peas"), ]

ggplot(data_Peas, aes(x = Species, y = TOR, fill= Species)) +
  geom_boxplot() +
  labs(title = "TOR gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(TOR~Species,data=data_Peas)
summary(mod)
anova(mod)

ggplot(data_Peas, aes(x = Species, y = WRKY, fill= Species)) +
  geom_boxplot() +
  labs(title = "WRKY gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(WRKY~Species,data=data_Peas)
summary(mod)
anova(mod)

ggplot(data_Peas, aes(x = Species, y = WRKY, fill= Species)) +
  geom_boxplot() +
  labs(title = "D14 gene",
       x = "Species",
       y = "Expression relative (ddCT)")


mod<-lm(D14~Species,data=data_Peas)
summary(mod)
anova(mod)

#Maize

data_maize <-ddCT %>% 
  mutate(Species="")

data_maize[data_maize$Group%in%c("Maize","A/M","P/M","P/A/M"),"Species"]<-"Maize" 
data_maize[data_maize$Group%in%c("Venezio"),"Species"]<-"Venezio"

data_Maize <- data_maize[data_maize$Species %in% c("Venezio", "Maize"), ]

ggplot(data_Maize, aes(x = Species, y = TOR, fill= Species)) +
  geom_boxplot() +
  labs(title = "TOR gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(TOR~Species,data=data_Maize)
summary(mod)
anova(mod)

ggplot(data_Maize, aes(x = Species, y =WRKY, fill= Species)) +
  geom_boxplot() +
  labs(title = "wRKY gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(WRKY~Species,data=data_Maize)
summary(mod)
anova(mod)

ggplot(data_Maize, aes(x = Species, y =D14, fill= Species)) +
  geom_boxplot() +
  labs(title = "D14 gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(WRKY~Species,data=data_Maize)
summary(mod)
anova(mod)

#Arabidopsis

data_arabi <-ddCT %>% 
  mutate(Species="")

data_arabi[data_arabi$Group%in%c("Arabidopsis","A/M","P/A","P/A/M"),"Species"]<-"Arabidopsis" 
data_arabi[data_arabi$Group%in%c("Venezio"),"Species"]<-"Venezio"

data_Arabi <- data_arabi[data_arabi$Species %in% c("Venezio", "Arabidopsis"), ]

ggplot(data_Arabi, aes(x = Species, y = TOR, fill= Species)) +
  geom_boxplot() +
  labs(title = "TOR gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(TOR~Species,data=data_Arabi)
summary(mod)
anova(mod)

ggplot(data_Arabi, aes(x = Species, y = WRKY, fill= Species)) +
  geom_boxplot() +
  labs(title = "WRKY gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(WRKY~Species,data=data_Arabi)
summary(mod)
anova(mod)


ggplot(data_Arabi, aes(x = Species, y = D14, fill= Species)) +
  geom_boxplot() +
  labs(title = "D14 gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(D14~Species,data=data_Arabi)
summary(mod)
anova(mod)

#Toute les espèces

data_species <- rbind(data_Peas,data_Arabi,data_Maize)

ddCT <-ddCT %>% 
  mutate(Species="")

ggplot(data_species, aes(x = Species, y = D14, fill= Species)) +
  geom_boxplot() +
  labs(title = "D14 gene",
       x = "Species",
       y = "Expression relative (ddCT)")

mod<-lm(D14~Species,data=data_species)
summary(mod)
anova(mod)

ddCT[ddCT$Group%in%c("A/M"),"Species"]<- "Arabidopsis,Maize"
ddCT[ddCT$Group%in%c("P/M"),"Species"]<- "Peas,Maize"
ddCT[ddCT$Group%in%c("P/A"),"Species"]<- "Peas,Arabidopsis"
ddCT[ddCT$Group%in%c("P/A/M"),"Species"]<- "Peas, Arabidopsis, Maize"
ddCT[ddCT$Group%in%c("Venezio"),"Species"]<-"Venezio"
ddCT[ddCT$Group%in%c("Arabidopsis"),"Species"]<-"Arabidopsis"
ddCT[ddCT$Group%in%c("Peas"),"Species"]<-"Peas"
ddCT[ddCT$Group%in%c("Maize"),"Species"]<-"Maize"


ggplot(ddCT, aes(x = Species, y = D14, fill= Species)) +
  geom_boxplot() +
  labs(title = "D14 gene",
       x = "Species",
       y = "Expression relative (ddCT)")


mod<-lm(D14~Species,data=ddCT)
summary(mod)
anova(mod)






#MODELE NON PARAMETRIQUE

library(ARTool)

   ##Modele effet des niveaux de diversité

#D14
ddCTD14 <- ddCT[, "D14", drop = FALSE]
ddCTD14 <- cbind(ddCTD14, diversity = ddCT$diversity)
ddCTD14 <- na.omit(ddCTD14)

mod <- art(D14~diversity,data=ddCTD14)
summary(mod)

#TOR
ddCTTOR <- ddCT[, "TOR", drop = FALSE]
ddCTTOR <- cbind(ddCTTOR, diversity = ddCT$diversity)
ddCTTOR <- na.omit(ddCTTOR)

mod <- art(TOR~diversity,data=ddCTTOR)
summary(mod)

#WRKY

ddCTWRKY <- ddCT[, "WRKY", drop = FALSE]
ddCTWRKY <- cbind(ddCTWRKY, diversity = ddCT$diversity)
ddCTWRKY <- na.omit(ddCTWRKY)

mod <- art(WRKY~diversity,data=ddCTWRKY)
summary(mod)

    #Effet de groupe

#D14
ddCTD14 <- ddCT[, "D14", drop = FALSE]
ddCTD14 <- cbind(ddCTD14, diversity = ddCT$Group)
ddCTD14 <- na.omit(ddCTD14)


mod <- art(D14~diversity,data=ddCTD14)
summary(mod)

#WRKY
ddCTWRKY <- ddCT[, "WRKY", drop = FALSE]
ddCTWRKY<- cbind(ddCTWRKY, diversity = ddCT$Group)
ddCTWRKY <- na.omit(ddCTWRKY)

mod <- art(WRKY~diversity,data=ddCTWRKY)
summary(mod)

#TOR

ddCTTOR <- ddCT[, "TOR", drop = FALSE]
ddCTTOR <- cbind(ddCTTOR, diversity = ddCT$Group)
ddCTTOR <- na.omit(ddCTTOR)

mod <- art(TOR~diversity,data=ddCTTOR)
summary(mod)


    ##Modèle 2 facteurs

ddCT$diversity <- factor(ddCT$diversity)
ddCT$Group <- factor(ddCT$Group)

#Virer les valeurs manquantes de D14 

ddCTD14 <- ddCT[, "D14", drop = FALSE]
ddCTD14 <- cbind(ddCTD14, Group = ddCT$Group)
ddCTD14 <- cbind(ddCTD14, diversity = ddCT$diversity)
ddCTD14 <- na.omit(ddCTD14)

#modèle deux facteur interaction(groupe, diveristy)

mod_art_2facteurs <- art(formula = D14 ~ interaction(diversity*Group), data = ddCTD14)
summary(mod_art_2facteurs)


