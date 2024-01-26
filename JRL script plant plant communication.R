
##Données qPCR

# Load CT data from the plate with the readxl package
library(readxl)

setwd("C:/Users/Elève/Desktop/D4/Plnt Plnt com/data")
dataset <- read_excel("dataCT.xlsx")

#Separate the targeted genes and the reference genes 

gene_reference <- dataset[, c(5,6,7)]
genes <- dataset[,c(2,3,4,8,9,10,11,12,13)]


#group conditions -> Create 8 groups

library(dplyr)

dataset$Condition <- cut(seq(nrow(dataset)),
                         breaks = 8, labels = 
                           c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))
               
colnames(dataset) <- c("Sample","TOR1","TOR2","TOR3","ABC-like1","ABC-like2","ABC-like3","WRKY1","WRKY2","WRKY3","D141","D142","D143","Conditions")

#Average of triplicats

dataset$TOR <- rowMeans(dataset[, c("TOR1", "TOR2", "TOR3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$`ABC-like`<- rowMeans(dataset[, c("ABC-like1", "ABC-like2", "ABC-like3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$WRKY<- rowMeans(dataset[, c("WRKY1", "WRKY2", "WRKY3")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

dataset$D14<- rowMeans(dataset[, c("D141", "D142", "D143")],na.rm = TRUE)
dataset <- dataset[, -c(2:4)]

#Overview of the dataset
  ##Mesure de la référence
ggplot(dataset, aes(x=Conditions, y= `ABC-like`, fill=Conditions)) +
  geom_boxplot()

mod<-lm(`ABC-like`~Conditions,data=dataset)
summary(mod)

#Test with the qPCRtools Package, and switch to ddCT

library(dplyr)
library(qPCRtools)

library(readxl)

#From ddCT of 3 genes

setwd("C:/Users/Elève/Desktop/D4")
ddCT <- read_excel("ddCT.xlsx")

ddCT$Group <- cut(seq(nrow(ddCT)),
                 breaks = 8, labels = 
                   c("Venezio", "Peas","Maize","Arabidopsis","A/M","P/M","P/A","P/A/M"))


#Plot of relative expression for each gene
library(ggplot2)

  #Plot Tor by group
    #Boxplot
ggplot(ddCT, aes(x = Group, y = TOR, fill=Group)) +
  geom_boxplot() +
  labs(title = "Tor gene",
       x = "Group",
       y = "Relative expression")

    #Violin
ggplot(ddCT, aes(x = Group, y = TOR, fill=Group)) +
  geom_violin() +
  labs(title = "gène TOR",
       x = "Groupe",
       y = "Expression relative (ddCT)")

  #Plot Wrky by group
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

  #Plot D14 by group
ggplot(ddCT, aes(x = Group, y = D14,fill=Group)) +
  geom_violin() +
  labs(title = "gene D14",
       x = "Groupe",
       y = "Expression relative (ddCT)")

#First fruits of analysis

#Linear regression and first ANOVA

mod<-lm(D14~Group,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(TOR~Group,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(WRKY~Group,data=ddCT)
summary(mod)
anova(mod)

#Plot of relative expression explained by the level of diversity

  #Regroup group by the level of diversity
ddCT<-ddCT %>% 
  mutate(diversity=1)

ddCT[ddCT$Group%in%c("A/M","P/M","P/A"),"diversity"]<-2
ddCT[ddCT$Group%in%c("P/A/M"),"diversity"]<-3
ddCT[ddCT$Group%in%c("Venezio"),"diversity"]<-0

ddCT$diversity<-as.factor(ddCT$diversity)

  #linear regression & anova
mod<-lm(D14~diversity,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(TOR~diversity,data=ddCT)
summary(mod)
anova(mod)

mod<-lm(WRKY~diversity,data=ddCT)
summary(mod)
anova(mod)

  #Plots
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

#Comparison between diversity levels for the D14 gene
library(stats)
library(multcomp)

tukey_results <- glht(mod, linfct = mcp(diversity = "Tukey"), alternative = "two.sided")

summary(tukey_results)
      # results
summary(tukey_results)
test <- TukeyHSD(mod)

##Are the normality of the residues respected?
shapiro.test(mod$residuals)

##Is homoscdecacity respected?

library(lmtest)
mod <- lm(D14 ~ diversity, data = ddCT)
bp_result <- bptest(mod)
print(bp_result)

##Robusty test

library(robustbase)
mod_robuste <- lmrob(D14~diversity, data = ddCT)
summary(mod_robuste)

plot(residuals(mod_robuste))

mod_non_robuste <- lm(D14 ~ diversity, data = ddCT)

    # Comparison of coefficients and standard errors
summary(mod_non_robuste)
summary(mod_robuste)

plot(residuals(mod_non_robuste))
plot(residuals(mod_robuste))

    ##Generalized linear model
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


#First fruit of analysis grouped by species 

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


#MODELE NON PARAMETRIQUE

#KRUSKAL WALLIS TEST 'retained for analysis"

#On d14 gene
mod<-lm(D14~qPCR_plate+levels,data=ddCT_M_3)
anova(mod)
plot(mod) 


kruskal.test(D14~levels,data=ddCT_M_3)

ggplot(ddCT_M_3, aes(x = levels, y = D14, fill= levels)) +
  geom_boxplot() +
  labs(title = "Relative expression of D14 gene by levels of diversity",
       x = "level of diversity",
       y = "Relative expression")

ggplot(ddCT_M_3, aes(x = levels, y = D14, fill = levels)) +
  geom_boxplot() +
  labs(title = "Relative expression of D14 gene by levels of diversity",
       x = "Level of diversity",
       y = "Relative expression") +
  theme_classic()

kruskal.test(D14~group,data=ddCT_M_3)
kruskal.test(D14~Peas, data=ddCT_M_3)
kruskal.test(D14~Arabidopsis, data=ddCT_M_3)
kruskal.test(D14~Maize, data=ddCT_M_3)

boxplot(D14 ~ levels, data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of D14 by Levels")

bxp_info <- boxplot(D14 ~ levels, data = ddCT_M_3, plot = FALSE)
medians <- bxp_info$stats[3, ]  # La troisième ligne correspond aux médianes

# Affichage des médianes
cat("Médianes par niveau:\n")
cat("Niveau 1:", medians[1], "\n")
cat("Niveau 2:", medians[2], "\n")
cat("Niveau 3:", medians[3], "\n")
cat("Niveau 4:", medians[4], "\n")


#On TOR gene
kruskal.test(TOR~levels,data=ddCT_M_3)

ggplot(ddCT_M_3, aes(x = levels, y = TOR, fill= levels)) +
  geom_violin() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

boxplot(TOR ~ levels, data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of TOR by Levels")
ggplot(ddCT_M_3, aes(x =groups, y = TOR, fill= levels)) +
  geom_violin() +
  labs(title = "D14  gene",
       x = "level of diversity",
       y = "Expression relative (ddCT)")

kruskal.test(TOR~group,data=ddCT_M_3)
boxplot(TOR ~ group, data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of TOR by Levels")

boxplot(TOR ~ group, data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of TOR by groups")
kruskal.test(TOR~Peas, data=ddCT_M_3)
kruskal.test(TOR~Arabidopsis, data=ddCT_M_3)
kruskal.test(TOR~Maize, data=ddCT_M_3)

alpha <- 0.05  # Significance level
n <- 14        # Sample size per group
num_simulations <- 1000  # Number of simulations

# Function to perform the Kruskal-Wallis test and return the result
run_kruskal_test <- function() {
  result <- kruskal.test(D14 ~ levels, data = ddCT_M_3)
  return(result$p.value < alpha)
}

# Perform simulations
set.seed(123)  # For reproducibility
sim_results <- replicate(num_simulations, run_kruskal_test())

# Calculate statistical power
power <- mean(sim_results)

# Display statistical power
cat("Statistical power of the Kruskal-Wallis test:", power, "\n")

boxplot(TOR ~ Maize, data = ddCT_M_3, col = "lightblue", main = "Relative expression of TOR by the presence of Maize", fill= "Maize")

bxp_info <- boxplot(TOR ~ Maize, data = ddCT_M_3, plot = FALSE)
medians <- bxp_info$stats[3, ]  # The third row corresponds to medians

# Display medians
cat("Medians by level:\n")
cat("Level 1:", medians[1], "\n")
cat("Level 2:", medians[2], "\n")


#On WRKY gene
kruskal.test(WRKY~levels,data=ddCT_M_3)

boxplot(WRKY ~ levels, data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of WRKY by Levels")

kruskal.test(WRKY~group,data=ddCT_M_3)
kruskal.test(WRKY~Peas, data=ddCT_M_3)
kruskal.test(WRKY~Arabidopsis, data=ddCT_M_3)
kruskal.test(WRKY~Maize, data=ddCT_M_3)

boxplot(WRKY~Maize,data = ddCT_M_3, col = "lightblue", main = "Boxplot of relative expression of WRKY by Maize"), y = TOR, fill= levels, col = "lightblue", main = "Boxplot of relative expression of TOR by Maize")) +
  geom_boxplot() +
  labs(y = "Relative expression of TOR")


num_groups <- 2
sample_size <- 4
effect_size <- 0.5

pwr.chisq.test(w=0.5,N=4, df = num_groups - 1)

power <- pwr.chisq.test(kruskal.test(WRKY~Maize, data=ddCT_M_3), df = num_groups - 1, N = sample_size * num_groups, sig.level = 0.05)$power
print(power)

result <- kruskal.test(Value ~ Group, data = data)

# Return TRUE if the null hypothesis is rejected (p-value < 0.05)
return(result$p.value < 0.05)

# Perform simulations
reject_null <- replicate(num_simulations, run_kruskal_test())

# Calculate statistical power
power <- mean(reject_null)
print(power)

bxp_info <- boxplot(WRKY ~ Maize, data = ddCT_M_3, plot = FALSE)
medians <- bxp_info$stats[3, ]  # The third row corresponds to medians

# Display medians
cat("Medians by level:\n")
cat("Level 1:", medians[1], "\n")
cat("Level 2:", medians[2], "\n")

num_groups <- 4
sample_size <- 14
effect_size <- 0.5

pwr.chisq.test(w=0.5, N=sample_size, df = num_groups - 1)



#By the package ART 

library(ARTool)

  ##Model: Effect of Diversity Levels.

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

    #Group effect

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


## Two-Factor Model

ddCT$diversity <- factor(ddCT$diversity)
ddCT$Group <- factor(ddCT$Group)

# Remove missing values from D14

ddCTD14 <- ddCT[, "D14", drop = FALSE]
ddCTD14 <- cbind(ddCTD14, Group = ddCT$Group)
ddCTD14 <- cbind(ddCTD14, diversity = ddCT$diversity)
ddCTD14 <- na.omit(ddCTD14)

# Two-factor interaction model (group, diversity)

mod_art_2factors <- art(formula = D14 ~ interaction(diversity*Group), data = ddCTD14)
summary(mod_art_2factors)

library(ggplot2)

## On a logarithmic scale
ddCTlog <- ddCT
ddCTlog$TOR <- log10(ddCT$TOR)
ddCTlog$WRKY <- log10(ddCT$WRKY)
ddCTlog$D14 <- log10(ddCT$D14)

# Homoscedasticity test

library(lmtest)
model <- lm(D14 ~ diversity, data = ddCTlog)
bptest(model)

model <- lm(TOR ~ diversity, data = ddCTlog)
bptest(model)

model <- lm(WRKY ~ diversity, data = ddCTlog)
bptest(model)

# Anova

model <- lm(D14 ~ diversity, data = ddCTlog)
summary(model)
anova(model)

model <- lm(TOR ~ diversity, data = ddCTlog)
summary(model)
anova(model)

model <- lm(WRKY ~ diversity, data = ddCTlog)
summary(model)
anova(model)

#TRIES OF MIXED MODEL 
#mixed model with lmer

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

