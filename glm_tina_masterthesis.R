# Tina Himmelsbach

inputDirectory = "xxx" 
setwd(inputDirectory)

library(psych)
library(corrplot)
library(tidyverse)
library(nnet) 
library(sandwich) 
library(lmtest)
library(ggplot2)
library(dplyr)

##### LOAD AND PREPARE DATA #####

## preparing the df with categories 0,1,2

#depredation conflicts = 1
df_confl_orig <- read.csv2("xxx.csv", sep=",")
df_confl_orig$confl <- 1 #define response variable
head(df_confl_orig)
df_confl_orig <- df_confl_orig[, c("confl", "SELS_mx", xxx)] #subset df

#trophy hunting = 2
df_trophy_orig <- read.csv2("xxx.csv", sep=",")
df_trophy_orig$confl <- 2 # define response variable
df_trophy_orig <- df_trophy_orig[, c("confl", "SELS_mx", xxx)] #subset df

#pseudo absences = 0 
library(sf)
ps1 <- read_sf("xxx.shp")
ps1$confl <- 0 # define response variable
ps1 <- ps1[, c("confl", "SELS_mx", xxx)] #subset df 
ps1 <- subset(ps1, !(SELS_mx %in% c("D2", "E2", "E3"))) #remove individual SELS

## stratified sampling (pseudo absences need to have the same amount of occurences per SELS)
set.seed(123) #so it is reproducible 

# Define the number of rows to select per value in SELS_mx
confltotal <- rbind(df_confl_orig, df_trophy_orig) #merge conflicts 
depr_per_SELS <- aggregate(confl ~ SELS_mx, data = confltotal, FUN = function(x) sum(grepl("1", x)))
hunt_per_SELS <- aggregate(confl ~ SELS_mx, data = confltotal, FUN = function(x) sum(grepl("2", x)))

rows_per_group <- c(A1=1, A2 = 7, A3 = 23, A4 = 7, B = 5, C1 = 36, C2 = 4, D1 = 2, D3 = 4, E1 = 1)
strat_ps1 <- data.frame()

# Loop through each group and sample the rows
for (group_value in names(rows_per_group)) {
  group_data <- ps1[ps1$SELS_mx == group_value, ]
  sampled_rows <- sample(nrow(group_data), size = rows_per_group[group_value], replace = TRUE)
  strat_ps1 <- rbind(strat_ps1, group_data[sampled_rows, ])
}

#save as shp for QGIS
#st_write(strat_ps1, "xxx.shp", append=TRUE) 

pseudoabsence <- st_drop_geometry(strat_ps1) #change data format

## combine all data 
df_total012 <- rbind(confltotal, pseudoabsence)
df_total012$confl <- as.factor(df_total012$confl)




##### EXPLORE CONFLICT DATA ###############

#Subset to relevant vars for LOGIT
df_plot <- dplyr::select(confltotal, 
                         SELS_mx, nrSELS, pavedkm, unpavedkm, PA_dist, 
                         vieh_sum, cattle_sum, livestock,
                         forest, shrub, grass, crop, woody, agri, 
                         set_area, set_num, pop_sum)

# historgrams 
par(mfrow = c(4, 4)) 
hist(df_plot$nrSELS, main="Number of SELS per bufer", xlab="count", ylab="Freq.")
hist(df_plot$pavedkm, main="Paved roads", xlab="km (sum)", ylab="Freq.")
hist(df_plot$unpavedkm, main="Unpaved roads", xlab="km (sum)", ylab="Freq.")
hist(df_plot$PA_dist, main="Distance to PA", xlab="km", ylab="Freq.")
#...


## Bivariate scatter plots
par(mfrow = c(1, 1)) 

ggplot(df_confl) +
  geom_point(aes(x = cattle_sum, y = vieh_sum))
p <- ggplot(df_confl, mapping = aes(x = cattle_sum, y = vieh_sum)) 
p+ geom_point()+
  geom_jitter() + # using geom_jitter to avoid overplotting of points 
  geom_smooth()
#...


## Correlation ##

#pearson matrix
par(mfrow = c(1, 1)) 
#version1 
sub_corr <- cor(df_plot[2:17], use = "na.or.complete") #pearson as default
#corrplot(sub_corr, method="shade",shade.col=NA, tl.col="black", tl.srt=45, tl.cex = 0.5)
corrplot(sub_corr, method="number", tl.col="black", number.cex = 0.80)#, main = "Pearson correlation matrix")

## correlation scatterplots
par(mfrow = c(2, 4)) 
plot(df_plot$pavedkm, df_plot$unpavedkm)
plot(df_plot$PA_dist, df_plot$set_area)
plot(df_plot$sml_livestock, df_plot$cattle_sum)
plot(df_plot$woody, df_plot$agri)
#...





##### PREPARE DATA FOR MODELING #####
#change variables to numeric/factor
df_total2 <- df_total012 %>%
  mutate(across(3:15, as.numeric))
df_total2$SELS_mx <- as.factor(df_total2$SELS_mx)

#calculate new variables
df_total2$woody <- df_total2$forest + df_total2$shrub
df_total2$agri <- df_total2$grass + df_total2$crop
df_total2$livestock <- df_total2$cattle_sum + df_total2$vieh_sum

#avoid zero values in PA_dist
df_total2$PA_dist <- replace(df_total2$PA_dist, df_total2$PA_dist == 0, 0.001)

#standardization of df 
#non-numeric columns in a separate df
non_numeric_columns <- df_total2[, c(1:2)]
#standardize only numeric columns
df_total_02 <- as.data.frame(scale(df_total2[3:18]))
#combine the standardized columns with the non-numeric columns
df_total_02 <- cbind(non_numeric_columns, df_total_02)


##### MODEL GLM ######
library(effectsize)
library(rcompanion)

glm.fit1 <- glm(confl ~ yourvariable, yourvariable, family= binomial, data = df_total_02)
glm.fit2 <- glm(confl ~ yourvariable, yourvariable, family= binomial, data = df_total_02)
glm.fit3 <- glm(confl ~ yourvariable, yourvariable, family= binomial, data = df_total_02)
glm.fit4 <- glm(confl ~ yourvariable, yourvariable, family= binomial, data = df_total_02)
glm.fit5 <- glm(confl ~ 1, data = df_total_02) #null model 

compareGLM(glm.fit1, glm.fit2s, glm.fit3, glm.fit4, glm.fit5) #check AIC, McFadden 

coef(glm.SELS.step1) #coefficients 



##### MULTINOMIAL LOGARITHMIC MODEL ######
library(nnet)
mlogit.final <- multinom(confl ~ PA_dist + unpavedkm + sml_livestock + set_num + woody, data = df_total_02) #
mlogit.null <- multinom(confl~ 1, data = df_total_02) #nullm-model 

AIC(mlogit.final) 
AIC(mlogit.null)
summary(mlogit.final) 
confint(mlogit.final, level = 0.90) #set conf. level to 90 (instead of default 95%)

#p-value calculation 
coefs.final <- coef(mlogit.final)
std_errors.final <- summary(mlogit.final)$standard.errors
z_values <- coefs.final / std_errors.final
p_values <- 2 * (1 - pnorm(abs(z_values)))
results <- data.frame(P_Value = p_values)
print(results)

# coefficient plot 
library(ggplot2)
library(dplyr)
library(broom)
coefs <- tidy(mlogit.final)

#90 % conf. int.
ggplot(coefs, aes(x = estimate, y = term, color = factor(y.level, levels = c("1", "2")))) +
  geom_point(position = position_dodge(width = -1/2), size = 3) + 
  geom_errorbarh(aes(xmin = estimate - 1.645 * std.error, xmax = estimate + 1.645 * std.error, y = term), 
                 height = 0.1, position = position_dodge(width = -0.5)) + 
  labs(x = "Coefficient", y = "Predictor Variable", color = "Response", fill = c("1", "2")) + #, title="90% conf. int.") +
  scale_color_manual(values = c("violetred2", "royalblue4"), labels = c("1" = "depr.", "2" = "hunt.")) +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) + # Adjust size as needed
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))#+

##### TESTS (here only for multinom. model) ######

# McFadden
library(pscl)
pR2(mlogit.final)

# chi-square
library(lmtest)
lrtest(mlogit.final)

# VIF test (1=no correlation, >5= high correlation)
library(car)
vif(mlogit.final)

# AUC --> Monte Carlo cross validations
library(pROC)
n_iterations <- 1000
auc_values <- numeric(n_iterations)
# Perform Monte Carlo cross-validation
for (i in 1:n_iterations) {
  sample <- sample(c(TRUE, FALSE), nrow(df_total_02), replace = TRUE, prob = c(0.70, 0.30))
  train <- df_total_02[sample, ]
  test <- df_total_02[!sample, ]
  
  predicted_probs <- predict(mlogit.final, newdata = test, type = "probs")
  roc_curve <- multiclass.roc(test$confl, predicted_probs)
  auc_values[i] <- auc(roc_curve)
}
# Calculate the average AUC
average_auc <- mean(auc_values)
print(average_auc)









