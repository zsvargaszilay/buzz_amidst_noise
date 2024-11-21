######## PACKAGES #########

library(dplyr)       
library(ggplot2)     
library(ggpubr)      
library(lubridate)   
library(scales)      
library(ggpirate)    
library(lmtest)
library(lme4)        
library(lmerTest)    
library(emmeans)     
library(multcomp)    
library(gridExtra)   
library(grid)        
library(png)         
library(patchwork)   

###############################

# Set working directory to the location of this script
actdir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(actdir)

# Load data
tomato_df <- read.csv("Tomato_work_data_Varga_Szilay_etal_2024.csv", header=T,
                      sep=';', encoding='UTF-8')

# Define colour scheme for treatments
# DN = direct noise, NN = non-noise, SF = self-pollination
treatment_colours <- c("DN" = "#888888", "NN" = "#DDCC77", "SF" = "#88CCEE")

############################################
########## BROWN PATCHES ANALYSIS ##########
############################################

tomato_df <- tomato_df[tomato_df$Brown_patches != "", ]
tomato_df <- tomato_df[tomato_df$Fertilization != " - ", ]

# Filter data to exclude Self-fertilisation (SF) treatment, 
# as it is not relevant for brown patches
brown_patches_df <- tomato_df[tomato_df$Treatment != "SF", ]

# Calculate presence of brown patches for DN and NN treatments
brown_patches_DN <- brown_patches_df[brown_patches_df$Treatment == "DN", "Brown_patches"]
brown_patches_NN <-  brown_patches_df[brown_patches_df$Treatment =="NN", "Brown_patches"]

yes_DN <- sum(brown_patches_DN == "Yes", na.rm = TRUE)
yes_NN <- sum(brown_patches_NN == "Yes", na.rm = TRUE)

table(brown_patches_df$Treatment, brown_patches_df$Brown_patches)
# No Yes
# DN  5  13
# NN  7  16

# Test for differences in fertilisation between DN and NN treatments
prop.test(x = c(yes_DN, yes_NN), 
          n = c(length(brown_patches_DN), length(brown_patches_NN)))

############## Result ######################
#2-sample test for equality of proportions with continuity correction
#data:  c(yes_DN, yes_NN) out of c(length(brown_patches_DN), length(brown_patches_NN))
#X-squared = 7.5946e-31, df = 1, p-value = 1
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.2795999  0.3327400
#sample estimates:
#  prop 1    prop 2 
#0.7222222 0.6956522 
############################################

############################################
## FERTILISATION (BERRY PRESENCE) ANALYSIS #
############################################

berry_table <- table(tomato_df$Treatment, tomato_df$Fertilization)
# No Yes
# DN  2  16
# NN  3  20
# SF  4  17

# Chi-négyzet teszt
chisq.test(berry_table)
# > chisq.test(berry_table)
# Pearson's Chi-squared test
# data:  berry_table
# X-squared = 0.55589, df = 2, p-value = 0.7573

############## Result ######################
#2-sample test for equality of proportions with continuity correction
#data:  c(yes_bogyo_DN, yes_bogyo_NN) out of c(length(bogyo_DN), length(bogyo_NN))
#X-squared = 1.6031e-31, df = 1, p-value = 1
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.2000541  0.2387014
#sample estimates:
#  prop 1    prop 2 
#0.8888889 0.8695652 

#Warning message:
# In prop.test(x = c(yes_berry_DN, yes_berry_NN), n = c(length(berry_DN),  :
# Chi-squared approximation may be incorrect

############################################

############################################

# Load the df again 
tomato_df <- read.csv("Tomato_work_data_Varga_Szilay_etal_2024.csv", header=T,
                      sep=';', encoding='UTF-8')

# Filter data to remove rows with missing or irrelevant values in key columns
tomato_df <- tomato_df[tomato_df$Fruit_marketing_value != " - ", ]
tomato_df <- tomato_df[tomato_df$Seed_count != "na", ] # Reduces dataset n = 50
tomato_df <- tomato_df[tomato_df$Marketing_value_score != "-", ] # Excludes unviable plants, n = 48

############################################
#### MARKETING VALUE OF FRUIT ANALYSIS #####
############################################

# Convert Marketing_value_score to numeric for analysis
tomato_df$Marketing_value_score <- as.numeric(tomato_df$Marketing_value_score)

# linear regression (lm) 
mylm2 <- lm(Marketing_value_score ~ Treatment, data = tomato_df)
summary(mylm2) #  p-value: 0.5993

residuals <- resid(mylm2)

gqtest(mylm2, order.by = ~Treatment, 
       data = tomato_df, 
       fraction = 10) # fraction = discard roughly 20% of the total observations
#GQ = 0.56109, df1 = 16, df2 = 16, p-value = 0.8708

# Shapiro-Wilk test
shapiro.test(residuals) # non normal distribution # W = 0.90439, p-value = 0.0008739
shapiro.test(tomato_df$Marketing_value_score)  # non normal distribution W = 0.85537, p-value = 3.029e-05

# We tried to transform the data 

# Generalised linear model with Poisson distribution due to non-normal distribution
mylm3 <- glm(Marketing_value_score ~ Treatment, data = tomato_df, family = "poisson")
coef(summary(mylm3))
anova(mylm3)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                         47     33.857         
# Treatment  2   0.5305        45     33.327    0.767

# Estimated marginal means for treatment levels
emmeans(mylm3, pairwise~Treatment, adjust = "none") # For Table 1

################ Result ###################
# $emmeans
# Treatment emmean    SE  df asymp.LCL asymp.UCL
# DN         0.847 0.169 Inf     0.516      1.18
# NN         1.012 0.151 Inf     0.716      1.31
# SF         0.928 0.152 Inf     0.629      1.23
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate    SE  df z.ratio p.value
# DN - NN   -0.1643 0.226 Inf  -0.725  0.4682
# DN - SF   -0.0807 0.228 Inf  -0.354  0.7230
# NN - SF    0.0836 0.214 Inf   0.390  0.6966
# 
# Results are given on the log (not the response) scale. 
############################################

############################################
Figure_5 <- ggplot(tomato_df, 
                   aes(x=factor(Treatment), 
                       y=Marketing_value_score, 
                       fill = Treatment), 
                   fill = "transparent") + 
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  geom_pirate(points = FALSE) +
  geom_jitter(size = 3, width = 0.3, height = 0.01, alpha = 0.6) +
  theme(rect = element_rect(fill = "transparent"))+
  scale_fill_manual(values = treatment_colours)+
  labs(x="Treatment", y = "Marketing value of the fruit")+
  theme_minimal()+
  theme(legend.position="top")+
  theme(text = element_text(size=20))

# ggsave(file = "Figure_5.pdf", Figure_5, width = 10, height = 8)
# ggsave(file = "Figure_5.png", Figure_5, width = 10, height = 8)

############################################
########### SEED NUMBER ANALYSIS ###########
############################################

# Filtering the rows where the seed count is 0 because "non developed fruit"
tomato_df <- tomato_df[tomato_df$Fruit_marketing_value != "non developed fruit", ]

# Convert Seed_count to numeric for analysis
tomato_df$Seed_count <- as.numeric(tomato_df$Seed_count)

# Perform linear regression 
mylm <- lm(Seed_count ~ Treatment, data = tomato_df)
summary(mylm) # Model summary
anova(mylm) # ANOVA table for model
# Df Sum Sq Mean Sq F value Pr(>F)
# Treatment  2   8444  4222.2  2.5715 0.08869 .

# Check for heteroscedasticity using Goldfeld-Quandt test
gqtest(mylm, order.by = ~Treatment, 
       data = tomato_df, 
       fraction = 9) # fraction = discard roughly 20% of the total observations
# GQ = 0.23034, df1 = 15, df2 = 14, p-value = 0.9961 #heteroscedasticity is not present, 

# Test normality of residuals and data distribution
residuals <- resid(mylm)
# Test for normal distribution of residuals
shapiro.test(residuals) 
# W = 0.96262, p-value = 0.163 # normal distribution
# Test for normal distribution of seed numbers
shapiro.test(tomato_df$Seed_count) 
# W = 0.98835, p-value = 0.9312 # normal distribution

# Generalised linear model with Poisson distribution due to non-normal distribution
mylm_seed <- glm(Seed_count ~ Treatment, data = tomato_df, family = "poisson")
coef(summary(mylm_seed))
anova(mylm_seed)

# Analysis of Deviance Table
# Model: poisson, link: log
# Response: Seed_count
# Terms added sequentially (first to last)
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         43     990.49              
# Treatment  2   92.877        41     897.61 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tomato_df_filtered <- tomato_df[, c("Treatment", "Seed_count")]

# Normality test 
by(tomato_df_filtered$Seed_count, tomato_df_filtered$Treatment, shapiro.test)

# Perform post-hoc pairwise comparisons without adjustment for multiple testing
emmeans(mylm_seed, pairwise~Treatment, adjust = "none") # For Table 1

############### Result ####################
# $emmeans
# Treatment emmean     SE  df asymp.LCL asymp.UCL
# DN          4.32 0.0309 Inf      4.26      4.38
# NN          4.68 0.0249 Inf      4.63      4.73
# SF          4.43 0.0282 Inf      4.38      4.49
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast estimate     SE  df z.ratio p.value
# DN - NN    -0.363 0.0396 Inf  -9.171  <.0001
# DN - SF    -0.114 0.0418 Inf  -2.732  0.0063
# NN - SF     0.249 0.0376 Inf   6.636  <.0001

# Results are given on the log (not the response) scale. 

############################################
Figure_6 <- ggplot(tomato_df, aes(x=factor(Treatment), 
                                  y=Seed_count, 
                                  fill = Treatment), 
                   fill = "transparent") + 
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  geom_pirate(points = FALSE) +
  geom_jitter(size = 3, width = 0.2, alpha = 0.6) +
  theme(rect = element_rect(fill = "transparent"))+
  scale_fill_manual(values = treatment_colours) +
  labs(x="Treatment", y = "Number of seeds")+
  theme_minimal()+
  theme(legend.position="top")+
  theme(text = element_text(size=20))

# ggsave("Figure_6_ver2.pdf", Figure_6, width = 10, height = 8)
# ggsave(Figure_6, file = "Figure_6_ver2.png", width = 10, height = 8)
############################################

# Try to use factors as random effects in the model
# tomato_df$Total_flowers_on_cluster <- as.numeric(tomato_df$Total_flowers_on_cluster)
# tomato_df$Cluster_location <- as.numeric(tomato_df$Cluster_location)
# 
# tomato_df$random_eff <- tomato_df$Total_flowers_on_cluster - tomato_df$Cluster_location
# 
# modell <- lmer(Seed_count ~ Treatment + (1|tomato_df$random_eff), data = tomato_df)
# 
# summary(modell)
# anova(modell)

