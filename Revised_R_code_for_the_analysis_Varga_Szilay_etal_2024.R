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
library(performance)
library(MuMIn)
library(car)

# Overdispersion vizsgálat # funkció #
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}

############################################
# Set working directory to the location of this script
actdir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(actdir)

# Load data
tomato_df <- read.csv("Tomato_work_data_Varga_Szilay_etal_2024.csv", header=T,
                      sep=';', encoding='UTF-8')

tomato_df$Treatment <- dplyr::recode(tomato_df$Treatment, 
                              "DN" = "DN-open", 
                              "NN" = "NN-open", 
                              "SF" = "NN-bagged")

############ FOR VISUALIZATION ##############
custom_colors_pontok = c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
  "#8c564b", "#e377c2", "#7f7f7f", "chartreuse2", "cornflowerblue", 
  "#f0e442", "#da3e52", "skyblue", "#d1e455", "darkorchid", 
  "deeppink", "darkslateblue", "#00fa9a", "#8b0000", "#5f9ea0"
)

treatment_colours <- c("DN-open" = "#888888", "NN-open" = "#DDCC77", "NN-bagged" = "#88CCEE")

############################################
################ BITE MARK  ################
############################################

tomato_df <- tomato_df[tomato_df$Brown_patches != "", ]
tomato_df <- tomato_df[tomato_df$Fertilization != " - ", ]

table(tomato_df$Treatment)
# Filter data to exclude Self-fertilisation (SF) treatment, 
# as it is not relevant for brown patches
brown_patches_df <- tomato_df[tomato_df$Treatment != "NN-bagged", ]

# Calculate presence of brown patches for DN and NN treatments
brown_patches_DN <- brown_patches_df[brown_patches_df$Treatment == "DN-open", "Brown_patches"]
brown_patches_NN <-  brown_patches_df[brown_patches_df$Treatment =="NN-open", "Brown_patches"]

yes_DN <- sum(brown_patches_DN == "Yes", na.rm = TRUE)
yes_NN <- sum(brown_patches_NN == "Yes", na.rm = TRUE)

table(brown_patches_df$Treatment, brown_patches_df$Brown_patches)

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

# Analysis for binary data - bite marks
brown_patches_df$Brown_patch_binary <- ifelse(brown_patches_df$Brown_patches 
                                              == "Yes", 1, 0)

# model with RANDOM EFFECTS
#brown_patches_df$Total_flowers_on_cluster <- as.numeric(brown_patches_df$Total_flowers_on_cluster)
#brown_patches_df$Cluster_location <- as.numeric(brown_patches_df$Cluster_location)
#brown_patches_df$random_eff <- brown_patches_df$Total_flowers_on_cluster - brown_patches_df$Cluster_location

brown_patches_model <- glm(Brown_patch_binary ~ Treatment,
                            data = brown_patches_df,
                            family = binomial)

overdisp_fun(brown_patches_model) # 0.3828682 
# with random effect(s) 'isSingular'

summary(brown_patches_model)
coef(summary(brown_patches_model))
anova(brown_patches_model)
AIC(brown_patches_model)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = brown_patches_model, plot = F)
residuals(simulationOutput)
windows()
plot(simulationOutput)

################ Result ###################
# > summary(brown_patches_model)
# 
# Call:
#   glm(formula = Brown_patch_binary ~ Treatment, family = binomial, 
#       data = brown_patches_df)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)        0.9555     0.5262   1.816   0.0694 .
# TreatmentNN-open  -0.1288     0.6945  -0.186   0.8528  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 49.572  on 40  degrees of freedom
# Residual deviance: 49.537  on 39  degrees of freedom
# AIC: 53.537
# 
# Number of Fisher Scoring iterations: 4
# 
# > coef(summary(brown_patches_model))
# Estimate Std. Error    z value   Pr(>|z|)
# (Intercept)       0.9555114  0.5262344  1.8157524 0.06940836
# TreatmentNN-open -0.1288329  0.6944636 -0.1855142 0.85282572
# > anova(brown_patches_model)
# Analysis of Deviance Table
# 
# Model: binomial, link: logit
# 
# Response: Brown_patch_binary
# 
# Terms added sequentially (first to last)
# 
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                         40     49.572         
# Treatment  1 0.034518        39     49.537   0.8526
# > AIC(brown_patches_model) 
# [1] 53.53747

############################################
############# INITIAL FRUIT SET ############
############################################

berry_table <- table(tomato_df$Treatment, tomato_df$Fertilization)

# Chi-négyzet teszt
chisq.test(berry_table)
# > chisq.test(berry_table)
# Pearson's Chi-squared test
# data:  berry_table
# X-squared = 0.55589, df = 2, p-value = 0.7573

# Binomial GLM
fertil_df <- tomato_df[tomato_df$Fertilization != "" & tomato_df$Fertilization != " - ", ]

fertil_df$Fertilized_binary <- ifelse(fertil_df$Fertilization == "Yes", 1, 0)

# model with RANDOM EFFECTS
fertil_df$Total_flowers_on_cluster <- as.numeric(fertil_df$Total_flowers_on_cluster)
fertil_df$Cluster_location <- as.numeric(fertil_df$Cluster_location)
fertil_df$random_eff <- fertil_df$Total_flowers_on_cluster - fertil_df$Cluster_location

fertil_model_binomial <- glmer(Fertilized_binary ~ Treatment 
                               + (1 | fertil_df$random_eff),
                          data = fertil_df,
                          family = binomial)

fertil_model <- glm(Fertilized_binary ~ Treatment,
                            data = fertil_df,
                            family = binomial)

summary(fertil_model_binomial)
coef(summary(fertil_model_binomial))
anova(fertil_model_binomial)
AIC(fertil_model_binomial) # 57.82734
AIC(fertil_model) # 56.82006 

# with + (1 | fertil_df$Set_number) 'isSingular'
# with + (1 | fertil_df$Set_number) and (1 | fertil_df$random_eff) 'isSingular'
# but with  (1 | fertil_df$random_eff) OK

overdisp_fun(fertil_model)

simulationOutput <- simulateResiduals(fittedModel = fertil_model, plot = F)
residuals(simulationOutput)
windows()
plot(simulationOutput)

################ Result ###################
# > summary(fertil_model_quasi)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) [glmerMod
# ]
# Family: binomial  ( logit )
# Formula: Fertilized_binary ~ Treatment + (1 | fertil_df$random_eff)
# Data: fertil_df
# 
# AIC      BIC   logLik deviance df.resid 
# 57.8     66.3    -24.9     49.8       58 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.9580  0.2495  0.3210  0.4185  0.6410 
# 
# Random effects:
#   Groups               Name        Variance Std.Dev.
# fertil_df$random_eff (Intercept) 0.7259   0.852   
# Number of obs: 62, groups:  fertil_df$random_eff, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          2.1448     0.8652   2.479   0.0132 *
#   TreatmentNN-bagged  -0.6256     0.9955  -0.628   0.5297  
# TreatmentNN-open     0.1037     1.0513   0.099   0.9214  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmntNN-b
# TrtmntNN-bg -0.719           
# TrtmntNN-pn -0.633  0.566    

# > coef(summary(fertil_model_quasi))
# Estimate Std. Error     z value   Pr(>|z|)
# (Intercept)         2.1448208  0.8651889  2.47902029 0.01317438
# TreatmentNN-bagged -0.6256325  0.9955054 -0.62845722 0.52970446
# TreatmentNN-open    0.1037024  1.0512814  0.09864384 0.92142106
# > anova(fertil_model_quasi)
# Analysis of Variance Table
# npar  Sum Sq Mean Sq F value
# Treatment    2 0.74075 0.37038  0.3704
# > AIC(fertil_model_quasi)
# [1] 57.82734

############################################

# Load the df again 
tomato_df <- read.csv("Tomato_work_data_Varga_Szilay_etal_2024.csv", header=T,
                      sep=';', encoding='UTF-8')

tomato_df$Treatment <- dplyr::recode(tomato_df$Treatment,
                              "DN" = "DN-open",
                              "NN" = "NN-open",
                              "SF" = "NN-bagged")

tomato_df$Treatment <- factor(tomato_df$Treatment, 
                              levels = c("DN-open", "NN-open", "NN-bagged"))

# Filter data to remove rows with missing or irrelevant values in key columns
tomato_df <- tomato_df[tomato_df$Fruit_marketing_value != " - ", ]
tomato_df <- tomato_df[tomato_df$Seed_count != "na", ] 
tomato_df <- tomato_df[tomato_df$Marketing_value_score != "-", ] 

#table(tomato_df$Treatment, tomato_df$Set_number)

############################################
#### MARKET VALUE OF FRUIT ANALYSIS #####
############################################

# Convert Marketing_value_score to numeric for analysis
tomato_df$Marketing_value_score <- as.numeric(tomato_df$Marketing_value_score)

# linear regression (lm) 
mylm2 <- lm(Marketing_value_score ~ Treatment, data = tomato_df)
summary(mylm2) #  p-value: 0.5993

residuals <- resid(mylm2)

bptest(mylm2) # BP = 0.36302, df = 2, p-value = 0.834

# Shapiro-Wilk test
shapiro.test(residuals) # non normal distribution # W = 0.90439, p-value = 0.0008739
shapiro.test(tomato_df$Marketing_value_score)  # non normal distribution W = 0.85537, p-value = 3.029e-05
# We tried to transform the data 

# Generalised linear model with Poisson distribution due to non-normal distribution
mylm3 <- glm(Marketing_value_score ~ Treatment, 
             data = tomato_df, family = "poisson")
coef(summary(mylm3))
anova(mylm3)
AIC(mylm3) # [1] 165.4106

# mylm3_summary <- summary(mylm3)
# 1-(mylm3_summary$deviance/mylm3_summary$null.deviance)

# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                         47     33.857         
# Treatment  2   0.5305        45     33.327    0.767

# # model with RANDOM EFFECTS
# tomato_df$Total_flowers_on_cluster <- as.numeric(tomato_df$Total_flowers_on_cluster)
# tomato_df$Cluster_location <- as.numeric(tomato_df$Cluster_location)
# tomato_df$random_eff <- tomato_df$Total_flowers_on_cluster - tomato_df$Cluster_location

# market_val_glmer <- glmer(Marketing_value_score ~ Treatment + (1|tomato_df$Set_number) +
#                            (1|tomato_df$random_eff),
#                          data = tomato_df,  family = "poisson")
# 
# summary(market_val_glmer) # boundary (singular) fit: see help('isSingular')

simulationOutput <- simulateResiduals(fittedModel = mylm3, plot = F)
residuals(simulationOutput)
windows()
plot(simulationOutput) 
testDispersion(simulationOutput)

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

Figure_5 <-  ggplot(tomato_df, 
                   aes(x=factor(Treatment), 
                       y=Marketing_value_score, 
                       fill = Treatment), 
                   fill = "transparent") + 
  geom_boxplot(width=0.1, outlier.shape = NA) +
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", exact = FALSE, 
  #                    label = "p.signif", label.y = 4) +
  geom_pirate(points = FALSE) +
  geom_jitter(size = 3, width = 0.3, height = 0.01, alpha = 0.6) +
  scale_color_manual(values = custom_colors_pontok) +
  #scale_color_manual()+
  theme(rect = element_rect(fill = "transparent"))+
  scale_fill_manual(values = treatment_colours)+
  labs(x="Treatment", y = "Market value of the fruit")+
  theme_minimal()+
  theme(legend.position="top")+
  theme(text = element_text(size=20))

ggsave(file = "Figure_5.pdf", Figure_5, width = 10, height = 8)
# ggsave(file = "Figure_5.png", Figure_5, width = 10, height = 8)

############################################
########### SEED NUMBER ANALYSIS ###########
############################################

# Filtering the rows where the seed count is 0 because "non developed fruit"
tomato_df <- tomato_df[tomato_df$Fruit_marketing_value != "non developed fruit", ]

#table(tomato_df$Treatment, tomato_df$Set_number)

# Convert Seed_count to numeric for analysis
tomato_df$Seed_count <- as.numeric(tomato_df$Seed_count)

# Perform linear regression 
mylm <- lm(Seed_count ~ Treatment, data = tomato_df)
summary(mylm) # Model summary
anova(mylm) # ANOVA table for model
# Df Sum Sq Mean Sq F value Pr(>F)
# Treatment  2   8444  4222.2  2.5715 0.08869 .

# Check for heteroscedasticity
bptest(mylm) #BP = 4.2555, df = 2, p-value = 0.1191 #heteroscedasticity is not present

# Test normality of residuals and data distribution
residuals <- resid(mylm)
shapiro.test(residuals) # W = 0.96262, p-value = 0.163 # normal distribution
shapiro.test(tomato_df$Seed_count) # W = 0.98835, p-value = 0.9312 # normal distribution

# Generalised linear model with Poisson distribution due to non-normal distribution
mylm_seed <- glm(Seed_count ~ Treatment, 
                 data = tomato_df, family = "poisson")
coef(summary(mylm_seed))
anova(mylm_seed)
AIC(mylm_seed) # 1170.244

# Analysis of Deviance Table
# Model: poisson, link: log
# Response: Seed_count
# Terms added sequentially (first to last)
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         43     990.49              
# Treatment  2   92.877        41     897.61 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# model with RANDOM EFFECTS
tomato_df$Total_flowers_on_cluster <- as.numeric(tomato_df$Total_flowers_on_cluster)
tomato_df$Cluster_location <- as.numeric(tomato_df$Cluster_location)
tomato_df$random_eff <- tomato_df$Total_flowers_on_cluster - tomato_df$Cluster_location

mylm_seed_glmer <- glmer(Seed_count ~ Treatment + (1|tomato_df$Set_number) +
                  (1|tomato_df$random_eff),
                data = tomato_df,  family = "poisson")
summary(mylm_seed_glmer)

coef(summary(mylm_seed_glmer))
anova(mylm_seed_glmer)
AIC(mylm_seed_glmer) # 754.123

Anova(mylm_seed_glmer)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: Seed_count
#           Chisq Df Pr(>Chisq)    
# Treatment 23.091  2  9.678e-06 ***

r.squaredGLMM(mylm_seed_glmer)

############### Result ####################
# > summary(mylm_seed_glmer)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: Seed_count ~ Treatment + (1 | tomato_df$Set_number) + (1 | tomato_df$random_eff)
# Data: tomato_df
# 
# AIC      BIC   logLik deviance df.resid 
# 754.1    763.0   -372.1    744.1       39 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -7.9032 -0.8697  0.0202  1.6776  6.2580 
# 
# Random effects:
#   Groups               Name        Variance Std.Dev.
# tomato_df$Set_number (Intercept) 0.1148   0.3389  
# tomato_df$random_eff (Intercept) 0.0959   0.3097  
# Number of obs: 44, groups:  tomato_df$Set_number, 20; tomato_df$random_eff, 8
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         4.40611    0.13952  31.580  < 2e-16 ***
#   TreatmentNN-open    0.24614    0.05164   4.767 1.87e-06 ***
#   TreatmentNN-bagged  0.10019    0.05411   1.851   0.0641 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmntNN-p
# TrtmntNN-pn -0.225           
# TrtmntNN-bg -0.211  0.499    
# > coef(summary(mylm_seed_glmer))
# Estimate Std. Error   z value      Pr(>|z|)
# (Intercept)        4.4061145 0.13952234 31.579992 6.950917e-219
# TreatmentNN-open   0.2461413 0.05163973  4.766511  1.874439e-06
# TreatmentNN-bagged 0.1001878 0.05411226  1.851481  6.410035e-02
# > anova(mylm_seed_glmer)
# Analysis of Variance Table
# npar Sum Sq Mean Sq F value
# Treatment    2 23.172  11.586  11.586
# > AIC(mylm_seed_glmer) # 754.123
# [1] 754.123

simulationOutput <- simulateResiduals(fittedModel = mylm_seed_glmer, plot = F)
residuals(simulationOutput)
windows()
plot(simulationOutput) # OK # dispersion = 0.66918, p-value = 0.632
testDispersion(simulationOutput)

###################################

tomato_df_filtered <- tomato_df[, c("Treatment", "Seed_count")]

# Normality test 
by(tomato_df_filtered$Seed_count, tomato_df_filtered$Treatment, shapiro.test)

# Perform post-hoc pairwise comparisons without adjustment for multiple testing
emmeans_results <- emmeans(mylm_seed_glmer, pairwise~Treatment, adjust = "none") # For Table 2

# Extract the contrasts from the results
contrasts_df <- as.data.frame(emmeans_results$contrasts)

p_values <- contrasts_df$p.value
############### Result ####################
# $emmeans
# Treatment emmean    SE  df asymp.LCL asymp.UCL
# DN-open     4.41 0.140 Inf      4.13      4.68
# NN-open     4.65 0.137 Inf      4.38      4.92
# NN-bagged   4.51 0.139 Inf      4.23      4.78
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate     SE  df z.ratio p.value
# (DN-open) - (NN-open)     -0.246 0.0516 Inf  -4.767  <.0001
# (DN-open) - (NN-bagged)   -0.100 0.0541 Inf  -1.851  0.0641
# (NN-open) - (NN-bagged)    0.146 0.0530 Inf   2.756  0.0059
# 
# Results are given on the log (not the response) scale. 

############################################
############## VISUALIZATION ###############
# Perform the ggplot
Figure_6 <- ggplot(tomato_df, aes(x = factor(Treatment), 
                                  y = Seed_count, 
                                  fill = Treatment)) + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_pirate(points = FALSE) +
  geom_jitter(size = 4, width = 0.2, alpha = 0.6,
              aes(color = Set_number)) +
  scale_color_manual(values = custom_colors_pontok) +
  theme(rect = element_rect(fill = "transparent")) +
  scale_fill_manual(values = treatment_colours) +
  labs(x = "Treatment", y = "Number of seeds") +
  guides(color = "none") +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 20)) +
  geom_signif(comparisons = list(c("DN-open", "NN-open"), 
                                 c("DN-open", "NN-bagged"), 
                                 c("NN-open", "NN-bagged")),
              annotations = c(sprintf("< 0.001"),  #p = %.3f", p_values[1]
                              sprintf("p = %.3f", p_values[2]),  
                              sprintf("p = %.3f", p_values[3])), 
              y_position = c(160, 170, 180),
              tip_length = 0.03)

ggsave("Figure_6.pdf", Figure_6, width = 10, height = 8)
# ggsave(Figure_6, file = "Figure_6_ver2.png", width = 10, height = 8)
