library(dplyr)
library(survival)
library(lmtest)

# CardioAG_corr: Age-bias corrected cardiovascular age gap (CardioAG)

# -----------------------------------------------
# Lifestyle Association Analysis for Non-CVD Cohort
# -----------------------------------------------

# Load dataset
df_lifestyle_noCVD <- read.csv(
  "/path_to_the_file/UKB_Non_CVD_Chohort_Lifestyles_Imputed.csv",
  check.names = FALSE
)

# Convert Sex to factor with labels
df_lifestyle_noCVD$Sex <- factor(
  df_lifestyle_noCVD$Sex,
  levels = c(0, 1),
  labels = c("_Female", "_Male")
)

# Example: Average weekly alcohol intake
var <- "Avg_weekly_alcohol_intake"

# Create a copy for categorical encoding
df_temp <- df_lifestyle_noCVD

# -------------------------------
# Recode Average weekly alcohol intake as categorical variable
# - < 6 drinks/week -> 0
# - 6–9 drinks/week -> 1
# - >= 10 drinks/week -> 2
# -------------------------------
df_temp <- df_temp %>%
  mutate(
    Avg_weekly_alcohol_intake = case_when(
      Avg_weekly_alcohol_intake < 6 ~ 0,
      Avg_weekly_alcohol_intake >= 10 ~ 2,
      TRUE ~ 1
    )
  )

# Convert to factor for regression
factor_list <- c(0, 1, 2)
df_temp[[var]] <- factor(df_temp[[var]], levels = factor_list)

# -----------------------------------------------
# Fit linear regression using categorical alcohol intake
# Dependent variable: CardioAG_corr, corrected cardiovascular age gap
# Covariates: age_imaging_derived, Sex
# -----------------------------------------------
formula_str <- as.formula(
  paste("CardioAG_corr ~", var, "+ age_imaging_derived + Sex")
)

linear_model_cat <- lm(formula_str, data = df_temp)

# -----------------------------------------------
# Optionally: Fit linear regression using alcohol intake as continuous variable
# -----------------------------------------------
linear_model_cont <- lm(formula_str, data = df_lifestyle_noCVD)

# -----------------------------------------------
# Extract regression coefficients, confidence intervals, and standard errors
# Example: categorical variable (factor)
# -----------------------------------------------
unique_len <- length(factor_list)

p_values <- summary(linear_model_cat)$coefficients[2:unique_len, 4]
estimates <- summary(linear_model_cat)$coefficients[2:unique_len, 1]
lowers <- confint(linear_model_cat)[2:unique_len, 1]
uppers <- confint(linear_model_cat)[2:unique_len, 2]
betas <- summary(linear_model_cat)$coefficients[2:unique_len, "Std. Error"]




# -----------------------------------------------
# Survival Analysis for Non-CVD Cohort, composite major adverse cardiovascular event (MACE) as example
# -----------------------------------------------

# Load dataset for survival analysis
df_cox_all_CVD <- read.csv("/path_to_the_file/UKB_Non_CVD_survival_Analysis.csv")

# Convert categorical covariates to factors with labels
df_cox_all_CVD$Diabetes <- factor(df_cox_all_CVD$Diabetes, levels = c(0,1), labels = c("_No", "_Yes"))
df_cox_all_CVD$Sex <- factor(df_cox_all_CVD$Sex, levels = c(0,1), labels = c("_Female", "_Male"))

# -----------------------------------------------
# Model 1: Unadjusted Cox proportional hazards model
# Outcome: Time to composite MACE
# Predictor: CardioAG_corr
# -----------------------------------------------
f1 <- coxph(Surv(Any_CVD_time, Any_CVD) ~ CardioAG_corr, data = df_cox_all_CVD)

# Test proportional hazards assumption
test.ph <- cox.zph(f1)
print(test.ph)  

# Summarize model
summary(f1)

# Compute confidence intervals (HR)
exp(confint(f1))

# Extract hazard ratio (HR), 95% CI, and p-value for CardioAG_corr
HR <- round(summary(f1)$coefficients["CardioAG_corr", "exp(coef)"], 2)  # HR
lower <- round(exp(confint(f1)["CardioAG_corr", 1]), 2)                  # Lower CI
upper <- round(exp(confint(f1)["CardioAG_corr", 2]), 2)                  # Upper CI
HR_CI <- paste0(HR, " (", lower, "-", upper, ")")
p_value <- round(summary(f1)$coefficients["CardioAG_corr", "Pr(>|z|)"], 5)

# -----------------------------------------------
# Model 2: Adjusted Cox proportional hazards model
# Adjusted for chronological age, sex, BMI, and diabetes
# -----------------------------------------------

f2 <- coxph(
  Surv(Any_CVD_time, Any_CVD) ~ CardioAG_corr + age_imaging_derived + Sex + BMI + Diabetes,
  data = df_cox_all_CVD
)

# Test proportional hazards assumption for adjusted model
test.ph <- cox.zph(f2)
print(test.ph)

# Summarize model
summary(f2)
exp(confint(f2))

# Extract HR, 95% CI, and p-value for CardioAG_corr in adjusted model
HR <- round(summary(f2)$coefficients["CardioAG_corr", "exp(coef)"], 2)
lower <- round(exp(confint(f2)["CardioAG_corr", 1]), 2)
upper <- round(exp(confint(f2)["CardioAG_corr", 2]), 2)
HR_CI <- paste0(HR, " (", lower, "-", upper, ")")
p_value <- round(summary(f2)$coefficients["CardioAG_corr", "Pr(>|z|)"], 5)


# -----------------------------------------------
# Likelihood Ratio Test for Cox Models
# -----------------------------------------------

# ---------------------------
# Base Cox proportional hazards model
# Covariates: chronological age, sex, BMI, diabetes
# Outcome: Time to composite MACE
# ---------------------------
cox_base <- coxph(
  Surv(Any_CVD_time, Any_CVD) ~ age_imaging_derived + Sex + BMI + Diabetes,
  data = df_cox_all_CVD
)

# ---------------------------
# Extended Cox model with CardioAG_corr
# ---------------------------
cox_cardioAG <- coxph(
  Surv(Any_CVD_time, Any_CVD) ~ age_imaging_derived + Sex + BMI + Diabetes + CardioAG_corr,
  data = df_cox_all_CVD
)

# ---------------------------
# Perform likelihood ratio test to compare nested models
# ---------------------------
likelihood_ratio_test <- lrtest(cox_base, cox_cardioAG)

# Print LRT results
print(likelihood_ratio_test)
