# -----------------------------------------------------------------------------#

# Created by Laura Lisa Khaukha
#Please Note: Load all the packages before running the downstream code
#Some Sections are labelled with the corresponding report sectioning

# -----------------------------------------------------------------------------#
#               #Packages installed & corresponding libraries used analysis 
# -----------------------------------------------------------------------------#

# Installing packages & loading the corresponding libraries 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(haven, tidyverse, data.table, broom, lmtest, sandwich, survival, foreach, doSNOW)
library(tableone) # Loading tableone seperately
library(doParallel) #for bootstrap later 

#-----------------------------------------------------------------------------#
#Section 1: Data Preparation & Initial EDA (Exploratory Data Analysis)
#-----------------------------------------------------------------------------#
#Setting up working directory and preparing data
setwd('~/Desktop')
dt_s <- read_dta("")

#Exploring the data & its structure and variable names
view(dt_s) 
str(dt_s)
names(dt_s)

#Checking missing for values in the data
colSums(is.na(dt_s)) #Observation: A_pr:2000; Confirms structural missingness in adherence variable, a expected finding

#-----------------------------------------------------------------------------#
#Section 1.1:  Descriptive Statistics & Table 1 Creation
#-----------------------------------------------------------------------------#
#Converting data-set
setDT(dt_s)

# Initial investigation of the variables (descriptive stats): generating (descriptive stats) across all variables
summary(dt_s)

############## Descriptive statistics by satin use at baseline ##################
str(dt_s$t) 

# Obtaining summary statistics for age and baseline BMI, at time 0 (baseline) for individuals not on statin
dt_s %>% filter(t==0, A_bas==0) %>%
  summarise(across(c(L1_bas, L4_bas), list(mean=mean, sd=sd, min=min, max=max)))

#Obtaining the summary statistics for baseline(month 0) age & BMI, for people on statins
dt_s %>% filter(t==0, A_bas==1) %>%
  summarise(across(c(L1_bas, L4_bas), list(mean=mean, sd=sd, min=min, max=max)))

#Tabulating sex(Gender) at (month 0) for statin users
dt_s %>% filter(t==0, A_bas==1) %>%
  count(L2_bas)

#Tabulating ethnicity at month (0) for people on statin (statin users)
dt_s %>% filter(t==0, A_bas==1) %>%
  count(L3_bas)

#Tabulating baseline sex(gender) at time (month 0) for statin non-users
dt_s %>% filter(t==0, A_bas==0) %>%
  count(L2_bas)

#Tabulating ethnicity at time (month 0) for statin non-users
dt_s %>% filter(t==0, A_bas==0) %>%
  count(L3_bas)

#Creating Table1 for export in the report, stratified by statin use at baseline ("A_bas")
#specifying the variables to include in the table
bvars <- c("L1_bas", "L2_bas", "L3_bas", "L4_bas")  # Age, Sex, Ethnicity, BMI
dt_baseline <- dt_s[dt_s$t==0, ] #Sub-setting the data to baseline time (month 0)
table1 <- CreateTableOne(vars = bvars, strata="A_bas", data=dt_baseline, factorVars=c("L2_bas", "L3_bas"))

#Printing the resulting table1, preparing for export and saving the table for report display
print(table1, showAllLevels=TRUE, quote=FALSE, noSpaces=TRUE)
table1_d <- print(table1, showAllLevels =TRUE, printToggle=FALSE)
write.csv(table1_d, "Table2_Baseline_Characteristics.csv", row.names =TRUE)

#######################  Distribution Visualization #############################

#Histogram for continuous variable (Age) at baseline (month 0)
ggplot(dt_s[dt_s$t == 0, ], aes(x = L1_bas)) +
  geom_histogram(fill= "skyblue", color ="black", bins = 30) +
  labs(title = "Distribution of Baseline Age", x = "Age", y = "Frequencies") +
  theme_minimal()

#Histograms for continuous variable (BMI) at baseline (month 0)
ggplot(dt_s[dt_s$t == 0, ], aes(x = L4_bas)) +
  geom_histogram(fill= "skyblue", color="black", bins=30) +
  labs(title = "Distribution of Baseline BMI", x = "BMI(kg/m²)", y="Frequencies") +
  theme_minimal()

#Barplots for categorical variable, sex(gender)
ggplot(dt_s[dt_s$t == 0, ], aes(x =factor(L2_bas))) +
  geom_bar(fill= "skyblue", color="black") +
  labs(title = "Sex Distribution at Baseline", x="Sex (0 = Female, 1 = Male)", y="Counts") +
  theme_minimal()

#Barplots for categorical variable,  ethnicity
ggplot(dt_s[dt_s$t == 0, ], aes(x = factor(L3_bas))) +
  geom_bar(fill ="skyblue", color ="black") +
  labs(title = "Ethnicity Distribution at Baseline", x = "Ethnicity (0 = White, 1= Other)", y="Counts") +
  theme_minimal()

#Box-plots for BMI, stratified by statin group
ggplot(dt_s[dt_s$t == 0, ], aes(x = factor(A_bas), y = L4_bas)) +
  geom_boxplot(fill="skyblue",color="black") +
  labs(title = "BMI Distribution by statin Initiation Status", x = "Statin(0 = No, 1 = Yes)", y = "BMI (kg/m²)") +
  theme_minimal()

#Box-plots for Age, stratified via statin group
ggplot(dt_s[dt_s$t == 0, ], aes(x=factor(A_bas), y =L1_bas)) +
  geom_boxplot(fill="skyblue",color="black") +
  labs(title = "Age Distribution by statin Initiation Status", x = "Statin(0 = No, 1 = Yes)", y = "Age") +
  theme_minimal()

#-----------------------------------------------------------------------------#
# Section 4.0 (same as in report): Estimation of the ITT Effect  
# Employment of pooled logistic regression and predicted risk curves
#-----------------------------------------------------------------------------#

#Fitting a pooled logistic regression model to estimate the ITT, hereby adjusting for time and baseline covarians
#Baseline statin use is the exposure, while recurrent CVD is the outcome
pol_itt <- glm(Y_outcome ~ A_bas + t + I(t^2) + I(t^3) + L1_bas + L2_bas + L3_bas + L4_bas,
               data = dt_s,
               family = binomial(link = "logit"))

#estimating variance-covarinace matrix, clustering by individual
cov_r <- sandwich::vcovCL(pol_itt, cluster = ~id) 
lmtest::coeftest(pol_itt,vcov.=cov_r) # Displaying the co-efficients with robust SEs and 95% CIs, and p-values

#Computing the Odds Ratios (OR) and 95% CIs, a tidy version 
tidy(pol_itt, exponentiate = TRUE, conf.int = TRUE)

#Predicting CVD outcome probabilities calculated from ITT model above, for each time point
dt_s$pred_prob <- predict(pol_itt, newdata = dt_s, type = "response")

#Calculating the average predicted risk over time by baseline statin group, in aim to obtain the cumulative risk curves below
r_curve <- dt_s %>%
  group_by(t, A_bas) %>% summarise(mean_risk = mean(pred_prob, na.rm = TRUE), .groups = "drop")

#Plotting the cumulative risk curves
ggplot(r_curve, aes(x = t, y = mean_risk, colour = factor(A_bas))) +
  geom_line(size = 1.2) +
  labs(title = "Unadjusted Risk of Recurrent CVD Over Time by Statin Initiation",
       x = "Time (Months) since Baseline",
       y = "Predicted Cumulative Risk (%)",
       colour = "Statin Initiation (A_bas)") +
  scale_colour_manual(values = c("0" = "#F1948A", "1" = "#85C1E9"),
                      labels = c("0" = "No Statin", "1" = "Statin")) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

#Saving the risk curve
ggsave("ITTriskcurve.png", width = 9.0, height = 5.5, dpi = 320, bg = "white")

# Evaluating the potential rare outcome assumption, for the justification of the logistic model over hazards
# calculating monthly CVD rate for each time point
dt_s %>% group_by(t) %>%
  summarise(monthly_rate = mean(Y_outcome == 1),
            events = sum(Y_outcome == 1),
            total = n()
  ) %>%
  arrange(t)

#-----------------------------------------------------------------------------#
#      Section 5.0 (as seen in the report): Estimation of PPI effect
#-----------------------------------------------------------------------------#

# Reloading the data-set just in cases
dt_s <- read_dta("Assignment_2025.dta")
setDT(dt_s) #Converting the data-set into a table for efficiency

#-----------------------------------------------------------------------------#
# Section 5.1 (In the report): IPW & Adherence Logistic Regression Models
#-----------------------------------------------------------------------------#
# Filtering data to only including post-baseline (t>0), where individuals have not been censored 
dt_sp <- dt_s[t >0 & !is.na(A_pr)]

#Checking the adherence switching overtime, stratified by statin status (treatment)
table(dt_sp$A_pr, dt_sp$A)

#Fitting a Logistic regression, modelling and computing adherence probabilities those previously not on statins (A_pr == 0)
#Hereby, adjusting for baseline co-variants and time
l_model <- glm(A ~ L1_bas + L2_bas + L3_bas + L4_bas + t + I(t^2) + I(t^3),
               data = dt_sp[A_pr == 0], family = binomial())

# Predicting the adherence probabilities for those not on statin (A_pr == 0)
dt_s[, p_A0 := predict(l_model, newdata = dt_s, type = "response")]

#Fitting another logistic model for adherence probabilities for those previously on statins (A_pr == 1)
#Hereby, similarly adjusting for baseline co-variants and time
l_model_2 <- glm(A ~ L1_bas + L2_bas + L3_bas + L4 + t + I(t^2) + I(t^3),
                 data = dt_sp[A_pr == 1], family = binomial())

# Predicting the adherence probabilities for those on statin (A_pr == 1)
dt_s[, p_A1 := predict(l_model_2, newdata = dt_s, type = "response")]

#Displaying the logistic model summaries
summary(l_model)
summary(l_model_2)

#--------------- Generating adherence weigths (IPW)for each individual ------------------# 

# Assigning the probabilities of adherence, using available adherence observations
dt_s[A==0, p_A := 1-p_A0] #In the case of non adherence, assign 1- predicted adherence probabilty
dt_s[A==1, p_A := p_A1]  #In the case of adherence, keep the predicted probabilty

# Setting adherence to 1 (p_A = 1) at  baseline (month 0)
dt_s[t == 0, p_A := 1]

#Computing the resulting cumulative adherence probability for each individual
#This is essential for establishing the subsequent IPW
dt_s[, cump_A := cumprod(p_A), by=id]

# Calculating the unstabilised inverse probability weights for adherence, 
# using the  cumulative adherence probability from above
dt_s[, IPW_ := 1 /cump_A]

#Viewing the summary statistics of these predicted adherence probabilities and weights
summary(dt_s$p_A)
summary(dt_s$p_A0)
summary(dt_s$p_A1)
summary(dt_s$IPW_) 

#-------- Plotting the resulting unstabilised IPW Distribution --------#

#Preparing export for display in the report
png("unstabiliseIPW.png", width = 2200, height = 1600, res = 300)

#Plotting the unstabilisised inverse probability weight
hist(dt_s$IPW_, breaks = 30, col = "lightpink",
     main = "Distribution of Unstabilised Inverse Probability Weights",
     xlab = "Inverse Probability Weights (IPW)",
     ylab = "Frequency")

dev.off()

#----- Truncating unstabilised weights: Plotting Truncated IPW Distribution --------#

# Calculating the 99th percentile of the unstabilised IPW, to reduce the effect of extreme weights
th_99 <- quantile(dt_s$IPW_, 0.99, na.rm = TRUE)

#Defining the truncation threshold, capping weights at minimum of 25 or99th percentile
tc_limi <- min(th_99, 25)
print(paste("truncation threshold:", tc_limi)) #viewing the truncation threshold

#Applying this truncation to weights
dt_s[, IPWtr := pmin(IPW_, tc_limi)]

#Preparing plot (truncated IPW distribution) for export and display in the report
png("t_IPW.png", width = 2200, height = 1600, res = 300)

#Plotting the distribution of the truncated weights
hist(dt_s$IPWtr,
     main = "Truncation of adherence weights (min of 25 & 99th percentile)",
     xlab = "Inverse probability weights (IPW)", 
     ylab = "Frequency",
     col = "lightblue", breaks = 30)

dev.off()

#Viewing the summary of truncation weights for reporting 
summary(dt_s$IPWtr)

#-----------------------------------------------------------------------------#
# Section 5.2 (as seen in the report): Weighted pooled logistic regression for PP effect
#-----------------------------------------------------------------------------#

#Sub-setting/Filtering followup data where individuals are uncensored (cens_time == 0) with available weights
dt_s_pp <- dt_s %>% filter(cens_time == 0, !is.na(IPWtr))

#Fitting a weighted pooled logistic regression model for PP effect, adjusting for time and baseline covariants and applying weights
pp_log  <- glm(Y_outcome ~ A_bas + t + I(t^2) + I(t^3) + L1_bas + L2_bas + L3_bas + L4_bas,
               data = dt_s_pp,
               weights = IPWtr,
               family = binomial(link = "logit"))

#estimating robust standard errors, clustering by individual
cov_rob <- sandwich::vcovCL(pp_log, cluster = ~id)
lmtest::coeftest(pp_log, vcov. = cov_rob) # Displaying the co-efficients with robust SEs and 95% CIs, and p-values

# Extracting the odds ratio OR, 95% confidence interval and p-value, in tidy manner by exponentiating log-odds
pp_re <-tidy(pp_log, exponentiate = TRUE, conf.int = TRUE)

#Filtering the PP treatment effect of baseline statin use on recurrent CVD t for reporting
pp_re %>%
  filter(term =="A_bas") %>%
  select(term, estimate, conf.low, conf.high, p.value)
print(pp_re, n = Inf)

#-----------------------------------------------------------------------------#
# Section 5.3 (In the report): IPW & Adherence Logistic Regression Models
#-----------------------------------------------------------------------------#
# Fitting a discrete-time hazards model which interaction between time and treatment
# Weights (IPW) are applied to adjust for adherence 
ppplog_rc <- glm(Y_outcome ~ (t + I(t^2) + I(t^3)) * A_bas + L1_bas + L2_bas + L3_bas + L4,
                 data = dt_s[cens_time == 0],
                 weights = IPWtr,
                 family = binomial(link = "logit"))

#estimating robust standard errors, clustering by individual
covob <- sandwich::vcovCL(ppplog_rc, cluster = ~id)
lmtest::coeftest(ppplog_rc, vcov. = covob)

#-----------------------------------------------------------------------------#
# Creating a followup data-set (hypothetical) for both interventions (A_bas=0: 1), for each individual across the 30 months
rc_dt <- dt_s[t == 0] #Filtering baseline (month 0) data
rc_dt <- rc_dt[, .SD[rep(1, 30)], by = id]
rc_dt[, t := seq_len(.N) - 1L, by = id] #Creating a time variable
rc_dt <- rc_dt[rep(seq_len(.N), each = 2)] #Creating duplicates of each row for each treatment group
rc_dt[, intervention := rep(0:1, times = .N/2)]
rc_dt[, A_bas := intervention] #assigning these hypothetical treatments

#-----------------------------------------------------------------------------#
# Generating outcome probabilities using the discrete-time hazards model fitted with IPW

#Predicting CVD event probabilities at each time point,from the discrete-time hazards model
rc_dt[, p_event := predict(ppplog_rc, newdata = rc_dt, type = "response")]
rc_dt[, p_healthy_:= 1 - p_event]

#Computing the cumulative survival probability and cumulative risk over time
setorder(rc_dt, id, intervention, t)
rc_dt[, p_healthy := cumprod(p_healthy_), by = .(id, intervention)]
rc_dt[, prisk := (1 - p_healthy) * 100]
#-----------------------------------------------------------------------------#
#Summarising & Printing the risk estimates at the end of follow up (month 29)
rc_summa <- rc_dt[t == 29, .(mean_risk = mean(prisk)), by = intervention]
print(rc_summa)
#-----------------------------------------------------------------------------#
#Preparing the data for the plotting of the risk curves
#Computing the average risk for each group at each follow up 
rcpldata <- rc_dt[, .(mean_prisk = mean(prisk)), by = .(intervention, ct = t + 1L) ]

#Adding baseline t month 0 with 0 risk
rcpldata_base <- data.table(intervention = c(0,1), ct = 0, mean_prisk = 0)
rcpldata <- rbind(rcpldata, rcpldata_base)[order(intervention, ct)] #Merging both datasets

#Converting data to wide format for subsequent plotting
rc_plot <- dcast(rcpldata, ct ~ intervention, value.var = "mean_prisk")
setnames(rc_plot, c("0", "1"), c("mean_prisk0", "mean_prisk1"))

#Computing absolute risk differences
rc_plot[, mean_priskdiff := mean_prisk1 - mean_prisk0]

#-----------------------------------------------------------------------------#
#Saving the estimates for reporting
write.csv(rc_plot, "rc_pointest.csv", row.names = FALSE)

#Plotting the standardised PP l risk curves
ggplot(rc_plot, aes(x = ct)) +
  geom_line(aes(y = mean_prisk0, color = "Non Statin User:A_bas = 0")) +
  geom_line(aes(y = mean_prisk1, color = "Statin User: A_bas = 1")) +
  labs(title = "Standardised Per-Protocol (PP) Risk Curves",
       x = "Time Since Baseline (Months)",
       y = "Predicted Cumulative Risk (%)",
       color = "Treatment Group") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  theme_minimal()

#saving the rc into the local directory
ggsave("5_3_ppriskcurve.png", width = 8, height = 5.5, dpi = 300)


#-----------------------------------------------------------------------------#
#       Section 5.3 (Part 2): Bootstrap-based Risk Curves with 95% CI
#-----------------------------------------------------------------------------#
#Filtering data for the retention of only uncensored individuals for PP effect
dt_s <- dt_s[cens_time == 0]

# Calculating  adherence cumulative probabillites and applying weights and truncating weights
dt_s[, cum_p_A := cumprod(p_A), by = id]
dt_s[, IPW_ := 1 / cum_p_A]
dt_s[, IPWtr := pmin(IPW_, 25)]

#Fitting a pooled logistic regression with interaction terms between time and outcome and truncated weights
plogipw <- glm(Y_outcome ~ (t + I(t^2) + I(t^3)) * A_bas,
               data = dt_s, weights = IPWtr,
               family = binomial(link = "logit"))

#-----------------------------------------------------------------------------#
# Creating a followup data-set (hypothetical) for both interventions scenarios (A_bas=0: 1), for each individual across the 30 months
risk_dt <- dt_s[t == 0]
risk_dt <- risk_dt[, .SD[rep(1, 30)], by = id]  
risk_dt[, t := seq_len(.N) - 1L, by = id]
risk_dt <- risk_dt[rep(seq_len(.N), each = 2)] # Creating duplicate rows for both interventions
risk_dt[, intervention := rep(0:1, times = .N/2)]
risk_dt[, A_bas := intervention]

#Predicting CVD event probabilities at each time point,using the plogipw regression model 
risk_dt[, p_event_k := predict(plogipw, newdata = risk_dt, type = "response")]
risk_dt[, p_healthy_k := 1 - p_event_k]

#Computing the cumulative survival probability and cumulative risk
setorder(risk_dt, id, intervention, t)
risk_dt[, p_healthy := cumprod(p_healthy_k), by = .(id, intervention)]
risk_dt[, p_risk := (1 - p_healthy) * 100]

#Summarising mean cumulative risk at each time-points
mean_risk <- risk_dt[, .(mean_risk = mean(p_risk)), by = .(intervention, t)]

#Coverting the datset into wide format for subsequent plotting and calculating the absolute risks
risk_curve <- dcast(mean_risk, t ~ intervention, value.var = "mean_risk")
setnames(risk_curve, c("0", "1"), c("risk0", "risk1"))
risk_curve[, riskdiff := risk1 - risk0]

#Printing the results
print(risk_curve)
#saving the estimates, just in case 
write.csv(risk_curve, "Results/rcestimates.csv")
#-----------------------------------------------------------------------------#
#Plotting the estimates of the standardised PP risk curves, just in case
ggplot(risk_curve, aes(x = t)) +
  geom_line(aes(y = risk0, color = "No Statin"), size = 1.2) +
  geom_line(aes(y = risk1, color = "Statin"), size = 1.2) +
  ylab("Cumulative Risk (%)") +
  xlab("Follow-up Month") +
  labs(title = "Standardised Risk Curves under PP Adherence") +
  theme_minimal() +
  theme(legend.position = "bottom")

#------------------- Bootstrap of 95% CIs  risk curves----------------------------------------------#
dir.create("Results", showWarnings = FALSE) #ensuring the existence of the results directiory

#Setting up the parallel  cluster enviromnet
n_cors <- parallel::detectCores() - 1 
cl <- makeCluster(n_cors) 
registerDoSNOW(cl)

# Setting up the bootstrap function for the re-estimation of risk curves across populations
bootstrap_risk_curves <- function(data, R = 200, seed = 123) {
  set.seed(seed)
  ids <- unique(data$id)
  
  # Start of the parallel bootstrap loop
  boot_results <- rbindlist(
    foreach(i = 1:R, .packages = c("data.table", "broom")) %dopar% {
      sampled_ids <- sample(ids, length(ids), replace = TRUE)  #Resampling individuals based on their ids with replacement
      boot_data <- rbindlist(lapply(sampled_ids, function(x) data[id == x]))
      
      # Recalculating cumulative adherence probabilities and IP weights, with truncation for each bootstrap sample
      boot_data[, cum_p_A := cumprod(p_A), by = id]
      boot_data[, IPW_ := 1 / cum_p_A]
      boot_data[, IPWtr := pmin(IPW_, 25)]
      
      #Refitting a pooled logistic regression model for each bootstrap sample
      model <- glm(Y_outcome ~ (t + I(t^2) + I(t^3)) * A_bas,
                   data = boot_data,
                   weights = IPWtr,
                   family = binomial("logit"))
      
      #Creating a  dataset for hypothetical followups accross the 30 months under both interventions
      risk_dt <- boot_data[t == 0]
      risk_dt <- risk_dt[, .SD[rep(1, 30)], by = id]
      risk_dt[, t := seq_len(.N) - 1L, by = id]
      risk_dt <- risk_dt[rep(seq_len(.N), each = 2)]
      risk_dt[, intervention := rep(0:1, times = .N/2)]
      risk_dt[, A_bas := intervention]
      
      #duplicating/creating two rows per individuals for each month for both intervention arms
      risk_dt[, p_event_k := predict(model, newdata = risk_dt, type = "response")]
      risk_dt[, p_healthy_k := 1 - p_event_k]
      setorder(risk_dt, id, intervention, t)
      risk_dt[, p_healthy := cumprod(p_healthy_k), by = .(id, intervention)] # Predicting outcome risk using the fitted model
      risk_dt[, p_risk := (1 - p_healthy) * 100]
      
      # Computing average risk and risk differences under both interventions for each time-point
      curve_summary <- risk_dt[, .(
        risk0 = mean(p_risk[intervention == 0], na.rm = TRUE),
        risk1 = mean(p_risk[intervention == 1], na.rm = TRUE),
        riskdiff = mean(p_risk[intervention == 1], na.rm = TRUE) - 
          mean(p_risk[intervention == 0], na.rm = TRUE)
      ), by = t]
      
      # Creating time and bootstrap labels within the loop 
      curve_summary[, `:=`(ct = t, bootstrap_iter = i)]
      curve_summary
    }
  )
  
  # Aggregating boot results to get 95% CI
  boot_summary <- boot_results[, list(
    cilb_risk0 = quantile(risk0, 0.025, na.rm = TRUE),
    ciub_risk0 = quantile(risk0, 0.975, na.rm = TRUE),
    cilb_risk1 = quantile(risk1, 0.025, na.rm = TRUE),
    ciub_risk1 = quantile(risk1, 0.975, na.rm = TRUE),
    cilb_riskdiff = quantile(riskdiff, 0.025, na.rm = TRUE),
    ciub_riskdiff = quantile(riskdiff, 0.975, na.rm = TRUE)
  ), by = ct]
  
  #saving the results
  fwrite(boot_summary, "Results/boot_results.csv")
  return(boot_summary)
}

#reassuring no conflict with cluster above
registerDoParallel(cl)

#Running the bootstrap with 200 iterations
boot_summary <- bootstrap_risk_curves(dt_s, R = 200)
stopCluster(cl) # stopping the cluster 

#-----------------------------------------------------------------------------#
# Plotting risk curves with bootstrap-based 95% CIs
# for the standardised per-protocol risk curves.

#Loading the estimates and 95% bootstrap CIs from the respective files
rc_pointest <- read.csv("rc_pointest.csv")
rc_95ci <- read.csv("Results/boot_results.csv")
rc_allest <- left_join(rc_pointest, rc_95ci, by = "ct")

#Selecting relevant CI columns merging with point estimates for merging
rc_95ci <- rc_95ci %>% dplyr::select(ct,
                                     cilb_risk0, ciub_risk0,
                                     cilb_risk1, ciub_risk1,
                                     cilb_riskdiff, ciub_riskdiff)

# Merging the point estimates with corresponding 95% CIs
rc_allest <- left_join(rc_pointest, rc_95ci, by = "ct") %>%
  relocate(ct,
           mean_prisk0, cilb_risk0, ciub_risk0,
           mean_prisk1, cilb_risk1, ciub_risk1,
           mean_priskdiff, cilb_riskdiff, ciub_riskdiff)
#Saving the final,combined dataset for reporting and rreproducibility
write.csv(rc_allest, "Results/rc_bootfinal.csv", row.names = FALSE)

#-----------------------------------------------------------------------------#
# Plotting the Standardised Per-Protocol Risk Curves with 95% CIs

fig3a <-ggplot(rc_allest, aes(x = ct)) + 
  #Shading the CI bands for the visulisation
  geom_ribbon(aes(ymin = cilb_risk0, ymax = ciub_risk0, fill = "No Statin"), alpha = 0.2) +
  geom_ribbon(aes(ymin = cilb_risk1, ymax = ciub_risk1, fill = "Statin"), alpha = 0.2) +
  
  #Point estimated risk curves
  geom_line(aes(y = mean_prisk0, colour = "No Statin"), size = 1.5) + 
  geom_line(aes(y = mean_prisk1, colour = "Statin"), size = 1.5) +
  
  #Creating the plots labels 
  xlab("Follow-up Time (Months)") +
  ylab("Cumulative Risk (%)") +
  labs(
    title = "Standardised Risk Curves with 95% CI",
    colour = "Group",
    fill = "95% Confidence Interval"
  ) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  #Increasing the font sizes of the text for display in the report 
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(fig3a)

#Savig the RC for reporting 
ggsave("Figure3A_RC.png", plot = fig3a, width = 8, height = 6, dpi = 320, bg = "white")

#-----------------------------------------------------------------------------#
# Plotting the Risk Difference Between Groups with 95% CI
fig3b <-ggplot(rc_allest, aes(x = ct)) +
  geom_ribbon(aes(ymin = cilb_riskdiff, ymax = ciub_riskdiff), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean_priskdiff), color = "black", size = 1.2) +
  labs(title = "Risk Difference Between Groups (Statin vs No Statin)",
       y = "Risk Difference (%)", x = "Follow-up Time (Months)") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 13),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 13),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(fig3b)


# Saving the figure for reporting 
ggsave("Figure3Briskdff.png", plot = fig3b, width = 8, height = 6, dpi = 320, bg = "white")

