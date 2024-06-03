#%% import libraries
library("parallel")
library("survival")
library("rstpm2")
library("coxme")
library("flexsurv")
library("simsurv")
library("ggplot2")
library("ggpubr")
library("survminer")
library("gsDesign")

#%% Load and process myeloid data
head(myeloid)
myeloid$trt[myeloid$trt == "A"] <- 0
myeloid$trt[myeloid$trt == "B"] <- 1
myeloid$trt <- as.numeric(myeloid$trt)

# Set seed for reproducibility
set.seed(4020)

#%% plot KM overall survival
fit <- survfit(Surv(futime, death) ~ trt, data = myeloid)
fit_cox <- coxph(Surv(futime, death) ~ trt, data = myeloid)

ggsurvplot(
  fit,
  data = myeloid,
  size = 1,                 # change line size
  conf.int = TRUE,          # Add confidence interval
  censor.shape="|",
  censor.size = 4,
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  legend.labs =
    c("A", "B"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  xlab = "Time [days]",
  ylab = "Overall survival"
)
#%%
# Fit a Weibull model to the data
mod_weib <- flexsurvspline(Surv(futime, death) ~ trt, data = myeloid, k = 0)

# Fit a FPM model to the data
mod_flex <- flexsurvspline(Surv(futime, death) ~ trt, data = myeloid, k = 2)

# Plot the survival curve for each model, overlaid with Kaplan-Meier
par(mfrow = c(1, 2), cex = 1.0) # graphics parameters
plot(mod_weib,
  main = "Weibull model",
  ylab = "Survival probability",
  xlab = "Time"
)
plot(mod_flex,
  main = "Flexible parametric model",
  ylab = "Survival probability",
  xlab = "Time"
)
#%% sample size est
trt1_rmst <- summary(fit)$table[9]
trt2_rmst <- summary(fit)$table[10]

lambda1 <- 1 / trt1_rmst
lambda2 <- 1 / trt2_rmst

eta <- 1 / summary(survfit(Surv(futime, !death) ~ 1, data = myeloid[myeloid$death == 1, ]))$table[5]

ratio <- 1.0

sided <- 2
hr <- exp(fit_cox$coefficients)

SSE <- nSurvival(
  lambda1 = lambda1,
  lambda2 = lambda2,
  Ts = 2500,
  Tr = 90,
  eta = eta,
  ratio = ratio,
  alpha = 0.05,
  beta = 0.1,
  approx = TRUE
)

#%% custom log cumulative harazard

# Define a function that returns the log cumulative hazard for an FPM model
logcumhaz_weib <- function(t, x, betas, knots) {
  basis <- flexsurv::basis(knots, log(t))
  res <-
    betas[["gamma0"]] * basis[[1]] +
    betas[["gamma1"]] * basis[[2]] +
    betas[["trt"]] * x[["trt"]]
  res
}

logcumhaz_flex <- function(t, x, betas, knots) {
  basis <- flexsurv::basis(knots, log(t))
  res <-
    betas[["gamma0"]] * basis[[1]] +
    betas[["gamma1"]] * basis[[2]] +
    betas[["gamma2"]] * basis[[3]] +
    betas[["gamma3"]] * basis[[4]] +
    betas[["trt"]] * x[["trt"]]
  res
}

#%%
cov <- data.frame(id = 1:ceiling(SSE$n), trt = rbinom(ceiling(SSE$n), 1, 0.5))
# Simulate event times using simsurv
dat_weib <- simsurv(
  betas = mod_weib$coefficients,
  x = cov,
  knots = mod_weib$knots,
  logcumhazard = logcumhaz_weib,
  maxt = 3000,
  interval = c(1E-8, 3001)
)

dat_flex <- simsurv(
  betas = mod_flex$coefficients,
  x = cov,
  knots = mod_flex$knots,
  logcumhazard = logcumhaz_flex,
  maxt = 3000,
  interval = c(1E-8, 3001)
)

# Merge covariate data and event times
dat_weib <- merge(cov, dat_weib)
# dat_weib$eventtime <- as.numeric(dat_weib$eventtime)
# dat_weib$status <- as.numeric(dat_weib$status)
# dat_weib$trt <- as.numeric(dat_weib$trt)
dat_flex <- merge(cov, dat_flex)
# dat_flex$eventtime <- as.numeric(dat_flex$eventtime)
# dat_flex$status <- as.numeric(dat_flex$status)
# dat_flex$trt <- as.numeric(dat_flex$trt)

#%%

fit_weib <- survfit(Surv(time = eventtime, event = status) ~ trt, data = dat_weib)
fit_flex <- survfit(Surv(time = eventtime, event = status) ~ trt, data = dat_flex)
ggsurvplot(
  fit_weib,
  data = dat_weib,
  size = 1,                 # change line size
  conf.int = TRUE,          # Add confidence interval
  censor.shape = "|",
  censor.size = 4,
  mark = 1,
  mark.time = dat_weib$eventtime[dat_weib$status == 1],
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  legend.labs =
    c("A", "B"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  xlab = "Time [days]",
  ylab = "Overall survival"
)
ggsurvplot(
  fit_flex,
  data = dat_flex,
  size = 1,                 # change line size
  conf.int = TRUE,          # Add confidence interval
  censor.shape="|",
  censor.size = 4,
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  legend.labs =
    c("A", "B"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  xlab = "Time [days]",
  ylab = "Overall survival"
)
ggsurvplot_combine(
  list(data = fit, weib = fit_weib, flex = fit_flex),
  data = NULL,
  size = 0.75,                 # change line size
  conf.int = FALSE,          # Add confidence interval
  censor.shape="|",
  censor.size = 0,
  pval = FALSE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  linetype = c("solid", "dashed", "solid", "dashed","solid", "dashed"),
  # palette = c("#F8766D", "#F8766D", "#00BFC4", "#00BFC4", "#00BA38", "#00BA38"),
  # legend.labs =
  #   c("A", "B"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw(),      # Change ggplot2 theme
  xlab = "Time [days]",
  ylab = "Overall survival"
)
#%%
# Fit the FPM model to the data, to obtain the values that will be
# used as 'true' parameter values when simulating the event times
true_mod <- mod_flex

# Define a function for conducting the simulation study
sim_run <- function(true_mod) {
  # Covariate data
  cov <- data.frame(id = 1:ceiling(SSE$n), trt = rbinom(ceiling(SSE$n), 1, 0.5))
  # Simulate event times using simsurv
  dat <- simsurv(
    betas = true_mod$coefficients,
    x = cov,
    knots = true_mod$knots,
    logcumhazard = logcumhaz_flex,
    maxt = 3000,
    interval = c(1E-8, 3001)
  )

  # Merge covariate data and event times
  dat <- merge(cov, dat)

  # Fit the two analysis models to the simulated dataset
  weib_mod <- flexsurvspline(Surv(eventtime, status) ~ trt,
    data = dat, k = 0
  )

  flex_mod <- flexsurvspline(Surv(eventtime, status) ~ trt,
    data = dat, k = 2
  )

  # Return the log HR for each model
  true_loghr <- true_mod$coefficients[["trt"]]
  weib_loghr <- weib_mod$coefficients[["trt"]]
  flex_loghr <- flex_mod$coefficients[["trt"]]

  # Return the bias in the log HR under each model
  c(
    true_loghr = true_loghr,
    weib_bias = weib_loghr - true_loghr,
    flex_bias = flex_loghr - true_loghr
  )
}

# Run 1000 replicates in simulation study and report the
# estimated mean bias under the Weibull and FPM models
out <- rowMeans(replicate(200, sim_run(true_mod = true_mod)))
#  true_loghr    weib_bias    flex_bias
#-0.345150476 -0.027285289 -0.003014806