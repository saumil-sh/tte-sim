library("MASS")
library("flexsurv")
library("rstpm2")
library("coxme")
library("simsurv")
library("ggplot2")
library("ggpubr")
#%%
# Load and examine myeloid data
data(myeloid, package = "survival")
head(myeloid)
myeloid$trt[myeloid$trt == "A"] <- 0
myeloid$trt[myeloid$trt == "B"] <- 1
myeloid$trt <- as.numeric(myeloid$trt)

fit <- survfit(Surv(futime, death) ~ trt, data = myeloid)

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
#%%

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
cov <- data.frame(id = 1:650, trt = rbinom(650, 1, 0.5))
  # Simulate event times using simsurv
dat_weib <- simsurv(
  betas = mod_weib$coefficients,
  x = cov,
  knots = mod_weib$knots,
  logcumhazard = logcumhaz_weib,
  maxt = 2500,
  interval = c(1E-8, 2501)
)

dat_flex <- simsurv(
  betas = mod_flex$coefficients,
  x = cov,
  knots = mod_flex$knots,
  logcumhazard = logcumhaz_flex,
  maxt = 2500,
  interval = c(1E-8, 2501)
)

# Merge covariate data and event times
dat_weib <- merge(cov, dat_weib)
dat_flex <- merge(cov, dat_flex)

#%%

fit_weib <- survfit(Surv(eventtime, status) ~ trt, data = dat_weib)
fit_flex <- survfit(Surv(eventtime, status) ~ trt, data = dat_flex)
ggsurvplot(
  fit_weib,
  data = dat_weib,
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
  pval = TRUE,              # Add p-value
  risk.table = FALSE,        # Add risk table
  linetype = c("solid", "dashed","solid", "dashed","solid", "dashed"),
  # palette = c("#f93414", "#06889b", "#fe6c31", "#006b8f", "#feae6c", "#011f51"),
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
  cov <- data.frame(id = 1:1000, trt = rbinom(1000, 1, 0.5))
  # Simulate event times using simsurv
  dat <- simsurv(
    betas = true_mod$coefficients,
    x = cov,
    knots = true_mod$knots,
    logcumhazard = logcumhaz,
    maxt = 2500,
    interval = c(1E-8, 2501)
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

# Set seed for reproducibility
set.seed(543543)

# Run 100 replicates in simulation study and report the
# estimated mean bias under the Weibull and FPM models
out <- rowMeans(replicate(20, sim_run(true_mod = true_mod)))
