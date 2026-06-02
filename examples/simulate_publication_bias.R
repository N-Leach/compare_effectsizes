library(tidyverse)





simulate_publication_bias <- function(
    target_published = 100,
    n_options = c(20, 50, 100, 200),
    alpha = 0.05,
    beta0 = 0,
    beta2 = 0,
    sig = 1,
    seed = 1234,
    max_iter = 10000
) {
  
  set.seed(seed)
  
  raw_data <- list()
  study_summaries <- list()
  
  n_published <- 0
  study_id <- 0
  
  while (n_published < target_published) {
    
    study_id <- study_id + 1
    
    ## safety stop
    if (study_id > max_iter) {
      warning("Maximum iterations reached.")
      break
    }
    
    n <- sample(n_options, 1)
    treatment <- rbinom(n, 1, 0.5)
    
    ## ground truth model 
    y <- beta0 + beta2 * treatment +
      rnorm(n, 0, sig)
    
    # fit on the simulated data 
    fit <- lm(y ~ treatment)
    
    coefs <- summary(fit)$coefficients
    
    est <- coefs["treatment", "Estimate"]
    se  <- coefs["treatment", "Std. Error"]
    p   <- coefs["treatment", "Pr(>|t|)"]
    
    ## publication rule
    if (p < alpha) {
      n_published <- n_published + 1
    }
    
    ## store raw data (study_id, y, treatment)
    raw_data[[study_id]] <- data.frame(
      study_id = study_id,
      y = y,
      treatment = treatment,
      published = p < alpha
    )
    
    ## store study-level summary
    study_summaries[[study_id]] <- data.frame(
      study_id = study_id,
      n = n,
      estimate = est,
      se = se,
      p = p
    )
  }
  
  list(
    raw    = do.call(rbind, raw_data),
    studies = do.call(rbind, study_summaries)
  )
}

t<- t.test(y ~ treatment, res$raw, var.equal = FALSE)
t$p.value
t$

res <- simulate_publication_bias(target_published = 100)
set.seed(123)

sig_n20<- subset(res$studies, subset = n == 20 & p<.05)
n_200 <- subset(res$studies, subset = n ==200)

ind1 <- sample(x = seq(1:nrow(sig_n20)), size = 1)
ind2 <- sample(x = seq(1:nrow(n_200)), size = 1)


two_studies <- rbind(sig_n20[ind1, ], n_200[ind2, ])

two_studies_raw_data <- subset(res$raw, subset = study_id== two_studies$study_id[1]|study_id== two_studies$study_id[2])

two_studies_raw_data$study_id <- as.factor(two_studies_raw_data$study_id)


dat <- two_studies_raw_data |>
  group_by(study_id, treatment) |>
  summarise(n = n(), m = mean(y), sd = sd(y), se = sd / sqrt(n),
            .groups = "drop" ) |>
  pivot_wider(
    names_from = treatment,
    values_from = c(n, m, sd, se),
    names_prefix = "")|>
  mutate(
    D = m_1 - m_0,
    s_p = sqrt(((n_0 - 1) * sd_0^2 + 
                  (n_1 - 1) * sd_1^2) / 
                 (n_0 + n_1 - 2))
  ) |>
  mutate(
    d = D/ s_p
  )



mod <- lm(y ~ study_id * treatment, data = two_studies_raw_data)

library(car)
Anova(mod, type = 2)

library(lme4)

mod_mixed <- lmer(y ~ treatment + (1 | study_id), data = two_studies_raw_data)

# error: boundary (singular) fit: see help('isSingular')
## because there are only two studies there is not enough information to estimate a variance 
## random intercept variance is being estimated at (or near) zero because two groups 
## can’t support a meaningful variance component.


library(metafor)
summary

raw_es <- escalc(measure = "MD",
                 n1i = n_0,
                 n2i = n_1,
                 m1i = m_0,
                 m2i = m_1, 
                 sd1i = sd_0,
                 sd2i = sd_1, 
                 data = dat)

fe.fit <- rma(yi = yi, 
              vi = vi, 
              method = "FE",
              data = raw_es)

fe.fit

smd_es <- escalc(measure = "SMD",
                 correct = FALSE, # {#eq-B_4.20} otherwise Hedges'g correction applied
                 n1i = n_0,
                 n2i = n_1,
                 di = d, 
                 data = dat)
fma.smd_es <- rma(yi, vi, data = smd_es, method = "FE")
fma.smd_es


