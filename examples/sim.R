library(tidyverse)
library(kableExtra)
library(metafor)

set.seed(123)

# simulating 
#y_{ij} = \beta_0 + \beta_1 \text{study}_{ij} + \beta_2 \text{treatment}_{ij}  + \beta_3(\text{study}_{ij} \cdot \text{treatment}_{ij}) + \epsilon \quad \epsilon_{ij} \sim N(0, v_{ij})

run_sim <- function(n1 = 50, n2 = 50, n_sims, beta0 = 0, beta1 = 0, beta2 = 0,
                    sigma = 1, seed = 137) {
  
  set.seed(seed)
  seeds <- sample.int(1e6, n_sims)
  n <- n1 + n2
  
  # (1) Raw data
  df_raw <- do.call(rbind, lapply(seq_len(n_sims), function(i) {
    set.seed(seeds[i])
    study     <- rep(0:1, c(n1, n2))
    treatment <- rep(rep(0:1, each = n1 / 2), 2)
    y         <- beta0 + beta1 * study + beta2 * treatment + rnorm(n, 0, sigma)
    data.frame(seed = seeds[i], study = study, treatment = treatment, y = y)
  }))
  
  # # (2) Summary grouped by seed and study
  s <- df_raw |> group_by(seed, study, treatment) |>
    summarise(n = n(),
              m = mean(y),
              sd = sd(y),
              .groups = "drop") |>
    pivot_wider(names_from = treatment,
                values_from = c(n, m, sd)) |>
    mutate(
      D = m_1 - m_0,
      s_p = sqrt(((n_0 - 1) * sd_0^2 + (n_1 - 1) * sd_1^2) /
                   (n_0 + n_1 - 2)),
      d = D / s_p
    )
  
  # (3) P-values: fit lm on df_raw per seed
  
  pvals <- do.call(rbind, lapply(seeds, function(s) {
    m  <- lm(y ~ study * treatment, data = df_raw[df_raw$seed == s, ])
    sm <- summary(m)$coefficients
    data.frame(
      seed          = s,
      p_treatment   = sm["treatment", "Pr(>|t|)"],
      p_interaction = sm["study:treatment", "Pr(>|t|)"]
    )
  }))
  
  list(raw_dat = df_raw, summary = s, p_values = pvals
  )
}


res <- run_sim(n_sims = 100)

s <- subset(res$raw_dat, seed == 65241)
fit <- lm(y ~ study * treatment, data = s)



example1 <- subset(res$summary,seed == 65241)

ex1 <- data.frame(X1 = example1$m_0,
                   X2 = example1$m_1, 
                   S1 = example1$sd_0,
                   S2 = example1$sd_1,
                   n1 = example1$n_0,
                   n2 = example1$n_1)



tab_summary <- example1|>
  select(study,n_0, m_0,sd_0, n_1, m_1, sd_1, D, d) |> 
  mutate(across(where(is.numeric), round, 3)) |>
  kable(
    align = "lcccccccc",
    col.names = c(
      "Context",
      "n", "Mean", "SD",
      "n", "Mean", "SD",
      "D", "d"
    )
  ) |>
  add_header_above(c(
    "     " = 1,
    "Control" = 3,
    "Treatment" = 3,
    " " = 2
  ))|>
  kable_styling(latex_options = "hold_position")


# for example 1 narrative 
ttest_1 <- t.test(y~treatment, var.equal= TRUE,
                  data = subset(s, subset = study == 0))

ttest_2  <- t.test(y~treatment, var.equal= TRUE,
                   data = subset(s, subset = study == 1))

# functions from Ellydee

apa_t <- function(t_test_obj) {
  sprintf("t(%d) = %.2f, p %s",
          t_test_obj$parameter,
          round(t_test_obj$statistic, 2),
          ifelse(t_test_obj$p.value < .001, "< .001",
                 paste0("= ", sub("^0", "", format(round(t_test_obj$p.value, 3))))))
}
apa_anova <- function(aov_obj, rownames = NULL) {
  if (!is.null(rownames)) {
    rownames(aov_obj) <- rownames
  }
  
  # Extract F-values from the rows that have them (not Residuals)
  terms <- rownames(aov_obj)
  terms <- terms[terms != "Residuals"]
  
  # Get residual df from the Residuals row
  residual_df <- aov_obj["Residuals", "Df"]
  
  results <- sapply(terms, function(term) {
    F_val  <- aov_obj[term, "F value"]
    df1    <- aov_obj[term, "Df"]
    p_val  <- aov_obj[term, "Pr(>F)"]
    
    p_str <- ifelse(p_val < .001, "< .001",
                    paste0("= ", sub("^0", "", format(round(p_val, 3)))))
    
    sprintf("F(%d, %d) = %.2f, p %s", df1, residual_df, round(F_val, 2), p_str)
  })
  
  names(results) <- terms
  return(results)
}


