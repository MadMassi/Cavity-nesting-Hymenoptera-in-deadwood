# clearing workspace
rm(list = ls())
graphics.off()
cat("\014") 

# ------------------------------------------------------------------------------
## Header ##
# Running title: Moisture reduces saproxylic bee, wasp, and parasitoid diversity in lying and standing deadwood
# Author: Massimo Martini
# Date: 10th April 2026

# ------------------------------------------------------------------------------

## Comments ##
# For simplicity, I leave out all diagnostic and plotting procedures

# ------------------------------------------------------------------------------
# 1. Load libraries and set working directory
# ------------------------------------------------------------------------------

# Set working directory
# setwd() 

library(vegan)
library(janitor)
library(dplyr)
library(tidyverse)
library(glmmTMB)
library(performance)
library(DHARMa)
library(emmeans)
library(car)
library(ggeffects)
library(iNEXT)
library(betapart)
library(ecodist)
library(lme4)
library(mediation)
library(piecewiseSEM)
library(indicspecies)
library(broom.mixed)

source("deadwd_functions.R")

# ------------------------------------------------------------------------------
# 2. Read and prepare data for analysis
# ------------------------------------------------------------------------------
log_data <- readr::read_csv("log_data.csv")
insect_data <- readr::read_csv("deadwd_insect_data.csv")

# calculate and transform variables
log_data$log_tr = log2(log_data$tree_r) #log-transformation of tree richness prior to scaling

log_data <- log_data %>%
  mutate(host_ab_tot = (par_abund + host_abund)) #calculate total host abundance

log_data <- log_data %>%
  mutate(p_rate = (par_abund / (par_abund + host_abund))) #calculate parasitism rate

log_data$log_hatot <- log1p(log_data$host_ab_tot) #log-transform host abundance

log_data$log_pa <- log1p(log_data$par_abund) #log-transform parasitoid abundance

log_data$rt_cwd <- sqrt(log_data$cwd_vol_imp) #square-root transform coarse-woody-debris volume

log_data$mc_new <- pmin(pmax(log_data$mc_mass, 0.0001), 1 - 0.0001) #reducing one value above 1 to slightly smaller than 1


#checking for correlations between variables
vars <- log_data[, c("log_tr", "cnpy_prop", "wantoc_imp", "host_ab_tot", "host_rich","par_abund", 
                     "par_rich", "ant_pres", "mc_imp", "cwd_vol_imp")]

cortable <- cor(vars, use = "complete.obs", method = "pearson")
cortable

cortable %>% 
  as.data.frame() %>% 
  rownames_to_column("var1") %>% 
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>% 
  filter(var1 < var2, abs(cor) > 0.50) %>% 
  arrange(desc(abs(cor))) # strong correlations between same-trophic-level abundance and richness

# change variables in factors with reference level
log_data$position <- as.factor(log_data$position)
log_data$treatment <- as.factor(log_data$treatment)
log_data$position <- relevel(log_data$position, ref = "S") # standing deadwood reference
log_data$treatment <- relevel(log_data$treatment, ref = "nA") # standing deadwood with ant exclusion reference

log_data$host_pres <- ifelse(log_data$host_ab_tot > 0, 1, 0) # calculate host presence/absence

log_data$ant_pres <- as.factor(log_data$ant_pres)
log_data$host_pres <- as.factor(log_data$host_pres)
log_data$geo_cat <- as.factor(log_data$geo_cat)
log_data$geo_cat <- relevel(log_data$geo_cat, ref = 1)

# ------------------------------------------------------------------------------
# 3. Exploratory analyses and scaling variables ####
# ------------------------------------------------------------------------------

table(log_data$treatment, log_data$host_abund > 0) # number of logs with hosts per treatment
table(log_data$treatment, log_data$ant_trap > 0) # number of logs with ants per treatment


with(insect_data,
     tapply(morphotype, list(treatment, level), function(x) length(unique(x)))) #species per treatment


# Comparing gmc variance between treatments
# Levene / Brown-Forsythe 
car::leveneTest(mc_mass ~ position, data = log_data, center = median)

tapply(log_data$mc_mass, log_data$position, IQR, na.rm = TRUE)

var_G  <- with(log_data, var(mc_mass[position=="G"],  na.rm=TRUE))
var_A  <- with(log_data, var(mc_mass[position=="S"],  na.rm=TRUE))

var_G / var_A # 1.8x
# lying deadwood has about twice the mean and variance in moiture content compared to standing


#Scaling all variables for modeling
#choose variables to scale
cols_to_scale <- c("log_tr", "host_rich", "log_hatot", "cnpy_prop",
                   "log_pa", "par_rich", "rt_cwd", "wantoc_imp")

#set new abbreviated scaled column names
name_map <- c(
  log_tr       = "sc_logtr",
  cnpy_prop    = "sc_cnpy",
  log_hatot    = "sc_loghatot",
  host_rich    = "sc_hr",
  log_pa       = "sc_logpa",
  par_rich     = "sc_pr",
  rt_cwd       = "sc_cwd",
  wantoc_imp   = "sc_wantoc"
)

#scale all selected variables in both dataframes
for (old in names(name_map)) {
  new <- name_map[[old]]
  log_data[[new]] <- as.numeric(scale(log_data[[old]]))
}

#scaling moisture content within each deadwood treatment
cols_to_scale_tr <- c("mc_imp")

name_map_tr <- c(
  mc_imp       = "sc_mctr")

scale_within_treatment <- function(df, vars, name_map, tr_var = "treatment") {  # function to scale  variables within-treatment
  
  for (old in names(name_map)) {
    new <- name_map[[old]]
    df[[new]] <- NA_real_
    
    for (tr in unique(df[[tr_var]])) {
      idx <- df[[tr_var]] == tr
      x   <- df[[old]][idx]
      
      if (sd(x, na.rm = TRUE) == 0) {
        df[[new]][idx] <- 0
      } else {
        df[[new]][idx] <- as.numeric(scale(x))
      }
    }
  }
  
  df
}

log_data <- scale_within_treatment(
  df       = log_data,
  vars     = cols_to_scale_tr,
  name_map = name_map_tr
)

# ------------------------------------------------------------------------------
# 4. Sampling completeness ####
# ------------------------------------------------------------------------------


log_data <- log_data %>%
  mutate(
    site    = factor(site),
    plot_id = factor(plot_id),
    treatment = factor(treatment),
    sample_id = interaction(site, plot_id, treatment, drop = TRUE)
  )

insect_data <- insect_data %>%
  mutate(
    treatment = factor(treatment),
    plot_id   = factor(plot_id),
    site      = factor(site),
    level     = factor(level)
  )

hosts <- insect_data %>%
  filter(level == "host")

hosts <- hosts %>%
  mutate(sample_id = interaction(site, plot_id, treatment, drop = TRUE))

host_comm <- hosts %>%                    # Species-by-log trap matrix
  group_by(sample_id, morphotype) %>%
  summarise(abund = sum(number_individuals), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = morphotype,
    values_from = abund,
    values_fill = 0
  )


host_df <- host_comm %>%
  left_join(log_data %>% dplyr::select(sample_id, plot_id, site, treatment),
            by = "sample_id")

#overall pooled sampling completeness
host_abund_only <- host_df[, !(names(host_df) %in% c("sample_id","plot_id","site","treatment"))]
host_totals     <- colSums(host_abund_only)

out_host_all <- iNEXT(host_totals, datatype = "abundance")
out_host_all$DataInfo$SC   # global sample coverage for hosts = 99%

#treatment x site assemblages
host_by_trtsite <- host_df %>%
  mutate(trt_site = paste(treatment, site, sep = "_")) %>%
  group_by(trt_site) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Make species-by-assemblage matrix
host_comm_trtsite <- as.data.frame(t(host_by_trtsite[, -1]))
colnames(host_comm_trtsite) <- host_by_trtsite$trt_site

out_host_trtsite <- iNEXT(host_comm_trtsite, datatype = "abundance")

out_host_trtsite$DataInfo # pooled sample coverage for each treatment & site

#plot x species assemblages
host_plot <- host_df %>%
  group_by(plot_id) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Species-only matrix
host_abund_plot <- host_plot[, -1]
rownames(host_abund_plot) <- host_plot$plot_id

chao_sample_coverage <- function(x) {
  x <- x[x > 0]
  n  <- sum(x)
  if (n == 0) return(NA_real_)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f1 == 0) return(1)
  if (f2 > 0) {
    1 - (f1 / n) * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
  } else {
    1 - (f1 / n) * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
  }
}

host_SC_plot <- apply(host_abund_plot, 1, chao_sample_coverage)

SC_plot_df <- data.frame(
  plot_id   = rownames(host_abund_plot),
  host_SC   = host_SC_plot,
  row.names = NULL
)

mean(SC_plot_df$host_SC) # mean per plot is 94%

plot_env <- log_data %>%
  distinct(plot_id, site, sc_logtr, sc_cnpy)

plot_env <- plot_env %>%
  left_join(SC_plot_df, by = "plot_id")
dat_host <- plot_env %>% filter(!is.na(host_SC)) # 3 plots in which hosts were not sampled at all

n_plot   <- nrow(dat_host)

# Avoid exact 0/1 for beta family
dat_host <- dat_host %>%
  mutate(host_SC_beta = (host_SC * (n_plot - 1) + 0.5) / n_plot)

dat_host <- dat_host %>%
  mutate(uncoverage = 1 - host_SC) %>%
  filter(uncoverage > 0)  # uncoverage


#model
m_host_SC <- glmmTMB(
  uncoverage ~ sc_logtr + sc_cnpy,
  data   = dat_host,
  family = Gamma(link = "log")
) 

summary(m_host_SC)

# plot-level sample coverage not affected by forest structure

#log-level
host_mat <- host_comm %>%
  as.data.frame()

rownames(host_mat) <- host_mat$sample_id
host_mat$sample_id <- NULL
host_SC_log <- apply(host_mat, 1, chao_sample_coverage)

SC_log_df <- data.frame(
  sample_id   = rownames(host_mat),
  host_SC   = host_SC_log,
  row.names = NULL
)

mean(SC_log_df$host_SC) # mean per log is 95%

plot_env <- log_data %>%
  distinct(sample_id, plot_id, site, treatment, sc_logtr, sc_cnpy)

log_env <- plot_env %>%
  left_join(SC_log_df, by = "sample_id")

dat_host <- log_env %>% filter(!is.na(host_SC)) # 125 logs with hosts

n_plot   <- nrow(dat_host)

# Avoid exact 0/1 for beta family
dat_host <- dat_host %>%
  mutate(host_SC_beta = (host_SC * (n_plot - 1) + 0.5) / n_plot)

dat_host <- dat_host %>%
  mutate(uncoverage = 1 - host_SC)

dat_host <- dat_host %>%
  mutate(incomplete = as.integer(uncoverage > 0))

m_inc <- glmmTMB(
  incomplete ~ treatment + sc_logtr + sc_cnpy +
    (1|plot_id),
  data = dat_host,
  family = binomial()) # probability of a log having incomplete sampling

m_sc <- glmmTMB(
  uncoverage ~ treatment + sc_logtr + sc_cnpy +
    (1|plot_id),
  data   = dat_host,
  ziformula = ~ treatment,
  family = tweedie(link = "log")
) 

summary(m_inc) # no bias in the probability of a log having incomplete sampling
summary(m_sc) # treatment A, G, and higher tree richness have lower sample coverage (higher uncoverage)

table(dat_host$uncoverage == 0, dat_host$treatment)

n_host <- rowSums(host_mat, na.rm = TRUE)

cor(dat_host$host_SC, n_host, use = "complete.obs")
cor(dat_host$uncoverage, n_host, use = "complete.obs") # sample coverage and uncoverage not
# correlated with abundance, suggesting that patters are driven primarily by rare species
# (singletons & doubletons)

stopifnot(all(dat_host$sample_id %in% rownames(host_mat)))

rare_df <- host_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(
    n_host   = rowSums(across(where(is.numeric)), na.rm = TRUE),
    S_obs    = rowSums(across(where(is.numeric)) > 0, na.rm = TRUE),
    f1       = rowSums(across(where(is.numeric)) == 1, na.rm = TRUE),
    f2       = rowSums(across(where(is.numeric)) == 2, na.rm = TRUE),
    rare2    = rowSums(across(where(is.numeric)) > 0 & across(where(is.numeric)) <= 2, na.rm = TRUE),
    rare_prop = ifelse(S_obs > 0, f1 / S_obs, NA_real_)
  ) %>%
  dplyr::select(sample_id, n_host, S_obs, f1, f2, rare2, rare_prop)

dat_host <- dat_host %>%
  left_join(rare_df, by = "sample_id") %>%
  mutate(
    ## guard against logs with n_host == 0 (if any slipped through)
    log_n_host = ifelse(n_host > 0, log(n_host), NA_real_),
    sc_log_n_host = as.numeric(scale(log_n_host))
  )

m_f1 <- glmmTMB(
  f1 ~ treatment + sc_logtr + sc_cnpy +
    offset(log_n_host) +
    (1 | plot_id),
  # ziformula = ~ treatment,
  family = nbinom2,
  data = dat_host
)

summary(m_f1)

m_rare2 <- glmmTMB(
  rare2 ~ treatment + sc_logtr + sc_cnpy +
    offset(log_n_host) +
    (1 | plot_id),
  ziformula = ~ treatment,
  family = nbinom2,
  data = dat_host
)

summary(m_rare2)

m_sc_effort <- glmmTMB(
  uncoverage ~ treatment + sc_logtr + sc_cnpy + log_n_host +
    (1|plot_id),
  # ziformula = ~ treatment,
  family = tweedie(link = "log"),
  data = dat_host
)

summary(m_sc_effort)

# Ground logs and tree-rich plots have fewer nests (see GLMM below of host abundance)
# Fewer nests = lower sample coverage.
# Additionally, the species that do occur there tend to be rarer (more singletons and doubletons).

# --------------------------------
## calculating estimated species richness based on 90% sample coverage
# for sensitivity check

# Species-only matrix with sample_id as rownames

# Turn each row into an abundance vector
host_list <- apply(host_mat, 1, function(x) x[x > 0])  # keep positive counts only
host_list <- host_list[lengths(host_list) > 0]         # drop logs with n = 0

# Estimate diversity at fixed sample coverage (SC = 0.90)
est_list <- lapply(names(host_list), function(id) {
  x <- host_list[[id]]
  
  # estimateD returns a data.frame; we keep q=0 (richness)
  ed <- estimateD(
    x,
    datatype = "abundance",
    base     = "coverage",
    level    = 0.90,
    conf     = 0.95
  )
  
  # Keep richness row (Order.q == 0)
  ed0 <- ed %>% filter(Order.q == 0)
  
  ed0$sample_id <- id
  ed0
})

host_rich_SC90 <- bind_rows(est_list)

# Inspect columns to see what iNEXT returned (names may vary slightly by version)
names(host_rich_SC90)

# observed SC per log
host_SC_log <- sapply(host_list, chao_sample_coverage)

SC_log_df <- data.frame(
  sample_id = names(host_SC_log),
  host_SC   = host_SC_log,
  row.names = NULL
)

mean(SC_log_df$host_SC) # mean SC per log is 95%

host_rich_SC90 <- host_rich_SC90 %>%
  left_join(SC_log_df, by = "sample_id")


# applying filters for species estimation, based on sample coverage cutoff and single species
dat_SC90_ok <- host_rich_SC90 %>%
  filter(host_SC >= 0.6)
dat_SC90_ok <- dat_SC90_ok %>%
  mutate(single_species = qD == 1)

log_data <- log_data %>% 
  left_join(dat_SC90_ok %>% dplyr::select(sample_id, qD),
            by = "sample_id")

hist(log1p(log_data$qD-1))

qD_mod <- glmmTMB(log1p(qD) ~ treatment + sc_mctr + ant_pres +
                    sc_logtr + sc_cnpy + 
                    sc_loghatot +
                    (1|plot_id),
                  data = log_data,
                  family = gaussian()) # only host abundance affects qD, same as for raw species richness

summary(qD_mod)

# -----------------------------------------------------------------------------------------
# 5. GLMMs ####
# -----------------------------------------------------------------------------------------

###
# ------ Moisture content
###
#log gravimetric moisture content
mc_mod <- glmmTMB(mc_new ~ treatment + sc_logtr + sc_cnpy +
                    geo_cat +
                    (1 | site/plot_id),
                  data = log_data,
                  dispformula = ~ site + treatment,
                  family = beta_family())

###
# ------ Ants
###
ap_mod <- glmmTMB(ant_pres ~ treatment + sc_logtr + 
                    sc_wantoc +  sc_cnpy + geo_cat +
                    (1| site/plot_id),
                  data = log_data,
                  family = binomial())

###
# ------ Hosts
###
#presence
hp_mod <- glmmTMB(host_pres ~ treatment + site + ant_pres +
                    (1| plot_id),
                  data = log_data,
                  family = binomial())

# abundance
ha_mod <- glmmTMB(host_ab_tot ~ treatment + sc_logtr + sc_mctr + ant_pres +
                    sc_cnpy + sc_cwd + 
                    (1 | site/plot_id),
                  data = log_data,
                  ziformula = ~ site + treatment,
                  family = nbinom2())

#abundance - design
ha_mod_exp <- glmmTMB(host_ab_tot ~ treatment + sc_logtr +
                        (1 | site/plot_id),
                      data = log_data,
                      ziformula = ~ site + treatment,
                      family = nbinom2())

#species richness
hr_mod <- glmmTMB(host_rich ~ treatment + sc_mctr +
                    sc_logtr + sc_cnpy + sc_loghatot + sc_cwd +
                    (1 | plot_id),   
                  data = log_data,
                  family = poisson())

hr_mod2 <- glmmTMB(host_rich ~ treatment + sc_mctr + ant_pres +
                     sc_logtr + sc_cnpy + sc_cwd +
                     (1 | site/plot_id),
                   data = log_data,
                   family = poisson()) # model without host abundance

#species richness - design
hr_mod_exp <- glmmTMB(host_rich ~ treatment + sc_logtr +
                        (1 | site/plot_id),
                      data = log_data,
                      ziformula = ~ site + treatment,
                      family = genpois())

###
# ------ Parasitoids
###

log_data_p <- subset(log_data, host_abund > 0) # only logs with at least one nest are considered

#abundance
pa_mod <- glmmTMB(par_abund ~ treatment + sc_mctr + ant_pres +
                    sc_logtr + sc_cnpy + sc_loghatot + sc_cwd +
                    (1 | site/plot_id),
                  data = log_data_p,
                  family = nbinom2())

pa_mod2 <- glmmTMB(par_abund ~ treatment + ant_pres + sc_mctr +
                     sc_logtr + sc_cnpy + sc_cwd +
                     (1 | site/plot_id),
                   data = log_data_p,
                   family = nbinom2()) # model without host mediation

#abundance - design
pa_mod_exp <- glmmTMB(par_abund ~ treatment + sc_logtr +
                        (1 | site/plot_id),
                      data = log_data_p,
                      family = nbinom1())

#species richness
pr_mod <- glmmTMB(par_rich ~ treatment + ant_pres + sc_mctr +
                    sc_logtr + sc_cnpy + sc_logpa +
                    (1 | plot_id),
                  data = log_data_p,
                  family = poisson()) # underdispersed model


pr_mod2 <- glmmTMB(par_rich ~ treatment + sc_hr + ant_pres + sc_mctr +
                     sc_logtr + sc_cnpy +
                     (1 | site/plot_id),
                   data = log_data_p,
                   family = poisson()) # parasiotid richness model without abundance, with host mediation

#species richness - design
pr_mod_exp <- glmmTMB(par_rich ~ treatment + sc_logtr +
                        (1 | site/plot_id),
                      data = log_data_p,
                      family = poisson())

###
# ------ Parasitism
###

p_modd <- glmmTMB(cbind(par_abund, host_abund) ~ treatment + sc_mctr + ant_pres +
                    sc_pr + sc_hr + 
                    sc_logtr + sc_cnpy +
                    (1 | site/plot_id),
                  data = log_data_p,
                  family = betabinomial())

p_modd2 <- glmmTMB(cbind(par_abund, host_abund) ~ treatment + sc_mctr + ant_pres +
                     sc_logtr + sc_cnpy +
                     (1 | plot_id),
                   data = log_data_p,
                   family = betabinomial()) # model without host and parasitoid variables

#parasitism - design
p_modd_exp <- glmmTMB(cbind(par_abund, host_abund) ~ treatment + sc_logtr +
                        (1 | plot_id),
                      data = log_data_p,
                      family = betabinomial())

models <- list(
  Moisture              = mc_mod,
  Ant_occurrence        = ap_mod,
  Host_occurrence       = hp_mod,
  Host_abundance        = ha_mod,
  Host_abundance_design = ha_mod_exp,
  Host_richness         = hr_mod,
  Host_richness_2       = hr_mod2,
  Host_richness_design  = hr_mod_exp,
  Parasitoid_abundance  = pa_mod,
  Parasitoid_abundance_2 = pa_mod2,
  Parasitoid_abundance_design = pa_mod_exp,
  Parasitoid_richness   = pr_mod,
  Parasitoid_richness2  = pr_mod2,
  Parasitoid_richness_design = pr_mod_exp,
  Parasitism            = p_modd,
  Parasitism_2          = p_modd2,
  Parasitism_design     = p_modd_exp
)

rm(list = names(models))

#export model summaries as markdown, PDF, and csv
model_summaries <- pretty_glmmTMB_summaries(models = models,
                                            file   = "model_summaries.md",
                                            digits = 3)
class(model_summaries)

combined_summaries <- bind_rows(
  lapply(names(models), function(m) {
    tidy(models[[m]], effects = "fixed", conf.int = TRUE) %>%
      mutate(model = m, .before = 1)
  })
)

write.csv(combined_summaries, "model_summaries.csv", row.names = FALSE)
rmarkdown::render("model_summaries.md", output_format = "pdf_document")


# -----------------------------------------------------------------------------------------
# 6. Community composition ####
# -----------------------------------------------------------------------------------------

log_data <- log_data %>%
  mutate(
    site    = factor(site),
    plot_id = factor(plot_id),
    treatment = factor(treatment),
    sample_id = interaction(site, plot_id, treatment, drop = TRUE)
  )

insect_data <- insect_data %>%
  mutate(
    treatment = factor(treatment),
    plot_id   = factor(plot_id),
    site      = factor(site),
    level     = factor(level)
  )

hosts <- insect_data %>%
  filter(level == "host")

hosts <- hosts %>%
  mutate(sample_id = interaction(site, plot_id, treatment, drop = TRUE))

host_comm <- hosts %>%     # Species-by-deadwood trap matrix
  group_by(sample_id, morphotype) %>%
  summarise(abund = sum(number_individuals), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = morphotype,
    values_from = abund,
    values_fill = 0
  )

host_comm_mat <- host_comm %>%
  dplyr::select(-sample_id) %>%
  as.data.frame()

host_meta <- hosts %>%
  distinct(sample_id, site, plot_id, treatment) %>%
  arrange(sample_id) # deadwood-level metadata for permanova

host_comm_tr <- host_comm[match(host_meta$sample_id, host_comm$sample_id), ]

host_comm_tr <- as.data.frame(host_comm_tr)
rownames(host_comm_tr) <- host_comm_tr$sample_id
host_comm_tr$sample_id <- NULL

summary(rowSums(host_comm_tr))         # total abundance per log
apply(host_comm_tr, 2, max)            # max per species
sort(colSums(host_comm_tr), decreasing = TRUE)[1:10]

colSums(host_comm_tr)  # Total abundance per species
plot(density(colSums(host_comm_tr), main = "Species Abundance Distribution"))

host_comm_tr <- sqrt(host_comm_tr) #square root recommended due to a few highly dominant species

env_vars <- log_data %>%
  dplyr::select(sample_id, mc_imp, ant_pres)  

host_meta <- host_meta %>%
  left_join(env_vars, by = c("sample_id"))

host_meta <- host_meta %>%
  group_by(treatment) %>%
  mutate(
    sc_mc = scale(mc_imp)[,1],     # GMC scaled within each treatment
    ant_pres = ant_pres  
  ) %>%
  ungroup()

all(rownames(host_comm_tr) == host_meta$sample_id) # should be TRUE if you set rownames

nrow(host_comm_tr); nrow(host_meta)   

#now creating matrix and metadata for plot-level permanova

host_comm2 <- hosts %>%     # Species-by-plot matrix
  group_by(plot_id, morphotype) %>%
  summarise(abund = sum(number_individuals), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = morphotype,
    values_from = abund,
    values_fill = 0
  )

host_comm_mat2 <- host_comm2 %>%
  dplyr::select(-plot_id) %>%
  as.data.frame()

host_meta2 <- hosts %>%
  distinct(site, plot_id) %>%
  arrange(plot_id) # plot-level metadata for permanova

host_comm_tr2 <- host_comm2[match(host_meta2$plot_id, host_comm2$plot_id), ]

host_comm_tr2 <- as.data.frame(host_comm_tr2)
rownames(host_comm_tr2) <- host_comm_tr2$plot_id
host_comm_tr2$plot_id <- NULL

host_comm_tr2 <- sqrt(host_comm_tr2) #square root recommended due to a few highly dominant species

env_vars <- log_data %>%
  dplyr::select(site, plot_id, cnpy_prop, log_tr, wantoc_imp, rt_cwd) %>%
  dplyr::distinct(plot_id, .keep_all = TRUE)

env_vars <- env_vars %>%
  semi_join(host_meta2, by = c("site", "plot_id")) %>%
  group_by(plot_id) %>%
  slice_head(n = 1) %>%
  ungroup() # filter env_vars to only plots with hosts

# Now left_join
host_meta2 <- host_meta2 %>%
  left_join(env_vars, by = c("site", "plot_id"))

host_meta2 <- host_meta2 %>%
  mutate(
    sc_cnpy = scale(cnpy_prop)[, 1],
    sc_wantoc  = scale(wantoc_imp)[, 1],
    sc_logtr = scale(log_tr)[, 1],
    sc_cwd = scale(rt_cwd)[, 1])  # Scale environmental predictors


all(rownames(host_comm_tr2) == host_meta2$plot_id) # should be TRUE if you set rownames
nrow(host_comm_tr2); nrow(host_meta2)   


###
# PERMANOVA
###
set.seed(123)

adon_treat <- adonis2(
  host_comm_tr ~ treatment + sc_mc + ant_pres,
  data         = host_meta,
  method       = "bray",
  permutations = 1000,
  strata       = host_meta$plot_id,  # restrict permutations within plots
  by           = "margin"   #or "term" to perform Type I
)

adon_treat # treatment significantly affects community composition

bc_dist <- vegdist(host_comm_tr, method = "bray")

anova(betadisper(bc_dist, host_meta$treatment)) # dispersion does not differ among treatments, results of permanova are safe

bd <- betadisper(bc_dist, host_meta$treatment)


anova(bd)              # parametric test
permutest(bd, permutations = 999)  # permutation-based test
plot(bd) # visualize


# Install pairwise Adonis function from Github
if (!require(pairwiseAdonis)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  library(pairwiseAdonis)
}

host_meta <- as.data.frame(host_meta)  # force data.frame class

adon_treat2 <- pairwise.adonis2(
  host_comm_tr ~ treatment,
  data         = host_meta,
  nperm        = 1000,     
  strata       = "plot_id"  # string name of column (quotes matter!)
)

adon_treat2 # G communities significantly different from SA and SnA


# Plot-level PERMANOVA
adon_env <- adonis2(
  host_comm_tr2 ~ sc_logtr + sc_cnpy + sc_wantoc,
  data         = host_meta2,
  method = "bray",
  binary = TRUE,          # use for turnover
  permutations = 1000,
  by           = "margin" 
)

adon_env 
# tree richness has a weak effect on community composition coming due to species turnover only
# see B-diversity analyses with MRM

bc_dist <- vegdist(host_comm_tr2, method = "bray")
bd <- betadisper(bc_dist, host_meta2$sc_logtr)

anova(bd)              # parametric test
permutest(bd, permutations = 999)  # permutation-based test
plot(bd) # visualize


###
# dbRDA
###

host_meta3 <- hosts %>%
  distinct(sample_id, site, plot_id, treatment) %>%
  arrange(sample_id)

env_vars <- log_data %>%
  dplyr::select(sample_id, site, plot_id, treatment,
                cnpy_prop, rt_cwd, log_tr, ant_pres, mc_imp)  

host_meta3 <- host_meta3 %>%
  left_join(env_vars, by = c("sample_id", "site", "plot_id", "treatment"))

host_meta3 <- host_meta3 %>%
  mutate(
    sc_logtr = scale(log_tr)[, 1]
  )  # Scale environmental predictors

host_cap <- capscale(
  host_comm_tr ~ treatment + sc_logtr,
  data     = host_meta3,
  distance = "bray"
)

anova(host_cap, by = "margin", permutations = 999)
scores(host_cap, display = "bp")

RsquareAdj(host_cap)

#extract information from dbRDA for ordination
#extract plot and species scores
site_sc <- as.data.frame(scores(host_cap, display = "sites", choices = 1:2))
bp_sc   <- scores(host_cap, display = "bp")

site_sc$sample_id <- rownames(site_sc)

host_sites_plot_env <- host_meta3 %>%
  inner_join(site_sc, by = "sample_id")

host_sites_plot_env <- host_sites_plot_env %>% # adding host abundance to plot point size
  left_join(
    log_data %>% dplyr::select(host_ab_tot, sample_id),
    by = "sample_id"
  )

env_fit <- envfit(
  host_cap ~ treatment + sc_logtr,
  data = host_meta3,
  permutations = 1000
)

eigenvals(host_cap)

eig <- eigenvals(host_cap, constrained = TRUE)
eig
eig / sum(eig) * 100
perc <- eig / sum(eig) * 100

spec_scores <- as.data.frame(
  scores(host_cap, display = "species", choices = 1:2)
)
spec_scores$morphotype <- rownames(spec_scores)

spec_scores$abs_score <- sqrt(spec_scores$CAP1^2 + spec_scores$CAP2^2)
top_species <- spec_scores %>%
  arrange(desc(abs_score)) %>%
  slice(1:3) # pick top 3 species

species_scaler <- 1.4   # or 6, or use ordiArrowMul

top_species$CAP1_scaled <- top_species$CAP1 * species_scaler
top_species$CAP2_scaled <- top_species$CAP2 * species_scaler


##
# Indicator species analysis 
##

# 1) metadata: one row per block, including uncolonized
meta_blocks <- log_data %>%
  distinct(log_id, treatment, position)  # adjust names

# 2) host counts per block x species (colonized only in raw, but we'll complete with zeros)
host_comm <- insect_data %>%
  filter(level == "host") %>%
  group_by(log_id, morphotype) %>%
  summarise(abund = sum(number_individuals), .groups = "drop") %>%
  pivot_wider(names_from = morphotype, values_from = abund, values_fill = 0)

# 3) add missing blocks as all-zero rows
host_comm_full <- meta_blocks %>%
  left_join(host_comm, by = "log_id") %>%
  mutate(across(-c(log_id, treatment, position), ~tidyr::replace_na(.x, 0)))

X <- host_comm_full %>%
  dplyr::select(-log_id, -treatment, -position) %>%
  as.matrix()

group_placement <- factor(host_comm_full$position)   # S vs G
group_treatment  <- factor(host_comm_full$treatment) # A, nA, G

ind_pos <- multipatt(X, group_placement, func="IndVal.g", control=how(nperm=999))
summary(ind_pos, alpha=0.05, indvalcomp = TRUE)

ind_trt <- multipatt(X, group_treatment, func="IndVal.g", control=how(nperm=999))
summary(ind_trt, alpha=0.05, indvalcomp = TRUE)


# -----------------------------------------------------------------------------------------
# 7. Beta-diversity ####
# -----------------------------------------------------------------------------------------

# First at deadwood treatment level (single pariwise comparison)
# presence–absence matrix S vs G (hosts)
host_pa_mat <- insect_data %>%
  filter(level == "host",
         !is.na(location),
         !is.na(morphotype)) %>%
  distinct(location, morphotype) %>%
  mutate(pres = 1) %>%
  pivot_wider(
    names_from  = morphotype,
    values_from = pres,
    values_fill = 0
  ) %>%
  mutate(location = factor(location, levels = c("S", "G"))) %>%
  arrange(location) %>%
  {
    mat <- as.matrix(dplyr::select(., -location))
    rownames(mat) <- .$location
    mat
  }

# Baselga decomposition (Sørensen family)
host_beta_pair <- beta.pair(host_pa_mat, index.family = "sorensen")

sor_mat  <- as.matrix(host_beta_pair$beta.sor)
sim_mat  <- as.matrix(host_beta_pair$beta.sim)
sne_mat  <- as.matrix(host_beta_pair$beta.sne)

beta_summary <- tibble(
  pair      = "S–G",
  beta_sor  = sor_mat["S", "G"],
  beta_turn = sim_mat["S", "G"],
  beta_nest = sne_mat["S", "G"]
) %>%
  mutate(
    prop_turn = beta_turn / beta_sor,
    prop_nest = beta_nest / beta_sor
  )

beta_summary # results

beta_summary <- tibble(
  pair      = "S–G",
  beta_sor  = sor_mat["S", "G"],
  beta_turn = sim_mat["S", "G"],
  beta_nest = sne_mat["S", "G"]
) %>%
  mutate(
    prop_turn = beta_turn / beta_sor,
    prop_nest = beta_nest / beta_sor
  )

beta_summary # results


#parasitoids
host_pa_mat <- insect_data %>%
  filter(level == "parasite",
         !is.na(location),
         !is.na(morphotype)) %>%
  distinct(location, morphotype) %>%
  mutate(pres = 1) %>%
  pivot_wider(
    names_from  = morphotype,
    values_from = pres,
    values_fill = 0
  ) %>%
  mutate(location = factor(location, levels = c("S", "G"))) %>%
  arrange(location) %>%
  {
    mat <- as.matrix(dplyr::select(., -location))
    rownames(mat) <- .$location
    mat
  }

host_beta_pair <- beta.pair(host_pa_mat, index.family = "sorensen")

sor_mat  <- as.matrix(host_beta_pair$beta.sor)
sim_mat  <- as.matrix(host_beta_pair$beta.sim)
sne_mat  <- as.matrix(host_beta_pair$beta.sne)

beta_summary <- tibble(
  pair      = "S–G",
  beta_sor  = sor_mat["S", "G"],
  beta_turn = sim_mat["S", "G"],
  beta_nest = sne_mat["S", "G"]
) %>%
  mutate(
    prop_turn = beta_turn / beta_sor,
    prop_nest = beta_nest / beta_sor
  )

beta_summary # results


Reduce(
  intersect,
  split(
    insect_data$morphotype[insect_data$level == "parasite"],
    insect_data$location[insect_data$level == "parasite"]
  )
)



# Second at plot-level for tree richness gradient, analyzed with MRM
insect_data2 <- insect_data %>%
  left_join(
    log_data %>% distinct(log_id, tree_r),
    by = "log_id"
  )

rich_levels <- c(1, 2, 4, 8, 16, 24)

insect_data2 <- insect_data2 %>%
  mutate(tree_lvl = factor(tree_r, levels = rich_levels))

host_pa_plot <- insect_data2 %>%
  filter(level == "host", !is.na(plot_id), !is.na(morphotype)) %>%
  distinct(plot_id, morphotype) %>%
  mutate(pres = 1) %>%
  pivot_wider(
    names_from  = morphotype,
    values_from = pres,
    values_fill = 0
  ) %>%
  arrange(plot_id) %>%
  {
    mat <- as.matrix(dplyr::select(., -plot_id))
    rownames(mat) <- as.character(.$plot_id)
    mat
  }

plot_meta <- insect_data2 %>%
  distinct(plot_id, tree_r, site)

# ensure consistent plot order
plots_all <- plot_meta %>%
  distinct(plot_id, tree_r, site) %>%
  mutate(plot_id = as.character(plot_id)) %>%
  arrange(plot_id)

keep <- plots_all$plot_id[plots_all$plot_id %in% rownames(host_pa_plot)]

mat_all <- host_pa_plot[keep, , drop = FALSE]

# Betapart decomposition (presence–absence)
beta_all <- beta.pair(mat_all, index.family = "sorensen")

beta_all

sne <- beta_all$beta.sne   # nestedness (dist)
sim <- beta_all$beta.sim   # turnover (dist)
sor <- beta_all$beta.sor   # total beta (dist)

# Tree richness distance
tr <- setNames(plots_all$tree_r, plots_all$plot_id)
tr <- tr[keep]
Dtree <- dist(log2(tr))
# Dtree <- dist(tr)    # for additive instead of multiplicative


# Site distance matrix (0 = same site, 1 = different site)
site_vec <- setNames(plots_all$site, plots_all$plot_id)[keep]
Dsite <- as.dist(outer(site_vec, site_vec, FUN = "!=") * 1)

stopifnot(length(Dtree) == length(Dsite))
stopifnot(length(Dtree) == length(sne))


mrm_nest <- MRM(sne ~ Dtree + Dsite, nperm = 999)
mrm_nest

mrm_turn <- MRM(sim ~ Dtree + Dsite, nperm = 999)
mrm_turn

mrm_total <- MRM(sor ~ Dtree + Dsite, nperm = 999)
mrm_total


# Convert distance objects to pairwise dataframe
n_plots <- length(keep)
plot_names <- keep

pairs_df <- data.frame()

for (i in 1:(n_plots - 1)) {
  for (j in (i + 1):n_plots) {
    # Get pairwise index for distance objects
    idx <- n_plots * (i - 1) - i * (i - 1) / 2 + j - i
    
    pairs_df <- rbind(pairs_df, data.frame(
      plot1 = plot_names[i],
      plot2 = plot_names[j],
      tree_richness_dist = as.numeric(Dtree)[idx],
      total_beta = as.numeric(sor)[idx],
      turnover = as.numeric(sim)[idx],
      nestedness = as.numeric(sne)[idx],
      site_dist = as.numeric(Dsite)[idx]
    ))
  }
}

# Add indicator for same vs different site
pairs_df <- pairs_df %>%
  mutate(same_site = ifelse(site_dist == 0, "Same site", "Different site"))

# Reshape for plotting
pairs_long <- pairs_df %>%
  pivot_longer(cols = c(turnover, nestedness),
               names_to = "component",
               values_to = "dissimilarity") %>%
  mutate(component = factor(component, 
                            levels = c("turnover", "nestedness"),
                            labels = c("Turnover", "Nestedness")))


# Summary statistics
cat("\nSummary statistics:\n")
cat("Total pairwise comparisons:", nrow(pairs_df), "\n")
cat("Mean Beta:", round(mean(pairs_df$total_beta), 3), "\n")
cat("Mean turnover:", round(mean(pairs_df$turnover), 3), "\n")
cat("Mean nestedness:", round(mean(pairs_df$nestedness), 3), "\n")

pairs_nonzero <- pairs_df %>% # calculate proportions only for pairs with non-zero total beta
  filter(total_beta > 0)

cat("Pairs with non-zero beta-diversity:", nrow(pairs_nonzero), "\n")
cat("Mean proportion turnover:", round(mean(pairs_nonzero$turnover / pairs_nonzero$total_beta), 3), "\n")
cat("Mean proportion nestedness:", round(mean(pairs_nonzero$nestedness / pairs_nonzero$total_beta), 3), "\n")
cat("\nPairs with zero beta-diversity (identical communities):", 
    nrow(pairs_df) - nrow(pairs_nonzero), "\n")
cat("Correlation between turnover and tree richness distance:", 
    round(cor(pairs_df$turnover, pairs_df$tree_richness_dist), 3), "\n")
cat("Correlation between nestedness and tree richness distance:", 
    round(cor(pairs_df$nestedness, pairs_df$tree_richness_dist), 3), "\n")


# -----------------------------------------------------------------------------------------
# 8. Mediation analyses ####
# -----------------------------------------------------------------------------------------

log_data$log_mc <- log(log_data$mc_imp) # transform
log_data$log_hatot <- log1p(log_data$host_ab_tot) # transform
log_data$sqrt_hr <- sqrt(log_data$host_rich)

# Host abundance
# mediation of moisture content
med_mc <- lmer(log_mc ~ position + (1 | plot_id), data = log_data)
out_mc <- lmer(log_hatot ~ log_mc + position + (1 | plot_id), data = log_data)

set.seed(22)

med_out_mc <- mediate(
  model.m      = med_mc,
  model.y      = out_mc,
  treat        = "position",
  mediator     = "log_mc",
  treat.value  = "G",
  control.value = "S",
  sims         = 1000
)

summary(med_out_mc) # approximately 34% of effect of ground contact is mediated by moisture

# mediation of ant presence
med_ant <- glmer(ant_pres ~ treatment + (1 | plot_id), data   = log_data,
                 family = binomial(link = "logit")) 

out_ant <- lmer(
  log_hatot ~ treatment + ant_pres + (1 | plot_id),
  data = log_data)

set.seed(22)

med_out_ant_pres <- mediate(
  model.m      = med_ant,
  model.y      = out_ant,
  treat        = "treatment",
  mediator     = "ant_pres",
  treat.value  = "G", # swap with A for comparison of A vs nA
  control.value = "nA",
  sims         = 1000
)

summary(med_out_ant_pres) # ant presence does not mediate effect of treatment


# Host richness
med_mc <- lmer(log_hatot ~ position + (1 | plot_id), data = log_data)
out_mc <- lmer(sqrt_hr ~ log_hatot + position + (1 | plot_id), data = log_data)

set.seed(22)

med_out_trt <- mediate(
  model.m      = med_mc,
  model.y      = out_mc,
  treat        = "position",
  mediator     = "log_hatot",
  treat.value  = "G",
  control.value = "S",
  sims         = 1000
)

summary(med_out_trt) # host abundance completely mediates effect of ground contact on host richness


# -----------------------------------------------------------------------------------------
# 9. SEM ####
# -----------------------------------------------------------------------------------------

# transform and scale variables
log_data$log_mc <- log(log_data$mc_imp) 
log_data$log_hatot <- log1p(log_data$host_ab_tot)
log_data$sqrt_hr <- sqrt(log_data$host_rich)

emp_logit <- function(y, m) {
  qlogis((y + 0.5) / (m + 1))
}

log_data$sem_cprate <- emp_logit(log_data$par_abund, log_data$host_ab_tot)

log_data$sem_mc <- as.numeric(scale(log_data$log_mc))
log_data$sem_ha <- as.numeric(scale(log_data$log_hatot))
log_data$sem_hr <- as.numeric(scale(log_data$sqrt_hr))
log_data$sem_pa <- as.numeric(scale(sqrt(log_data$par_rich)))
log_data$sem_pr <- as.numeric(scale(log1p(log_data$par_abund)))

log_data$ant_pres_num <- as.numeric(as.character(log_data$ant_pres))

#glmm sub-models for path analysis
prt_sem <- glmmTMB(sem_cprate ~ sem_pr + sem_hr +
                     (1 | site/plot_id),
                   data = log_data,
                   family = gaussian())

pr_sem <- glmmTMB(sem_pr ~ sem_pa + sem_hr +
                    (1 | plot_id),
                  data = log_data,
                  family = gaussian())

pa_sem <- glmmTMB(sem_pa ~ sem_ha + sc_cnpy + sc_cwd +
                    (1 | plot_id),
                  data = log_data,
                  family = gaussian())

hr_sem <- glmmTMB(sem_hr ~ sem_ha + sc_logtr +
                    (1 | plot_id),
                  data = log_data,
                  family = gaussian())

ha_sem <- glmmTMB(sem_ha ~ treatment + sem_mc + ant_pres_num +
                    (1 | site/plot_id),
                  data = log_data,
                  family = gaussian())

mc_sem <- glmmTMB(sem_mc ~ treatment +
                    (1 | site/plot_id),
                  data = log_data,
                  family = gaussian())

ant_sem <- glmmTMB(ant_pres_num ~ treatment +
                     (1 | plot_id),
                   data = log_data,
                   family = binomial())

cnpy_sem <- glmmTMB(sc_cnpy ~ sc_logtr +
                      (1 | site),
                    data = log_data,
                    family = gaussian())


sem_mods <- list(
  prt_sem = prt_sem,
  pr_sem  = pr_sem,
  pa_sem  = pa_sem,
  hr_sem  = hr_sem,
  ha_sem  = ha_sem,
  mc_sem  = mc_sem,
  ant_sem  = ant_sem,
  cnpy_sem = cnpy_sem,
  sem_cprate %~~% sem_ha,
  sem_cprate %~~% treatment,
  sem_cprate %~~% sem_pa
)


sem_all <- do.call(piecewiseSEM::psem, c(sem_mods, list(data = log_data)))
summary(sem_all, standardize = "scale", conserve = TRUE)

sem_coef <- as.data.frame(coefs(sem_all))
write.csv(sem_coef, "sem_coef.csv", row.names = FALSE) #export results

#calculating contrasts - 1 ground contact vs 2 standing & 1 ant exclusion vs 2 ant-accessible 
# ground contact vs standing on host abundance
emm_ha <- emmeans(ha_sem, ~ treatment)
contrast(
  emm_ha,
  method = list(G_vs_standing = c(-0.5, -0.5, 1))
)

# ground contact vs standing on moisture
emm_mc <- emmeans(mc_sem, ~ treatment)
contrast(
  emm_mc,
  method = list(G_vs_standing = c(-0.5, -0.5, 1))
)

emm_ant <- emmeans(ant_sem, ~ treatment)

# ant exclusion vs access on ant presence
contrast(
  emm_ant,
  method = list(nA_vs_access = c(1, -0.5, -0.5)),
  adjust = "none"
)

# ground contact vs standing on ant presence
contrast(
  emm_ant,
  method = list(G_vs_A_only = c(0, -1, 1)),
  adjust = "none"
)

message("Script ran successfully ✔️")
