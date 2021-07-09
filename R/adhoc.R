library(dplyr)
relaxed <- readRDS("processed_stanfits/relaxed_priors/nf_overall_table2.rds")
s3s4mix <- readRDS("processed_stanfits/s3s4mix/nf_overall_table2.rds")
discrete <- readRDS("processed_stanfits/discrete_pairs/nf_overall_table2.rds")
x <- list(
  relaxed = relaxed, s3s4mix = s3s4mix,
  discrete_pairs = discrete
)
x <- bind_rows(x, .id = "meta-model")

## same underlying data
same_data <- list(
  relaxed = relaxed, s3s4mix = s3s4mix
) %>% bind_rows(.id = "meta-model")
same_data <- arrange(same_data, DIC)


## S3 fit for illustration has a problem.
## Fit again.
sim_si <- readRDS('data/illustration/sim_data_manuscript.rds')


### For each analyses, get the best model
### nf_priors
nf_discrete <- readRDS("processed_stanfits/maxshed21_nfpriors/discrete_pairs/nf_overall_table2.rds")
nf_discrete <- arrange(nf_discrete, DIC)
nf_discrete <- nf_discrete[1, ]

nf_s3s4 <- readRDS("processed_stanfits/maxshed21_nfpriors/s3s4mix/nf_overall_table2.rds")
nf_s3s4 <- arrange(nf_s3s4, DIC)
nf_s3s4 <- nf_s3s4[1, ]


nf_full <- readRDS("processed_stanfits/maxshed21_nfpriors/release/nf_overall_table2.rds")
nf_full <- arrange(nf_full, DIC)
nf_full <- nf_full[1, ]


relaxed_discrete <- readRDS("processed_stanfits/discrete_pairs/nf_overall_table2.rds")
relaxed_discrete <- arrange(relaxed_discrete, DIC)
relaxed_discrete <- relaxed_discrete[1, ]


relaxed_all <- readRDS("processed_stanfits/relaxed_priors/nf_overall_table2.rds")
relaxed_all <- arrange(relaxed_all, DIC)
relaxed_all <- relaxed_all[1, ]


relaxed_s3s4 <- readRDS("processed_stanfits/s3s4mix/nf_overall_table2.rds")
relaxed_s3s4 <- arrange(relaxed_s3s4, DIC)
relaxed_s3s4 <- relaxed_s3s4[1, ]

best_summary <- data.frame(
  priors = c("relaxed", "relaxed", "relaxed", "nf", "nf", "nf"),
  analysis = c("full", "discrete", "s3s4mix", "full", "discrete", "s3s4mix"),
  best_model = c(relaxed_all$model, relaxed_discrete$model, relaxed_s3s4$model,
                 nf_full$model, nf_discrete$model, nf_s3s4$model),
  mean_inf = c(relaxed_all$mean_inf, relaxed_discrete$mean_inf, relaxed_s3s4$mean_inf,
               nf_full$mean_inf, nf_discrete$mean_inf, nf_s3s4$mean_inf),
  sd_inf = c(relaxed_all$sd_inf, relaxed_discrete$sd_inf, relaxed_s3s4$sd_inf,
               nf_full$sd_inf, nf_discrete$sd_inf, nf_s3s4$sd_inf),
  mean_si = c(relaxed_all$mean_si, relaxed_discrete$mean_si, relaxed_s3s4$mean_si,
               nf_full$mean_si, nf_discrete$mean_si, nf_s3s4$mean_si),
  sd_si = c(relaxed_all$sd_si, relaxed_discrete$sd_si, relaxed_s3s4$sd_si,
               nf_full$sd_si, nf_discrete$sd_si, nf_s3s4$sd_si)
)


## Similar summary for parameters
nf_discrete <- readRDS("processed_stanfits/maxshed21_nfpriors/discrete_pairs/scenario4a_nf_tab1.rds")
nf_discrete$priors <- "nf"
nf_discrete$analysis <- "discrete"
nf_discrete <- tibble::rownames_to_column(nf_discrete)

nf_s3s4 <- readRDS("processed_stanfits/maxshed21_nfpriors/s3s4mix/scenario4a_nf_tab1.rds")
nf_s3s4$priors <- "nf"
nf_s3s4$analysis <- "s3s4"
nf_s3s4 <- tibble::rownames_to_column(nf_s3s4)


nf_full <- readRDS("processed_stanfits/maxshed21_nfpriors/release/scenario4a_nf_tab1.rds")
nf_full$priors <- "nf"
nf_full$analysis <- "full"
nf_full <- tibble::rownames_to_column(nf_full)

rel_discrete <- readRDS("processed_stanfits/discrete_pairs/scenario4a_nf_tab1.rds")
rel_discrete$priors <- "relaxed"
rel_discrete$analysis <- "discrete"
rel_discrete <- tibble::rownames_to_column(rel_discrete)

rel_s3s4 <- readRDS("processed_stanfits/s3s4mix/scenario4a_nf_tab1.rds")
rel_s3s4$priors <- "relaxed"
rel_s3s4$analysis <- "s3s4"
rel_s3s4 <- tibble::rownames_to_column(rel_s3s4)

rel_full <- readRDS("processed_stanfits/relaxed_priors/scenario4a_nf_tab1.rds")
rel_full$priors <- "relaxed"
rel_full$analysis <- "full"
rel_full <- tibble::rownames_to_column(rel_full)

best_summary <- rbind(nf_discrete, nf_s3s4, nf_full,
                      rel_discrete, rel_s3s4, rel_full)


best_summary <- select(
  best_summary, rowname, priors, analysis, best
)
best_summary <- spread(best_summary, rowname, best)
