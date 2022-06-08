x <- readRDS("processed_stanfits/maxshed21_nfpriors/release/best_si_nf.rds")
s4 <- x[[1]]
baseline <- x[[2]]

s3s4 <- readRDS("processed_stanfits/maxshed21_nfpriors/s3s4mix/best_si_nf.rds")[[1]]
obs_data <- cowling_data
## baseline - -15 to 34
## s3s4 - -12 to 42
## s4 -17 to 34
p1 <- SI_figure(baseline[["unconditional"]], obs_data) +
  expand_limits(x = c(-17, 43), y = 0.12) +
  theme(text = element_text(size = 24)) +
  coord_cartesian(clip = "off")


p2 <- plot_both_SIs(
  SI1 = s3s4[[1]], SI2 = s3s4[[2]],
  data = obs_data
) + expand_limits(x = c(-17, 43), y = 0.12) +
  theme(text = element_text(size = 24)) +
    coord_cartesian(clip = "off")

p3 <- plot_both_SIs(
  SI1 = s4[[1]], SI2 = s4[[2]],
  data = obs_data
) + expand_limits(x = c(-17, 43), y = 0.12) +
  theme(text = element_text(size = 24)) +
    coord_cartesian(clip = "off")


ggsave("baseline.png", p1)
ggsave("s4.png", p3)
ggsave("s3s4.png", p2)

baseline_t2 <- readRDS("processed_stanfits/maxshed21_nfpriors/release/scenario3a_nf_tab2.rds")
s4_t2 <- readRDS("processed_stanfits/maxshed21_nfpriors/release/scenario4a_nf_tab2.rds")
s3s4_t2 <- readRDS("processed_stanfits/maxshed21_nfpriors/s3s4mix/scenario4a_nf_tab2.rds")
