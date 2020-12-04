# looking at importation events
data <- readRDS("data/cowling_data_clean.rds")

data_c <- data%>%
  mutate(si = as.numeric(si))%>%
  dplyr::rename(nu = onset_first_iso)%>%
  dplyr::filter(!is.na(nu))%>%
  mutate(delay_import = infector_returnDate_fromOtherCity - infector_onsetDate)

# most people developed symptoms after they returned but some very shortly after
ggplot(data = data_c%>%filter(!is.na(delay_import)))+
  geom_histogram(aes(x = as.numeric(delay_import)), binwidth = 1)

# index cases significantly likely to be local later in the pandemic cf earlier
ggplot(data = data_c, aes(group = infector_returned_fromOtherCity))+
  geom_boxplot(aes(x = infector_returned_fromOtherCity, y = infector_onsetDate))
