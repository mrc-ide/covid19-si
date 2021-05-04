
data <- readRDS("data/cowling_data_clean.rds")

data <- data%>%
  mutate(delay_import = infector_returnDate_fromOtherCity - infector_onsetDate)

# relationship between SI and importation (Y/N)
dplyr::count(data, infector_returned_fromOtherCity)

ggplot(data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity, y = si))

# relationship between SI and importation (days)
ggplot(data%>%filter(infector_returned_fromOtherCity == "imported"), aes(x = delay_import, y = si))+
  geom_point()+
  geom_smooth(method='lm')

# relationship between importation (Y/N) and time
ggplot(data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity, y = infector_onsetDate))

# relationship between importation (days) and time
ggplot(data%>%filter(infector_returned_fromOtherCity == "imported"))+
  geom_point(aes(x = delay_import, y = infector_onsetDate))+
  geom_abline(slope = 1, intercept = 0)

#comparing delay to isolation in imported versus local index cases
data %>%
  group_by(infector_returned_fromOtherCity) %>%
  dplyr::summarize(Mean = mean(onset_first_iso, na.rm=TRUE))

ggplot(data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity, y = onset_first_iso))

#plot date of isolation and date of importation

data <- data%>%
  mutate(first_iso = pmin(infector_isolateDate_beforeSymptom, infector_isolateDate_afterSymptom, na.rm = TRUE))

ggplot(data)+
  geom_point(aes(x = infector_returnDate_fromOtherCity, y = first_iso))+
  geom_abline(slope = 1, intercept = 0)

sum(data$first_iso == data$infector_returnDate_fromOtherCity, na.rm = TRUE)
