
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
  #+
  #geom_abline(slope = 1, intercept = 0)

# relationship between importation (Y/N) and time
ggplot(data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity, y = infector_onsetDate))

# relationship between importation (days) and time
ggplot(data%>%filter(infector_returned_fromOtherCity == "imported"))+
  geom_point(aes(x = delay_import, y = infector_onsetDate))+
  geom_abline(slope = 1, intercept = 0)
