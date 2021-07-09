## Read and clean data here so all scripts use the
## same data. Model specific filters to be applied
## in model specific files
cowling_data <- readRDS("data/cowling_data_clean.rds") %>%
  mutate(si = as.numeric(si))%>%
   rename(nu = onset_first_iso)%>%
   filter(!is.na(nu))

cowling_data <- arrange(cowling_data, nu)


# discrete pairs
data_discrete_pairs <- filter(cowling_data, cluster_size == 2)
data_discrete_pairs <- arrange(data_discrete_pairs, nu)

# whole data set with nu replacement for s3 type pairs
data_s3_s4mix <- cowling_data %>%
  mutate(
    infector_first_isolateDate = pmap_dbl(
      list(infector_isolateDate_beforeSymptom,
           infector_isolateDate_afterSymptom), min, na.rm = TRUE
    ),
    infectee_first_isolateDate = pmap_dbl(
      list(infectee_isolateDate_beforeSymptom,
           infectee_isolateDate_afterSymptom), min, na.rm = TRUE
    )
  )

## when a value is NA it results in Inf
data_s3_s4mix$infector_first_isolateDate[! is.finite(data_s3_s4mix$infector_first_isolateDate)] <- NA
data_s3_s4mix$infector_first_isolateDate <- as.Date(data_s3_s4mix$infector_first_isolateDate, origin = "1970-01-01")
## when both values are NA it results in Inf
data_s3_s4mix$infectee_first_isolateDate[! is.finite(data_s3_s4mix$infectee_first_isolateDate)] <- NA
data_s3_s4mix$infectee_first_isolateDate <- as.Date(data_s3_s4mix$infectee_first_isolateDate, origin = "1970-01-01")

for(i in 1:(length(data_s3_s4mix$infector_first_isolateDate))){
  if(!is.na(data_s3_s4mix$infectee_first_isolateDate[i])){
    if(data_s3_s4mix$infector_first_isolateDate[i] >= data_s3_s4mix$infectee_onsetDate[i] && data_s3_s4mix$infector_first_isolateDate[i] >= data_s3_s4mix$infectee_first_isolateDate[i]){
      if(data_s3_s4mix$infectee_first_isolateDate[i] >= data_s3_s4mix$infectee_onsetDate[i]){
        data_s3_s4mix$nu[i] <- 41
      }
    }
  }
}
data_s3_s4mix <- arrange(data_s3_s4mix, nu)

## discrete pairs with nu replacement for s3 type pairs
data_discrete_pairs_s3_s4mix <- filter(data_s3_s4mix, cluster_size == 2)
data_discrete_pairs_s3_s4mix <- arrange(data_discrete_pairs_s3_s4mix, nu)
