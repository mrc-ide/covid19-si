## Read and clean data here so all scripts use the
## same data. Model specific filters to be applied
## in model specific files
cowling_data <- readRDS("data/cowling_data_clean.rds") %>%
  mutate(si = as.numeric(si))%>%
   rename(nu = onset_first_iso)%>%
   filter(!is.na(nu))


# offset - must be a negative number!
offset <- -20
# sub-set to only incude those SIs that are possible under our assumed offset
data_offset <- filter(cowling_data, si > offset, nu > offset)


# discrete pairs
data_discrete_pairs <- filter(data_offset, cluster_size ==2)

# whole data set with nu replacement for s3 type pairs
data_s3_s4mix <- data_offset %>%
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

# discrete pairs with nu replacement for s3 type pairs
data_discrete_pairs_s3_s4mix <- filter(data_s3_s4mix, cluster_size == 2)

data_offset <- arrange(data_offset, nu)
