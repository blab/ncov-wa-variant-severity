# Script to prepare data for analysis: severity and vaccine effectiveness by variant
# for sequenced case data in Washington state.
#
# Data curated by WA-DOH from public health surveillance records and GISAID sequence data.
# Variant calls using pangolin v1.2.13. 
#

library(tidyverse)
library(lubridate)

# load data
raw <- read.delim('data_pull_2021-09-02_subset.csv',sep=',',header = TRUE)

names(raw)



######################################################
### set up workspace and initialize helper functions
######################################################


# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
# https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/variant-surveillance/variant-info.html
WHO_names <- list(B.1.1.7 = 'Alpha (B.1.1.7)',
                  Q.1 = 'Alpha (B.1.1.7)',
                  Q.2 = "Alpha (B.1.1.7)",
                  Q.3 = "Alpha (B.1.1.7)",
                  Q.4 = "Alpha (B.1.1.7)",
                  Q.5 = "Alpha (B.1.1.7)",
                  Q.6 = "Alpha (B.1.1.7)",
                  Q.7 = "Alpha (B.1.1.7)",
                  Q.8 = "Alpha (B.1.1.7)",
                  B.1.351 = "Beta (B.1.351)",
                  B.1.351.1 = "Beta (B.1.351)",
                  B.1.351.2 = "Beta (B.1.351)",
                  B.1.351.3 = "Beta (B.1.351)",
                  B.1.351.4 = "Beta (B.1.351)",
                  B.1.351.5 = "Beta (B.1.351)",
                  AY.1 = "Delta (B.1.617.2)",
                  AY.10 = "Delta (B.1.617.2)",
                  AY.100 = "Delta (B.1.617.2)",
                  AY.101 = "Delta (B.1.617.2)",
                  AY.102 = "Delta (B.1.617.2)",
                  AY.102.1 = "Delta (B.1.617.2)",
                  AY.102.2 = "Delta (B.1.617.2)",
                  AY.103 = "Delta (B.1.617.2)",
                  AY.104 = "Delta (B.1.617.2)",
                  AY.105 = "Delta (B.1.617.2)",
                  AY.106 = "Delta (B.1.617.2)",
                  AY.107 = "Delta (B.1.617.2)",
                  AY.108 = "Delta (B.1.617.2)",
                  AY.109 = "Delta (B.1.617.2)",
                  AY.11 = "Delta (B.1.617.2)",
                  AY.110 = "Delta (B.1.617.2)",
                  AY.111 = "Delta (B.1.617.2)",
                  AY.112 = "Delta (B.1.617.2)",
                  AY.112.1 = "Delta (B.1.617.2)",
                  AY.113 = "Delta (B.1.617.2)",
                  AY.114 = "Delta (B.1.617.2)",
                  AY.115 = "Delta (B.1.617.2)",
                  AY.116 = "Delta (B.1.617.2)",
                  AY.116.1 = "Delta (B.1.617.2)",
                  AY.117 = "Delta (B.1.617.2)",
                  AY.118 = "Delta (B.1.617.2)",
                  AY.119 = "Delta (B.1.617.2)",
                  AY.119.1 = "Delta (B.1.617.2)",
                  AY.119.2 = "Delta (B.1.617.2)",
                  AY.12 = "Delta (B.1.617.2)",
                  AY.120 = "Delta (B.1.617.2)",
                  AY.120.1 = "Delta (B.1.617.2)",
                  AY.120.2 = "Delta (B.1.617.2)",
                  AY.120.2.1 = "Delta (B.1.617.2)",
                  AY.121 = "Delta (B.1.617.2)",
                  AY.121.1 = "Delta (B.1.617.2)",
                  AY.122 = "Delta (B.1.617.2)",
                  AY.122.1 = "Delta (B.1.617.2)",
                  AY.122.2 = "Delta (B.1.617.2)",
                  AY.122.3 = "Delta (B.1.617.2)",
                  AY.123 = "Delta (B.1.617.2)",
                  AY.123.1 = "Delta (B.1.617.2)",
                  AY.124 = "Delta (B.1.617.2)",
                  AY.124.1 = "Delta (B.1.617.2)",
                  AY.125 = "Delta (B.1.617.2)",
                  AY.126 = "Delta (B.1.617.2)",
                  AY.127 = "Delta (B.1.617.2)",
                  AY.127.1 = "Delta (B.1.617.2)",
                  AY.128 = "Delta (B.1.617.2)",
                  AY.129 = "Delta (B.1.617.2)",
                  AY.13 = "Delta (B.1.617.2)",
                  AY.130 = "Delta (B.1.617.2)",
                  AY.131 = "Delta (B.1.617.2)",
                  AY.132 = "Delta (B.1.617.2)",
                  AY.133 = "Delta (B.1.617.2)",
                  AY.14 = "Delta (B.1.617.2)",
                  AY.15 = "Delta (B.1.617.2)",
                  AY.16 = "Delta (B.1.617.2)",
                  AY.16.1 = "Delta (B.1.617.2)",
                  AY.17 = "Delta (B.1.617.2)",
                  AY.18 = "Delta (B.1.617.2)",
                  AY.19 = "Delta (B.1.617.2)",
                  AY.2 = "Delta (B.1.617.2)",
                  AY.20 = "Delta (B.1.617.2)",
                  AY.20.1 = "Delta (B.1.617.2)",
                  AY.21 = "Delta (B.1.617.2)",
                  AY.22 = "Delta (B.1.617.2)",
                  AY.23 = "Delta (B.1.617.2)",
                  AY.23.1 = "Delta (B.1.617.2)",
                  AY.23.2 = "Delta (B.1.617.2)",
                  AY.24 = "Delta (B.1.617.2)",
                  AY.25 = "Delta (B.1.617.2)",
                  AY.25.1 = "Delta (B.1.617.2)",
                  AY.25.1.1 = "Delta (B.1.617.2)",
                  AY.26 = "Delta (B.1.617.2)",
                  AY.26.1 = "Delta (B.1.617.2)",
                  AY.27 = "Delta (B.1.617.2)",
                  AY.28 = "Delta (B.1.617.2)",
                  AY.29 = "Delta (B.1.617.2)",
                  AY.29.1 = "Delta (B.1.617.2)",
                  AY.3 = "Delta (B.1.617.2)",
                  AY.3.1 = "Delta (B.1.617.2)",
                  AY.3.2 = "Delta (B.1.617.2)",
                  AY.3.3 = "Delta (B.1.617.2)",
                  AY.30 = "Delta (B.1.617.2)",
                  AY.31 = "Delta (B.1.617.2)",
                  AY.32 = "Delta (B.1.617.2)",
                  AY.33 = "Delta (B.1.617.2)",
                  AY.33.1 = "Delta (B.1.617.2)",
                  AY.34 = "Delta (B.1.617.2)",
                  AY.34.1 = "Delta (B.1.617.2)",
                  AY.34.1.1 = "Delta (B.1.617.2)",
                  AY.34.2 = "Delta (B.1.617.2)",
                  AY.35 = "Delta (B.1.617.2)",
                  AY.36 = "Delta (B.1.617.2)",
                  AY.37 = "Delta (B.1.617.2)",
                  AY.38 = "Delta (B.1.617.2)",
                  AY.39 = "Delta (B.1.617.2)",
                  AY.39.1 = "Delta (B.1.617.2)",
                  AY.39.1.1 = "Delta (B.1.617.2)",
                  AY.39.1.2 = "Delta (B.1.617.2)",
                  AY.39.1.3 = "Delta (B.1.617.2)",
                  AY.39.2 = "Delta (B.1.617.2)",
                  AY.4 = "Delta (B.1.617.2)",
                  AY.4.1 = "Delta (B.1.617.2)",
                  AY.4.10 = "Delta (B.1.617.2)",
                  AY.4.2 = "Delta (B.1.617.2)",
                  AY.4.2.1 = "Delta (B.1.617.2)",
                  AY.4.2.2 = "Delta (B.1.617.2)",
                  AY.4.2.3 = "Delta (B.1.617.2)",
                  AY.4.3 = "Delta (B.1.617.2)",
                  AY.4.4 = "Delta (B.1.617.2)",
                  AY.4.5 = "Delta (B.1.617.2)",
                  AY.4.6 = "Delta (B.1.617.2)",
                  AY.4.7 = "Delta (B.1.617.2)",
                  AY.4.8 = "Delta (B.1.617.2)",
                  AY.4.9 = "Delta (B.1.617.2)",
                  AY.40 = "Delta (B.1.617.2)",
                  AY.41 = "Delta (B.1.617.2)",
                  AY.42 = "Delta (B.1.617.2)",
                  AY.42.1 = "Delta (B.1.617.2)",
                  AY.43 = "Delta (B.1.617.2)",
                  AY.43.1 = "Delta (B.1.617.2)",
                  AY.43.2 = "Delta (B.1.617.2)",
                  AY.43.3 = "Delta (B.1.617.2)",
                  AY.43.4 = "Delta (B.1.617.2)",
                  AY.43.5 = "Delta (B.1.617.2)",
                  AY.43.6 = "Delta (B.1.617.2)",
                  AY.43.7 = "Delta (B.1.617.2)",
                  AY.44 = "Delta (B.1.617.2)",
                  AY.45 = "Delta (B.1.617.2)",
                  AY.46 = "Delta (B.1.617.2)",
                  AY.46.1 = "Delta (B.1.617.2)",
                  AY.46.2 = "Delta (B.1.617.2)",
                  AY.46.3 = "Delta (B.1.617.2)",
                  AY.46.4 = "Delta (B.1.617.2)",
                  AY.46.5 = "Delta (B.1.617.2)",
                  AY.46.6 = "Delta (B.1.617.2)",
                  AY.46.6.1 = "Delta (B.1.617.2)",
                  AY.47 = "Delta (B.1.617.2)",
                  AY.48 = "Delta (B.1.617.2)",
                  AY.49 = "Delta (B.1.617.2)",
                  AY.5 = "Delta (B.1.617.2)",
                  AY.5.1 = "Delta (B.1.617.2)",
                  AY.5.2 = "Delta (B.1.617.2)",
                  AY.5.3 = "Delta (B.1.617.2)",
                  AY.5.4 = "Delta (B.1.617.2)",
                  AY.5.5 = "Delta (B.1.617.2)",
                  AY.50 = "Delta (B.1.617.2)",
                  AY.51 = "Delta (B.1.617.2)",
                  AY.52 = "Delta (B.1.617.2)",
                  AY.53 = "Delta (B.1.617.2)",
                  AY.54 = "Delta (B.1.617.2)",
                  AY.55 = "Delta (B.1.617.2)",
                  AY.56 = "Delta (B.1.617.2)",
                  AY.57 = "Delta (B.1.617.2)",
                  AY.58 = "Delta (B.1.617.2)",
                  AY.59 = "Delta (B.1.617.2)",
                  AY.6 = "Delta (B.1.617.2)",
                  AY.60 = "Delta (B.1.617.2)",
                  AY.61 = "Delta (B.1.617.2)",
                  AY.62 = "Delta (B.1.617.2)",
                  AY.63 = "Delta (B.1.617.2)",
                  AY.64 = "Delta (B.1.617.2)",
                  AY.65 = "Delta (B.1.617.2)",
                  AY.66 = "Delta (B.1.617.2)",
                  AY.67 = "Delta (B.1.617.2)",
                  AY.68 = "Delta (B.1.617.2)",
                  AY.69 = "Delta (B.1.617.2)",
                  AY.7 = "Delta (B.1.617.2)",
                  AY.7.1 = "Delta (B.1.617.2)",
                  AY.7.2 = "Delta (B.1.617.2)",
                  AY.70 = "Delta (B.1.617.2)",
                  AY.71 = "Delta (B.1.617.2)",
                  AY.72 = "Delta (B.1.617.2)",
                  AY.73 = "Delta (B.1.617.2)",
                  AY.74 = "Delta (B.1.617.2)",
                  AY.75 = "Delta (B.1.617.2)",
                  AY.75.1 = "Delta (B.1.617.2)",
                  AY.75.2 = "Delta (B.1.617.2)",
                  AY.75.3 = "Delta (B.1.617.2)",
                  AY.76 = "Delta (B.1.617.2)",
                  AY.77 = "Delta (B.1.617.2)",
                  AY.78 = "Delta (B.1.617.2)",
                  AY.79 = "Delta (B.1.617.2)",
                  AY.8 = "Delta (B.1.617.2)",
                  AY.80 = "Delta (B.1.617.2)",
                  AY.81 = "Delta (B.1.617.2)",
                  AY.82 = "Delta (B.1.617.2)",
                  AY.83 = "Delta (B.1.617.2)",
                  AY.84 = "Delta (B.1.617.2)",
                  AY.85 = "Delta (B.1.617.2)",
                  AY.86 = "Delta (B.1.617.2)",
                  AY.87 = "Delta (B.1.617.2)",
                  AY.88 = "Delta (B.1.617.2)",
                  AY.89 = "Delta (B.1.617.2)",
                  AY.9 = "Delta (B.1.617.2)",
                  AY.9.1 = "Delta (B.1.617.2)",
                  AY.9.2 = "Delta (B.1.617.2)",
                  AY.9.2.1 = "Delta (B.1.617.2)",
                  AY.9.2.2 = "Delta (B.1.617.2)",
                  AY.90 = "Delta (B.1.617.2)",
                  AY.91 = "Delta (B.1.617.2)",
                  AY.91.1 = "Delta (B.1.617.2)",
                  AY.92 = "Delta (B.1.617.2)",
                  AY.93 = "Delta (B.1.617.2)",
                  AY.94 = "Delta (B.1.617.2)",
                  AY.95 = "Delta (B.1.617.2)",
                  AY.96 = "Delta (B.1.617.2)",
                  AY.97 = "Delta (B.1.617.2)",
                  AY.98 = "Delta (B.1.617.2)",
                  AY.98.1 = "Delta (B.1.617.2)",
                  AY.99 = "Delta (B.1.617.2)",
                  AY.99.1 = "Delta (B.1.617.2)",
                  AY.99.2 = "Delta (B.1.617.2)",
                  B.1.617.2 = "Delta (B.1.617.2)",
                  B.1.427 = "Epsilon (B.1.427/B.1.429)",
                  B.1.429 = "Epsilon (B.1.427/B.1.429)",
                  B.1.525 = "Eta (B.1.525)",
                  P.1 = "Gamma (P.1)",
                  P.1.1 = "Gamma (P.1)",
                  P.1.10 = "Gamma (P.1)",
                  P.1.10.1 = "Gamma (P.1)",
                  P.1.10.2 = "Gamma (P.1)",
                  P.1.11 = "Gamma (P.1)",
                  P.1.12 = "Gamma (P.1)",
                  P.1.12.1 = "Gamma (P.1)",
                  P.1.13 = "Gamma (P.1)",
                  P.1.14 = "Gamma (P.1)",
                  P.1.15 = "Gamma (P.1)",
                  P.1.16 = "Gamma (P.1)",
                  P.1.17 = "Gamma (P.1)",
                  P.1.17.1 = "Gamma (P.1)",
                  P.1.2 = "Gamma (P.1)",
                  P.1.3 = "Gamma (P.1)",
                  P.1.4 = "Gamma (P.1)",
                  P.1.5 = "Gamma (P.1)",
                  P.1.6 = "Gamma (P.1)",
                  P.1.7 = "Gamma (P.1)",
                  P.1.7.1 = "Gamma (P.1)",
                  P.1.8 = "Gamma (P.1)",
                  P.1.9 = "Gamma (P.1)",
                  B.1.617.1 = "Kappa (B.1.617.1)",
                  B.1.621 = "Mu (B.1.621)",
                  B.1.621.1 = "Mu (B.1.621)",
                  B.1.526='Iota (B.1.526)',
                  B.1.525='Eta (B.1.525)',
                  C.37='Lambda (C.37)',
                  B.1.1.529 = "Omicron (B.1.1.529)",
                  BA.1 = "Omicron (B.1.1.529)",
                  BA.2 = "Omicron (B.1.1.529)",
                  BA.3 = "Omicron (B.1.1.529)",
                  P.2 = "Zeta (P.2)",
                  P.3='other', # P.3 assigned to non-variant due to lack of VOI/VOC classification by WHO or lack of detection in WA state
                  other_voc_voi='other_voc_voi', # this is here for the colormap
                  other='other') 


### Helper function

compareNA <- function(v1,v2,test='equal') { # http://www.cookbook-r.com/Manipulating_data/Comparing_vectors_or_factors_with_NA/ via https://stackoverflow.com/a/16822770/2184306
  # This function returns TRUE wherever elements are the same, including NA's,
  # and false everywhere else.
  if (test=='equal')
    same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  else if(test=='<')
    same <- (v1 < v2)  |  (is.na(v1) & is.na(v2))
  else if(test=='>')
    same <- (v1 > v2)  |  (is.na(v1) & is.na(v2))
  else if(test=='<=')
    same <- (v1 <= v2)  |  (is.na(v1) & is.na(v2))
  else if(test=='>=')
    same <- (v1 >= v2)  |  (is.na(v1) & is.na(v2))
  
  same[is.na(same)] <- FALSE
  return(same)
}




#########################################
#### prepare data for analysis
#########################################


#### cleaning data frame
d<-raw

# data.frame for keeping track of exclusions
exclusions <- data.frame(data_view=c(),n_kept=c(),reason=c())
exclusions <- exclusions %>% rbind(data.frame(data_view='all',reason='raw data',n_kept=nrow(d)))

sum((d$lineage=='None'),na.rm=TRUE)
sum(is.na(d$lineage=='None'),na.rm=TRUE)

# Excluding poor quality seqs
exclude_seqs <- read.delim('exclusion_ids_lf.csv',sep=',',header = TRUE)
#exclude_seqs <- read.delim('Y:/Confidential/DCHS/CDE/Zoonoses/Molecular/COVID/Analysis Projects/Variant severity analysis/exclusion_ids_lf.csv',sep=',',header = TRUE)

sum((d$lineage=='None'),na.rm=TRUE)
sum(is.na(d$lineage=='None'),na.rm=TRUE)

## drop those with lineage=='None' or no lineage (n=223)
d <- d %>% filter( !(lineage=='None' | is.na(lineage)))
exclusions <- exclusions %>% rbind(data.frame(data_view='with_known_lineage',reason='filtered out lineage==None',n_kept=nrow(d)))


## exclude those with bad sequence quality
d <- d %>%  filter(!(d$CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER %in% exclude_seqs$ID))
exclusions <- exclusions %>% rbind(data.frame(data_view='good sequence quality',reason='>10% ambiguity in sequencing',n_kept=nrow(d)))


## WHO names
d$who_lineage <- as.character(d$lineage)
for(n in names(WHO_names)){
  idx <- d$who_lineage==n
  d$who_lineage[idx] <- WHO_names[[n]]
}

## relabel all non-VOC as "other"
d$who_lineage[!(d$who_lineage %in% WHO_names)] <- 'other'
d$who_lineage <- factor(d$who_lineage,levels=unique(WHO_names))
d$who_lineage <- relevel(d$who_lineage,ref='other')

#format AGE_YEARS to numeric
d$AGE_YEARS <- as.numeric(d$AGE_YEARS)

# age bin
# include >99 in the 90-99 age bracket 
d$age_bin = factor(5+floor(pmin(99,d$AGE_YEARS)/10)*10) 


# drop miscoded age
d <- d %>% filter(AGE_YEARS!=141) %>% droplevels()
exclusions <- exclusions %>% rbind(data.frame(data_view='valid ages',reason='miscoded age',n_kept=nrow(d)))


# dates
d$collection_date <- as.Date(d$collection_date, format = "%m/%d/%Y")
####
#ASK STEPH TO CHECK ADMITDATE 
####
d$admitdate <- as.Date(d$admitdate, format = "%m/%d/%Y")
d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 <- as.Date(d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1, format = "%m/%d/%Y")
d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 <- as.Date(d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2, format = "%m/%d/%Y")
d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3 <- as.Date(d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3, format = "%m/%d/%Y")
d$earliest_positive_test_date <- as.Date(d$earliest_positive_test_date, format = "%m/%d/%Y")

d$week_collection <- paste(year(d$collection_date),sprintf('%02d',week(d$collection_date)),sep='-')
d$week_collection <- factor(d$week_collection, levels=sort(unique(d$week_collection)))
d$week_collection_number <- as.numeric(d$week_collection)

# outcomes
d$mhosp <- as.character(d$mhosp) 
d$mhosp[is.na(d$mhosp)] <- 'No'
d$mhosp <- factor(d$mhosp, levels=c('No','Yes'))

#exluding those with admission but no admit date
d <- d %>% filter(!((mhosp == "Yes")&is.na(admitdate)))
exclusions <- exclusions %>% rbind(data.frame(data_view='valid admit date',reason='no admitdate for hosp admission',n_kept=nrow(d)))

# make county shorter
d$COUNTY <- sub(' County','',d$COUNTY)
d$COUNTY[is.na(d$COUNTY)] <- "Unknown"
d$COUNTY <- relevel(factor(d$COUNTY),ref='King')

# race as factor
race_cats <- names(d)[grepl('RACE',names(d))]
d$race <- NA
for (r in 1:length(race_cats)){
  
  race_idx <- d[,race_cats[r]]=='Y'
  multiple_idx <- !is.na(d$race)
  
  d$race[race_idx & !multiple_idx] <- race_cats[r]
  d$race[race_idx & multiple_idx] <- 'RACE_MULTIPLE' # define multiple race category
}
d$race<-sub('RACE_','',d$race)
d$race[is.na(d$race)] <- 'UNKNOWN' # code NA as UNKNOWN

d$race<-relevel(factor(str_to_sentence(d$race)),ref='White')

# sex unknown
d$SEX_AT_BIRTH <- as.character(d$SEX_AT_BIRTH)
d$SEX_AT_BIRTH[is.na(d$SEX_AT_BIRTH)] <- 'Unknown'
d$SEX_AT_BIRTH <- relevel(factor(d$SEX_AT_BIRTH),ref='Female')

# clean up test date vs collection date
d$days_between_collection_and_earliest_test <- as.numeric(d$collection_date - d$earliest_positive_test_date)

  # collection date before earliest test date. 
  tmp <- d %>% filter(d$days_between_collection_and_earliest_test>=0) %>%
    select(CASE_ID,days_between_collection_and_earliest_test,collection_date,earliest_positive_test_date,CDC_N_COV_2019_SEQUENCE_REASON) %>% 
    arrange(days_between_collection_and_earliest_test) %>% distinct()
  
  # collection date well after earliest date, mostly reinfections
  ggplot(tmp)+stat_ecdf(aes(x=days_between_collection_and_earliest_test)) +
    geom_vline(aes(xintercept=21),linetype='dashed')+
    scale_y_continuous(trans='log10') + ylab('cumulative distribution')
  
  d %>% filter(d$days_between_collection_and_earliest_test>21) %>%
    select(CASE_ID,days_between_collection_and_earliest_test,collection_date,earliest_positive_test_date,CDC_N_COV_2019_SEQUENCE_REASON) %>% 
    arrange(days_between_collection_and_earliest_test) %>% distinct() %>%
    write.table(file='collection_21+_days_after_earliest_test.csv',sep=',',row.names = FALSE)


# define infection_type to track possible reinfection found during sentinel surveillance
d$infection_type <- NA

# define best available infection event date
  
  # first, assume it's the earliest positive test date
  d$best_infection_event_date <- d$earliest_positive_test_date
  
  # if earliest test is way after collection, use collection date and label sample reinfection
  idx<-compareNA(d$days_between_collection_and_earliest_test, -21, '<')
  d$best_infection_event_date[idx] <- d$collection_date[idx]
  d$days_between_collection_and_earliest_test[idx] <- 0
  d$infection_type[idx] <- 'suspected reinfection'
  
  # if earliest test is a little after collection, assume collection date
  idx<-compareNA(d$days_between_collection_and_earliest_test, -21, '>=') &
    compareNA(d$days_between_collection_and_earliest_test, 0, '<')
  d$best_infection_event_date[idx] <- d$collection_date[idx]
  d$days_between_collection_and_earliest_test[idx] <- 0

  # if earliest test is way before collection, use collection date and assume suspected reinfection
  idx<-compareNA(d$days_between_collection_and_earliest_test, 21, '>')
  d$best_infection_event_date[idx] <- d$collection_date[idx]
  d$days_between_collection_and_earliest_test[idx] <- 0
  d$infection_type[idx] <- 'suspected reinfection'
  
  # for NA earliest dates, assume collection date
  idx <- is.na(d$best_infection_event_date)
  d$best_infection_event_date[idx] <- d$collection_date[idx]
  d$days_between_collection_and_earliest_test[idx] <- 0
  

# vaccine type per shot
d$first_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1==207] <- 'Moderna'
d$first_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1==208] <- 'Pfizer/BioNTech'
d$first_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1==217] <- 'Pfizer/BioNTech'
d$first_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1==212] <- 'J&J'

d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==207] <- 'Moderna'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==208] <- 'Pfizer/BioNTech'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==217] <- 'Pfizer/BioNTech'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==212] <- 'J&J'

d$third_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_3==207] <- 'Moderna'
d$third_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_3==208] <- 'Pfizer/BioNTech'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==217] <- 'Pfizer/BioNTech'
d$third_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_3==212] <- 'J&J'

# drop rows with unknown vaccine type 
### Miguel's thoughts: if we're just doing vaccinated/unvaccinated, maybe we can include these people back in?  
d <- d %>% filter(!compareNA(IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1,213)) %>%
  filter(!compareNA(IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2,213))
exclusions <- exclusions %>% rbind(data.frame(data_view='known vaccine',reason='unknown vaccine type',n_kept=nrow(d)))



  
# mixed vaccine brands
  d %>% filter(first_shot!=second_shot) %>% nrow()
  d %>% filter(first_shot!=second_shot)
  
  d$vaccine_brand=d$first_shot
  d$vaccine_brand[d$first_shot!=d$second_shot] <- 'Mixed'
  d$vaccine_brand[(d$first_shot == "J&J") & (d$second_shot %in% c('Pfizer/BioNTech','Moderna'))] <- "J&J_mRNA_booster"
  
  d$vaccine_brand[is.na(d$vaccine_brand)] <- 'None'
  d$vaccine_brand <- factor(d$vaccine_brand,levels=c('Moderna','Pfizer/BioNTech','J&J','J&J_mRNA_booster', 'Mixed','None'))
    

# vaccine "active" date
# for analyses of vaccine effectiveness and breakthru, we classify vaccine as 
# "active" 21 days after administration assuming 14 days for protection against symptoms and 7 more days 
#for protection against hospitalization
d$date_first_shot_active_date <- d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 + 21
d$date_second_shot_active_date <- d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 + 21
d$date_third_shot_active_date <- d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_3 + 21

  
# vaccine active binaries. was vaccination active at date of collection?
d$first_shot_active <- c('No','Yes')[1+as.numeric((d$date_first_shot_active_date - d$best_infection_event_date)<=0)]
d$second_shot_active <- c('No','Yes')[1+as.numeric((d$date_second_shot_active_date - d$best_infection_event_date)<=0)]
d$third_shot_active <- c('No','Yes')[1+as.numeric((d$date_third_shot_active_date - d$best_infection_event_date)<=0)]

# is this a first shot breakthru or a second shot breakthru?
d$first_shot_breakthru <- c('No','Yes')[1+as.numeric((d$first_shot_active=='Yes') & !compareNA(d$second_shot_active,'Yes'))]
d$second_shot_breakthru <- c('No','Yes')[1+as.numeric((d$second_shot_active=='Yes') )]
d$third_shot_breakthru <- c('No','Yes')[1+as.numeric((d$third_shot_active=='Yes') )]

# doses active at best_infection_event_date
d$doses_active <- 0
d$doses_active[compareNA(d$first_shot_active,'Yes')&!compareNA(d$second_shot_active,'Yes')] <-1
d$doses_active[compareNA(d$second_shot_active,'Yes')&!compareNA(d$second_shot_active,'Yes')] <-2
d$doses_active[compareNA(d$third_shot_active,'Yes')] <-3

# vaccine brand and dose active at time of collection

  # define brand-dose interaction
  d <- d %>% mutate(active_vaccine_brand_dose = as.character(interaction(vaccine_brand,doses_active,sep=' - dose ')))
  # set all 0-dose to none
  d$active_vaccine_brand_dose[grepl(' - dose 0',d$active_vaccine_brand_dose)] <- 'None'
  # set mixed:1 to whatever the first dose was
  idx <- grepl('Mixed -dose 1',d$active_vaccine_brand_dose)
  d$active_vaccine_brand_dose[idx] <- paste(d$first_shot[idx]," - dose 1",sep='')
  
  # order factors
  d$active_vaccine_brand_dose <- factor(d$active_vaccine_brand_dose,
                                        levels=c('None',"Pfizer/BioNTech - dose 1","Moderna - dose 1",
                                                 "J&J - dose 1","Pfizer/BioNTech - dose 2","Moderna - dose 2","Mixed - dose 2", "J&J - dose 2", "J&J_mRNA_booster - dose 2","Pfizer/BioNTech - dose 3","Moderna - dose 3","Mixed - dose 3" ))  
  
  
d$vaccination_active <- NA
d$vaccination_active <- "No Vaccination to \n <21 days post dose one"
d$vaccination_active[compareNA(d$first_shot_active,'Yes')&!compareNA(d$third_shot_active,'Yes')] <- "≥21 days post dose one to \n <21 days post booster"
d$vaccination_active[compareNA(d$third_shot_active,'Yes')] <- "≥21 days post booster"
d$vaccination_active[(d$vaccine_brand == "J&J_mRNA_booster")&compareNA(d$second_shot_active,'Yes')] <- "≥21 days post booster"
d$vaccination_active <- factor(d$vaccination_active, levels = c("No Vaccination to \n <21 days post dose one", "≥21 days post dose one to \n <21 days post booster", "≥21 days post booster")) 

with(d, table(vaccination_active,active_vaccine_brand_dose, useNA = "ifany")) #check recoding

# exclude active mixed vaccinations since there are few and results can be misleading
##actually let's keep these in for now based on a reviewer comment that these shouldnt be different if we're only doing collapsed categories
# d <- d %>% filter(!grepl('Mixed',active_vaccine_brand_dose)) %>% droplevels()
# exclusions <- exclusions %>% rbind(data.frame(data_view='defined vaccinations active at time of sample (excluding mixed vaccine)',
#                                               reason='exclude mixed vaccinations since there are very few and shrinkage estimators may be pooled misleading close to the most-no-vaccine mean',
#                                               n_kept=nrow(d)))
  
# figure out which sample(s) to use for multiples
# SPECIMEN__COLLECTION__DTTM is most recent
  # multiples
  tmp <- d %>% group_by(CASE_ID) %>% summarize(n=n()) %>%
    filter(n>1)
  
  # among multiples, which ones have the same lineage?
  tmp <- d %>% group_by(CASE_ID) %>% mutate(n=n()) %>% ungroup() %>% 
    filter(n>1) %>%
    group_by(CASE_ID) %>% 
    mutate(count_lineages = length(unique(who_lineage))) %>%
    group_by(CASE_ID) %>% 
    mutate(days_apart = max(best_infection_event_date)-min(best_infection_event_date))
  
  # almost all have same lineage
  mean(tmp$count_lineages==1)
  
  # hand label
  multiple_lineage_cases = tmp %>% filter(count_lineages>1) %>% select(CASE_ID) %>% 
    distinct()
  
  # exclude suspected coinfections labeled earlier
  multiple_lineage_cases <- multiple_lineage_cases %>% 
    filter(!(CASE_ID %in% d$CASE_ID[!is.na(d$infection_type)]))
 
  idx <- d$CASE_ID %in% multiple_lineage_cases$CASE_ID 
  d$infection_type[idx] <- 'possible coinfection or miscalled lineage'
  
  # included suspected reinfections
  idx <- d$CDC_N_COV_2019_SEQUENCE_REASON == 'SUSPECTED REINFECTION'
  d$infection_type[idx] <- 'suspected reinfection'
  
  d$infection_type[is.na(d$infection_type)] <- 'monoinfection'
  
# filter distinct, keeping first one
d <- d %>% filter(!(CASE_ID %in% multiple_lineage_cases$CASE_ID))

exclusions <- exclusions %>% rbind(data.frame(data_view='only single lineage cases',reason='a small number of cases have multiple lineages for various reasons. But since number is small we just exclude',n_kept=nrow(d)))

##### intervals and censoring for survival analysis
# interval from sample to outcome
d$hosp_days_at_risk <- as.numeric(d$admitdate-d$best_infection_event_date)

# most data is clustered within reasonable infection dynamics timescales of the 
# sample time, which is good
hist(d$hosp_days_at_risk)
mean(d$hosp_days_at_risk>=0,na.rm=TRUE)
mean(d$hosp_days_at_risk>=-14,na.rm=TRUE)

# exclude implausibly early admission dates, where this sequence may not represent
# the infection that lead to hospitalization
d <- d %>% filter(!compareNA(hosp_days_at_risk,-14,test='<'))


exclusions <- exclusions %>% rbind(data.frame(data_view='sequence collection date within 14 days of admission date if hospitalized',
                                              reason='exclude hospitalizations where it is unlikely that this sequence is responsible for the hospitalization',
                                              n_kept=nrow(d)))

## get earliest sample for cases with multiple sequences
d <- d %>%
  group_by(CASE_ID) %>%
  arrange(hosp_days_at_risk) %>%
  slice_max(hosp_days_at_risk,with_ties=FALSE) %>%
  ungroup()

hist(d$hosp_days_at_risk)
mean(d$hosp_days_at_risk>=0,na.rm=TRUE)
mean(d$hosp_days_at_risk>=-14,na.rm=TRUE)

exclusions <- exclusions %>% rbind(data.frame(data_view='first detection in cases with multiple sequences of a single lineage',
                                              reason='no multiple counting and no value in modeling multiple sequences per case',
                                              n_kept=nrow(d)))


#creating datasets for sens analysis where we limit time between collection and hosp to 14,21,30 days
hist(d$hosp_days_at_risk)
mean(d$hosp_days_at_risk>=0,na.rm=TRUE)
mean(d$hosp_days_at_risk>=-14,na.rm=TRUE)

d %>% filter(!is.na(hosp_days_at_risk)) %>% summarize(range = range(hosp_days_at_risk))
d %>% filter(!is.na(hosp_days_at_risk)) %>% count()



d_14 <- d %>% filter(!compareNA(hosp_days_at_risk,14,test='>'))
hist(d_14$hosp_days_at_risk)
d_14 %>% filter(!is.na(hosp_days_at_risk)) %>% summarize(range = range(hosp_days_at_risk))
d_14_count <- d_14 %>% filter(!is.na(hosp_days_at_risk)) %>% count()

d_21 <- d %>% filter(!compareNA(hosp_days_at_risk,21,test='>'))
hist(d_21$hosp_days_at_risk)
d_21 %>% filter(!is.na(hosp_days_at_risk)) %>% summarize(range = range(hosp_days_at_risk))
d_21_count <- d_21 %>% filter(!is.na(hosp_days_at_risk)) %>% count()


d_30 <- d %>% filter(!compareNA(hosp_days_at_risk,30,test='>'))
hist(d_30$hosp_days_at_risk)
d_30 %>% filter(!is.na(hosp_days_at_risk)) %>% summarize(range = range(hosp_days_at_risk))
d_30_count <- d_30 %>% filter(!is.na(hosp_days_at_risk)) %>% count()


##### NEED TO UPDATE
# format days at risk for people who haven't been hospitalized 
no_hosp_idx <- is.na(d$hosp_days_at_risk)
d$hosp_days_at_risk[no_hosp_idx] <- as.Date('01/05/2022', format = "%m/%d/%Y")-d$collection_date[no_hosp_idx]


# start time at risk in the 14 days preceding hospitalization or sample collection
d$hosp_days_at_risk <- d$hosp_days_at_risk + 14
hist(d$hosp_days_at_risk)
hist(d$hosp_days_at_risk[d$mhosp=='Yes'])

weird_dates <- d %>% filter(compareNA(hosp_days_at_risk,0,test='<'))


###### add analysis type label field
exclusions$analysis <- 'all'


## write out data process
write.table(exclusions,'output/sample_size_and_exclusions_summary.csv',sep=',',row.names=FALSE)



#doing the same for the 14,21,30 day datasets

# format days at risk for people who haven't been hospitalized 
no_hosp_idx_14 <- is.na(d_14$hosp_days_at_risk)
d_14$hosp_days_at_risk[no_hosp_idx_14] <- as.Date('01/05/2022', format = "%m/%d/%Y")-d_14$collection_date[no_hosp_idx_14]


# start time at risk in the 14 days preceding hospitalization or sample collection
d_14$hosp_days_at_risk <- d_14$hosp_days_at_risk + 14
hist(d_14$hosp_days_at_risk)
hist(d_14$hosp_days_at_risk[d_14$mhosp=='Yes'])

# format days at risk for people who haven't been hospitalized 
no_hosp_idx_21 <- is.na(d_21$hosp_days_at_risk)
d_21$hosp_days_at_risk[no_hosp_idx_21] <- as.Date('01/05/2022', format = "%m/%d/%Y")-d_21$collection_date[no_hosp_idx_21]


# start time at risk in the 14 days preceding hospitalization or sample collection
d_21$hosp_days_at_risk <- d_21$hosp_days_at_risk + 14
hist(d_21$hosp_days_at_risk)
hist(d_21$hosp_days_at_risk[d_21$mhosp=='Yes'])

# format days at risk for people who haven't been hospitalized 
no_hosp_idx_30 <- is.na(d_30$hosp_days_at_risk)
d_30$hosp_days_at_risk[no_hosp_idx_30] <- as.Date('01/05/2022', format = "%m/%d/%Y")-d_30$collection_date[no_hosp_idx_30]


# start time at risk in the 14 days preceding hospitalization or sample collection
d_30$hosp_days_at_risk <- d_30$hosp_days_at_risk + 14
hist(d_30$hosp_days_at_risk)
hist(d_30$hosp_days_at_risk[d_30$mhosp=='Yes'])
