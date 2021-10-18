# Script to prepare data for analysis: severity and vaccine effectiveness by variant
# for sequenced case data in Washington state.
#
# Data curated by WA-DOH from public health surveillance records and GISAID sequence data.
# Variant calls using pangolin v1.2.13. 
#

library(tidyverse)
library(lubridate)

# load data
raw <- read.delim('data_pull_2021-09-02_prepped.csv',sep=',',header = TRUE)
# raw <-readRDS("data_pull_2021-09-02_prepped.RDS")
names(raw)



######################################################
### set up workspace and initialize helper functions
######################################################


# https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
# https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/variant-surveillance/variant-info.html

WHO_names <- list(B.1.1.7 = 'Alpha (B.1.1.7)',
                  Q.1 = 'Alpha (B.1.1.7)',
                  Q.2 = 'Alpha (B.1.1.7)',
                  Q.3 = 'Alpha (B.1.1.7)',
                  Q.4 = 'Alpha (B.1.1.7)',
                  Q.8 = 'Alpha (B.1.1.7)',
                  B.1.351 = 'Beta (B.1.351)',
                  B.1.351.1='Beta (B.1.351)',
                  B.1.351.2='Beta (B.1.351)',
                  B.1.351.3='Beta (B.1.351)', #all 1.351 sublineages assigned to parent based on CDC 
                  P.1='Gamma (P.1)',
                  P.1.1='Gamma (P.1)', # assigned to P.1 since there are few that show up in sensitivity analysis   
                  P.1.2='Gamma (P.1)',
                  P.1.3='Gamma (P.1)',
                  P.1.10='Gamma (P.1)',
                  P.1.11='Gamma (P.1)',
                  B.1.617.2='Delta (B.1.617.2)',
                  AY.1='Delta (B.1.617.2)',  
                  AY.2='Delta (B.1.617.2)',
                  AY.3='Delta (B.1.617.2)',
                  AY.3.1='Delta (B.1.617.2)',
                  AY.4='Delta (B.1.617.2)',
                  AY.5='Delta (B.1.617.2)',
                  AY.6='Delta (B.1.617.2)',
                  AY.7='Delta (B.1.617.2)',
                  AY.7.1='Delta (B.1.617.2)',
                  AY.8='Delta (B.1.617.2)',
                  AY.9='Delta (B.1.617.2)',
                  AY.10='Delta (B.1.617.2)',
                  AY.11='Delta (B.1.617.2)',
                  AY.12='Delta (B.1.617.2)',
                  AY.13='Delta (B.1.617.2)',
                  AY.14='Delta (B.1.617.2)',
                  AY.15='Delta (B.1.617.2)',
                  AY.16='Delta (B.1.617.2)',
                  AY.17='Delta (B.1.617.2)',
                  AY.18='Delta (B.1.617.2)',
                  AY.19='Delta (B.1.617.2)',
                  AY.20='Delta (B.1.617.2)',
                  AY.21='Delta (B.1.617.2)',
                  AY.22='Delta (B.1.617.2)',
                  AY.23='Delta (B.1.617.2)',
                  AY.24='Delta (B.1.617.2)',
                  AY.25='Delta (B.1.617.2)',
                  B.1.617.1='Kappa (B.1.617.1)',
                  B.1.427='Epsilon (B.1.427/B.1.429)',
                  B.1.429='Epsilon (B.1.427/B.1.429)',
                  `B.1.427/429`='Epsilon (B.1.427/B.1.429)',
                  
                  B.1.526='Iota (B.1.526)',
                  B.1.525='Eta (B.1.525)',
                  C.37='Lambda (C.37)',
                  
                  # if you want to group these sublineages with the parent,
                  # just change the renaming they are assigned to the parent name above
                  other_voc_voi='other_voc_voi', # this is here for the colormap
                  other='other',
                  P.2='other',
                  P.3='other') #P.2 and P.3 assigned to non-variant due to lack of VOI/VOC classification by WHO or lack of detection in WA state


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

## removing duplicates 
d <- d %>% group_by(CASE_ID, CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER) %>% unique() %>% ungroup()
exclusions <- exclusions %>% rbind(data.frame(data_view='No duplicates',reason='duplicates were created in the data pulling',n_kept=nrow(d)))

sum((d$lineage=='None'),na.rm=TRUE)
sum(is.na(d$lineage=='None'),na.rm=TRUE)

# Excluding poor quality seqs
exclude_seqs <- read.delim('exclusion_ids_lf.csv',sep=',',header = TRUE)

## which ones have bad sequence quality in raw data
list_bad_quality <- raw %>% filter((raw$CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER %in% exclude_seqs$ID))
sum((d$lineage=='None'),na.rm=TRUE)
sum(is.na(d$lineage=='None'),na.rm=TRUE)

## drop those with lineage=='None' or no lineage (n=223)
d <- d %>% filter( !(lineage=='None' | is.na(lineage)))
exclusions <- exclusions %>% rbind(data.frame(data_view='with_known_lineage',reason='filtered out lineage==None',n_kept=nrow(d)))

# which ones have bad sequence quality now (n=96)
list_bad_quality_afternone <- d %>% filter((d$CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER %in% exclude_seqs$ID))

## exclude those with bad sequence quality
d <- d %>%  filter(!(d$CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER %in% exclude_seqs$ID))
exclusions <- exclusions %>% rbind(data.frame(data_view='good sequence quality',reason='>10% ambiguity in sequencing',n_kept=nrow(d)))
list_bad_quality_bad_seq <- d %>% filter((d$CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER %in% exclude_seqs))

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


# age bin
# include >99 in the 90-99 age bracket 
d$age_bin = factor(5+floor(pmin(99,d$AGE_YEARS)/10)*10) 


# drop miscoded age
d <- d %>% filter(AGE_YEARS!=141) %>% droplevels()
exclusions <- exclusions %>% rbind(data.frame(data_view='valid ages',reason='miscoded age',n_kept=nrow(d)))


# dates
d$collection_date <- as.Date(d$collection_date)
d$admitdate <- as.Date(d$admitdate)
d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 <- as.Date(d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1)
d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 <- as.Date(d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2)
d$earliest_positive_test_date <- as.Date(d$earliest_positive_test_date)

d$Discharge_Date_Time[d$Discharge_Date_Time=='none'] <- NA
d$Discharge_Date_Time <- as.Date(d$Discharge_Date_Time)

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

# sex unknown
d$SEX_AT_BIRTH <- as.character(d$SEX_AT_BIRTH)
d$SEX_AT_BIRTH[is.na(d$SEX_AT_BIRTH)] <- 'Unknown'
d$SEX_AT_BIRTH <- relevel(factor(d$SEX_AT_BIRTH),ref='Female')

# clean up test date vs collection date
d$days_between_collection_and_earliest_test <- as.numeric(d$collection_date - d$earliest_positive_test_date)

  # collection date before earliest test date. This shouldn't happen....
  d %>% filter(d$days_between_collection_and_earliest_test<0) %>%
    select(CASE_ID,days_between_collection_and_earliest_test,collection_date,earliest_positive_test_date,sequence_reason_clean) %>% 
    arrange(days_between_collection_and_earliest_test) %>% distinct() %>%
    write.table(file='collection_before_earliest_test.csv',sep=',',row.names = FALSE)
  
  # collection date before earliest test date. This shouldn't happen....
  tmp <- d %>% filter(d$days_between_collection_and_earliest_test>=0) %>%
    select(CASE_ID,days_between_collection_and_earliest_test,collection_date,earliest_positive_test_date,sequence_reason_clean) %>% 
    arrange(days_between_collection_and_earliest_test) %>% distinct()
  
  # collection date well after earliest date, mostly reinfections
  ggplot(tmp)+stat_ecdf(aes(x=days_between_collection_and_earliest_test)) +
    geom_vline(aes(xintercept=21),linetype='dashed')+
    scale_y_continuous(trans='log10') + ylab('cumulative distribution')
  
  d %>% filter(d$days_between_collection_and_earliest_test>21) %>%
    select(CASE_ID,days_between_collection_and_earliest_test,collection_date,earliest_positive_test_date,sequence_reason_clean) %>% 
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
d$first_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1==212] <- 'J&J'

d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==207] <- 'Moderna'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==208] <- 'Pfizer/BioNTech'
d$second_shot[d$IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2==212] <- 'J&J'

# drop rows with unknown vaccine type
d <- d %>% filter(!compareNA(IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_1,213)) %>%
  filter(!compareNA(IIS_VACCINE_INFORMATION_AVAILABLE_ADMINISTERED_2,213))
exclusions <- exclusions %>% rbind(data.frame(data_view='known vaccine',reason='unknown vaccine type',n_kept=nrow(d)))


# clean up 7 people got third dose
  sum(d$IIS_EVER_RECEIVED_VACCINE_NUM_DOSES==3,na.rm=TRUE)

  #what vaccines did they get for the first two?
  d %>% filter(compareNA(d$IIS_EVER_RECEIVED_VACCINE_NUM_DOSES,3)) %>% 
    select(first_shot,second_shot)
  # all got mRNA
  
  # gonna lump these into 2 dose mRNA
  d$doses_received = d$IIS_EVER_RECEIVED_VACCINE_NUM_DOSES
  d$doses_received[compareNA(d$IIS_EVER_RECEIVED_VACCINE_NUM_DOSES,3)]=2
  d$doses_received[is.na(d$doses_received)] <- 0
  
# mixed vaccine brands
  d %>% filter(first_shot!=second_shot) %>% nrow()
  d %>% filter(first_shot!=second_shot)
  
  d$vaccine_brand=d$first_shot
  d$vaccine_brand[d$first_shot!=d$second_shot] <- 'Mixed'
  
  d$vaccine_brand[is.na(d$vaccine_brand)] <- 'None'
  d$vaccine_brand <- factor(d$vaccine_brand,levels=c('Moderna','Pfizer/BioNTech','J&J','Mixed','None'))
    
# collapsed vaccine type
  d$vaccine_type <- as.character(d$vaccine_brand)
  d$vaccine_type[d$vaccine_type %in% c('Pfizer/BioNTech','Moderna')] <- 'mRNA'
  d$vaccine_type[d$first_shot %in% c('Pfizer/BioNTech','Moderna') &
                   d$second_shot %in% c('Pfizer/BioNTech','Moderna')] <- 'mRNA'
  d$vaccine_type <- relevel(factor(d$vaccine_type),'mRNA')
  
# vaccine "active" date
# for analyses of vaccine effectiveness and breakthru, we classify vaccine as 
# "active" 21 days after administration assuming 14 days for protection against symptoms and 7 more days 
#for protection against hospitalization
d$date_first_shot_active_date <- d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_1 + 21
d$date_second_shot_active_date <- d$IIS_VACCINE_INFORMATION_AVAILABLE_DATE_2 + 21
  
# vaccine active binaries. was vaccination active at date of collection?
d$first_shot_active <- c('No','Yes')[1+as.numeric((d$date_first_shot_active_date - d$best_infection_event_date)<=0)]
d$second_shot_active <- c('No','Yes')[1+as.numeric((d$date_second_shot_active_date - d$best_infection_event_date)<=0)]

# is this a first shot breakthru or a second shot breakthru?
d$first_shot_breakthru <- c('No','Yes')[1+as.numeric((d$first_shot_active=='Yes') & !compareNA(d$second_shot_active,'Yes'))]
d$second_shot_breakthru <- c('No','Yes')[1+as.numeric((d$second_shot_active=='Yes') )]

# doses active at best_infection_event_date
d$doses_active <- 0
d$doses_active[compareNA(d$first_shot_active,'Yes')&!compareNA(d$second_shot_active,'Yes')] <-1
d$doses_active[compareNA(d$second_shot_active,'Yes')] <-2
 
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
                                                 "J&J - dose 1","Pfizer/BioNTech - dose 2","Moderna - dose 2","Mixed - dose 2"))  
  
  
d$vaccination_active <- NA
d$vaccination_active <- "No Vaccination to \n <21 days post dose one"
d$vaccination_active[compareNA(d$first_shot_active,'Yes')] <- "≥21 days post dose one"
d$vaccination_active <- factor(d$vaccination_active, levels = c("No Vaccination to \n <21 days post dose one", "≥21 days post dose one")) 

# exclude active mixed vaccinations since there are few and results can be misleading
d <- d %>% filter(!grepl('Mixed',active_vaccine_brand_dose)) %>% droplevels()
exclusions <- exclusions %>% rbind(data.frame(data_view='defined vaccinations active at time of sample (excluding mixed vaccine)',
                                              reason='exclude mixed vaccinations since there are very few and shrinkage estimators may be pooled misleading close to the most-no-vaccine mean',
                                              n_kept=nrow(d)))
  
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

  # for ones with same lineage, mostly during obviously the same infection
  # view(tmp %>% filter(count_lineages==1) %>%
  #        group_by(days_apart) %>% summarize(n=n()) )

  # also ones with multiple lineages look like the same infection for the most part
  # view(tmp %>% filter(count_lineages>1) %>%
  #        group_by(days_apart) %>% summarize(n=n()) )
  
  
  # what's going on with multiple lineages in short time?
  # tmp %>% filter(count_lineages>1) %>%
  #    select(CASE_ID,lineage,collection_date,sequence_reason_clean) %>%
  #    view()
  # these look like a mix of coinfections, reinfections, and pango calling related lineages differently
  # hand label
  multiple_lineage_cases = tmp %>% filter(count_lineages>1) %>% select(CASE_ID) %>% 
    distinct()
  
  # exclude suspected coinfections labeled earlier
  multiple_lineage_cases <- multiple_lineage_cases %>% 
    filter(!(CASE_ID %in% d$CASE_ID[!is.na(d$infection_type)]))
 
  idx <- d$CASE_ID %in% multiple_lineage_cases$CASE_ID 
  d$infection_type[idx] <- 'possible coinfection or miscalled lineage'
  
  # included suspected reinfections
  idx <- d$sequence_reason_clean == 'SUSPECTED REINFECTION'
  d$infection_type[idx] <- 'suspected reinfection'
  
  d$infection_type[is.na(d$infection_type)] <- 'monoinfection'
  
  

  
  # check that multiples have same metadata
  # looks good
  # tmp %>% filter(count_lineages>1) %>% 
  #   select(CASE_ID,mhosp,DEATH_DATE,ICU_STATUS,MECHANICAL_VENTILATION,age_bin,race) %>%
  #   view()
  
  
# since only a small amount look like possible reinfections or miscalled lineages, we just drop them

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

# format days at risk for people who haven't been hospitalized
no_hosp_idx <- is.na(d$hosp_days_at_risk)
d$hosp_days_at_risk[no_hosp_idx] <- as.Date('2021-07-23')-d$collection_date[no_hosp_idx]

# I can't easily do left censored analyses for hospitalizations that precede the collection date
# but I think it is reasonable to approximate the start of exposure as roughly when infection could start
# results are insensitive to this choice
# start time at risk in the 14 days preceding hospitalization or sample collection
d$hosp_days_at_risk <- d$hosp_days_at_risk + 14
hist(d$hosp_days_at_risk)
hist(d$hosp_days_at_risk[d$mhosp=='Yes'])


## let's only use data through July 23
d <- d %>% filter(best_infection_event_date <= '2021-07-23')

exclusions <- exclusions %>% rbind(data.frame(data_view='with best date on or before July 23',
                                              reason='very incomplete sequencing and outcomes after that',
                                              n_kept=nrow(d)))


###### add analysis type label field
exclusions$analysis <- 'all'


## write out data process
write.table(exclusions,'output/sample_size_and_exclusions_summary.csv',sep=',',row.names=FALSE)



