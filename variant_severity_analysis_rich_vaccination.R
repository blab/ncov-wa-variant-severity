# Script to analyze measures of severity and vaccine effectiveness by variant
# for sequenced case data in Washington state.
# 
# Methods: hierarchical Cox proportional hazards regression
# 
# PRIMARY ANALYSIS. Case-Hospitalization ratio in sentinel surveillance
# looking at vaccine type, dose, and vaccine-variant interactions
#

# MUST RUN variant_severity_data_prep.R first to set up workspace

library(coxme)
library(survival)
library(lme4)
library(kableExtra)
library(survminer)
library(lmtest)


# package:merTools is also required in a function, but not loaded into global namespace since it conflicts with dplyr::Select



# if you are working with this but don't have the raw data, you can load
# saved cox models, and saved who_lineage hazard ratio parameter data frames here.
# To use, comment out sections that don't run (data organization and model fitting),
# and otherwise regenerate figures, or make new ones. 
# For the cox models, because the model object is available, it is also possible to consider
# new contrasts, or make figures or tables describing control variables. 
# load(file='variant_severity_simple_vaccination_cox_models.Rdata')




########################
### helper functions ###
########################

coxme_random_params <- function(mod,data,group='who_lineage'){
  # formats coefficient estimates from coxme_object (hierarchical model)
  
  res <- data.frame(logRR = ranef(mod)[[group]])
  res$logRR <- res$logRR-res$logRR[1] # set reference to zero
  
  res[group]=sub(group,'',row.names(res))
  row.names(res)<-NULL
  
  # find relevant block of variance-covariance matrix
  # random effects are first, then fixed https://cran.r-project.org/web/packages/coxme/coxme.pdf
  # order should be order they appear in the formula
  mod_vcov <- mod$variance
  idx_type <- rep(names(mod$frail),lengths(mod$frail))
  group_idx <- which(idx_type==group)
  mod_vcov <- mod_vcov[group_idx,group_idx]
  
  res$sd <- diag(mod_vcov)
  # propagate reference variance, variance of difference of two random variable
  res$sd[2:nrow(res)]<- res$sd[2:nrow(res)] + 
    res$sd[1] - 2*c(as.matrix(mod_vcov)[(row(mod_vcov) == (col(mod_vcov) - 1))])[1:nrow(res)-1]
  res$sd[1]<-0
  res$sd <- sqrt(res$sd)
  res$lower95 <- res$logRR -2*res$sd
  res$upper95 <- res$logRR +2*res$sd
  
  res <- res %>% 
    filter(!!as.symbol(group) != levels(data[[group]])[1]) %>%
    arrange(logRR) 
  res[[group]] <- factor(res[[group]],levels=res[[group]])
  
  return(res)
}

coxph_params <- function(mod,ref,group='who_lineage'){
  # formats coefficient estimates from coxph_object (fixed effects)
  
  res <- data.frame(confint(mod))
  res$logRR <- coef(mod)
  
  res <- res %>% filter(grepl(group,rownames(res)))
  res[[group]]=factor(sub(group,'',row.names(res)))
  row.names(res)<-NULL
  names(res) <- c('lower95','upper95','logRR',group)
  
  res[[group]] <- factor(res[[group]],levels=levels(ref[[group]]))
  
  return(res)
}


glmer_random_params <- function(mod,ref, group='who_lineage'){
  # format coefficients from poisson regression (hierarchical model)
  
  if(group=='who_lineage'){
    ref_set =  c('other',levels(ref$who_lineage))
  } else if (group == 'vaccination_active'){
    ref_set =  c('No',levels(ref$vaccination_active))
  }
  res<-merTools::REextract(mod)
  res <- res %>% filter(rownames(res) %in% ref_set)
  res$`(Intercept)`<-res$`(Intercept)`-res$`(Intercept)`[1] # relevel 'other'
  res$logRR <- res$`(Intercept)`
  
  res$sd <- res$`(Intercept)_se`
  # propagate reference level variance
  # covariance is missing so conservatively, sum
  res$sd[2:nrow(res)]<- res$sd[2:nrow(res)]^2 + res$sd[1]^2 
  res$sd[1]<-0
  res$sd <- sqrt(res$sd)
  res$lower95 <- res$logRR -2*res$sd
  res$upper95 <- res$logRR +2*res$sd
  
  res <- res %>% select(groupID,logRR,sd,lower95,upper95)
  names(res)[1]<-group
  
  res <- res %>% 
    filter(!!as.symbol(group) != ref_set[1]) 
  res[[group]] <- factor(res[[group]],levels=ref_set[2:length(ref_set)])
  
  return(res)
}


######################################
#### hospitalization analyses ########
######################################

names(d)

############# PRIMARY ANALYSIS. Case-Hospitalization ratio in sentinel surveillance
############# hierarchical cox model
############# vaccination with rich model of vaccine type and variant

# hospital sentinel only cox hierarchical model
## added in excludsion for new reinfection flag since that only came into play post setp 2021
cox_dat <- d %>% 
  filter(sequence_reason_clean=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>% 
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()


# track which data went into this analysis
exclusions <- exclusions %>% rbind(data.frame(data_view='all sentinel surveillance and lineages with at least 1 hospitalization but excluding suspected reinfections',
                                              reason='primary analysis',
                                              n_kept=nrow(cox_dat),
                                              analysis='cox hierarchical sentinel'))


hosp_surv <- Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes'))

cox_sentinel <- coxme(hosp_surv ~ (1|who_lineage) + 
                        age_bin + (1|SEX_AT_BIRTH) +
                         (1|vaccination_active) + week_collection_number, 
                      data=cox_dat,
                      x=FALSE,y=FALSE)
summary(cox_sentinel)

cox_sentinel_lineage_params <- coxme_random_params(cox_sentinel,cox_dat,group='who_lineage')


ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95,color=who_lineage)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16)),
                     limits=log(c(1/2,9))) +
  scale_color_manual(values=cmap,guide=FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk.svg',units='in',width=5,height=3,device='svg')

# risk reduction of hospitalization from vaccine
cox_sentinel_vaccine_params <- coxme_random_params(cox_sentinel,cox_dat,group='vaccination_active')

ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=vaccination_active,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/16,1/8,1/4,1/2,1,2,4,8,16)),
                     labels=(c(1/16, 1/8,1/4,1/2,1,2,4,8,16)),
                     limits=log(c(1/16,2))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk.svg',units='in',width=5,height=3,device='svg')


cox_sentinel_vax_rel_risk=data.frame(vaccine_type=cox_sentinel_vaccine_params$vaccination_active,
                                     rel_risk_hosp_given_case=round(exp(cox_sentinel_vaccine_params$logRR),2),
                                     lower95 = round(exp(cox_sentinel_vaccine_params$lower95),2),
                                     upper95 = round(exp(cox_sentinel_vaccine_params$upper95),2))

cox_sentinel_vax_rel_risk %>%
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  save_kable(file = "output/rich_vaccination/variant_severity_simple_vaccine_relative_risk_hospitalization.png",
             density=600,zoom=3) 
write.table(cox_sentinel_vax_rel_risk,'output/rich_vaccination/variant_severity_simple_vaccine_relative_risk_hospitalization.csv',sep=',',row.names = FALSE)

##Likelihood ratio test for global effect:

#we already have out main model: cox_sentienl

#now create a reduced model with df -1 by removing the effect that we're interested in, which is variants
cox_sentinel_test <- coxme(hosp_surv ~  age_bin + (1|SEX_AT_BIRTH) +
                             (1|vaccination_active) + week_collection_number, 
                           data=cox_dat,
                           x=FALSE,y=FALSE)

#now conduct a LR test

lrtest(cox_sentinel, cox_sentinel_test)

## adding in Kaplan meier curves 
#### NEED TO HARDCODE OMICRON IN WHEN WE UPDATE DATA
survdiff(Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes')) ~ who_lineage, data= cox_dat)

ggsurv <- ggsurvplot(
  fit = survfit(Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes')) ~ who_lineage, data= cox_dat), 
  ylim = c(0.925, 1),
  risk.table = TRUE, risk.table.col = "who_lineage", risk.table.height = 0.25,
  
  legend.labs = c("Ancestral", "Alpha", 'Beta', "Gamma", "Delta", "Kappa", "Epsilon", "Iota", "Eta", "Lambda"))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 10, face = "bold"),legend.key.height = unit(1, 'cm'), legend.key.width = unit(1, 'cm'))
ggsurv

##### Control variables
# age
cox_sentinel_age_params <- data.frame(age_bin=sub('age_bin','',names(coef(cox_sentinel))),
                                      logRR = coef(cox_sentinel))
cox_sentinel_age_params$lower95<-confint(cox_sentinel)[,1]
cox_sentinel_age_params$upper95<-confint(cox_sentinel)[,2]

cox_sentinel_age_params <- cox_sentinel_age_params %>% filter(!(age_bin %in% c('vaccination_active','week_collection_number')))

ggplot() +
  geom_pointrange(data=cox_sentinel_age_params,aes(y=age_bin,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/4,1,4,16,64,256,1024,4096)),
                     labels=(c(1/4,1,4,16,64,256,1024,4096))) +
  scale_y_discrete( breaks=seq(15,95, by=10),
                    labels=c('10-19','20-29','30-39','40-49','50-59','60-69','70-79',
                             '80-89','90+'))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('risk of hospitalization compared to 0-9 year olds') 

ggsave('output/rich_vaccination/case_hospitalization_age_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_age_relRisk.svg',units='in',width=5,height=3,device='svg')

# race ADDING THIS BACK IN FOR REVIEWERS
#### For it to work, we need to add Race back into the model. If not, this will throw an error
cox_sentinel_race_params <- coxme_random_params(cox_sentinel,cox_dat,group='race')

ggplot() +
  geom_pointrange(data=cox_sentinel_race_params,aes(y=race,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/2,1,2,4)),
                     labels=(c(1/2,1,2,4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('risk of hospitalization compared to White') 

ggsave('output/rich_vaccination/case_hospitalization_race_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_race_relRisk.svg',units='in',width=5,height=3,device='svg')

# sex
cox_sentinel_sex_params <- coxme_random_params(cox_sentinel,cox_dat,group='SEX_AT_BIRTH')

ggplot() +
  geom_pointrange(data=cox_sentinel_sex_params,aes(y=SEX_AT_BIRTH,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/2,1,2,4)),
                     labels=(c(1/2,1,2,4))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('risk of hospitalization compared to female sex at birth') 

ggsave('output/rich_vaccination/case_hospitalization_sex_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_sex_relRisk.svg',units='in',width=5,height=3,device='svg')


######################################
### SENSITIVITY ANALYSES #############
######################################

########## temporal trend identifiability

# time only
cox_sentinel_time_only <- coxme(hosp_surv ~ #(1|who_lineage) + 
                                  age_bin + (1|SEX_AT_BIRTH) +
                                  (1|vaccination_active) + week_collection_number, 
                                data=cox_dat,
                                x=FALSE,y=FALSE)
summary(cox_sentinel_time_only)

coef(cox_sentinel_time_only)['week_collection_number']
vcov(cox_sentinel_time_only)['week_collection_number','week_collection_number']

cox_sentinel_time_only_vaccine_params <- coxme_random_params(cox_sentinel_time_only,cox_dat,group='vaccination_active')

# variant and NO time
cox_sentinel_voc_and_no_time <- coxme(hosp_surv ~ (1|who_lineage) + 
                                     age_bin + (1|SEX_AT_BIRTH) +
                                     (1|vaccination_active) , 
                                   data=cox_dat,
                                   x=FALSE,y=FALSE)
summary(cox_sentinel_voc_and_no_time)

coef(cox_sentinel_voc_and_no_time)['week_collection_number']
vcov(cox_sentinel_voc_and_no_time)['week_collection_number','week_collection_number']


cox_sentinel_voc_and_no_time_vaccine_params <- coxme_random_params(cox_sentinel_voc_and_no_time,cox_dat,group='vaccination_active')
cox_sentinel_voc_and_no_time_lineage_params <- coxme_random_params(cox_sentinel_voc_and_no_time,cox_dat,group='who_lineage')


# risk reduction of hospitalization from vaccine


ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=as.numeric(who_lineage),x=logRR,xmin=lower95,xmax=upper95, color='Cox Sentinel (VOC/VOI)')) +
  geom_pointrange(data=cox_sentinel_voc_and_no_time_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (VOC/VOI & time)')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','gray','cornflowerblue'),
                     breaks=c('Cox Sentinel (VOC/VOI)','Cox Sentinel (VOC/VOI & time)','Cox Sentinel (time only)'),
                     name='Model') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16,32))) +
  scale_y_continuous(breaks=1:length(cox_sentinel_lineage_params$who_lineage),
                     labels=levels(cox_sentinel_lineage_params$who_lineage))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_voc_and_time_sensitivity.png',units='in',width=7,height=4,device='png')
ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_voc_and_time_sensitivity.svg',units='in',width=7,height=4,device='svg')


ggplot() +
  geom_pointrange(data=,aes(y=as.numeric(vaccination_active)+0.1,x=logRR,xmin=lower95,xmax=upper95, color='Cox Sentinel (VOC/VOI)')) +
  geom_pointrange(data=cox_sentinel_voc_and_no_time_vaccine_params,aes(y=as.numeric(vaccination_active),x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (VOC/VOI & time)')) +
  geom_pointrange(data=cox_sentinel_time_only_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (time only)')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','gray','cornflowerblue'),
                     breaks=c('Cox Sentinel (VOC/VOI)','Cox Sentinel (VOC/VOI & time)','Cox Sentinel (time only)'),
                     name='Model') +
  scale_x_continuous(breaks=log(c(1/16, 1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     limits=log(c(1/14,3))) +
  scale_y_continuous(breaks=1:length(cox_sentinel_vaccine_params$vaccination_active),
                     labels=levels(cox_sentinel_vaccine_params$vaccination_active))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk_time_sensitivity.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk_time_sensitivity.svg',units='in',width=5,height=3,device='svg')


# mediation of time by variant
# around 80% of the time effect disappears when you include variant, as expected. 

time_effect_size <- data.frame(model = c('Cox Sentinel VOC/VOI','Cox Sentinel VOC/VOI & time','Cox Sentinel time only'),
                               coef=c(NA,coef(cox_sentinel)['week_collection_number'],coef(cox_sentinel_time_only)['week_collection_number']),
                               se=c(NA,vcov(cox_sentinel)['week_collection_number','week_collection_number'],vcov(cox_sentinel_time_only)['week_collection_number','week_collection_number'])) %>%
  mutate(z = coef/se) %>%
  mutate(coeff_fraction_after_mediation = coef/coef[3]) %>%
  mutate(z_fraction_after_mediation = z/z[3]) %>%
  mutate(effect_size = coef*mean(cox_dat$week_collection_number))

time_effect_size


#######################################
### time from collection to hospitalization 

## 14 day cutoff

cox_dat_14 <- d_14 %>% 
  filter(sequence_reason_clean=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>% 
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

hosp_surv_14 <- Surv(time=cox_dat_14$hosp_days_at_risk,event=as.numeric(cox_dat_14$mhosp=='Yes'))


cox_sentinel_14 <- coxme(hosp_surv_14 ~ (1|who_lineage) + 
                        age_bin + (1|SEX_AT_BIRTH) +
                        (1|vaccination_active) + week_collection_number, 
                      data=cox_dat_14,
                      x=FALSE,y=FALSE)
summary(cox_sentinel_14)

cox_sentinel_lineage_params_14 <- coxme_random_params(cox_sentinel_14,cox_dat_14,group='who_lineage')

## 21 day cutoff
cox_dat_21 <- d_21 %>% 
  filter(sequence_reason_clean=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>% 
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

hosp_surv_21 <- Surv(time=cox_dat_21$hosp_days_at_risk,event=as.numeric(cox_dat_21$mhosp=='Yes'))

cox_sentinel_21 <- coxme(hosp_surv_21 ~ (1|who_lineage) + 
                           age_bin + (1|SEX_AT_BIRTH) +
                           (1|vaccination_active) + week_collection_number, 
                         data=cox_dat_21,
                         x=FALSE,y=FALSE)
summary(cox_sentinel_21)

cox_sentinel_lineage_params_21 <- coxme_random_params(cox_sentinel_21,cox_dat_21,group='who_lineage')


cox_dat_30 <- d_30 %>% 
  filter(sequence_reason_clean=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>% 
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

hosp_surv_30 <- Surv(time=cox_dat_30$hosp_days_at_risk,event=as.numeric(cox_dat_30$mhosp=='Yes'))


cox_sentinel_30 <- coxme(hosp_surv_30 ~ (1|who_lineage) + 
                           age_bin + (1|SEX_AT_BIRTH) +
                           (1|vaccination_active) + week_collection_number, 
                         data=cox_dat_30,
                         x=FALSE,y=FALSE)
summary(cox_sentinel_30)

cox_sentinel_lineage_params_30 <- coxme_random_params(cox_sentinel_30,cox_dat_30,group='who_lineage')


ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95, color = 'Cox Sentinel (no filter)')) +
  geom_pointrange(data=cox_sentinel_lineage_params_14,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='14 day')) +
  geom_pointrange(data=cox_sentinel_lineage_params_21,aes(y=as.numeric(who_lineage)-0.2,x=logRR,xmin=lower95,xmax=upper95,color='21 day')) +
  geom_pointrange(data=cox_sentinel_lineage_params_30,aes(y=as.numeric(who_lineage)-0.3,x=logRR,xmin=lower95,xmax=upper95,color='30 day')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','gray','cornflowerblue','forestgreen'),
                     breaks=c('Cox Sentinel (no filter)','14 day','21 day','30 day'),
                     name=' Max Time from Collection to Hosp') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16,32))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_hosp_time_cutoff_sensitivity.png',units='in',width=7,height=4,device='png')
ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_hosp_time_cutoff_sensitivity.svg',units='in',width=7,height=4,device='svg')

########################################
### fixed vs hierarchical vs poisson


# fixed effect model sentinel (typical thing people do)
exclusions <- exclusions %>% rbind(data.frame(data_view='all sentinel surveillance and lineages with at least 1 hospitalization except suspected reinfections',
                                              reason='sensitivity analysis',
                                              n_kept=nrow(cox_dat),
                                              analysis='cox fixed effect sentinel'))


cox_fixed_sentinel <- coxph(hosp_surv ~ who_lineage + 
                              age_bin + SEX_AT_BIRTH +
                              vaccination_active,
                            data=cox_dat,
                            x=FALSE,y=FALSE)
summary(cox_fixed_sentinel)

cox_fixed_sentinel_lineage_params <- coxph_params(cox_fixed_sentinel,cox_sentinel_lineage_params)

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=cox_fixed_sentinel_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='blue') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16,32))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 


cox_fixed_sentinel_vaccine_params <- coxph_params(cox_fixed_sentinel,cox_sentinel_vaccine_params,group='vaccination_active')

ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=vaccination_active,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=cox_fixed_sentinel_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='blue') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     limits=log(c(1/16,2))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 



# hospitalization all samples
cox_dat <- d %>% 
  filter(infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes" ) %>% # can toggle to look at reinfection
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk,
         vaccination_active) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

exclusions <- exclusions %>% rbind(data.frame(data_view='all sequences and lineages with at least 1 hospitalization except suspected reinfections',
                                              reason='sensitivity analysis',
                                              n_kept=nrow(cox_dat),
                                              analysis='cox hierarchical all samples'))


hosp_surv <- Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes'))

cox_all <- coxme(hosp_surv ~ (1|who_lineage) + 
                   age_bin + (1|SEX_AT_BIRTH) +
                   (1|vaccination_active),
                 data=cox_dat,
                 x=FALSE,y=FALSE)
summary(cox_all)

cox_all_lineage_params <- coxme_random_params(cox_all,cox_dat)

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=cox_all_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='blue') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8)),
                     limits=log(c(1/2,9))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 


# risk reduction of hospitalization from vaccine

# using all significantly attenuates vaccine risk effects, because of 
# including people tested for possible vaccine breakthrough
cox_all_vaccine_params <- coxme_random_params(cox_all,cox_dat,group='vaccination_active')

ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=vaccination_active,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=cox_all_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='blue') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16)),
                     limits=log(c(1/16,2))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 



# SENSITIVITY POISSON REGRESSION
# hospitalization sentinel samples
pois_dat <- d %>% 
  filter(sequence_reason_clean=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>%
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,
         vaccination_active, 
         week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage,SEX_AT_BIRTH,age_bin,vaccination_active,week_collection_number) %>%
  summarize(n_hosp = sum(mhosp=='Yes'),
            cases=n()) %>%
  filter(who_lineage %in% c('other',levels(cox_sentinel_lineage_params$who_lineage))) %>%
  droplevels()

sum(pois_dat$cases)

exclusions <- exclusions %>% rbind(data.frame(data_view='all sentinel surveillance and lineages with at least 1 hospitalization except suspected reinfections',
                                              reason='sensitivity analysis',
                                              n_kept=sum(pois_dat$cases),
                                              analysis='poisson hierarchical sentinel samples'))


pois_sentinel <- glmer(n_hosp ~ (1|who_lineage) + log(offset(cases)) + 
                         age_bin + (1|SEX_AT_BIRTH) +
                         (1|vaccination_active),
                       data=pois_dat,
                       family='poisson')
summary(pois_sentinel)

pois_sentinel_lineage_params <- glmer_random_params(pois_sentinel,cox_sentinel_lineage_params)

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=pois_sentinel_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='red') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  geom_jitter() +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8)),
                     limits=log(c(1/2,9))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 


pois_sentinel_vaccine_params <- glmer_random_params(pois_sentinel,cox_sentinel_vaccine_params,group='vaccination_active')

ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=vaccination_active,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=pois_sentinel_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='red') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  geom_jitter() +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8)),
                     limits=log(c(1/16,2))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 



# SENSITIVITY POISSON REGRESSION
# hospitalization all samples
pois_dat <- d %>% 
  # filter(sequence_reason_clean=='SENTINEL SURVEILLANCE') %>%
  filter(infection_type != 'suspected reinfection' & REINFECTION_FLAG != "Yes") %>% # can toggle to include reinfections
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp, week_collection_number,
         vaccination_active) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage,SEX_AT_BIRTH,age_bin,vaccination_active,week_collection_number) %>%
  summarize(n_hosp = sum(mhosp=='Yes'),
            cases=n()) %>%
  filter(who_lineage %in% c('other',levels(cox_sentinel_lineage_params$who_lineage))) %>%
  droplevels()

exclusions <- exclusions %>% rbind(data.frame(data_view='all sequences and lineages with at least 1 hospitalization except suspected reinfections',
                                              reason='sensitivity analysis',
                                              n_kept=sum(pois_dat$cases),
                                              analysis='poisson hierarchical all samples'))


pois_all <- glmer(n_hosp ~ (1|who_lineage) + log(offset(cases)) + 
                    age_bin + (1|SEX_AT_BIRTH) + 
                    (1|vaccination_active),
                  data=pois_dat,
                  family='poisson')
summary(pois_all)

pois_all_lineage_params <- glmer_random_params(pois_all,cox_sentinel_lineage_params)

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel')) +
  geom_pointrange(data=pois_all_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Poisson All')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','cornflowerblue','orangered','forestgreen'),
                     breaks=c('Cox Sentinel','Cox All','Poisson Sentinel','Poisson All'),
                     name='Model') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8)),
                     limits=log(c(1/2,9))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

pois_all_vaccine_params <- glmer_random_params(pois_all,cox_sentinel_vaccine_params,group='vaccination_active')

ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=vaccination_active,x=logRR,xmin=lower95,xmax=upper95)) +
  geom_pointrange(data=pois_all_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95),color='red') +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  geom_jitter() +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8)),
                     limits=log(c(1/16,2))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 




## all sensitivity analyses together now

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=as.numeric(who_lineage)+0.2,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (hierarchical)')) +
  geom_pointrange(data=cox_fixed_sentinel_lineage_params,aes(y=as.numeric(who_lineage)+0.1,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (fixed effects)')) +
  geom_pointrange(data=cox_all_lineage_params,aes(y=as.numeric(who_lineage),x=logRR,xmin=lower95,xmax=upper95,color='Cox All (hierarchical)')) +
  geom_pointrange(data=pois_sentinel_lineage_params,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Poisson Sentinel (hierarchical)')) +
  geom_pointrange(data=pois_all_lineage_params,aes(y=as.numeric(who_lineage)-0.2,x=logRR,xmin=lower95,xmax=upper95,color='Poisson All (hierarchical)')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','gray','cornflowerblue','orangered','forestgreen'),
                     breaks=c('Cox Sentinel (hierarchical)','Cox Sentinel (fixed effects)','Cox All (hierarchical)','Poisson Sentinel (hierarchical)','Poisson All (hierarchical)'),
                     name='Model') +
  scale_x_continuous(breaks=log(c(1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/8,1/4,1/2,1,2,4,8,16,32))) +
  scale_y_continuous(breaks=1:length(cox_sentinel_lineage_params$who_lineage),
                     labels=levels(cox_sentinel_lineage_params$who_lineage))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_sensitivity.png',units='in',width=7,height=4,device='png')
ggsave('output/rich_vaccination/case_hospitalization_variant_relRisk_sensitivity.svg',units='in',width=7,height=4,device='svg')


ggplot() +
  geom_pointrange(data=cox_sentinel_vaccine_params,aes(y=as.numeric(vaccination_active)+0.2,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (mixed)')) +
  geom_pointrange(data=cox_fixed_sentinel_vaccine_params,aes(y=as.numeric(vaccination_active)+0.1,x=logRR,xmin=lower95,xmax=upper95,color='Cox Sentinel (fixed effects)')) +
  geom_pointrange(data=cox_all_vaccine_params,aes(y=as.numeric(vaccination_active),x=logRR,xmin=lower95,xmax=upper95,color='Cox All (mixed)')) +
  geom_pointrange(data=pois_sentinel_vaccine_params,aes(y=as.numeric(vaccination_active)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Poisson Sentinel (mixed)')) +
  geom_pointrange(data=pois_all_vaccine_params,aes(y=as.numeric(vaccination_active)-0.2,x=logRR,xmin=lower95,xmax=upper95,color='Poisson All (mixed)')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','gray','cornflowerblue','orangered','forestgreen'),
                     breaks=c('Cox Sentinel (mixed)','Cox Sentinel (fixed effects)','Cox All (mixed)','Poisson Sentinel (mixed)','Poisson All (mixed)'),
                     name='Model') +
  scale_x_continuous(breaks=log(c(1/16, 1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     limits=log(c(1/14,3))) +
  scale_y_continuous(breaks=1:length(cox_sentinel_vaccine_params$vaccination_active),
                     labels=levels(cox_sentinel_vaccine_params$vaccination_active))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk_sensitivity.png',units='in',width=7,height=4,device='png')
ggsave('output/rich_vaccination/case_hospitalization_vaccine_relRisk_sensitivity.svg',units='in',width=7,height=4,device='svg')




# write out exclusions table
write.table(exclusions,'output/rich_vaccination/sample_size_and_exclusions_summary.csv',sep=',',row.names=FALSE)


# save model objects where DUA allows (so any that do not contain data matrix)
save(cox_sentinel,
     cox_sentinel_lineage_params,cox_sentinel_race_params,cox_sentinel_age_params,cox_sentinel_sex_params,cox_sentinel_vaccine_params,cox_sentinel_vax_rel_risk,
     cox_fixed_sentinel,cox_all,
     cox_fixed_sentinel_lineage_params,cox_all_lineage_params,pois_all_lineage_params,pois_sentinel_lineage_params,
     cox_fixed_sentinel_vaccine_params,cox_all_vaccine_params,pois_all_vaccine_params,pois_sentinel_vaccine_params,
     file='output/rich_vaccination/cached_models.Rdata')



