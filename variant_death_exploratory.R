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
    ref_set =  c('No Vaccination to \n <21 days post dose one',levels(ref$vaccination_active))
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
d$REINFECTION_FLAG <- as.character(d$REINFECTION_FLAG)
############# PRIMARY ANALYSIS. Case-Hospitalization ratio in sentinel surveillance
############# hierarchical cox model
############# vaccination with rich model of vaccine type and variant

# hospital sentinel only cox hierarchical model
## added in excludsion for new reinfection flag since that only came into play post setp 2021
cox_dat <- d %>% 
  filter(CDC_N_COV_2019_SEQUENCE_REASON=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & is.na(REINFECTION_FLAG)) %>% 
  select(who_lineage,SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number, race, COUNTY, dead_days_at_risk, dead) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_dead = sum(dead=='Yes')) %>%
  filter(n_dead>0) %>%
  ungroup() %>%
  droplevels() %>%
  filter(!(who_lineage %in% c("Kappa (B.1.617.1)", "Mu (B.1.621)", 'Eta (B.1.525)','Lambda (C.37)')))


# track which data went into this analysis
exclusions <- exclusions %>% rbind(data.frame(data_view='all sentinel surveillance and lineages with at least 1 hospitalization but excluding suspected reinfections',
                                              reason='primary analysis',
                                              n_kept=nrow(cox_dat),
                                              analysis='cox hierarchical sentinel'))


dead_surv <- Surv(time=cox_dat$dead_days_at_risk,event=as.numeric(cox_dat$dead=='Yes'))

cox_sentinel <- coxme(dead_surv ~ (1|who_lineage) + 
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
                     limits=log(c(1/4,9))) +
  scale_color_manual(values=cmap,guide=FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for deaths') 

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
cox_sentinel_test <- coxme(dead_surv ~  age_bin + (1|SEX_AT_BIRTH) +
                             (1|vaccination_active) + week_collection_number, 
                           data=cox_dat,
                           x=FALSE,y=FALSE)

#now conduct a LR test

lrtest(cox_sentinel, cox_sentinel_test)

## adding in Kaplan meier curves 
#### always make sure to check labels since the legend is hardcoded in
dead_counts <- survdiff(Surv(time=cox_dat$dead_days_at_risk,event=as.numeric(cox_dat$dead=='Yes')) ~ who_lineage, data= cox_dat)
dead_counts %>%
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  save_kable(file = "output/rich_vaccination/cases_and_deads.png",
             density=600,zoom=3) 
write.table(dead_counts,'output/rich_vaccination/variant_severity_cases_and_dead.csv',sep=',',row.names = FALSE)

ggsurv <- ggsurvplot(
  fit = survfit(Surv(time=cox_dat$dead_days_at_risk,event=as.numeric(cox_dat$dead=='Yes')) ~ who_lineage, data= cox_dat), 
  ylim = c(0.925, 1),
  # risk.table = TRUE, risk.table.col = "who_lineage", risk.table.height = 0.25,
  
  legend.labs = c("Ancestral", "Alpha", 'Beta',"Delta", "Epsilon", "Gamma", "Iota", "Omicron"))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 7, face = "bold"),legend.key.height = unit(0.75, 'cm'), legend.key.width = unit(0.5, 'cm'))
ggsurv

ggsave('output/rich_vaccination/km_curves.png',units='in',width=5,height=3,device='png')



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
