library(coxme)
library(survival)
library(lme4)
library(kableExtra)
library(survminer)
library(lmtest)


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
### this adds reinfections back in 
cox_dat <- d %>% 
  filter(CDC_N_COV_2019_SEQUENCE_REASON=='SENTINEL SURVEILLANCE' ) %>% 
  select(who_lineage, collection_date, SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

cox_dat <- cox_dat %>% filter(collection_date >= "2021-09-01")
cox_dat <- cox_dat %>% filter(who_lineage %in% c('Delta (B.1.617.2)', "Omicron (B.1.1.529)"))
cox_dat$who_lineage <- factor(cox_dat$who_lineage,levels=c('Delta (B.1.617.2)', "Omicron (B.1.1.529)"))
cox_dat$who_lineage <- relevel(cox_dat$who_lineage,ref='Delta (B.1.617.2)')

hosp_surv <- Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes'))

cox_sentinel <- coxph(hosp_surv ~ who_lineage + 
                        age_bin + SEX_AT_BIRTH +
                        vaccination_active + week_collection_number , 
                      data=cox_dat,
                      x=FALSE,y=FALSE)
summary(cox_sentinel)

cox_sentinel_lineage_params <- coxph_params(cox_sentinel,cox_dat,group='who_lineage')


ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95,color=who_lineage)) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_x_continuous(breaks=log(c(1/64,1/32,1/16,1/8,1/4,1/2,1,2,4,8,16)),
                     labels=(c(1/64,1/32,1/16,1/8,1/4,1/2,1,2,4,8,16)),
                     limits=log(c(1/32,9))) +
  scale_color_manual(values=cmap,guide=FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_variant_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_variant_relRisk.svg',units='in',width=5,height=3,device='svg')

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

ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_vaccine_relRisk.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_vaccine_relRisk.svg',units='in',width=5,height=3,device='svg')


cox_sentinel_vax_rel_risk=data.frame(vaccine_type=cox_sentinel_vaccine_params$vaccination_active,
                                     rel_risk_hosp_given_case=round(exp(cox_sentinel_vaccine_params$logRR),2),
                                     lower95 = round(exp(cox_sentinel_vaccine_params$lower95),2),
                                     upper95 = round(exp(cox_sentinel_vaccine_params$upper95),2))

cox_sentinel_vax_rel_risk %>%
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  save_kable(file = "output/rich_vaccination/variant_severity_simple_vaccine_relative_risk_hospitalization.png",
             density=600,zoom=3) 
write.table(cox_sentinel_vax_rel_risk,'output/rich_vaccination/fixed_variant_severity_simple_vaccine_relative_risk_hospitalization.csv',sep=',',row.names = FALSE)

##Likelihood ratio test for global effect:

#we already have out main model: cox_sentienl

#now create a reduced model with df -1 by removing the effect that we're interested in, which is variants
cox_sentinel_test <- coxph(hosp_surv ~  age_bin + SEX_AT_BIRTH +
                             vaccination_active + week_collection_number, 
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
  ## risk.table = TRUE, risk.table.col = "who_lineage", risk.table.height = 0.25,
  
  legend.labs = c( "Delta", "Omicron"))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 10, face = "bold"),legend.key.height = unit(1, 'cm'), legend.key.width = unit(1, 'cm'))
ggsurv

ggsave('output/rich_vaccination/fixed_omicron_delta_km_curves.png',units='in',width=5,height=3,device='png')


#########
### Sensitivity for reinfections

cox_dat_no_re <- d %>% 
  filter(CDC_N_COV_2019_SEQUENCE_REASON=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' & is.na(REINFECTION_FLAG)) %>% 
  select(who_lineage, collection_date, SEX_AT_BIRTH,age_bin,
         mhosp,hosp_days_at_risk, vaccination_active,
         vaccination_active,week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup() %>%
  droplevels()

cox_dat_no_re <- cox_dat_no_re %>% filter(collection_date >= "2021-09-01")
cox_dat_no_re <- cox_dat_no_re %>% filter(who_lineage %in% c('Delta (B.1.617.2)', "Omicron (B.1.1.529)"))
cox_dat_no_re$who_lineage <- factor(cox_dat_no_re$who_lineage,levels=c('Delta (B.1.617.2)', "Omicron (B.1.1.529)"))
cox_dat_no_re$who_lineage <- relevel(cox_dat_no_re$who_lineage,ref='Delta (B.1.617.2)')

hosp_surv <- Surv(time=cox_dat_no_re$hosp_days_at_risk,event=as.numeric(cox_dat_no_re$mhosp=='Yes'))

cox_sentinel_no_re <- coxph(hosp_surv ~ who_lineage + 
                              age_bin + SEX_AT_BIRTH +
                              vaccination_active + week_collection_number , 
                            data=cox_dat_no_re,
                            x=FALSE,y=FALSE)
summary(cox_sentinel_no_re)

cox_sentinel_lineage_params_no_re <- coxph_params(cox_sentinel_no_re,cox_dat_no_re,group='who_lineage')


ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=who_lineage,x=logRR,xmin=lower95,xmax=upper95,color='Reinfections included')) +
  geom_pointrange(data=cox_sentinel_lineage_params_no_re,aes(y=as.numeric(who_lineage)-0.1,x=logRR,xmin=lower95,xmax=upper95,color='Reinfections excluded')) +
  geom_vline(aes(xintercept=0),linetype='dashed') +
  scale_color_manual(values=c('black','cornflowerblue'),
                     breaks=c('Reinfections included','Reinfections excluded'),
                     name='reinfections sens') +
  scale_x_continuous(breaks=log(c(1/64,1/32,1/16,1/8,1/4,1/2,1,2,4,8,16)),
                     labels=(c(1/64,1/32,1/16,1/8,1/4,1/2,1,2,4,8,16)),
                     limits=log(c(1/32,9))) +
  scale_color_manual(values=cmap,guide=FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank()) +
  ylab('') +
  xlab('hazard ratio for hospitalization') 

ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_variant_relRisk_reinfection_sens.png',units='in',width=5,height=3,device='png')
ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_variant_relRisk_reinfection_sens.svg',units='in',width=5,height=3,device='svg')


save(cox_sentinel, cox_sentinel_lineage_params,
     file='output/omicron_delta_cached_variant_models.Rdata') 


###############
#vaccine facet for delta/omicron
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

###think about whether to add or exclude reinfections
cox_dat <- d %>%
  filter(CDC_N_COV_2019_SEQUENCE_REASON=='SENTINEL SURVEILLANCE' ) %>%
  select(who_lineage,SEX_AT_BIRTH,age_bin, collection_date, admitdate,
         mhosp,hosp_days_at_risk, vaccination_active, week_collection_number) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup()  %>%
  mutate(active_vaccine_type_dose_lineage = interaction(vaccination_active,
                                                        who_lineage,
                                                        sep=' : ')
  ) %>%
  filter(who_lineage %in% c("Omicron (B.1.1.529)",'Delta (B.1.617.2)')) %>% #Will need to hardcode omicron in 
  droplevels()

cox_dat$who_lineage <- factor(cox_dat$who_lineage,levels=c('Delta (B.1.617.2)', "Omicron (B.1.1.529)"))
cox_dat$who_lineage <- relevel(cox_dat$who_lineage,ref='Delta (B.1.617.2)')

# track which data went into this analysis
exclusions <- exclusions %>% rbind(data.frame(data_view='all sentinel surveillance and lineages with at least 1 hospitalization',
                                              reason='primary analysis',
                                              n_kept=nrow(cox_dat),
                                              analysis='cox hierarchical sentinel'))


hosp_surv <- Surv(time=cox_dat$hosp_days_at_risk,event=as.numeric(cox_dat$mhosp=='Yes'))

cox_sentinel <- coxph(hosp_surv ~
                        active_vaccine_type_dose_lineage +
                        age_bin + SEX_AT_BIRTH + week_collection_number,
                      data=cox_dat,
                      x=FALSE,y=FALSE)
summary(cox_sentinel)

cox_sentinel_lineage_params <- coxph_params(cox_sentinel,cox_dat,group='active_vaccine_type_dose_lineage', by_lineage=FALSE) # by_lineage flag make "none" the reference for each variant




ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=active_vaccine_type_dose,x=logRR,xmin=lower95,xmax=upper95,color=who_lineage)) +
  geom_vline(data=data.frame(xint=0),mapping=aes(xintercept=xint),linetype='dashed') +
  facet_grid(rows=vars(fct_rev(who_lineage)))+#, scales = "free_y", labeller=as_labeller(lineage_names)) +
  scale_color_manual(values=cmap, guide=FALSE)+
  scale_x_continuous(breaks=log(c(1/32, 1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/32,1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     limits=log(c(1/32,16))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7)) +
  ylab('') +
  xlab('hospitalization hazard ratio')


ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_vaccine_variant_interaction_relRisk.png',units='in',width=5,height=5,device='png')
ggsave('output/rich_vaccination/fixed_omicron_delta_case_hospitalization_vaccine_variant_interaction__relRisk.svg',units='in',width=5,height=5,device='svg')

save(cox_sentinel, cox_sentinel_lineage_params,
     file='output/omicron_delta_cached_variant_models_interaction.Rdata')

