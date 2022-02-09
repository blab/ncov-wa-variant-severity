# MUST RUN variant_severity_data_prep.R first to set up workspace

library(coxme)
library(survival)
library(lme4)
library(kableExtra)
# package:merTools is also required in a function, but not loaded into global namespace since it conflicts with dplyr::Select



# if you are working with this but don't have the raw data, you can load
# saved cox models, and saved who_lineage hazard ratio parameter data frames here.
# To use, comment out sections that don't run (data organization and model fitting),
# and otherwise regenerate figures, or make new ones.
# For the cox models, because the model object is available, it is also possible to consider
# new contrasts, or make figures or tables describing control variables.
# load(file='variant_severity_simple_vaccination_cox_models.Rdata')

#compared to the rich vaccination script, the coxme_random_params function was updated to use 
#the None category of each variant as the reference for that variant. 
#it's not the cleanest script but it was hard to change without being able to see the dataset
#this script however, requires you to change the size if the number of categories changes
#it's currently set up for excluding J&J but please change accordingly if J&J is to be added back in


########################
### helper functions ###
########################


coxme_random_params <- function(mod,data,group='who_lineage',by_lineage=FALSE){
  # formats coefficient estimates from coxme_object (hierarchical model)
  
  res <- data.frame(logRR = ranef(mod)[[group]])
  res$logRR <- res$logRR-res$logRR[1] # set reference to zero
  
  res[group]=sub(group,'',row.names(res))
  row.names(res)<-NULL
  
  if (group =='active_vaccine_type_dose_lineage'){
    res$who_lineage <- factor(sub(".*: ", "", res$active_vaccine_type_dose_lineage),
                              levels=levels(data$who_lineage))
    res$active_vaccine_type_dose <- factor(sub(" :.*", "", res$active_vaccine_type_dose_lineage),
                                           levels=levels(data$vaccination_active))
  }
  
  # find relevant block of variance-covariance matrix
  # random effects are first, then fixed https://cran.r-project.org/web/packages/coxme/coxme.pdf
  # order should be order they appear in the formula
  if (by_lineage){
    
    mod_vcov <- mod$variance
    idx_type <- rep(names(mod$frail),lengths(mod$frail))
    group_idx <- which(idx_type==group)
    mod_vcov <- mod_vcov[group_idx,group_idx]
    
    lineage_type <- res$who_lineage
    
    res$sd <- diag(mod_vcov)
    
    for (L in levels(res$who_lineage)){
      lineage_idx <- which(lineage_type==L)
      
      tmp_vcov <- mod_vcov[lineage_idx,lineage_idx]
      
      # propagate reference variance, variance of difference of two random variable
      tmp_idx <- lineage_idx[2:length(lineage_idx)]
      res$sd[tmp_idx]<- res$sd[tmp_idx] + 
        res$sd[lineage_idx[1]] - 2*c(as.matrix(tmp_vcov)[(row(tmp_vcov) == (col(tmp_vcov) - 1))])[1:length(tmp_idx)]
      res$sd[lineage_idx[1]]<-0
      
      # propagate mean
      res$logRR[tmp_idx] <- res$logRR[tmp_idx] - res$logRR[lineage_idx[1]] 
      res$logRR[lineage_idx[1]] <- 0 # set reference to zero
    }
    
  } else {
    mod_vcov <- mod$variance
    idx_type <- rep(names(mod$frail),lengths(mod$frail))
    group_idx <- which(idx_type==group)
    mod_vcov <- mod_vcov[group_idx,group_idx]
    
    res$sd <- diag(mod_vcov)
    # propagate reference variance, variance of difference of two random variable
    res$sd[2:nrow(res)]<- res$sd[2:nrow(res)] + 
      res$sd[1] - 2*c(as.matrix(mod_vcov)[(row(mod_vcov) == (col(mod_vcov) - 1))])[1:nrow(res)-1]
    res$sd[1]<-0
  }
  
  res$sd <- sqrt(res$sd)
  res$lower95 <- res$logRR -2*res$sd
  res$upper95 <- res$logRR +2*res$sd
  
  res <- res %>% 
    # filter(!!as.symbol(group) != levels(data[[group]])[1]) %>%
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



############# PRIMARY ANALYSIS. Case-Hospitalization ratio in sentinel surveillance
############# hierarchical cox model
############# vaccination with rich model of vaccine type and variant

# hospital sentinel only cox hierarchical model 
cox_dat <- d %>%
  filter(CDC_N_COV_2019_SEQUENCE_REASON=='SENTINEL SURVEILLANCE' &
           infection_type != 'suspected reinfection' &is.na(REINFECTION_FLAG)) %>%
  select(who_lineage,SEX_AT_BIRTH,age_bin, collection_date, admitdate,
         mhosp,hosp_days_at_risk, vaccination_active, week_collection_number, CDC_N_COV_2019_SEQUENCE_ACCESSION_NUMBER) %>%
  # drop lineages with no hospitalization outcomes
  group_by(who_lineage) %>%
  mutate(n_hosp = sum(mhosp=='Yes')) %>%
  filter(n_hosp>0) %>%
  ungroup()  %>%
  mutate(active_vaccine_type_dose_lineage = interaction(vaccination_active,
                                                        who_lineage,
                                                        sep=' : ')
  ) %>%
  filter(who_lineage %in% c('other', 'Alpha (B.1.1.7)', 'Gamma (P.1)','Delta (B.1.617.2)', "Omicron (B.1.1.529)")) %>% #Will need to hardcode omicron in 
  droplevels()

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

cox_sentinel_lineage_params$active_vaccine_type_dose <- as.character(cox_sentinel_lineage_params$active_vaccine_type_dose)
cox_sentinel_lineage_params <- cox_sentinel_lineage_params %>% filter(!(cox_sentinel_lineage_params$active_vaccine_type_dose_lineage %in% c("≥21 days post dose one to \n <21 days post booster : other", "≥21 days post booster : Gamma (P.1)","≥21 days post booster : Alpha (B.1.1.7)", "≥21 days post booster : other" )))
cox_sentinel_lineage_params$active_vaccine_type_dose <- factor(cox_sentinel_lineage_params$active_vaccine_type_dose, levels=c("No Vaccination to \n <21 days post dose one", "≥21 days post dose one to \n <21 days post booster", "≥21 days post booster"))

lineage_names <- c(
  `Omicron (B.1.1.529)` = "Omicron (B.1.1.529)",
  `Gamma (P.1)`="Gamma (P.1)",
  `Delta (B.1.617.2)`="Delta (B.1.617.2)",
  `Alpha (B.1.1.7)`="Alpha (B.1.1.7)",
  `other`="REF:Ancestral"
)

ggplot() +
  geom_pointrange(data=cox_sentinel_lineage_params,aes(y=active_vaccine_type_dose,x=logRR,xmin=lower95,xmax=upper95,color=who_lineage)) +
  geom_vline(data=data.frame(xint=0),mapping=aes(xintercept=xint),linetype='dashed') +
  facet_grid(rows=vars(fct_rev(who_lineage)), scales = "free_y", labeller=as_labeller(lineage_names)) +
  scale_color_manual(values=cmap, guide=FALSE)+
  scale_x_continuous(breaks=log(c(1/32, 1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     labels=(c(1/32,1/16,1/8,1/4,1/2,1,2,4,8,16,32)),
                     limits=log(c(1/24,16))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7)) +
  ylab('') +
  xlab('hospitalization hazard ratio')

ggsave('output/rich_vaccination/case_hospitalization_vaccine_variant_interaction_relRisk.png',units='in',width=5,height=5,device='png')
ggsave('output/rich_vaccination/case_hospitalization_vaccine_variant_interaction__relRisk.svg',units='in',width=5,height=5,device='svg')

save(cox_sentinel, cox_sentinel_lineage_params,
     file='output/cached_variant_models.Rdata')


hosp_by_variant_vaccine <- with(cox_dat, table(mhosp, active_vaccine_type_dose_lineage))

hosp_by_variant_vaccine %>%
  kbl() %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
  save_kable(file = "output/rich_vaccination/hosp_by_variant_and_vaccine.png",
             density=600,zoom=3) 
write.table(hosp_by_variant_vaccine,'output/rich_vaccination/hosp_by_variant_and_vaccine.csv',sep=',',row.names = FALSE)


ancestral_booster <- cox_dat %>% filter(active_vaccine_type_dose == "≥21 days post booster : other")
view(ancestral_booster)