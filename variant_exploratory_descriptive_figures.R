# script to generate some exploratory and descriptive figures


# MUST RUN variant_severity_data_prep.R first to set up workspace
# a couple plots also require the IDM line list
# IDM_linelist_EXPANDED_AsOf_2021-07-16_exported_2021-07-16.csv

library(zoo)
library(data.table)
library(lme4)


# make colormap from mpn65
# https://github.com/google/palette.js/blob/79a703df344e3b24380ce1a211a2df7f2d90ca22/palette.js#L802
# reordered to better match old colormap

cmap <- toupper(paste('#',c(
  'ff0029', '377eb8', '66a61e', '00d2d5', 'b3e900', '984ea3', 'ff7f00',
  'af8d00', '7f80cd', 'c42e60', 'a65628', 'f781bf', '8dd3c7', 'bebada',
  'fb8072', '80b1d3', 'fdb462', 'fccde5', 'bc80bd', 'ffed6f', 'c4eaff',
  'cf8c00', '1b9e77', 'd95f02', 'e7298a', 'e6ab02', 'a6761d', '0097ff',
  '00d067', '000000', '252525', '525252', '737373', '969696', 'bdbdbd',
  'f43600', '4ba93b', '5779bb', '927acc', '97ee3f', 'bf3947', '9f5b00',
  'f48758', '8caed6', 'f2b94f', 'eff26e', 'e43872', 'd9b100', '9d7a00',
  '698cff', 'd9d9d9', '00d27e', 'd06800', '009f82', 'c49200', 'cbe8ff',
  'fecddf', 'c27eb6', '8cd2ce', 'c4b8d9', 'f883b0', 'a49100', 'f48800',
  '27d0df', 'a04a9b'),sep=''))
# reduce saturation
cmap<-rgb(t(round(col2rgb(cmap)*0.9)), maxColorValue=255)

cmap <- cmap[1:length(levels(d$who_lineage))]
names(cmap) <- levels(d$who_lineage)[c(2:length(cmap),1)]
cmap['other']='#999999'


#########################
#### outcome counts  ####
#########################


# fraction of cases sequenced
# lineages over time
sum(d$mhosp=="Yes")  # only enough hospitalizations to be worth analyzing
sum(d$dead=="Yes",na.rm=TRUE)
sum(d$ICU_STATUS=="Yes",na.rm=TRUE)
sum(d$MECHANICAL_VENTILATION=="Yes",na.rm=TRUE)


########################################
## descriptive variant over time plots
########################################

## variant fraction running average

plot_dat <- d %>% 
  select(best_infection_event_date,who_lineage) %>%
  group_by(best_infection_event_date,who_lineage) %>%
  summarize(variant_count=n()) %>%
  group_by(best_infection_event_date) %>%
  mutate(n=sum(variant_count)) %>% droplevels()

plot_dat <- plot_dat %>% right_join(
    expand.grid(best_infection_event_date=seq.Date(min(plot_dat$best_infection_event_date),max(plot_dat$best_infection_event_date),by=1),
                who_lineage=levels(plot_dat$who_lineage))
  ) %>% arrange(best_infection_event_date) %>%
  replace_na(list(variant_count=0,n=0))


tau<-21
plot_dat<-plot_dat %>% group_by(who_lineage) %>%
  mutate(running_average =na.approx(rollapply(variant_count,tau,sum,align='center',fill=NA),na.rm=FALSE)/
           na.approx(rollapply(n,tau,sum,align='center',fill=NA),na.rm=FALSE)) %>%
  drop_na() %>%
  group_by(best_infection_event_date) %>%
  mutate(running_average=running_average/sum(running_average)) %>% 
  mutate(who_lineage=factor(who_lineage,levels=levels(d$who_lineage))) %>%
  droplevels() 

# area plot

ggplot(plot_dat) + 
  geom_area(aes(x=best_infection_event_date, y=running_average, group=who_lineage,fill=who_lineage)) +
  # geom_line(aes(x=best_infection_event_date, y=running_average, group=who_lineage,color=who_lineage)) +
  scale_fill_manual(values=cmap) +
  scale_color_manual(values=cmap) +
  guides(fill=guide_legend(title="Variant")) +
  scale_x_date(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('observed variant fraction') + xlab('') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black'),
        legend.key.size = unit(0.14, 'in'))

ggsave('output/data_description/WA_voc_fraction_21_day_running_average.png',units='in',width=6,height=3,device='png')
ggsave('output/data_description/WA_voc_fraction_21_day_running_average.svg',units='in',width=6,height=3,device = 'svg')

# logit facet plot
ggplot(plot_dat) + 
  geom_point(aes(x=best_infection_event_date, y=variant_count/n, group=who_lineage,color=who_lineage,size=n)) +
  geom_line(aes(x=best_infection_event_date, y=running_average, group=who_lineage,color=who_lineage)) +
  facet_wrap('who_lineage') +
  scale_fill_manual(values=cmap) +
  scale_color_manual(values=cmap) +
  guides(color=guide_legend(title="Variant")) +
  scale_x_date(expand = c(0, 0)) + 
  scale_y_continuous(limits=c(1e-3,0.99),trans='logit',breaks=c(0.001,0.003,0.01,0.03,0.1,0.3,0.7,0.9,0.97)) +
  scale_size_continuous(range = c(0.25,2),guide='none') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('observed variant fraction') + xlab('') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
    axis.line = element_line(), 
    axis.ticks = element_line(colour='black'),
    legend.key.size = unit(0.14, 'in'))

ggsave('output/data_description/WA_voc_fraction_logit.png',units='in',width=8,height=6,device='png')
ggsave('output/data_description/WA_voc_fraction_logit.svg',units='in',width=8,height=6,device = 'svg')


# code to append estimated variant fractions to smoothed case data from IDM linelist

## case data from IDM linelist
# using fread from data.table because it's faster
cdat<-fread('IDM_linelist_EXPANDED_AsOf_2021-07-16_exported_2021-07-16.csv',
            select = c('SPECIMEN__COLLECTION__DTTM',"LabTestResult"),
            sep=',',header=TRUE,stringsAsFactors = TRUE)

# format date and select data cutoff
cdat$best_infection_event_date <- as.Date(cdat$SPECIMEN__COLLECTION__DTTM)
cdat <- cdat %>% filter(best_infection_event_date>=min(d$best_infection_event_date) &
                          best_infection_event_date <=max(d$best_infection_event_date))


# create smooth cases per day

# summarize cases per day
cases <- cdat %>% filter(LabTestResult=='Positive') %>% 
  group_by(best_infection_event_date) %>% summarize(n=n()) 

# smooth cases
cases<-cases %>% 
  mutate(cases_running_average =na.approx(rollapply(n,tau,mean,align='right',fill=NA),na.rm=FALSE)) 


# join cases with fractions
cases <- cases %>% right_join(plot_dat %>% select(best_infection_event_date,who_lineage,running_average),by='best_infection_event_date')

# construct estimate of absolute case count running average by variant
cases <- cases %>% mutate(variant_count=running_average*cases_running_average)


# area plot
ggplot(cases) +
  geom_area(aes(x=best_infection_event_date,y=variant_count,group=who_lineage,fill=who_lineage)) +
  scale_fill_manual(values=cmap) +
  scale_x_date(date_breaks='2 months',date_labels="%b") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ylab('estimated cases per day') + xlab('') +
  guides(fill=guide_legend(title="Variant")) +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.13, 'in'))

ggsave('output/data_description/WA_estimated_variant_cases_21_day_running_average.png',units='in',width=6,height=3,device='png')
ggsave('output/data_description/WA_estimated_variant_cases_21_day_running_average.svg',units='in',width=6,height=3,device = 'svg')

# facet plot
ggplot(cases) +
  geom_area(aes(x=best_infection_event_date,y=variant_count,group=who_lineage,fill=who_lineage)) +
  facet_wrap('who_lineage',scale='free_y')+
  scale_fill_manual(values=cmap) +
  scale_x_date(date_breaks='2 months',date_labels="%b") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  guides(color=guide_legend(title="Variant")) +
  ylab('estimated cases per day') + xlab('') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
        axis.line = element_line(), 
        axis.ticks = element_line(colour='black'),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.key.size = unit(0.13, 'in')) 

ggsave('output/data_description/WA_estimated_variant_cases_facet.png',units='in',width=8,height=6,device='png')
ggsave('output/data_description/WA_estimated_variant_cases_facet.svg',units='in',width=8,height=6,device = 'svg')




###########################
#### exploratory stats ####
###########################

CHR <- d %>% group_by(who_lineage,age_bin) %>%
  summarize(n=n(),
            n_hosp=sum(mhosp=='Yes'))

ggplot(CHR) + 
  geom_line(aes(x=age_bin,y=n_hosp/n,group=who_lineage,color=who_lineage))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values=cmap) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  guides(color=guide_legend(title="Variant")) +
  ylab('hospitalizations per case') + xlab('age bin (center)') +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border=element_blank(), 
    axis.line = element_line(), 
    axis.ticks = element_line(colour='black'),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
    legend.key.size = unit(0.13, 'in')) 

# slopes pretty parallel, no obvious interaction with age and lineage

# regression intereaction?
CHR <- d %>% group_by(who_lineage,AGE_YEARS) %>%
  summarize(n=n(),
            n_hosp=sum(mhosp=='Yes'))

# no significant age lineage interactions in a fixed effect model
summary(mod <- glm(n_hosp~AGE_YEARS*who_lineage,
                   family='poisson',
                   data=CHR))

# hierarchical/smoothed model?
summary(mod <- glmer(n_hosp~ AGE_YEARS + (1+AGE_YEARS|who_lineage),
                     family='poisson',
                     data=CHR))


# get effect and look at Wald CI
tmp <- merTools::REextract(mod)

# recenter 'other'
tmp$AGE_YEARS <- tmp$AGE_YEARS - tmp$AGE_YEARS[1]
# propagate uncertainty
tmp$AGE_YEARS_se <- tmp$AGE_YEARS_se + tmp$AGE_YEARS_se[1]
tmp$AGE_YEARS_se[1]<-0


tmp$age_lineage_lower95 <- tmp$AGE_YEARS-2*tmp$AGE_YEARS_se
tmp$age_lineage_upper95 <- tmp$AGE_YEARS+2*tmp$AGE_YEARS_se

tmp$z_lower95=tmp$age_lineage_lower95/abs(tmp$AGE_YEARS)
tmp$z_upper95=tmp$age_lineage_upper95/abs(tmp$AGE_YEARS)

# CI scatter around zero, some overlapping, some not that far away. 
# we could think about this harder but I don't think I'm seeing much. 
# If I had more time, I woudl try harder.
tmp %>% select(groupID,age_lineage_lower95,age_lineage_upper95,z_lower95,z_upper95)

