# Script for Table One, general characteristics of study population
# for sequenced case data in Washington state.
#
# Data curated by WA-DOH from public health surveillance records and GISAID sequence data.
# Variant calls using pangolin v1.2.13. 
#

# run `variant_severity_data_prep.R` and `variant_severity_analysis_rich_vaccination.R` scripts prior to running this script

library(tidyverse)

# run variant_severity_data.prep and then variant_severity_analysis_rich_vaccination
with(cox_dat, table(who_lineage, useNA = "ifany"))
with(cox_dat, prop.table(table(who_lineage))) %>% round(3)

with(cox_dat, table(SEX_AT_BIRTH, useNA = "ifany"))
with(cox_dat, table(SEX_AT_BIRTH, who_lineage, useNA = "ifany"))
with(cox_dat, prop.table(table(SEX_AT_BIRTH, who_lineage),2)) %>% round(3)
with(cox_dat, chisq.test(table(SEX_AT_BIRTH, who_lineage), correct = FALSE))

with(cox_dat, table(age_bin, useNA = "ifany"))
with(cox_dat, table(age_bin, who_lineage, useNA = "ifany"))
with(cox_dat, prop.table(table(age_bin, who_lineage),2)) %>% round(3)
with(cox_dat, chisq.test(table(age_bin, who_lineage), correct = FALSE))


with(cox_dat, table(mhosp, useNA = "ifany"))
with(cox_dat, table(mhosp, who_lineage, useNA = "ifany"))
with(cox_dat, prop.table(table(mhosp, who_lineage),2)) %>% round(3)
with(cox_dat, chisq.test(table(mhosp, who_lineage), correct = FALSE))


with(cox_dat, table(vaccination_active, useNA = "ifany"))
with(cox_dat, table(vaccination_active, who_lineage, useNA = "ifany"))
with(cox_dat, prop.table(table(vaccination_active, who_lineage),2)) %>% round(2)
with(cox_dat, chisq.test(table(vaccination_active, who_lineage), correct = FALSE))
with(cox_dat, table(mhosp, vaccination_active))

#median time to hospialization with IQR by variant. 
cox_dat %>%  filter(hosp_days_at_risk >0) %>% group_by(who_lineage) %>% summarize(median = median(hosp_days_at_risk), IQR = IQR(hosp_days_at_risk))

     
