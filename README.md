# Associations between SARS-CoV-2 variants and risk of COVID-19 hospitalization among confirmed cases in Washington State: a retrospective cohort study

Miguel I. Paredes1,2✝*, Stephanie M. Lunn3✝, Michael Famulare4, Lauren A. Frisbie3, Ian Painter, 3, Roy Burstein4, Pavitra Roychoudhury2,5,  Hong Xie5, Shah A. Mohamed Bakhash5, Ricardo Perez5, Maria Lukes5, Sean Ellis5, Saraswathi Sathees5, Patrick C. Mathias,5, Alexander Greninger,2,5, Lea M. Starita6,7, Chris D. Frazar6, Erica Ryke6, Weizhi Zhong7, Luis Gamboa7, Machiko Threlkeld6, Jover Lee2, Deborah A. Nickerson6,7, Daniel L. Bates8, Matthew E. Hartman8,9, Eric Haugen8, Truong N. Nguyen8, Joshua D. Richards8, Jacob L. Rodriguez8, John A. Stamatoyannopoulos8, Eric Thorland8, Geoff Melly3, Philip E. Dykema3, Drew C. MacKellar3, Hannah K. Gray3, Avi Singh3, JohnAric MoonDance Peterson3, Denny Russell3,  Laura Marcela Torres3, Scott Lindquist3, Trevor Bedford1,2,6, Krisandra J. Allen3, Hanna N. Oltean3*

1 Department of Epidemiology, University of Washington, Seattle, WA, USA
2 Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, Washington, USA
3 Washington State Department of Health, Shoreline, WA USA
4 Institute for Disease Modeling, Bill and Melinda Gates Foundation, Seattle, WA USA
5 Department of Laboratory Medicine and Pathology, University of Washington, Seattle, WA, USA
6 Department of Genome Sciences, University of Washington, Seattle, WA, USA
7 Brotman Baty Institute for Precision Medicine, Seattle, WA USA
8 Altius Institute for Biomedical Sciences, Seattle, WA USA
9 Department of Cardiovascular Services, Swedish Medical Center, Seattle, WA USA
✝These authors contributed equally to this work

*Abstract* 

_Background_: The COVID-19 pandemic is now dominated by variant lineages; the resulting impact on disease severity remains unclear. Using a retrospective cohort study, we assessed the risk of hospitalization following infection with nine variants of concern or interest (VOC/VOI).

_Methods_: Our study includes individuals with positive SARS-CoV-2 RT-PCR in the Washington Disease Reporting System and with available viral genome data, from December 1, 2020 to July 30, 2021. The main analysis was restricted to cases with specimens collected through sentinel surveillance. Using a Cox proportional hazards model with mixed effects, we estimated hazard ratios (HR) for the risk of hospitalization following infection with a VOC/VOI, adjusting for age, sex, and vaccination status.

_Findings_: Of the 27,814 cases, 23,170 (83.3%) were sequenced through sentinel surveillance, of which 726 (3.1%) were hospitalized due to COVID-19. Higher hospitalization risk was found for infections with Gamma (HR 3.17, 95% CI 2.15-4.67), Beta (HR: 2.97, 95% CI 1.65–5.35), Delta (HR: 2.30, 95% CI 1.69-3.15), and Alpha (HR 1.59, 95% CI 1.26–1.99) compared to infections with an ancestral lineage. Following VOC infection, unvaccinated patients show a similar higher hospitalization risk, while vaccinated patients show no significant difference in risk, both when compared to unvaccinated, ancestral lineage cases.   

_Interpretation_: Infection with a VOC results in a higher hospitalization risk, with an active vaccination attenuating that risk. Our findings support promoting hospital preparedness, vaccination, and robust genomic surveillance.  

Code:
While data is not included due to patient privacy concerns, the code contained in this repository represent the analytic code for the above manuscript. If running with data from the Washington Department of Health, please begin by first running `variant_severity_data_prep.R` followed by `variant_exploratory_descriptive_figures.R` before continuing with the other analytic scripts. 
