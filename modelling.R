library(tidyverse)
library(glmmTMB)
library(easystats)
library(performance)
library(lme4)

read_csv("tax_and_spec_div_values.csv") %>%
  filter(image_type=="masked")->masked_df

m2<-glmmTMB(species_richness~scale(CV)+(1|site/subplot_id),data=masked_df)
r2_nakagawa(m2)

library(lme4)
m_masked2<-glmmTMB(species_richness~scale(SV)+(1|site/subplot_id),data=masked_df)
r2_nakagawa(m_masked2)


a %>%
  filter(image_type=="unmasked")->unmasked_df

m_unmasked1<-glmmTMB(species_richness~scale(CV)+(1|site/subplot_id),data=unmasked_df)
r2_nakagawa(m_unmasked1)

m_unmasked2<-glmmTMB(species_richness~SV+(1|site/subplot_id),data=unmasked_df)
r2_nakagawa(m_unmasked2)



