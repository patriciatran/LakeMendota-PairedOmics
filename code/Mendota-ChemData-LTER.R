
library(lubridate)
library(tidyverse)

chem.data <- read.csv("~/Downloads/chemphys (5).csv")

chem.data <- chem.data %>% gather("nutrient","value",cl:mn)

chem.data <- chem.data %>% filter(!is.na(value))
chem.data$sampledate <- ymd(chem.data$sampledate)
chem.data$month <- month(chem.data$sampledate)

may_to_sept <- chem.data %>% filter(between(month, 5, 9))

ggplot(may_to_sept, aes(x=value, y=depth, 
                        col=sampledate))+
  geom_point(alpha=0.8)+
  facet_wrap(nutrient~.,scales="free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        text = element_text(size=12),
        panel.grid = element_blank()
        )+
  xlab("concentration (mg/L)")+
  ylab("Depth (m)")

ggsave("Figures/Mendota_Chem_Data_LTER.pdf", width = 6, height = 6)
