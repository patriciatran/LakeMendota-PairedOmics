#2020 Buoy data see email from Mark Gahler:
library(tidyverse)
library(lubridate)

buoy.2020 <- read.csv("Data/ntl130_2020_v1.csv")

buoy.2020$sampledate <- ymd(buoy.2020$sampledate)

plot.buoy <- buoy.2020 %>% group_by(sampledate, depth) %>% mutate(averagewtemp = mean(wtemp),
                                                     minwtemp = min(wtemp),
                                                     maxwtemp = max(wtemp)) %>%
  select(-sampletime) %>%
  distinct() %>%
  ggplot(aes(x=sampledate, y=-depth, fill=averagewtemp))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient(low="navy blue", high="yellow")+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90))+
  xlab("2020")+
  geom_vline(xintercept = ymd("2020-10-18"))+
  ylab("Depth (m)")+
  scale_x_date(date_labels = "%m-%d",date_breaks="1 week")

ggsave(plot = plot.buoy , 
       filename = "Figures/Mendota-Buoy2020.pdf",
       width=11, 
       height = 7)
  
