ysi_2021<-read.csv("Data/cat_YSI_MO_2021.csv")
library(tidyverse)
library(lubridate)

ysi_2021$Timestamp<-mdy_hm(ysi_2021$Timestamp)

ysi_2021$Date<-date(ysi_2021$Timestamp)

ggplot(ysi_2021, aes(x=Temperature..C., y=-Folder, col=as.factor(Date)))+
  geom_point()+
  facet_wrap(Date~.,)+
  geom_line(orientation ="y")+
  theme_bw()+
  ylab("Depth(m)")+
  xlab("Temperature (C)")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank())

ggsave("Figures/YSI_MO_Temp_facet.pdf", width = 11, height = 11)

ggplot(ysi_2021 %>% filter(Dissolved.Oxygen..mg.L.>0), 
       aes(x=Dissolved.Oxygen..mg.L., y=-Folder, col=as.factor(Date)))+
  facet_wrap(Date~.,)+
  geom_point()+
  geom_line(orientation="y")+
  theme_bw()+
  ylab("Depth(m)")+
  xlab("Dissolved Oxygen (mg/L)")+
  theme(legend.position = "none",
        strip.background = element_blank(),
        panel.grid = element_blank())

ggsave("Figures/YSI_MO_DO_mgL_facet.pdf", width = 11, height = 11)

# It's kinda rare to  get exactly 0 mg/L
anoxic_depths_dates <- ysi_2021 %>% filter(Dissolved.Oxygen..mg.L.>=0) %>% 
  filter(Dissolved.Oxygen..mg.L. < 0.5)


ysi_sub <- ysi_2021 %>% filter(Dissolved.Oxygen..mg.L.>0)

ysi_sub <- ysi_sub %>% dplyr::select(Date, Folder, Dissolved.Oxygen..mg.L.,Temperature..C.)
ysi_sub <- ysi_sub %>% gather("toplot","value", Dissolved.Oxygen..mg.L.,Temperature..C.)

ggplot(ysi_sub, aes(x=value, y=-Folder, col=toplot))+
  facet_wrap(Date~.,)+
  geom_point()+
  geom_line(orientation="y")+
  theme_bw()+
  ylab("Depth(m)")+
  xlab("Value - Temp (C) or DO (mg/L)")+
  scale_color_manual(values=c("red","blue"))+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave("Figures/YSI_MO_DO_mgL_Temp_C_facet_2021.pdf", width = 11, height = 11)


