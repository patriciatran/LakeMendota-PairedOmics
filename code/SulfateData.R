# Sulfate  and Sulfide Data:
library(tidyverse)
library(lubridate)
library(tidyverse)
library(scales)

# Load tables:
brock <- read.csv("Data/Brockdata-sulfate-sulfide.csv")

ggplot(brock, aes(x=Concentration.in.uM, y=Depth..m., col=Dataset))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values=c("red","blue"))+
  ylab("Depth (m)")+
  xlab("Concentration (uM)")+
  theme(panel.grid = element_blank())+
  ggtitle("September 02, 1979 (Brock 1982)")+
  geom_hline(yintercept = -12)+
  annotate("text", x=20, y=-11, label="Thermocline")+
  geom_hline(yintercept = -15, linetype="dashed")+
  annotate("text", x=20, y=-14, label="Oxycline")

ggsave("Figures/Sulfur-profile-uM-Brock.pdf", width = 6, height = 7)

sulfur.2019 <- read.csv("Data/sulfur-2019-09-20.csv")

sulfur.2019 %>% group_by(Date, Depth, Dataset) %>% mutate(average=mean(Concentration_uM)) %>%
ggplot()+
  geom_point(aes(x=average, y=-Depth, col=Dataset), size=2)+
  geom_point(aes(x=Concentration_uM, y=-Depth, col=Dataset), alpha=.5)+
  theme_bw()+
  scale_color_manual(values=c("red","blue"))+
  ylab("Depth (m)")+
  xlab("Concentration (uM)")+
  theme(panel.grid = element_blank())+
  ggtitle("September 20, 2019")

ggsave("Figures/Sulfur-profile-uM-Tran-2019.pdf", width = 6, height = 7)

delfino_lee <- read.csv("Data/delfino_Lee_sulfur.csv")

delfino_lee %>%
  ggplot()+
  geom_point(aes(x=Sulfide_uM, y=Depth), size=2, color="blue")+
  facet_wrap(Date~.)+
  theme_bw()+
  ylab("Depth (m)")+
  xlab("Concentration (uM)")+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ggtitle("Delfino & Lee 1971")

ggsave("Figures/Sulfur_Delfino_Lee_1971.pdf", width = 6, height = 7)


# Sulfide and sulfate plots forr 2021 Lake Mendota
# Data collected by Samantha Bachand
# Sulfide data analysed by Patricia Tran and Samantha Bachand 
# Sulfate data analysed by Samantha Bachand

## Sulfate 2021:

sulfate.df <- read.csv("~/Downloads/Summary-Sulfate-2021 - Sheet1 (1).csv")


sulfate.df$Date_collected <- mdy(sulfate.df$Date_collected)

ggplot(sulfate.df, 
       aes(x=Sulfate_ppm, y=-Depth, col=Date_collected, group=Date_collected))+
  geom_point()+
  #geom_line()+
  #  facet_wrap(Date_collected~.)+
  theme_bw()+
  #  my_theme+
  ylab("Depth (m)")

## Sulfide

sulfide <- read.csv("Data/SulfideMeasurement_DataSheet_Template_2021-09-02_LR.csv")

sulfide <- sulfide %>% filter(!is.na(Sample.Depth...m.))
sulfide$Sample.Date <- mdy(sulfide$Sample.Date)

sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression <- replace(sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression, 
                                                                                                 sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression<0,0) 

head(sulfide)

ggplot(sulfide, aes(x=Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression,
                    y=-Sample.Depth...m.,
                    group_by(Sample.Date)))+
  geom_point(aes(col=as.factor(Sample.Date), shape=Replicate),
             size=10, alpha=0.70)+
  geom_line(aes(col=as.factor(Sample.Date)))+
  xlim(0,60)+
  theme_bw()+
  ylim(-25,0)+
  xlab("Sulfide uM")+
  ylab("Depth (m)")+
  # ggtitle("Sulfide in Lake Mendota 2021",
  #        subtitle="Measured 2021-09-02")+
  theme(text=element_text(size=20))+
  labs(color='Sample date') 

ggsave("Figures/sulfide2021.pdf", width = 11, height = 8.5, dpi=300)

ggplot(sulfide, aes(x=Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression,
                    y=-Sample.Depth...m.))+
  geom_point(aes(col=(as.factor(Sample.Date)), shape=Replicate),
             size=10, alpha=0.70)+
  xlim(0,40)+
  theme_bw()+
  ylim(-25,0)+
  xlab("Sulfide uM")+
  ylab("Depth (m)")+
  #ggtitle("Sulfide in Lake Mendota 2021",
  #        subtitle="Measured 2021-09-02")+
  facet_wrap(Sample.Date ~ ., scales="free")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(vjust =, hjust=0),
        text= element_text(size=15, colour="black"))+
  guides(fill=guide_legend(title="Date"))

ggsave("Figures/sulfide2021.by.date.pdf", width =14.5, height = 8, dpi=300)

## Both in the same figure:
sulfide.to.join <- sulfide %>% dplyr::select(Sample.Depth...m., Sample.Date, Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression, Replicate)
colnames(sulfide.to.join) <- c("Depth","Date_collected","Sulfide_uM", "Replicate")

sulfate.df.to.join <- sulfate.df %>% dplyr::select(Depth, Date_collected, Sulfate_ppm)
sulfate.df.to.join$Sulfate_uM <- (sulfate.df.to.join$Sulfate_ppm * 1000)/96.06

sulfur.data <- full_join(sulfide.to.join, sulfate.df.to.join)

sulfur.data <- sulfur.data %>% gather("Sulfur","Concentration",c(Sulfide_uM, Sulfate_uM))



ggplot(sulfur.data %>% 
         filter(!is.na(Concentration), 
                !is.na(Date_collected)), 
       aes(x=Concentration, 
           y=-Depth, 
           col=Sulfur,
           group=interaction(Sulfur, Date_collected)))+
  geom_point()+
  facet_wrap(Date_collected~., scales="free")+
  theme_bw()+
  #  my_theme+
  ylab("Depth (m)")+
  #ggtitle("Sulfur concentrations in Lake Mendota in 2021")+
  scale_color_manual(values=c("red","blue"))+
  #  geom_line(orientation="y")+
  geom_smooth(orientation="y")+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position="bottom")+
  xlab("Concentration in uM")

# Supp Figure.
ggsave("Figures/Sulfate_Sulfide_uM_Mendota_2021.pdf", width = 5, height = 4)






