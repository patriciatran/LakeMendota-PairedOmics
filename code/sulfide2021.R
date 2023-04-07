# Load packages
library(lubridate)
library(tidyverse)
library(scales)

# Load Data:
sulfide <- read.csv("Data/SulfideMeasurement_DataSheet_Template_2021-09-02_LR.csv")


sulfide <- sulfide %>% filter(!is.na(Sample.Depth...m.))
sulfide$Sample.Date <- mdy(sulfide$Sample.Date)
sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression <- replace(sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression, 
                sulfide$Concentration.of.H2S.in.uM.if.not.taking.the.50uM.point.for.linear.regression<0,0) 


# Plot
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


# Save
ggsave("Figures/sulfide2021.pdf", width = 11, height = 8.5, dpi=300)


# Plot
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

# Save
ggsave("Figures/sulfide2021.by.date.pdf", width =14.5, height = 8, dpi=300)
