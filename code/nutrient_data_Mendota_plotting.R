# Plot the most recent LTER MO nutrient data 2013 - 2018:
library(tidyverse)
library(lubridate)

nutrient.data <- read.csv("Data/robin_MEMO_nutrients_2013-2018.tsv", 
                          sep="\t", header=TRUE)


nutrient.data$Date <- str_replace(nutrient.data$Sample.Name, "ME","")
unique(nutrient.data$Date)
nutrient.data$Date <- str_replace(nutrient.data$Date, "D.*","")
unique(nutrient.data$Date)
nutrient.data$Date <- str_replace(nutrient.data$Date, "s.*","")

nutrient.data$Date <- dmy(nutrient.data$Date)
#2013, 2018, 2014 have some depth discrete samples for nutrients.

nutrient.data$Date[342]<-dmy("01Dec2017")
nutrient.data$Date[343] <- dmy("01Dec2017")

my_theme <- theme_bw()+
  theme(panel.border = element_rect(color="black"),
        text = element_text(size=15, colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(color="black"))


ggplot(nutrient.data %>% filter(Biol.Rep != "S") %>% gather("name","value",3:10), 
       aes(x=Date, y=value)) +
  #geom_jitter(width = 0.2) +
  #geom_line(aes(group = TP.ug.L), stat = "summary", fun = mean) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, 
  width = 0.1, color="grey")+
  geom_point(aes(fill=Biol.Rep), alpha=0.5, pch=21, size=2)+
  my_theme+
  facet_grid(name~., scales="free")+
  xlab("Date")+
  ylab("Concentration")

ggsave("Figures/Field_Nutrients_2014-2018.pdf", width = 8.5, height = 11)
