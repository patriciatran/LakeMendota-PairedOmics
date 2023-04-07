#  Figure 1 :  Environmental Data Plots

library(dplyr)
library(readr)
library(ggplot2)
library(akima)
library(lubridate)
library(fields)
#install.packages("viridis")
library(viridis)
library(tidyverse)
library(ggrepel)

df <- list.files(path="Data/2020/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

df.2 <- df %>%filter(!is.na(sampledate))


df.2$sampledate <- mdy(df.2$sampledate)

df.2 <- df.2 %>% mutate(sampledata_depth = paste(sampledate, depth))

sample.metadata <- read.csv("Data/Samples-mendota.csv")
unique(sample.metadata$Date.YYYY.MM.DD)

list.of.sequenced.data <- unique(sample.metadata$Date.YYYY.MM.DD)
sequenced.dates <- paste0(list.of.sequenced.data,"|", collapse="")
nchar(sequenced.dates)

sequenced.dates <- substring(sequenced.dates, 1, nchar(sequenced.dates)-1)

df.2$sequenced <- ifelse(grepl(sequenced.dates, 
                               df.2$sampledate),
                         "sequenced date", "not sequenced")

unique(df.2$sequenced)



## Quick plot of all variables:

df.3 <- df.2 %>% gather("Variable","Value",wtemp:turb_fnu)
df.3$mixing <- ""


for (i in 1:nrow(df.3)){
  if (as.Date(df.3$sampledate[i]) > ymd("2020-10-18")){
    df.3$mixing[i] <- "Mixed"
  }
  else{
    df.3$mixing[i] <- "Stratified"
  }
print(i)
  }

df.3$mixing <- factor(df.3$mixing, levels = c("Stratified", "Mixed"))

my_theme <- theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(size=12),
        strip.background = element_blank())

ggplot(df.3)+
  geom_point(aes(x=Value, y=-depth, color=mixing), size=2, alpha=0.7)+
  facet_wrap(Variable~., scales="free", ncol=4)+
  my_theme+
  ylab("Depth (m)")+
  scale_color_manual(values=c("#547D6D","#55DDE0"))

  
ggsave("Figures/2020_profiles.plots.together.pdf", width = 8.5, height = 11, units="in")

below.20m <- df.3 %>% filter(depth >= 20)

df.3


ggplot(df.3 %>% filter(Variable %in% c("chlor_rfu", "do_raw", "ph", "turb_fnu")))+
  geom_point(aes(x=Value, y=-depth, color=sequenced, shape=mixing), size=1, alpha=0.9)+
  geom_line(aes(x=Value, y=-depth, 
                group=interaction(sampledate,mixing), 
                color=sequenced), orientation="y")+
  facet_wrap(Variable~., scales="free", ncol=4)+
  my_theme+
  ylab("Depth (m)")+
  scale_color_manual(values=c("grey","red"))+
  theme(legend.position="bottom")

# Figure 1C:
ggsave("Figures/Env_data_profile_subset.pdf", 
       height = 3.5, width = 5)

## Table for paper:
head(df.2)
TableSamples_metadata <- df.2 %>% select(sampledate, wtemp, do_raw, ph, depth) %>%
  group_by(sampledate) %>% mutate(avg_tmp = mean(wtemp),
                                         min_tmp = min(wtemp),
                                         max_temp = max(wtemp),
                                         avg_DO = mean(do_raw),
                                         min_DO = min(do_raw),
                                         max_DO = max(do_raw),
                                         avg_ph = mean(ph)) %>%
  select(-depth, -wtemp,-do_raw, -ph) %>%
  distinct()

TableSamples_metadata

write.table(x = TableSamples_metadata, file = "../../Data/Table_Samples_Envr_metadata.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

#2020 and 2021:
head(df.2)
df.2 <- df.2 %>% filter(sampledate >= as.Date("2020-01-01") & sampledate <= as.Date("2020-12-31"))
df.2 <- df.2 %>% select(sampledate, depth, wtemp, do_sat, do_raw)



# Make the heatmap
heatmap.data <- interp(x = as.Date(df.2$sampledate), y = -(df.2$depth), z = df.2$wtemp, duplicate = "strip")

## Sample points ## 
## Samples
samples.list <- read.csv("Data/2020-Lake-Mendota-Samples-PQT.csv")
samples.list

samples.list$Date <- mdy(samples.list$Date)


samples.list$labelling <- paste0(samples.list$Sample.Type,"_",samples.list$Replicate)

ggplot(samples.list, aes(x=Date, y=-Depth..m.))+
  ylim(c(-25,0))+
  geom_point(color = "black", shape=4, size=1)+
  #geom_text_repel()+
  theme_bw()+
  ylab("Depth (m)")+
  xlab("2020, 2021")+
  ggtitle("Collection dates and depths in Lake Mendota Summer 2020 & 2021")


# viruses 
samples.virus <- samples.list %>% 
  filter(Sample.Type == "Virus") %>% 
  select(Date, Depth..m.) %>% unique()


## LTER MO:


# Plotting the Secchi Data 2020 and 2021
secchi <- read.csv("Data/Secchi-2020-2021.csv")
avg_secchi_df <- secchi %>% group_by(Date) %>% mutate(avg_secchi = mean(Secchi.Depth.meters)) %>%
  select(Date, avg_secchi)
avg_secchi_df$Date <- mdy(avg_secchi_df$Date)
avg_secchi_df$Date

## Temp plot
par(mar = c(4.5,3,2,0.5))
unique.dates <- list.of.sequenced.data
#heat.colors()


list.of.sequenced.data

# FIGURE 1:
pdf("Figures/Temp2020.pdf", width = 11, height = 8.5)
image.plot(heatmap.data, axes = F, col = viridis(20),zlim=c(0,30))
axis(side = 1, at= ymd(unique.dates), col = "black", labels=unique.dates, las=2)
axis(side = 2, at = seq(from = -25, to = 0, by = 5), labels = T, las = 1)
title("Temperature in Lake Mendota 2020")
points(x=df.2$sampledate, y=-(df.2$depth)) # YSI profile
points(x=samples.list$Date, y=-(samples.list$Depth..m.), pch=19, col="red") # Patricia's samples
points(x=samples.virus$Date, y=-(samples.virus$Depth..m.), pch=0, cex=2, col="red") # Patricia's samples
#points(x=MGMT$Sample.date, y=-(MGMT$Depth..m.), pch = 1, lwd= 2, cex=3, col = "white") #MG and MT

#abline(h=-12, col="grey") # LTER MO
#abline(v=as.Date("2020-09-02"), col="blue" ) # Ben's start incub
#abline(v=as.Date("2020-10-10"), col="blue") # Ben end incub
points(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=18, col="black") #Secchi
lines(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=1, col="black") #Secchi
dev.off()

## DO sat:
df.2$do_raw[df.2$do_raw<0] <- 0

for (i in 1:nrow(df.2)){
  if(is.na(df.2$do_raw[i])){
    df.2$do_raw[i] <- df.2$do_raw_adj[i]
  }
}


heatmap.data2 <- interp(x = as.Date(df.2$sampledate), y = -(df.2$depth), z = df.2$do_raw, duplicate = "strip")

pdf("Figures/DO2020.pdf", width = 11, height = 8.5)
par(mar = c(4.5,3,2,0.5))
image.plot(heatmap.data2, axes = F, col = viridis(20),zlim=c(0,15))
axis(side = 1, at= ymd(unique.dates), col = "black", labels=unique.dates, las=2)
axis(side = 2, at = seq(from = -25, to = 0, by = 5), labels = T, las=1)
title("Dissolved dissolved oxygen (mg/L) in Lake Mendota 2020")
points(x=df.2$sampledate, y=-(df.2$depth), col="black")
points(x=samples.list$Date, y=-(samples.list$Depth), pch=19, col="red") # Patricia's samples
points(x=samples.virus$Date, y=-(samples.virus$Depth), pch=0, cex=2, col="red") # Patricia's samples
#points(x=MGMT$Sample.date, y=-(MGMT$Depth..m.), pch = 1, lwd= 2, cex=3, col = "white") #MG and MT

abline(h=-12, col="grey") # 12 meter for the LTER MO
#abline(v=as.Date("2020-09-02"), col="blue" ) # Ben start incub
#abline(v=as.Date("2020-10-10"), col="blue") # Ben stop incub
points(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=18, col="black") #Secchi
lines(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=1, col="black") #Secchi
# Legend
# Red dots: Samples Bacterial (Patricia)
# Red squares: Samples viral (Patricia)
# Black dots: YSI profile
# Green line: 12m integrated samples
# Blue vertical line: Ben's start of Incubation
dev.off()

# Temp 2021
ysi2021 <- read.csv("Data/2021/cat_YSI_MO_2021.csv")
ysi2021$sampledate <- mdy_hm(ysi2021$Timestamp)
ysi2021$sampledate <- date(ysi2021$sampledate)
ysi2021 <- ysi2021 %>% select(sampledate, Folder, Temperature..C., Dissolved.Oxygen..mg.L.)

colnames(ysi2021) <- c("sampledate", "depth","wtemp","DO_mgL")

# Temp heatmap 2021

heatmap.data <- interp(x = as.Date(ysi2021$sampledate), y = -(ysi2021$depth), z = ysi2021$wtemp, duplicate = "strip")

pdf("Figures/Temp2021.pdf", width = 11, height = 8.5)
image.plot(heatmap.data, axes = F, col = viridis(20),zlim=c(0,30))
axis(side = 1, at= ymd(unique.dates), col = "black", labels=unique.dates, las=2)
axis(side = 2, at = seq(from = -25, to = 0, by = 5), labels = T, las = 1)
title("Temperature in Lake Mendota 2021")
points(x=ysi2021$sampledate, y=-(ysi2021$depth)) # YSI profile
points(x=samples.list$Date, y=-(samples.list$Depth..m.), pch=19, col="red") # Patricia's samples
points(x=samples.virus$Date, y=-(samples.virus$Depth..m.), pch=0, cex=2, col="red") # Patricia's samples
#points(x=MGMT$Sample.date, y=-(MGMT$Depth..m.), pch = 1, lwd= 2, cex=3, col = "white") #MG and MT
abline(h=-12, col="grey") # LTER MO
points(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=18, col="black") #Secchi
lines(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=1, col="black") #Secchi
dev.off()

# DO 2021 
ysi2021$DO_mgL[ysi2021$DO_mgL<0] <- 0
heatmap.data2 <- interp(x = as.Date(ysi2021$sampledate), y = -(ysi2021$depth), z = ysi2021$DO_mgL, duplicate = "strip")

unique.dates <- unique(as.character(ysi2021$sampledate))

pdf("Figures/DO2021.pdf", width = 11, height = 8.5)
par(mar = c(4.5,3,2,0.5))
image.plot(heatmap.data2, axes = F, col = viridis(20),zlim=c(0,15))
axis(side = 1, at= ymd(unique.dates), col = "black", labels=unique.dates, las=2)
axis(side = 2, at = seq(from = -25, to = 0, by = 5), labels = T, las=1)
title("Dissolved dissolved oxygen (mg/L) in Lake Mendota 2021")
points(x=ysi2021$sampledate, y=-(ysi2021$depth), col="black")
points(x=samples.list$Date, y=-(samples.list$Depth), pch=19, col="red") # Patricia's samples
#points(x=samples.virus$Date, y=-(samples.virus$Depth), pch=0, cex=2, col="red") # Patricia's samples
#points(x=MGMT$Sample.date, y=-(MGMT$Depth..m.), pch = 1, lwd= 2, cex=3, col = "white") #MG and MT

abline(h=-12, col="grey") # 12 meter for the LTER MO
points(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=18, col="black") #Secchi
lines(x=avg_secchi_df$Date, y=-(avg_secchi_df$avg_secchi), pch=1, col="black") #Secchi

dev.off()

# Plotting the PAR data
PAR <- read.csv("Data/PAR2021.csv")

ggplot(PAR, aes(x=Date, y=-Depth, fill=PAR, group_by(Date)))+
  geom_tile(aes(size=PAR), shape=21, colour="white", size=0.5)+
  #geom_line(orientation = "y")+
  theme_bw()+
  ylim(-20,0)+
  scale_fill_gradient(low="red",high="yellow")+
  theme(axis.text.x = element_text(angle=90))+
  theme(panel.grid = element_blank())
  
ggsave("Figures/PAR2021.pdf")

