metabolic_to_plot <- read.csv("Data/metabolic_mt_mg_to_plot.csv", header=TRUE)

metabolic_to_plot$Date <- mdy(metabolic_to_plot$Date)
metabolic_to_plot <- metabolic_to_plot %>% mutate(mixing = ifelse(Date == "2020-10-19","Mixed","Stratified"))

unique(metabolic_to_plot$Reaction)

library(ggpubr)



for (i in 1:length(unique(metabolic_to_plot$Reaction))){
  
  
  subset_to_plot <- metabolic_to_plot %>%
    filter(Reaction == unique(metabolic_to_plot$Reaction)[i])
  
  plot1 <- ggplot(subset_to_plot, aes(x=paste(Date,Depth), 
                                      y=Percent.of.Community,
                                      col=Type))+
    geom_point()+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90))+
    geom_vline(xintercept =13.5)+
    ggtitle(unique(metabolic_to_plot$Reaction)[i],
            subtitle=paste0(subset_to_plot$Number.of.Genomes[1], " genomes"))+
    xlab("Sample Date and Depth (m)")
  
  plot2 <- ggplot(subset_to_plot, aes(x=mixing, y=Percent.of.Community, col=Type))+
    geom_violin()+
    ggtitle(unique(metabolic_to_plot$Reaction)[i],
            paste0(subset_to_plot$Number.of.Genomes[1], " genomes"))+
    theme_bw()
  
  ggarrange(plot1, plot2, labels=c("A","B"), common.legend= TRUE)
  
  ggsave(paste0("Figures/",unique(metabolic_to_plot$Reaction)[i],"panel_abund_expr.pdf"),
         width = 11,
         height = 6)
  
  print(i)
  
}

