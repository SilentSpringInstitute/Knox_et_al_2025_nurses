#Application of a non-targeted biomonitoring method to characterize occupational chemical exposures of women nurses relative to office workers 
#Kristin E. Knox, Dimitri Abrahamsson, Jessica Trowbridge, June-Soo Park, Miaomiao Wang, Erin Carrera, Lisa Hartmayer, Rachel Morello-Frosch, R. A. Rudel

# This script makes the figures for the paper
# AUTHOR: K Knox
# WRITTEN IN R VERSION: R version 4.4.0 (2024-04-24)


library(tidyverse)
library(patchwork)


# Set the working directory
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)
setwd('..') 


options(scipen=999)


########################################

# Figure 3

# load in the positive and negative data for the volcano plots
volcano_neg <- read_rds("data/volcano_neg.rds")
volcano_pos <- read_rds("data/volcano_pos.rds")

negplot <-
ggplot(data=volcano_neg, aes(x=log2fold, y=neg_log10_padj)) +
  geom_point() +
  geom_hline(yintercept=1.30103, colour="red", linetype="dashed", size=.6) + #because -log10(.05) #1.30103
  scale_y_continuous(expand=c(0,0), limits=c(0,8)) +
  scale_x_continuous(limits=c(-4,4), breaks=c(-4, -3, -2, -1, 0, 1, 2, 3, 4)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=10)) +
  theme(axis.text.x = element_text(size=10)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.title.y = element_text(size=10)) +
  theme(axis.title.x = element_text(size=10))+
  ylab("-log10 adj. p-value") +
  xlab("log 2 fold change") + 
  ggtitle("ESI-")


posplot <-
ggplot(data=volcano_pos, aes(x=log2fold, y=neg_log10_padj)) +
  geom_point() +
  geom_hline(yintercept=1.30103, colour="red", linetype="dashed", size=.6) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8)) +
  scale_x_continuous(limits=c(-4,4), breaks=c(-4, -3, -2, -1, 0, 1, 2, 3, 4)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=10)) +
  theme(axis.text.x = element_text(size=10)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.title.y = element_text(size=10)) +
  theme(axis.title.x = element_text(size=10))+
  ylab("-log10 adj. p-value") +
  xlab("log 2 fold change") + 
  ggtitle("ESI+")


(negplot | posplot) + plot_annotation(title="Figure 3") & theme(plot.title = element_text(hjust = 0.5, size=10))
#ggsave(filename='figures/final paper figures/figure3.svg', width=7, height=7, units="in", dpi=300)


########################################

# Figure 4

#read in the negative data for the plots
plot_neg <-read_rds("data/data_to_plot_neg.rds")

# read in the positive file for the plots
plot_pos <-read_rds("data/data_to_plot_pos.rds")


disinfect_4hydroxy<-
  ggplot(filter(plot_neg, feature_ID=="feature_452"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
    geom_point(position=position_jitter(h=NULL, w=0.075)) +
    scale_shape_manual(values=c(16,1)) +
    scale_y_continuous(trans="log10") +
    scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
    theme(plot.title = element_text(hjust=0.5, size=8)) +
    theme(axis.text.x = element_text(size=7)) +
    theme(axis.text.y = element_text(size=7)) +
    theme(axis.title.y = element_text(size=7)) +
    ylab("mean intensity") +
    xlab("") + 
    theme(legend.position = "none")+
    ggtitle("4-hydroxyquinoline")


disinfect_salicylic <-
  ggplot(filter(plot_neg, feature_ID=="feature_319"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
    geom_point(position=position_jitter(h=NULL, w=0.075)) +
    scale_shape_manual(values=c(16,1)) +
    scale_y_continuous(trans="log10") +
    scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
    theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
    theme(plot.title = element_text(hjust=0.5, size=8)) +
    theme(axis.text.x = element_text(size=7)) +
    theme(axis.text.y = element_text(size=7)) +
    theme(axis.title.y = element_text(size=7)) +
    ylab("mean intensity") +
    xlab("") + 
    theme(legend.position = "none")+
    ggtitle("salicylic acid")


theophylline <-
ggplot(filter(plot_neg, feature_ID=="feature_1112"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("theophylline")


acetominophen <-
ggplot(filter(plot_pos, feature_ID=="feature_514"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("acetominophen")


tridecanedioic <-
ggplot(filter(plot_neg, feature_ID=="feature_2688"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("tridecanedioic acid")


phthalate<-
ggplot(filter(plot_neg, feature_ID=="feature_5148"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("dicyclohexyl phthalate")


pfas<-
  ggplot(filter(plot_neg, feature_ID=="feature_7403"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("6:2 fluorotelomer sulfonic acid")


level1 <- (pfas | disinfect_salicylic | theophylline | acetominophen)  + plot_annotation(title="Level 1") & theme(plot.title = element_text(hjust = 0.5, size=8))

level2 <- (tridecanedioic + plot_spacer())  + plot_annotation(title="Level 2") & theme(plot.title = element_text(hjust = 0.35, size=8))

atc <- (phthalate | disinfect_4hydroxy) + plot_annotation(title="          Annotation Not Confirmed") & theme(plot.title = element_text(hjust = 0.5, size=8))

(wrap_elements(level1)) / (wrap_elements(level2) | wrap_elements(atc)) + plot_annotation(title="Figure 4") & theme(plot.title = element_text(hjust = 0.5, size=8))
#ggsave(filename='figures/final paper figures/figure4.svg', width=7, height=8.499, units="in", dpi=300)


#################################################

# Make Figure 5

#################################################

tris<-
  ggplot(filter(plot_pos, feature_ID=="feature_6069"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("\n\ntris(2-butoxyethyl) phosphate")


phthal<-
  ggplot(filter(plot_pos, feature_ID=="feature_1587"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("multiple phthalates")


quin8<-
  ggplot(filter(plot_neg, feature_ID=="feature_452"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("\n\n8-Hydroxyquinoline")


methpred<-
  ggplot(filter(plot_pos, feature_ID=="feature_5652"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("methylprednisolone")


cet<-
  ggplot(filter(plot_pos, feature_ID=="feature_5890"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("cetirizine dihydrochloride")


octyl_gem<-
  ggplot(filter(plot_neg, feature_ID=="feature_2824"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("octylparaben / gemfibrozil")

hydroxy_methylpropan <-
  ggplot(filter(plot_pos, feature_ID=="feature_1634"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  #  ggtitle("2-Hydroxy-1-[4-(2-hydroxyethoxy)phenyl]-2-methylpropan-1-one")
  ggtitle("2-hydroxy-1-\n[4-(2-hydroxyethoxy)phenyl]\n-2-methylpropan-1-one")


hexamethyltetralin <-
  ggplot(filter(plot_pos, feature_ID=="feature_2386"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("\n6-Acetyl-1,1,2,4,4,7\n-hexamethyltetralin")


(tris | quin8 | hydroxy_methylpropan | hexamethyltetralin) / (methpred | cet | octyl_gem | phthal)+ plot_annotation(title="Figure 5 - Level 2 Without Spectra Tentative Identifications") & theme(plot.title = element_text(hjust = 0.5, size=8))

#ggsave(filename='figures/final paper figures/figure5.svg', width=7, height=5, units="in", dpi=300)


#################################################

# Make supplemental figures

#################################################

# Figure S3 - diet


ah_acid<-
  ggplot(filter(plot_neg, feature_ID=="feature_1434"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("alpha-hydroxyhippuric acid")


eic_acid<-
  ggplot(filter(plot_neg, feature_ID=="feature_4434"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("eicosapentaenoic acid")


hepe<-
  ggplot(filter(plot_neg, feature_ID=="feature_4884"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("5-HEPE")


caffeine<-
  ggplot(filter(plot_pos, feature_ID=="feature_1137"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("caffeine*")


parax<-
  ggplot(filter(plot_pos, feature_ID=="feature_945"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("paraxanthine")


trypt<-
  ggplot(filter(plot_pos, feature_ID=="feature_1302"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("tryptophan")


pip<-
  ggplot(filter(plot_pos, feature_ID=="feature_3090"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("piperine*")


(caffeine | parax | trypt | pip) / (ah_acid | eic_acid | hepe) + plot_annotation(title="Figure S3 - Confirmed and Tentative Identifications from Dietary Sources") & theme(plot.title = element_text(hjust = 0.5, size=8))

#ggsave(filename='figures/final paper figures/figureS3.svg', width=7, height=8.666, units="in", dpi=300)


# Figure S4 - endogenous
hode<-
  ggplot(filter(plot_neg, feature_ID=="feature_4241"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("9-HODE")


dihete<-
  ggplot(filter(plot_neg, feature_ID=="feature_5295"), aes(x=worker, y=mean_intensity, color=worker, shape=detect)) +
  geom_point(position=position_jitter(h=NULL, w=0.075)) +
  scale_shape_manual(values=c(16,1)) +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(limits=c("nurse","office"), labels=c("nurses","office\nworkers")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank()) +
  theme(panel.border=element_rect(size=0.5, linetype="solid", colour="black", fill="transparent")) +
  theme(plot.title = element_text(hjust=0.5, size=8)) +
  theme(axis.text.x = element_text(size=7)) +
  theme(axis.text.y = element_text(size=7)) +
  theme(axis.title.y = element_text(size=7)) +
  ylab("mean intensity") +
  xlab("") + 
  theme(legend.position = "none")+
  ggtitle("5,12-DiHETE")


(hode |dihete)+ plot_annotation(title="Figure S4 - Tentative Endogenous Molecules") & theme(plot.title = element_text(hjust = 0.5, size=8))

#ggsave(filename='figures/final paper figures/figureS4.svg', width=7, height=5, units="in", dpi=300)


