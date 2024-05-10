library('LSD')
library(ggplot2)
library(plyr)
library(ggridges)
library(RColorBrewer)
library(dplyr)
library(ggprism)
library(ggpubr)
library(cowplot)
library(MASS)
library(viridis)
theme_set(theme_bw(base_size = 16))
library(grid)
library(gridExtra)
library(scales)
library(ggthemes)
library(ggridges)
library(pheatmap)

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

outdir = "output/main"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

## IM
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/RNA-seq_data/Challenge_2_RNA-Seq_Data.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)
colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(15,17,16,3,0);



x1 <- ggplot(dat,aes(x = Tools, y = Irreproducibility))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('IM') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 1.5), breaks = seq(0,1.5,0.5))

x2 <- ggplot(dat,aes(x = Protocols_Platforms, y = Irreproducibility))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('IM') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 1.5), breaks = seq(0,1.5,0.5)) 

p <- ggarrange(x1,x2,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_RNA-Seq_Data_IM",".pdf"), 
           plot = p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )




## ACVC
x3 <- ggplot(dat,aes(x = Tools, y = ACVC))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('') + ylab('ACVC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 6), breaks = seq(0,6,2)) 


x4 <- ggplot(dat,aes(x = Protocols_Platforms, y = ACVC))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('ACVC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 6), breaks = seq(0,6,2)) 


p <- ggarrange(x3,x4,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_RNA-Seq_Data_ACVC",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )


## CM
x5 <- ggplot(dat,aes(x = Tools, y = Consistency))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('') + ylab('CM') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0.4, 1), breaks = seq(0.4,1,0.2)) 

x6 <- ggplot(dat,aes(x = Protocols_Platforms, y = Consistency))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('CM') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0.4, 1), breaks = seq(0.4,1,0.2))

p <- ggarrange(x5,x6,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_RNA-Seq_Data_CM",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )






# ACC
x7 <- ggplot(dat,aes(x = Tools, y = ACC))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('') + ylab('ACC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(7, 10), breaks = seq(7,10,1))  

x8 <- ggplot(dat,aes(x = Protocols_Platforms, y = ACC))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('ACC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(7, 10), breaks = seq(7,10,1)) 


p <- ggarrange(x7,x8,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_RNA-Seq_Data_ACC",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )



## RE

x9 <- ggplot(dat,aes(x = Tools, y = Average.Resolution.Entropy))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('RE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 0.015), breaks = seq(0,0.015,0.005)) 

x10 <- ggplot(dat,aes(x = Protocols_Platforms, y = Average.Resolution.Entropy))  +
geom_boxplot(width = 0.8,outlier.shape = NA,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('RE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
scale_y_continuous(limits=c(0, 0.015), breaks = seq(0,0.015,0.005)) 

p <- ggarrange(x9,x10,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_RNA-Seq_Data_RE",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )





## Cell mixing experiment
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Cell_mixing_experiment/Challenge_2_Cell_mixing_experiment.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)
colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(16);


x11 <- ggplot(dat,aes(x = Tools, y = SCC))  + scale_color_viridis() + 
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1), breaks = seq(0.3,1,0.35), oob = rescale_none) 


x12 <- ggplot(dat,aes(x = Protocols_Platforms, y = SCC))  +
    stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1), breaks = seq(0.3,1,0.35), oob = rescale_none)


p <- ggarrange(x11,x12,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Cell_mixing_experiment_SCC",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )





x13 <- ggplot(dat,aes(x = Tools, y = MRD))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.05, 0.5), breaks = seq(0.05,0.5,0.15), oob = rescale_none) 


x14 <- ggplot(dat,aes(x = Protocols_Platforms, y = MRD))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.05, 0.5), breaks = seq(0.05,0.5,0.15), oob = rescale_none) 

p <- ggarrange(x13,x14,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Cell_mixing_experiment_MRD",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )





x15 <- ggplot(dat,aes(x = Tools, y = NRMSE))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('NRMSE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1.5), breaks = seq(0.3,1.5,0.4), oob = rescale_none) 


x16<- ggplot(dat,aes(x = Protocols_Platforms, y = NRMSE))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1.5), breaks = seq(0.3,1.5,0.4), oob = rescale_none) 

p <- ggarrange(x15,x16,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Cell_mixing_experiment_NRMSE",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )











## SIRV data
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/SIRV-set4_data/Challenge_2_SIRV_MRD_SCC_NRMSE_metrics.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)

colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(15,17,16,3);



x17 <- ggplot(dat,aes(x = Tools, y = SCC))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 

x18 <- ggplot(dat,aes(x = Protocols_Platforms, y = SCC))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
 ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 


p <- ggarrange(x17,x18,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_SIRV_SCC",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )






x19 <- ggplot(dat,aes(x = Tools, y = MRD))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=0.5) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 

x20 <- ggplot(dat,aes(x = Protocols_Platforms, y = MRD))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 


p <- ggarrange(x19,x20,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_SIRV_MRD",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )



x21 <- ggplot(dat,aes(x = Tools, y = NRMSE))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=0.5) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('NRMSE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 2.5), breaks = seq(0.3,2.5,1.1), oob = rescale_none) 


x22 <- ggplot(dat,aes(x = Protocols_Platforms, y = NRMSE))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 2.5), breaks = seq(0.3,2.5,1.1), oob = rescale_none) 


p <- ggarrange(x21,x22,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_SIRV_NRMSE",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )





## ERCC_PET_metric
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/SIRV-set4_data/Challenge_2_ERCC_PET_metric.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)

colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(15,17,16,3);



x23 <- ggplot(dat,aes(x = Tools, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=0.5) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('PET (ERCC)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none)

x24 <- ggplot(dat,aes(x = Protocols_Platforms, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('PET (ERCC)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 

p <- ggarrange(x23,x24,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_ERCC_PET",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )



## Long_SIRV_PET_metric
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/SIRV-set4_data/Challenge_2_Long_SIRV_PET_metric.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)

colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(15,17,16,3);



x25 <- ggplot(dat,aes(x = Tools, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=0.5) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('PET (Long SIRV)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none)

x26 <- ggplot(dat,aes(x = Protocols_Platforms, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('PET (Long SIRV)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 

p <- ggarrange(x25,x26,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Long_SIRV_PET",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )





## SIRV_PET_metric
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/SIRV-set4_data/Challenge_2_SIRV_PET_metric.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","CapTrap_PacBio","CapTrap_ONT","R2C2_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)

colorset = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(15,17,16,3);



x27 <- ggplot(dat,aes(x = Tools, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=0.5) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('PET (SIRV)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none)

x28 <- ggplot(dat,aes(x = Protocols_Platforms, y = PET))  +
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('PET (SIRV)') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.5), oob = rescale_none) 

p <- ggarrange(x27,x28,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_SIRV_PET",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )










## Simulation data SCC
dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Simulation_data/Challenge_2_Simulation_Metrics.txt',header = T, sep = "\t");
dat$Protocols_Platforms <- factor(dat$Protocols_Platforms,levels=c("cDNA_PacBio","cDNA_ONT","dRNA_ONT","cDNA_Illumina"))
dat$Tools <- factor(dat$Tools,levels=c("Bambu","FLAIR","FLAMES","IsoQuant","IsoTools","TALON","NanoSim","RSEM","kallisto-lr"));
head(dat)

#colorset1 = c("#C06636","#646E3B","#2B5851","#E8B960");
colorset1 = c("#C06636","#646E3B","#2B5851","#802417","#17486F","#508EA2","#E8B960");
shapeset = c(4,8);



x29 <- ggplot(dat,aes(x = Tools, y = SCC))  + scale_color_viridis() + 
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1), breaks = seq(0.3,1,0.35), oob = rescale_none) 



x30 <- ggplot(dat,aes(x = Protocols_Platforms, y = SCC))  +
    stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset1,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('SCC') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.3, 1), breaks = seq(0.3,1,0.35), oob = rescale_none) 


p <- ggarrange(x29,x30,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Simulation_Data_SCC",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )


## Simulation data MRD

x31 <- ggplot(dat,aes(x = Tools, y = MRD))  + 
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset1,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.04, 0.5), breaks = seq(0.04,0.5,0.23), oob = rescale_none) 


x32 <- ggplot(dat,aes(x = Protocols_Platforms, y = MRD))  +
    stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset1,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('MRD') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.04, 0.5), breaks = seq(0.04,0.5,0.23), oob = rescale_none) 
p <- ggarrange(x31,x32,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Simulation_Data_MRD",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )


## Simulation data NRMSE

x33 <- ggplot(dat,aes(x = Tools, y = NRMSE))  + 
  stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Protocols_Platforms),size=1) +
scale_colour_manual(values=alpha(colorset1,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
xlab('Tools') + ylab('NRMSE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.1, 1.2), breaks = seq(0.1,1.2,0.55), oob = rescale_none) 


x34 <- ggplot(dat,aes(x = Protocols_Platforms, y = NRMSE))  +
    stat_summary(fun.y = "mean", geom = "bar", position = "dodge", alpha = 0.8,fill="white",color="black",size=0.2) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width = 0.2,size=0.2) +
geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.5),aes(shape = Samples, colour = Tools),size=1) +
scale_colour_manual(values=alpha(colorset1,0.8))  + scale_shape_manual(values=shapeset)  +
theme_bw() + theme(panel.grid=element_blank()) +
ylab('NRMSE') +
theme(axis.text.x=element_text(angle = 90,hjust = 1)) + scale_y_continuous(limits=c(0.1, 1.2), breaks = seq(0.1,1.2,0.55), oob = rescale_none) 

p <- ggarrange(x33,x34,nrow = 1, align = "v",hjust=0)

ggsave(filename = paste0(outdir, "/fig_Simulation_Data_NRMSE",".pdf"), 
           plot =  p, 
           width=30, 
           height=15,
           units = "cm",
           dpi=300 )


## Ranking for each metrics

dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Ranking_for_each_metrics/tools.txt',header = FALSE, sep = "\t");
head(dat)

colorset = c("#873600","#BA4A00","#DC7633","#EDBB99");


x35 <- ggplot(dat,aes(x=V4,y=V1,color=V2)) + geom_point(size=4) +
#      scale_size(range=c(1,4),breaks=c(4,3,2,1),labels=c("1","2","3","4"),name="Ranking") +
    scale_colour_manual(values=colorset) +
    theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
    theme_prism(base_size = 5)+ 
   # theme(legend.position="none") +
    theme(axis.text.x=element_text(angle = 90,hjust = 1)) +
    scale_x_discrete(labels = c("IM_real_data","ACVC_real_data","CM_real_data","ACC_real_data","RE_real_data","MRD_cell_mixing","NRMSE_cell_mixing","SCC_cell_mixing","MRD_SIRV","NRMSE_SIRV","SCC_SIRV","PET_SIRV","PET_Long_SIRV","PET_ERCC","MRD_simulation","NRMSE_simulation","SCC_simulation"))


dat <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Ranking_for_each_metrics/protocols_platforms.txt',header = FALSE, sep = "\t");
head(dat)

x36 <- ggplot(dat,aes(x=V4,y=V1,color=V2)) + geom_point(size=4) +
#      scale_size(range=c(1,4),breaks=c(4,3,2,1),labels=c("1","2","3","4"),name="Ranking") +
    scale_colour_manual(values=colorset) +
    theme(axis.text.x=element_text(angle = 90,hjust = 1))  + 
    theme_prism(base_size = 5)+ 
   # theme(legend.position="none") +
    theme(axis.text.x=element_text(angle = 90,hjust = 1)) +
    scale_x_discrete(labels = c("IM_real_data","ACVC_real_data","CM_real_data","ACC_real_data","RE_real_data","MRD_cell_mixing","NRMSE_cell_mixing","SCC_cell_mixing","MRD_SIRV","NRMSE_SIRV","SCC_SIRV","PET_SIRV","PET_Long_SIRV","PET_ERCC","MRD_simulation","NRMSE_simulation","SCC_simulation"))



p1 <- ggarrange(x35,x36,ncol = 2, align = "v")

ggsave(filename = paste0(outdir, "/fig_Ranking",".pdf"), 
           plot =  p1, 
           width=32, 
           height=8,
           units = "cm",
           dpi=300 )


## Evaluation of quantification tools with respect to multiple transcript features

## isoform abundance
dat1 <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Evaluation_with_multiple_transcript_features/3_true_abund.tsv',header = F, sep = "\t");
head(dat1)

rownames(dat1) = c("RSEM","IsoTools","TALON","IsoQuant","Bambu","FLAIR")
colnames(dat1) = c("K1","K2","K3","K4");
head(dat1)

#color.key <- c("#CC0000","#FF3333","white","#3399FF","#3300CC")
pheatmap(dat1[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(5))


## Kvalue
dat2 <- read.table('Challenge2_Figures_Data/Evaluation_with_multiple_transcript_features/3_K_value.tsv',header = F, sep = "\t");
head(dat2)

rownames(dat2) = c("RSEM","IsoTools","TALON","IsoQuant","Bambu","FLAIR")
colnames(dat2) = c("K1","K2","K3","K4");
head(dat2)

color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(dat2[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE,  show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))



## Iso_num
dat3 <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Evaluation_with_multiple_transcript_features/3_num_isoforms.tsv',header = F, sep = "\t");
head(dat3)

rownames(dat3) = c("RSEM","IsoTools","TALON","IsoQuant","Bambu","FLAIR")
colnames(dat3) = c("K1","K2","K3","K4");
head(dat3)

color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(dat3[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F,show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))


## Exon_num
dat4 <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Evaluation_with_multiple_transcript_features/3_num_exons.tsv',header = F, sep = "\t");
head(dat4)

rownames(dat4) = c("RSEM","IsoTools","TALON","IsoQuant","Bambu","FLAIR")
colnames(dat4) = c("K1","K2","K3","K4");
head(dat4)

color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(dat4[,1:4], legend = T, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F,show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))


## Isoform length
dat5 <- read.table('/Users/rloving/Desktop/kallisto-lr/LRGASP/Challenge2_Figures_Data/Evaluation_with_multiple_transcript_features/3_isoform_length.tsv',header = F, sep = "\t");
head(dat5)

rownames(dat5) = c("RSEM","IsoTools","TALON","IsoQuant","Bambu","FLAIR")
colnames(dat5) = c("K1","K2","K3","K4");
head(dat5)

color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pheatmap(dat5[,1:4], legend = F, cluster_col = FALSE,  cluster_row = FALSE, show_colnames=F,show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))

x1<- pheatmap(dat1[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50));
x2 <- pheatmap(dat2[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))
x3 <- pheatmap(dat3[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE,  show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))
x4 <- pheatmap(dat4[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))
x5 <- pheatmap(dat5[,1:4], legend = F, cluster_col = FALSE, cluster_row = FALSE, show_colnames=F, show_rownames=T,scale = "row",color=colorRampPalette(c("#EAF2F8","#2980B9","#154360"))(50))

require(ggplotify)
x1 = as.ggplot(x1)
x2 = as.ggplot(x2)
x3 = as.ggplot(x3)
x4 = as.ggplot(x4)
x5 = as.ggplot(x5)



plot_tmp <- plot_grid(plot_grid(x1,x2,x3,x4,x5,ncol=5),nrow=1)

ggsave(filename = paste0(outdir, "/fig_heatmap_evaluation_multiple_transcript_features",".pdf"), 
           plot =  plot_tmp, 
           width=90, 
           height=18,
           units = "cm",
           dpi=300 )







