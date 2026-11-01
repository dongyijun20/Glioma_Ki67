#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root<-"/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/20261031_count_sensitivity"
out_deliv <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/delivery_figures_20261101"
dir.create(out_deliv, recursive=TRUE, showWarnings=FALSE)
g<-read.csv(file.path(root,"donor_global_thresholds.csv"));l<-read.csv(file.path(root,"donor_lineage_thresholds.csv"));e<-read.csv(file.path(root,"donor_ecotype_thresholds.csv"))
lev<-c("MKI67_ge1","MKI67_ge2","Phase_only","Intersection_ge1","Union_ge1","Union_ge2");labs<-c("MKI67 >=1","MKI67 >=2","Phase only","Intersection","Union >=1","Union >=2")
g$definition<-factor(g$definition,levels=lev,labels=labs);l$definition<-factor(l$definition,levels=lev,labels=labs);e$definition<-factor(e$definition,levels=lev,labels=labs);e$ecotype<-factor(e$ecotype,levels=c("Tumor","Immune","Glial","Stromal"))
pal<-c("Tumor"="#D96B4A","Immune"="#3C8DAD","Glial"="#64A876","Stromal"="#9B7EBD")
theme_ms <- theme_classic(base_size=9, base_family="Arial") +
  theme(plot.title=element_text(face="bold", size=10, family="Arial"),
        axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=9, family="Arial"),
        legend.text=element_text(size=8, family="Arial"),
        strip.text=element_text(size=8, family="Arial", face="bold"))
A<-ggplot(g,aes(definition,global_fraction*100))+
  geom_boxplot(width=.62,outlier.shape=NA,fill="#E9EEF2",color="#4B5563")+
  geom_jitter(width=.12,size=1.1,alpha=.5,color="#3C8DAD")+
  labs(title="A  Global donor-level fraction",x=NULL,y="Criterion-positive cells (%)")+
  theme_ms+theme(axis.text.x=element_text(angle=25,hjust=1))
B<-ggplot(l,aes(definition,within_lineage_fraction*100))+
  geom_boxplot(width=.62,outlier.shape=NA,fill="#F2EEE8",color="#4B5563")+
  geom_jitter(width=.12,size=.7,alpha=.4,color="#D96B4A")+
  facet_wrap(~lineage,nrow=1)+
  labs(title="B  Within-lineage positive fraction",x=NULL,y="Positive cells within lineage (%)")+
  theme_ms+theme(axis.text.x=element_text(angle=25,hjust=1))
C<-ggplot(e,aes(definition,fill=ecotype))+
  geom_bar()+scale_fill_manual(values=pal)+
  scale_y_continuous(breaks=seq(0,38,10),limits=c(0,38))+
  labs(title="C  Four-compartment donor assignment",x=NULL,y="Donors (n=38)",fill="Dominant compartment")+
  theme_ms+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="bottom")
top<-A+C+plot_layout(nrow=1,widths=c(1,1.15))
fig<-top/B+plot_layout(heights=c(1,1.15))
cairo_pdf(file.path(root,"FigureS7_Ki67_definition_comparison.pdf"), width=11.5, height=7.5, family="Arial", bg="white")
print(fig); dev.off()
cairo_pdf(file.path(out_deliv,"FigureS7_Ki67_definition_comparison.pdf"), width=11.5, height=7.5, family="Arial", bg="white")
print(fig); dev.off()
png(file.path(out_deliv,"FigureS7_Ki67_definition_comparison_600dpi.png"), width=11.5, height=7.5, units="in", res=600, type="cairo", bg="white", family="Arial")
print(fig); dev.off()
cat("FIG_S7_DONE\n")
