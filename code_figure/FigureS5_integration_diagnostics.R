#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/20261031_integration_metrics"
out_deliv <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/delivery_figures_20261101"
dir.create(out_deliv, recursive=TRUE, showWarnings=FALSE)
g <- read.csv(file.path(root,"integration_metrics_global.csv"))
w <- read.csv(file.path(root,"integration_metrics_by_lineage.csv"))
cols <- c(PCA="#3C8DAD",Harmony="#D96B4A")
g$embedding <- factor(g$embedding,levels=c("pca","harmony"),labels=c("PCA","Harmony"))
w$embedding <- factor(w$embedding,levels=c("pca","harmony"),labels=c("PCA","Harmony"))
theme_ms <- theme_classic(base_size=9, base_family="Arial") +
  theme(plot.title=element_text(face="bold", size=10, family="Arial"),
        axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=9, family="Arial"),
        legend.text=element_text(size=8, family="Arial"),
        legend.title=element_text(size=8, family="Arial"))
A <- ggplot(g,aes(embedding,batch_iLISI_normalized,group=1))+
  geom_line(color="grey60",linewidth=.6)+geom_point(aes(color=embedding),size=3)+
  scale_color_manual(values=cols,guide="none")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2))+
  labs(title="A  Global dataset mixing",x=NULL,y="Normalized inverse-Simpson index")+theme_ms
B <- ggplot(g,aes(embedding,batch_silhouette_mean,group=1))+
  geom_line(color="grey60",linewidth=.6)+geom_point(aes(color=embedding),size=3)+
  scale_color_manual(values=cols,guide="none")+
  labs(title="B  Dataset silhouette",x=NULL,y="Mean silhouette width")+theme_ms
C <- ggplot(w,aes(lineage,batch_iLISI_median,color=embedding,group=embedding))+
  geom_line(position=position_dodge(width=.25),linewidth=.6)+
  geom_point(position=position_dodge(width=.25),size=2.5)+
  scale_color_manual(values=cols)+
  labs(title="C  Within-lineage dataset mixing",x=NULL,y="Median inverse-Simpson index",color=NULL)+
  theme_ms+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="top")
fig <- A+B+C+plot_layout(nrow=1,widths=c(1,1,1.6))
cairo_pdf(file.path(root,"FigureS5_integration_diagnostics.pdf"), width=10.5, height=4.2, family="Arial", bg="white")
print(fig); dev.off()
cairo_pdf(file.path(out_deliv,"FigureS5_integration_diagnostics.pdf"), width=10.5, height=4.2, family="Arial", bg="white")
print(fig); dev.off()
png(file.path(out_deliv,"FigureS5_integration_diagnostics_600dpi.png"), width=10.5, height=4.2, units="in", res=600, type="cairo", bg="white", family="Arial")
print(fig); dev.off()
cat("FIG_S5_DONE\n")
