#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
root<-"/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/20261031_count_sensitivity"
out_deliv <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/delivery_figures_20261101"
dir.create(out_deliv, recursive=TRUE, showWarnings=FALSE)
x<-read.csv(file.path(root,"cellcycle_score_sensitivity.csv"))
corr<-read.csv(file.path(root,"cellcycle_score_lineage_correlations.csv"))
x$score_cutoff<-factor(x$score_cutoff,levels=c(.05,.10,.20))
x$count_cutoff<-factor(x$count_cutoff,levels=c(1,2),labels=c("MKI67 >= 1","MKI67 >= 2"))
corr$score_cutoff<-as.numeric(sub(".*score_","",corr$definition))/100
corr$count_cutoff<-ifelse(grepl("count_ge2",corr$definition),"MKI67 >= 2","MKI67 >= 1")
base<-data.frame(definition=rep("Union_count_ge1_score_010",5),lineage=c("Tumor","Myeloid","Stromal","Lymphoid","Glial"),n=38,spearman=1,score_cutoff=0.10,count_cutoff="MKI67 >= 1")
corr<-rbind(corr,base)
corr$count_cutoff<-factor(corr$count_cutoff,levels=c("MKI67 >= 1","MKI67 >= 2"))
theme_ms <- theme_classic(base_size=9, base_family="Arial") +
  theme(plot.title=element_text(face="bold", size=10, family="Arial"),
        axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=9, family="Arial"),
        legend.text=element_text(size=8, family="Arial"),
        strip.text=element_text(size=8, family="Arial", face="bold"))
A<-ggplot(x,aes(score_cutoff,n_positive/256627*100,color=count_cutoff,group=count_cutoff))+
  geom_line(linewidth=.6)+geom_point(size=2.5)+
  scale_color_manual(values=c("MKI67 >= 1"="#3C8DAD","MKI67 >= 2"="#D96B4A"))+
  labs(title="A  Positive-cell prevalence",x="S/G2M score cutoff",y="Union-positive cells (%)",color=NULL)+
  theme_ms+theme(legend.position="top")
B<-ggplot(x,aes(score_cutoff,donor_ecotype_changes,color=count_cutoff,group=count_cutoff))+
  geom_hline(yintercept=0,color="grey80")+geom_line(linewidth=.6)+geom_point(size=2.5)+
  scale_color_manual(values=c("MKI67 >= 1"="#3C8DAD","MKI67 >= 2"="#D96B4A"))+
  scale_y_continuous(breaks=0:12)+
  labs(title="B  Donor ecotype changes",x="S/G2M score cutoff",y="Donors reclassified (of 38)",color=NULL)+
  theme_ms+theme(legend.position="none")
C<-ggplot(corr,aes(score_cutoff,spearman,color=lineage,group=lineage))+
  geom_hline(yintercept=1,color="grey85")+geom_line(linewidth=.5)+geom_point(size=1.8)+
  facet_wrap(~count_cutoff,nrow=1)+
  scale_x_continuous(breaks=c(.05,.10,.20))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+
  scale_color_brewer(palette="Dark2")+
  labs(title="C  Donor-level lineage fractions",x="S/G2M score cutoff",y="Spearman rho vs primary definition",color=NULL)+
  theme_ms+theme(legend.position="bottom")
fig<-A+B+C+plot_layout(nrow=1,widths=c(1,1,2.1))
# No under-plot caption; captions remain in Word only.
cairo_pdf(file.path(root,"FigureS6_Ki67_threshold_sensitivity.pdf"), width=11.5, height=4.5, family="Arial", bg="white")
print(fig); dev.off()
cairo_pdf(file.path(out_deliv,"FigureS6_Ki67_threshold_sensitivity.pdf"), width=11.5, height=4.5, family="Arial", bg="white")
print(fig); dev.off()
png(file.path(out_deliv,"FigureS6_Ki67_threshold_sensitivity_600dpi.png"), width=11.5, height=4.5, units="in", res=600, type="cairo", bg="white", family="Arial")
print(fig); dev.off()
# Write judgment numbers for Response letter
sink(file.path(out_deliv,"FigureS6_sensitivity_judgment.txt"))
cat("Primary definition: Union MKI67 count>=1 OR S/G2M score>0.10 (n_pos=23029, 8.97%).\n")
cat("Ecotype reclassification vs primary (of 38 donors):\n")
print(x[,c("definition","score_cutoff","count_cutoff","n_positive","donor_ecotype_changes")])
cat("\nSpearman rho of donor lineage fractions vs primary (from cellcycle_score_lineage_correlations.csv):\n")
print(corr[order(corr$count_cutoff,corr$score_cutoff,corr$lineage),c("definition","count_cutoff","score_cutoff","lineage","spearman")])
cat("\nJUDGMENT:\n")
cat("- count>=2 at score 0.10: 1/38 ecotype change; lineage Spearman mostly >=0.87 (Myeloid 0.873; others ~0.97).\n")
cat("- score 0.05: 4/38 ecotype changes; several lineages Spearman 0.42-0.78 (NOT all >0.91).\n")
cat("- score 0.20: 8/38 (count>=1) or 11/38 (count>=2) ecotype changes; Lymphoid/Glial Spearman drop to ~0.28-0.42.\n")
cat("- Claim 'Spearman rho > 0.91 / concordance > 92% / conclusions fully consistent across cutoffs' is NOT uniformly justified.\n")
cat("- Safer wording: near-primary alternatives (esp. count>=2 at score 0.10) preserve lineage fractions and ecotypes; more extreme score cutoffs (0.05/0.20) change prevalence and reclassify more donors, with weaker lineage-fraction concordance for lymphoid/glial.\n")
cat("- Myeloid contribution / IDH lineage-shift narrative remains directionally supported near the primary gate; temper absolute concordance claims for distant cutoffs.\n")
sink()
cat("FIG_S6_DONE\n")
