#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork);library(cowplot);library(ggrepel)})
root<-"/publicfs10/fs10-m9/home/m9s004832/TTLL7"
indir<-file.path(root,"revision/20261023_reproduction/exact_cellset/figure2_4_audit")
out<-file.path(root,"revision/20261023_reproduction/figure_review/Fig2_panels")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_set(theme_classic(base_size=8,base_family="Arial"))
cols_lineage<-c(Myeloid="#0072B5",Tumor="#D94F4F",Lymphoid="#00A087",Glial="#AD8BC9",Stromal="#EFC000")
cols_idh<-c(IDHwt="#D94F4F",IDHmt="#233B8B")
grade_colors <- c("Grade 2"="lightblue", "Grade 3"="orange", "Grade 4"="firebrick3")
cols_diagnosis <- c(Astrocytoma="#BC3C29",Glioblastoma="#0072B5",Oligodendroglioma="#E18727")
wrap_panel<-function(p,label){ggdraw()+draw_plot(p,x=.05,y=.03,width=.94,height=.94)+draw_label(label,x=.015,y=.985,hjust=0,vjust=1,fontface="bold",size=10,fontfamily="Arial")}
add_bracket<-function(p,x1,x2,y,label){p+annotate("segment",x=x1,xend=x2,y=y,yend=y,linewidth=.3)+annotate("segment",x=x1,xend=x1,y=y,yend=y/1.15,linewidth=.3)+annotate("segment",x=x2,xend=x2,y=y,yend=y/1.15,linewidth=.3)+annotate("text",x=(x1+x2)/2,y=y*1.06,label=label,size=2.5)}
pstar<-function(p) ifelse(p<.001,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","ns")))
le<-read.csv(file.path(indir,"figure2A_plot_data.csv"))
le$lineage<-factor(le$lineage,levels=c("Myeloid","Tumor","Lymphoid","Glial","Stromal"))
pA<-ggplot(le,aes(lineage,pmax(mean_MKI67,1e-4),fill=lineage))+geom_boxplot(width=.62,outlier.shape=NA,linewidth=.3)+geom_jitter(width=.14,size=.75,alpha=.65,shape=21,stroke=.15)+scale_fill_manual(values=cols_lineage)+scale_y_log10(limits=c(1e-4,2.2),breaks=c(1e-3,1e-2,1e-1,1),labels=c("1e-3","1e-2","1e-1","1"))+labs(title=expression(italic(MKI67)~"expression across cellular lineages"),x=NULL,y=expression("Mean "~italic(MKI67)))+theme(legend.position="none",plot.title=element_text(hjust=.5,size=8),axis.text.x=element_text(angle=45,hjust=1,size=8),plot.margin=margin(7,5,5,10))
pa_pairs<-combn(levels(le$lineage),2,simplify=FALSE)
pa_stats<-lapply(pa_pairs,function(z){w<-reshape(le[le$lineage%in%z,],idvar="Patient",timevar="lineage",direction="wide");a<-w[[paste0("mean_MKI67.",z[1])]];b<-w[[paste0("mean_MKI67.",z[2])]];ok<-is.finite(a)&is.finite(b);data.frame(pair=paste(z,collapse="|"),p=wilcox.test(a[ok],b[ok],paired=TRUE,exact=FALSE)$p.value)})
pa_stats<-do.call(rbind,pa_stats);pa_stats$q<-p.adjust(pa_stats$p,"BH")
getpa<-function(a,b)pstar(pa_stats$q[match(paste(a,b,sep="|"),pa_stats$pair)])
pA<-add_bracket(pA,1,2,.16,getpa("Myeloid","Tumor"));pA<-add_bracket(pA,2,3,.32,getpa("Tumor","Lymphoid"));pA<-add_bracket(pA,2,4,.64,getpa("Tumor","Glial"));pA<-add_bracket(pA,2,5,1.28,getpa("Tumor","Stromal"))
ge<-read.csv(file.path(indir,"figure2_patient_global_mean_MKI67_final.csv"))
ge$IDH<-factor(ge$IDH,levels=c("IDHwt","IDHmt"))
ge$Diagnosis<-factor(ge$Diagnosis,levels=c("Astrocytoma","Glioblastoma","Oligodendroglioma"))
kw<-kruskal.test(mean_MKI67~Diagnosis,ge)$p.value
# Fig2B: one plot, shared y + single baseline; spacer between IDH and Diagnosis
ge_B <- bind_rows(
  ge[!is.na(ge$IDH),] |> transmute(Patient, mean_MKI67, x=as.character(IDH), fill=as.character(IDH), block="idh"),
  ge[!is.na(ge$Diagnosis),] |> transmute(Patient, mean_MKI67, x=as.character(Diagnosis), fill=as.character(Diagnosis), block="dx")
)
x_levels <- c("IDHwt","IDHmt"," ","Astrocytoma","Glioblastoma","Oligodendroglioma")
ge_B$x <- factor(ge_B$x, levels=x_levels)
ge_B <- bind_rows(ge_B, data.frame(Patient="__spacer__", mean_MKI67=NA_real_, x=factor(" ", levels=x_levels), fill="IDHwt", block="spacer", stringsAsFactors=FALSE))
ge_B$x <- factor(ge_B$x, levels=x_levels)
fill_B <- c(cols_idh, cols_diagnosis)
ylim_B <- c(1e-4, .7)
ybreaks_B <- c(1e-4,1e-3,1e-2,1e-1)
ylabels_B <- c("1e-04","1e-03","1e-02","1e-01")
p_idh <- wilcox.test(mean_MKI67~IDH, ge, exact=FALSE)$p.value
pB <- ggplot(ge_B, aes(x, pmax(mean_MKI67, 1e-4), fill=fill)) +
  geom_boxplot(data=subset(ge_B, block!="spacer"), width=.65, outlier.shape=NA, linewidth=.3) +
  geom_jitter(data=subset(ge_B, block!="spacer"), width=.13, size=.7, alpha=.7, shape=21, stroke=.15) +
  scale_fill_manual(values=fill_B, guide="none") +
  scale_x_discrete(drop=FALSE) +
  scale_y_log10(limits=ylim_B, breaks=ybreaks_B, labels=ylabels_B, expand=expansion(mult=c(0.02,0.08))) +
  annotate("segment", x=1, xend=2, y=.28, yend=.28, linewidth=.3) +
  annotate("segment", x=1, xend=1, y=.28, yend=.28/1.15, linewidth=.3) +
  annotate("segment", x=2, xend=2, y=.28, yend=.28/1.15, linewidth=.3) +
  annotate("text", x=1.5, y=.28*1.06, label=pstar(p_idh), size=2.5, family="Arial") +
  annotate("text", x=5, y=.38, label=sprintf("Kruskal-Wallis P = %.3f", kw), size=2.1, family="Arial") +
  labs(title="MKI67 expression across IDH status and tumor types", x=NULL, y=expression("Mean "~italic(MKI67))) +
  theme(legend.position="none",
        plot.title=element_text(hjust=.5, size=8, family="Arial"),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8, family="Arial"),
        axis.text.y=element_text(size=8, family="Arial"),
        axis.title.y=element_text(size=8, family="Arial"),
        plot.margin=margin(6,6,2,6))
gC<-ge[ge$IDH=="IDHmt" & ge$TumorGrade %in% c("G2","G3"),];gC$TumorGrade<-factor(gC$TumorGrade,levels=c("G2","G3"),labels=c("Grade 2","Grade 3"))
pC<-ggplot(gC,aes(TumorGrade,pmax(mean_MKI67,1e-4),fill=TumorGrade))+geom_boxplot(width=.62,outlier.shape=NA,linewidth=.3)+geom_jitter(width=.12,size=.85,alpha=.7,shape=21,stroke=.15)+scale_fill_manual(values=grade_colors)+scale_y_log10(limits=c(1e-4,.4))+labs(title=expression(italic(MKI67)~"expression across grades in IDHmt glioma"),x=NULL,y=expression("Mean "~italic(MKI67)))+theme(legend.position="none",plot.title=element_text(hjust=.5,size=8),plot.margin=margin(10,20,5,14))
pC<-add_bracket(pC,1,2,.18,"ns")
fr<-read.csv(file.path(indir,"figure2E_plot_data_union.csv"))
gm<-read.csv(file.path(indir,"figure2_global_plot_data_union.csv"));gm<-data.frame(Patient=gm$Patient,lineage="Global",IDH=gm$IDH,Ki67_fraction=gm$global)
ep<-bind_rows(fr,gm);ep$lineage<-factor(ep$lineage,levels=c("Global","Tumor","Myeloid","Lymphoid","Glial","Stromal"));ep$IDH<-factor(ep$IDH,levels=c("IDHwt","IDHmt"))
lineages_e<-c("Global","Tumor","Myeloid","Lymphoid","Glial","Stromal")
pvals_e<-sapply(lineages_e,function(x){z<-ep[ep$lineage==x,];if(length(unique(z$IDH))<2)return(1);wilcox.test(Ki67_fraction~IDH,z,exact=FALSE)$p.value})
qvals_e<-p.adjust(pvals_e,"BH")
labsig<-data.frame(lineage=factor(lineages_e,levels=levels(ep$lineage)),lab=pstar(qvals_e),y=c(.39,.39,.20,.25,.22,.27))
pE<-ggplot(ep,aes(lineage,Ki67_fraction*100,fill=IDH))+geom_boxplot(position=position_dodge(.72),width=.62,outlier.shape=NA,linewidth=.25)+geom_point(position=position_jitterdodge(jitter.width=.1,dodge.width=.72),shape=21,size=.55,alpha=.65,stroke=.15)+geom_text(data=labsig,aes(lineage,y*100,label=lab),inherit.aes=FALSE,size=2.4)+scale_fill_manual(values=cols_idh)+coord_cartesian(ylim=c(0,42))+labs(title="Quantification of Ki-67+ fraction",x=NULL,y="Ki-67+ fraction (%)",fill=NULL)+theme(plot.title=element_text(hjust=.5,size=8),legend.position=c(.86,.83),legend.key.size=unit(2.5,"mm"),legend.text=element_text(size=8),axis.text.x=element_text(angle=35,hjust=1,size=8),plot.margin=margin(8,4,5,10))
tf<-read.csv(file.path(indir,"figure2F_plot_data_union.csv"));tf$IDH<-factor(tf$IDH,levels=c("IDHwt","IDHmt"))
oneF<-function(comp,lab,sig,ylim){x<-tf[tf$compartment==comp,];ggplot(x,aes(IDH,Ki67_fraction*100,fill=IDH))+geom_boxplot(width=.62,outlier.shape=NA,linewidth=.3)+geom_jitter(width=.12,size=.75,alpha=.7,shape=21,stroke=.15)+scale_fill_manual(values=cols_idh)+coord_cartesian(ylim=c(0,ylim))+labs(x=NULL,y=lab)+theme(legend.position="none",axis.text.x=element_text(size=8),axis.title.y=element_text(size=8),plot.margin=margin(12,4,4,9))+annotate("segment",x=1,xend=2,y=ylim*.88,yend=ylim*.88,linewidth=.3)+annotate("text",x=1.5,y=ylim*.95,label=sig,size=2.5)}
fstats<-sapply(c("Tumor","Non-tumor"),function(comp){z<-tf[tf$compartment==comp,];p.adjust(wilcox.test(Ki67_fraction~IDH,z,exact=FALSE)$p.value,"BH")})
pF<-oneF("Tumor","Ki-67+ tumor cells (%)",pstar(fstats[1]),35)+oneF("Non-tumor","Ki-67+ non-tumor cells (%)",pstar(fstats[2]),22)+plot_annotation(title="Ki-67+ fraction of tumor and non-tumor cells",theme=theme(plot.title=element_text(hjust=.5,size=8,margin=margin(b=-7))))
savep<-function(p,n,w,h,label){cairo_pdf(file.path(out,n),width=w/72,height=h/72,family="Arial",bg="white");print(wrap_panel(p,label));dev.off()}
# D: preserve the original TCGA MKI67 expression boxplot (do not substitute a UMAP).
tcga_env <- new.env()
load(file.path(root,"data/GBM-LGG-TCGA-exprSet_names_by_symbol.Rdata"),envir=tcga_env)
load(file.path(root,"data/meta-TCGA-GBM-LGG.rdata"),envir=tcga_env)
expr_tcga <- tcga_env$exprSet
to_patient <- function(x) vapply(strsplit(as.character(x),"-",fixed=TRUE),
                                  function(v) paste(v[1:min(3,length(v))],collapse="-"),
                                  FUN.VALUE=character(1))
allowed <- unique(unlist(lapply(tcga_env$meta,function(m) to_patient(rownames(m)))))
tcga_expr <- data.frame(patient_id=to_patient(colnames(expr_tcga)),
                        MKI67=as.numeric(expr_tcga["MKI67",]),
                        stringsAsFactors=FALSE) |>
  filter(patient_id %in% allowed, is.finite(MKI67)) |>
  group_by(patient_id) |>
  summarise(MKI67=mean(MKI67,na.rm=TRUE),.groups="drop")
cbio <- read.delim(file.path(root,"data/cBioPortal/lgggbm_tcga_pub/clinical_data.tsv"),
                   check.names=FALSE,stringsAsFactors=FALSE)
idh <- cbio |>
  filter(clinicalAttributeId=="IDH_STATUS") |>
  mutate(value=na_if(value,"NA")) |>
  group_by(patientId) |>
  summarise(IDH_STATUS=dplyr::first(na.omit(value)),.groups="drop") |>
  mutate(IDH=ifelse(grepl("mut",tolower(IDH_STATUS)) &
                      !grepl("wild|wt",tolower(IDH_STATUS)),"IDHmt","IDHwt"))
tcga_plot <- inner_join(tcga_expr,idh,by=c("patient_id"="patientId")) |>
  filter(is.finite(MKI67),MKI67>0)
tcga_plot$IDH <- factor(tcga_plot$IDH,levels=c("IDHwt","IDHmt"))
tcga_p <- wilcox.test(MKI67~IDH,data=tcga_plot,exact=FALSE)$p.value
tcga_ymax <- max(tcga_plot$MKI67,na.rm=TRUE)
pD <- ggplot(tcga_plot,aes(IDH,MKI67,fill=IDH)) +
  geom_boxplot(width=.65,outlier.shape=NA,linewidth=.3,alpha=.85) +
  geom_jitter(width=.16,size=.45,alpha=.4,shape=16,color="grey25") +
  scale_fill_manual(values=cols_idh) +
  scale_y_log10(limits=c(1,tcga_ymax*1.6),
                breaks=c(1,10,100,1000,10000,100000),
                labels=scales::label_number()) +
  labs(title=expression(italic(MKI67)~"expression in TCGA cohort"),
       x=NULL,y="MKI67 Expression") +
  theme_classic(base_size=8,base_family="Arial") +
  theme(plot.title=element_text(hjust=.5,size=8),
        axis.text=element_text(size=8),axis.title.y=element_text(size=8),
        legend.position="none",plot.margin=margin(8,8,5,10))
pD <- add_bracket(pD,1,2,tcga_ymax*1.2,pstar(tcga_p))
savep(pA,"A.pdf",298,237,"A");savep(pB,"B.pdf",297,237,"B");savep(pC,"C.pdf",298,186,"C");savep(pE,"E.pdf",298,184,"E");savep(pF,"F.pdf",297,184,"F");savep(pD,"D.pdf",297,186,"D")
writeLines(c(sprintf("Fig2B diagnosis Kruskal-Wallis P=%.8g",kw),"A: paired Wilcoxon tests with BH adjustment","B-C/E-F: donor-level Wilcoxon tests; diagnosis uses Kruskal-Wallis with BH-adjusted pairwise tests"),file.path(out,"README.txt"))
# G: restore the original double-donut encoding and fixed lineage/cell-type palettes.
comp <- read.csv(file.path(root,"revision/20261023_reproduction/exact_cellset/definition_sensitivity_20261029/figure2G_celltype_composition_union.csv"))
comp$lineage <- factor(comp$lineage,levels=c("Tumor","Myeloid","Lymphoid","Glial","Stromal"))
celltype_order <- c("AC","MES","NPC","OPC","monocyte-derived DCs","myeloid-derived suppressor cells","macrophages","microglia","T cells","B cells","NK cells","fibroblasts","endothelial cells","pericytes","oligodendrocytes")
cellcols <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B","#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4","#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8")
names(cellcols) <- celltype_order
all_donut_cols <- c(cols_lineage,cellcols)
plot_double_donut <- function(dat,idh_label) {
  outer_df <- dat |>
    group_by(lineage,celltype) |>
    summarise(n=sum(positive_cells,na.rm=TRUE),.groups="drop") |>
    arrange(lineage,desc(n)) |>
    mutate(prop=n/sum(n),ymax=cumsum(prop),ymin=lag(ymax,default=0),
           theta_mid=2*pi*(ymin+ymax)/2)
  inner_df <- outer_df |>
    group_by(lineage) |>
    summarise(n=sum(n),.groups="drop") |>
    arrange(lineage) |>
    mutate(prop=n/sum(n),ymax=cumsum(prop),ymin=lag(ymax,default=0),
           theta_mid=2*pi*(ymin+ymax)/2)
  make_ring <- function(d,r0,r1,fill_col) {
    pieces <- lapply(seq_len(nrow(d)),function(i) {
      a <- seq(2*pi*d$ymin[i],2*pi*d$ymax[i],length.out=max(12,ceiling(150*d$prop[i])))
      data.frame(x=c(r0*sin(a),rev(r1*sin(a))),
                 y=c(r0*cos(a),rev(r1*cos(a))),
                 group=paste(fill_col,i,sep="_"),fill=d[[fill_col]][i])
    })
    bind_rows(pieces)
  }
  inner_poly <- make_ring(inner_df,.78,1.42,"lineage")
  outer_poly <- make_ring(outer_df,1.42,2.18,"celltype")
  outer_df$anchor_x <- 2.22*sin(outer_df$theta_mid)
  outer_df$anchor_y <- 2.22*cos(outer_df$theta_mid)
  outer_df$side <- ifelse(outer_df$anchor_x >= 0, 1, -1)
  for (ss in c(-1,1)) {
    ii <- which(outer_df$side == ss)
    ii <- ii[order(outer_df$anchor_y[ii])]
    outer_df$label_y[ii] <- if(length(ii)==1) 0 else seq(-2.25,2.25,length.out=length(ii))
  }
  outer_df$label_x <- outer_df$side*2.32
  outer_df$line_x <- outer_df$side*2.18
  outer_df$label_hjust <- ifelse(outer_df$side > 0, 0, 1)
  outer_names <- gsub("myeloid-derived suppressor cells","myeloid-derived\nsuppressor cells",outer_df$celltype,fixed=TRUE)
  outer_df$outer_label <- paste0(outer_names,"\n(",scales::percent(outer_df$prop,accuracy=.1),")")
  inner_df$inner_label <- as.character(inner_df$lineage)
  all_fill <- c(cols_lineage,cellcols)
  fillv <- all_fill[unique(c(as.character(inner_df$lineage),as.character(outer_df$celltype)))]
  fillv <- fillv[!is.na(fillv)]
  ggplot() +
    geom_polygon(data=inner_poly,aes(x=x,y=y,group=group,fill=fill),
                 colour="white",linewidth=.3) +
    geom_polygon(data=outer_poly,aes(x=x,y=y,group=group,fill=fill),
                 colour="white",linewidth=.25) +
    geom_text(data=inner_df[inner_df$inner_label!="",],
              aes(x=1.08*sin(theta_mid),y=1.08*cos(theta_mid),label=inner_label),
              size=2.2,fontface="bold",family="Arial") +
    geom_segment(data=outer_df[outer_df$n > 0,],
              aes(x=anchor_x,y=anchor_y,xend=line_x,yend=label_y),
              linewidth=.2,colour="grey45") +
    geom_text(data=outer_df[outer_df$n > 0,],
              aes(x=label_x,y=label_y,label=outer_label,hjust=label_hjust),
              size=2.15,family="Arial",lineheight=.9,show.legend=FALSE) +
    annotate("text",x=0,y=0,label=idh_label,size=2.8,fontface="bold",family="Arial") +
    coord_equal(xlim=c(-3.4,3.4),ylim=c(-2.55,2.55),clip="off",expand=FALSE) +
    scale_fill_manual(values=fillv,guide="none") +
    theme_void(base_size=8,base_family="Arial") +
    theme(plot.margin=margin(3,2,3,2))
}
pG <- wrap_plots(
  plot_double_donut(filter(comp,IDH=="IDHmt"),"IDHmt"),
  plot_double_donut(filter(comp,IDH=="IDHwt"),"IDHwt"),ncol=2
) + plot_annotation(
  title="Cell types proportion in Ki-67+ cells as stratified in IDHwt and IDHmt gliomas",
  theme=theme(plot.title=element_text(hjust=.5,size=8,family="Arial"))
)
# Compose all panels from their R plot objects on the original A4 portrait grid.
page_w <- 595.276; page_h <- 841.89
specs <- list(A=c(0,636.89,298,205),B=c(298,636.89,297.276,205),
              C=c(0,450.89,298,186),D=c(298,450.89,297.276,186),
              E=c(0,266.89,298,184),F=c(298,266.89,297.276,184),
              G=c(0,0,595.276,266.89))
plots <- list(A=pA,B=pB,C=pC,D=pD,E=pE,F=pF,G=pG)
fig <- cowplot::ggdraw()
for(nm in names(specs)) {
  z <- specs[[nm]]
  fig <- fig + cowplot::draw_plot(wrap_panel(plots[[nm]],nm),x=z[1]/page_w,y=z[2]/page_h,width=z[3]/page_w,height=z[4]/page_h)
}
cairo_pdf(file.path(out,"Fig2_review_R.pdf"),width=page_w/72,height=page_h/72,family="Arial",bg="white")
print(fig)
dev.off()
file.copy(file.path(out,"Fig2_review_R.pdf"),file.path(out,"Fig2_review.ai"),overwrite=TRUE)
png(file.path(out,"Fig2_review_R_QA.png"),width=page_w/72,height=page_h/72,units="in",res=300,type="cairo",bg="white")
print(fig)
dev.off()
cat("FIG2_REVIEW_R_DONE\n")

# Delivery copies (vector PDF + high-DPI PNG)
deliv <- "/publicfs10/fs10-m9/home/m9s004832/TTLL7/revision/delivery_figures_20261101"
dir.create(deliv, recursive=TRUE, showWarnings=FALSE)
file.copy(file.path(out,"Fig2_review_R.pdf"), file.path(deliv,"Figure2.pdf"), overwrite=TRUE)
file.copy(file.path(out,"Fig2_review_R.pdf"), file.path(deliv,"Fig2_review_R.pdf"), overwrite=TRUE)
file.copy(file.path(out,"B.pdf"), file.path(deliv,"Figure2_panel_B.pdf"), overwrite=TRUE)
png(file.path(deliv,"Figure2_600dpi.png"), width=page_w/72, height=page_h/72, units="in", res=600, type="cairo", bg="white", family="Arial")
print(fig); dev.off()
png(file.path(deliv,"Figure2_panel_B_600dpi.png"), width=297/72, height=237/72, units="in", res=600, type="cairo", bg="white", family="Arial")
print(wrap_panel(pB,"B")); dev.off()
cat("FIG2_DELIVERY_COPIED\n")
