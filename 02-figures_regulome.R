
#####figures####
out<-'outputs/figures_regulomes'
dir.create(out)
source('../utils/r_utils.R')
library(ggnet)
library(network)
library(sna)
library(pheatmap)
library(DESeq2)


#functions####
CompDEGsPathways<-function(res_gsea,
                           res_de,
                           top.n=NULL,
                           FC_col='log2FoldChange',
                           pval_col='padj',
                           col_range=c(-2.5,2.5),
                           transpose=FALSE,
                           show_rownames=FALSE,
                           show_pval=TRUE){
  
  
  #get leading edges
  degs_list<-LeadingEdges(res_gsea)
  
  #filter degs
  if(!is.null(top.n))
    degs_list<-lapply(degs_list, function(x)head(x,top.n))
  
  #trans in dataframe
  degs_pathways<-Reduce(rbind,lapply(names(degs_list),function(p)data.table(pathway=p,
                                                                            gene=degs_list[[p]])))
  
  
  #merge pathway by degs
  res_de_p<-merge(res_de,degs_pathways,by='gene')
  res_de_p<-unique(res_de_p,by=c('gene','pathway'))
  #create heatmaps
  dep_mat<-data.frame(dcast(res_de_p,gene~pathway,value.var =FC_col),row.names = 'gene')
  dep_mat[is.na(dep_mat)]<-0
  
  #add pvalue
  if(show_pval){
    
    res_de_p[,padjsig:=lapply(.SD,function(x)ifelse(x<0.001,'***',ifelse(x<0.01,'**',ifelse(x<0.05,'*','')))),.SDcols=pval_col]
    dep_matp<-data.frame(dcast(res_de_p,gene~pathway,value.var ='padjsig'),row.names = 'gene')
    dep_matp[is.na(dep_matp)]<-''
  }
  
  
  #plot heatmap
  if (transpose) {
    if(show_pval) dep_matp<-t(dep_matp)
    dep_mat<-t(dep_mat)
    
  }
  col_breaks<-c(((col_range[1]*10):(col_range[2]*10))/10)
  if(show_pval){
    print(pheatmap(dep_mat,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize= 7,
                   main='Top DEGs',
                   show_rownames = show_rownames,
                   display_numbers = dep_matp,
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }else{
    print(pheatmap(dep_mat,
                   breaks =col_breaks,
                   color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                         "RdBu")))(length(col_breaks)-1),
                   fontsize= 7,
                   main='Top DEGs',
                   show_rownames = show_rownames,
                   # display_numbers = dep_matp,
                   # cellwidth =20,
                   # cellheight =  8,
                   
                   fontsize_number = 8))
  }
  
  
  
}




#analysis####

#Regulons identified
regulons<-rbind(fread('outputs/01-tf-downstream_gene-inference/TF_regulons_up.csv')[,sens:='up'],
                fread('outputs/01-tf-downstream_gene-inference/TF_regulons_down.csv')[,sens:='down'])
regulons[,TF_regulation:=ifelse(sens=='up','positive regulation','negative regulation')]
linkList<-fread('outputs/01-tf-downstream_gene-inference/tf_gene_weight_and_spearman.csv.gz')
linkList[,sens:=ifelse(r>0,'up','down')]

regulons<-merge(regulons,linkList[,TF:=str_remove(regulatoryGene,'_prot')][,gene:=str_remove(targetGene,'_trs')][,.(TF,gene,weight,r)])
regulons[,avg_regul_score:=mean(sign(r)*weight),by=.(TF,sens)]
regulons[,TF:=factor(TF,levels =unique(regulons[order(abs(avg_regul_score))]$TF ))]
table(regulons$TF)
ggplot(regulons)+geom_bar(aes(x=TF,fill=avg_regul_score))+
  facet_grid(~TF_regulation,scales = 'free_x',space = 'free')+
  #scale_x_discrete(limits=))+
  scale_fill_gradient2(low = 'darkblue',high = 'darkred',mid = 'white',midpoint = 0)+
  theme_bw()+labs(x='Transcription Factor',y='Number of regulated genes')

ggsave(fp(out,'fig1-barplot_regulons.pdf'))


#HFD enrichemnent in regulons
#NES pval
res_tfgsea<-Reduce(rbind,lapply(c('up','down'),
                                function(s)fread(ps('outputs/01-tf-downstream_gene-inference/res_fgsea_degs_in_regulons_',s,'.csv'))[,regul_type:=s]
))
res_tfgsea[,TF:=pathway]

res_tfgsea[,regulon:=paste(pathway,regul_type,sep='_')]
res_tfgsea[,sex:=ifelse(sex=='femelle','female',sex)]
res_tfgsea[,correlation:=ifelse(regul_type=='up','positive','negative')]

res_tfgsea[,TF_regulation:=ifelse(regul_type=='up','positive regulation','negative regulation')]

#7M####
#enriched by sex
res_tfgsea[,sens_regulation:=ifelse(regul_type=='up',+1,-1)]

res_tfgsea[,regulon:=factor(regulon,levels = regulon[order(NES)]),by=.(sex,age)]

ggplot(res_tfgsea[padj<0.001&abs(NES)>2][age=='7M'])+
  geom_col(aes(x=regulon,y=NES,fill=-log10(pval)))+
  facet_wrap('sex',scales = 'free_x')+
  scale_fill_gradient2(high = 'darkred',mid='white',low = 'darkblue')+
  scale_x_discrete(guide = guide_axis(angle=60))+
  theme_minimal()
ggsave(fp(out,'figA-barplot_regulons_enriched_in_7M_HFD_DEGs_padj0.001_NES2.pdf'))


#GO enrichement of these regulons
reg_7Mm<-res_tfgsea[padj<0.001&abs(NES)>2][age=='7M'&sex=='male']$regulon
reg_7Mf<-res_tfgsea[padj<0.001&abs(NES)>2][age=='7M'&sex=='female']$regulon

res_enrdt<-fread('outputs/01-tf-downstream_gene-inference/res_annot_regulons.csv.gz')
#collapse go family: if 2 of the same family take the most general
# res_enrdt[p_value<0.001&term_size<2000,redundant:=any(sapply(strsplit(parents,'\\|')[[1]],function(x)x%in%term_id)),by='query']
# res_enrdt[redundant==T]
top5gosource<-res_enrdt[term_size<2000&str_detect(source,'GO|KEGG')&source!='GO:CC'&p_value<0.001][,.SD[head(order(p_value),5)],by=.(query,source)]

termsmale<-unique(top5gosource[query%in%reg_7Mm]$term_name)
termsfemale<-unique(top5gosource[query%in%reg_7Mf]$term_name)
res_enrdt
res_enrdt[,pct.query.overlap:=intersection_size/query_size]
res_enrdt[,pct.term.background:=term_size/effective_domain_size] 

res_enrdt[,fold.enrichment:=pct.query.overlap/pct.term.background]

res_enrdt[,mlog10p:=-log10(p_value)]

top5go_mat<-dcast(res_enrdt[term_id%in%top5gosource$term_id],
                  term_name~query,value.var = 'mlog10p')
top5go_mat<-data.frame(top5go_mat,row.names = 'term_name')

top5go_mat[is.na(top5go_mat)]<-0


res_enrdt[,pvalsig:=ifelse(p_value<0.0001,'***',ifelse(p_value<0.001,'**',ifelse(p_value<0.01,'*','')))]

top5go_matp<-data.frame(dcast(res_enrdt[term_id%in%top5gosource$term_id],term_name~query,value.var ='pvalsig'),row.names = 'term_name')
top5go_matp[is.na(top5go_matp)]<-''

summary(top5go_mat)
col_breaks<-c(0,1,2,3,4,5)


#col_breaks<-(0:200)/20
col_breaks<-c(0,0.2,0.4,0.8,1,1.5,2,3,4,6,8,10,12,15,18,20,25,30)

# col_breaks<-c(0,5:30)
# col_breaks<-sort(union(-col_breaks,col_breaks))
pdf(fp(out,'figB-heatmap_top5_GOKEGG_enrichments_regulons_sig_7M_male.pdf'),width = 5,height = 8.5)

print(pheatmap(top5go_mat[termsmale,reg_7Mm],
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(length(col_breaks)-1),
               fontsize=6,
               display_numbers = top5go_matp[termsmale,reg_7Mm],
               
               breaks =col_breaks ))
dev.off()


pdf(fp(out,'figB-heatmap_top5_GOKEGG_enrichments_regulons_sig_7M_female.pdf'),width = 5,height = 8.5)

print(pheatmap(top5go_mat[termsfemale,reg_7Mf],
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(length(col_breaks)-1),
               fontsize=6,
               display_numbers = top5go_matp[termsfemale,reg_7Mf],
               
               breaks =col_breaks ))
dev.off()


#Gene expression heatmap
#heatmap lead edges top50
#install.packages('pheatmap')
library(pheatmap)
source('../utils/visualisation.R')
res_tfgseasig<-res_tfgsea[padj<0.001&abs(NES)>2]
regulon_leading7Mf<-LeadingEdges(res_tfgseasig[sex=='female'&age=='7M'])
regulon_leading7Mf50<-lapply(regulon_leading7M, function(x)head(x,50))

col_breaks<-(-40:40)/20
lapply(names(regulon_leading7Mf50), function(reg){
  print(reg)
  reg50<-regulon_leading7Mf50[[reg]]
  expr_mat<-norm_countf[rna_matf[gene_name%in%reg50]$gene_id,]
  rownames(expr_mat)<-rna_matf[gene_name%in%reg50]$gene_name
  
  pdf(fp(out,ps('figC-heatmap_',reg,'regulon_top50DE_female7M.pdf')),width = 5,height = 6)
  print(pheatmap(expr_mat[,mtd[sex=='F'&age=='7M']$sample],scale = 'row',
                 annotation_col = data.frame(mtd,row.names = 'sample')[colnames(expr_matle50_7M),c('sex','age','condition')],
                 color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "RdYlBu")))(length(col_breaks)-1),fontsize=6,
                 breaks =col_breaks ))
  dev.off()
})

regulon_leading7Mm<-LeadingEdges(res_tfgseasig[sex=='male'&age=='7M'])
regulon_leading7Mm50<-lapply(regulon_leading7Mm, function(x)head(x,50))

lapply(names(regulon_leading7Mm50), function(reg){
  print(reg)
  reg50<-regulon_leading7Mm50[[reg]]
  expr_mat<-norm_countf[rna_matf[gene_name%in%reg50]$gene_id,]
  rownames(expr_mat)<-rna_matf[gene_name%in%reg50]$gene_name
  
  pdf(fp(out,ps('figC-heatmap_',reg,'regulon_top50DE_male7M.pdf')),width = 5,height = 6)
  print(pheatmap(expr_mat[,mtd[sex=='M'&age=='7M']$sample],scale = 'row',
                 annotation_col = data.frame(mtd,row.names = 'sample')[colnames(expr_matle50_7M),c('sex','age','condition')],
                 color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "RdYlBu")))(length(col_breaks)-1),fontsize=6,
                 breaks =col_breaks ))
  dev.off()
})



#DEGs heatmap
out<-'outputs/figures_regulomes/V2/'
dir.create(out)
#top20 by regulon/Sex [todo]
source('../utils/visualisation.R')
res_de<-Reduce(rbind,lapply(c('4M','7M'),function(m)Reduce(rbind,lapply(c('male','femelle'),function(s)fread(fp('ref-data',paste0('res_de_',m,'_',s,'.csv')),
                                                                                                             select = 1:5,
                                                                                                             col.names = c('gene_id','gene_name','log2FC','pvalue','padj'))[
                                                                                                               ,age:=m][
                                                                                                                 ,sex:=s][
                                                                                                                   ,compa:='HFDvsChow']))))

regulons<-rbind(fread('outputs/01-tf-downstream_gene-inference/TF_regulons_up.csv')[,sens:='up'],
                fread('outputs/01-tf-downstream_gene-inference/TF_regulons_down.csv')[,sens:='down'])

# res_de[gene_name%in%regulons$gene,padj2:=p.adjust(pvalue,method = 'BH'),by=.(sex,age)]
# table(res_de[padj2<0.05][,.(sex,age)])
# #          age
# # sex        4M  7M
# #   femelle   0   6
# #   male    109   4


CompDEGsPathways(res_tfgsea[sex=='male'&age=='7M'&regulon%in%reg_7Mm],
                 res_de[,gene:=gene_name][sex=='male'&age=='7M'],
                 show_pval = F,
                 show_rownames = T,
                 top.n=50,
                 FC_col = 'log2FC',col_range = c(-1,1),
                 pval_col = 'pvalue',transpose = T)


pdf(fp(out,'figD-heatmap_top20_DEGs_regulons_sig_7M_female.pdf'),width = 5,height = 8.5)
CompDEGsPathways(res_tfgsea[sex=='female'&age=='7M'&regulon%in%reg_7Mf],
                 res_de[,gene:=gene_name][sex=='femelle'&age=='7M'],
                 top.n=20,FC_col = 'log2FC',col_range = c(-1,1),pval_col = 'pvalue',transpose = F)
dev.off()


pdf(fp(out,'fig3V2-heatmap_top20_DEGs_regulons_sig_7M_male.pdf'),width = 5,height = 8.5)
CompDEGsPathways(res_tfgsea[sex=='male'&age=='7M'&regulon%in%reg_7Mm],
                 res_de[,gene:=gene_name][sex=='male'&age=='7M'],
                 top.n=20,FC_col = 'log2FC',col_range = c(-1,1),pval_col = 'pvalue',transpose = F)
dev.off()


#4M####
#enriched by sex
ggplot(res_tfgsea[padj<0.001&abs(NES)>2][age=='4M'])+
  geom_col(aes(x=regulon,y=NES,fill=-log10(pval)))+
  facet_wrap('sex',scales = 'free_x')+
  scale_fill_gradient2(high = 'darkred',mid='white',low = 'darkblue')+
  scale_x_discrete(guide = guide_axis(angle=60))+
  theme_minimal()
ggsave(fp(out,'fig3V2-barplot_regulons_enriched_in_4M_HFD_DEGs_padj0.001_NES2.pdf'))

ggplot(res_tfgsea[age=='4M'][padj<0.01][regulon%in%res_tfgsea[padj<0.001&abs(NES)>2][age=='7M']$regulon])+
  geom_col(aes(x=regulon,y=NES,fill=-log10(pval)))+
  facet_wrap('sex',scales = 'free_x')+
  scale_fill_gradient2(high = 'darkred',mid='white',low = 'darkblue',)+
  scale_x_discrete(guide = guide_axis(angle=60))+
  theme_minimal()
ggsave(fp(out,'fig3V2-barplot_7Mregulons_enrichment_in_4M_HFD_padj0.01.pdf'))


#GO enrichement of these regulons
reg_4Mm<-res_tfgsea[padj<0.001&abs(NES)>2][age=='4M'&sex=='male']$regulon
reg_4Mf<-res_tfgsea[padj<0.001&abs(NES)>2][age=='4M'&sex=='female']$regulon


termsmale<-unique(top5gosource[query%in%reg_4Mm]$term_name)
termsfemale<-unique(top5gosource[query%in%reg_4Mf]$term_name)


#col_breaks<-(0:200)/20
col_breaks<-c(0,0.2,0.4,0.8,1,1.5,2,3,4,6,8,10,12,15,18,20,25,30)

# col_breaks<-c(0,5:30)
# col_breaks<-sort(union(-col_breaks,col_breaks))
pdf(fp(out,'fig3V2-heatmap_top5_GOKEGG_enrichments_regulons_sig_4M_male.pdf'),width = 5,height = 8.5)

print(pheatmap(top5go_mat[termsmale,reg_4Mm],
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(length(col_breaks)-1),
               fontsize=6,
               display_numbers = top5go_matp[termsmale,reg_4Mm],
               
               breaks =col_breaks ))
dev.off()


pdf(fp(out,'fig3V2-heatmap_top5_GOKEGG_enrichments_regulons_sig_4M_female.pdf'),width = 5,height = 8.5)

print(pheatmap(top5go_mat[termsfemale,reg_4Mf],
               color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                     "RdYlBu")))(length(col_breaks)-1),
               fontsize=6,
               display_numbers = top5go_matp[termsfemale,reg_4Mf],
               
               breaks =col_breaks ))
dev.off()

#DEGs heatmap
pdf(fp(out,'fig3V2-heatmap_top20_DEGs_regulons_sig_4M_male.pdf'),width = 5,height = 8.5)
CompDEGsPathways(res_tfgsea[sex=='male'&age=='4M'&regulon%in%reg_4Mm],
                 res_de[,gene:=gene_name][sex=='male'&age=='4M'],
                 top.n=20,FC_col = 'log2FC',col_range = c(-1,1),pval_col = 'pvalue',transpose = F)
dev.off()


pdf(fp(out,'fig3V2-heatmap_top20_DEGs_regulons_sig_4M_female.pdf'),width = 5,height = 8.5)
CompDEGsPathways(res_tfgsea[sex=='female'&age=='4M'&regulon%in%reg_4Mf],
                 res_de[,gene:=gene_name][sex=='femelle'&age=='4M'],
                 top.n=20,FC_col = 'log2FC',col_range = c(-1,1),pval_col = 'pvalue',transpose = F)
dev.off()

#SUPP validation regulon####
#Mecp2 ChipSeq
expressed_genes_mecp2pbind<-fread('outputs/01-tf-downstream_gene-inference/validation_Mecp2/expressed_genes_ChIPscore_Mecp2.csv.gz')
expressed_genes_mecp2pbind[,`Mecp2 Binding score`:=score]
ggplot(unique(expressed_genes_mecp2pbind[order(gene_name,-score)],by='gene_name'))+
  geom_boxplot(aes(x=`in Mecp2 regulon`,y=`Mecp2 Binding score`,fill=`in Mecp2 regulon`))+
  coord_cartesian(ylim = c(0,1000))+
  theme_minimal()
ggsave(fp(out,'figsupp-Mecp2_Regulon_ChIPseqGSE71126_RAW_validation.pdf'))
