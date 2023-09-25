
#Regulome analysis using proteomic and transcriptomics
out<-'outputs/01-regulome_analysis'
dir.create(out,recursive = T)
#install.packages('BiocManager')
#BiocManager::install("GENIE3")
#install.packages('R.utils')

library(GENIE3)
library(fgsea)
library(edgeR)
library(data.table)
library(ggplot2)
fp<-function(...)file.path(...)


#answer the biolofical questions: Trscriptomic change induce by maternal HFD is  associated to a specific TF? 


#O) #data quality check####
# get the full matrix of transcripto /proteo

#rna####
#metadata RNA
mtd<-fread('ref-data/raw_expr/metadata_mice.csv',select = 1:4,col.names = c('sample','sex','age','condition'),skip = 1)
mtd[,file_key:=paste(ifelse(sex=='M','male','female'),age,sep='_')]

mat_male7m<-fread('ref-data/raw_expr/raw_rna_counts_male_7M.csv',select = 1:(nrow(mtd[sex=='M'&age=='7M'])+2),col.names = c('gene_name','gene_id',mtd[sex=='M'&age=='7M']$sample))
mat_male4m<-fread('ref-data/raw_expr/raw_rna_counts_male_4M.csv',select = c('gene_name','gene_id',mtd[sex=='M'&age=='4M']$sample))
mat_female7m<-fread('ref-data/raw_expr/raw_rna_counts_female_7M.csv',select = c('gene_name','gene_id',mtd[sex=='F'&age=='7M']$sample))
mat_female4m<-fread('ref-data/raw_expr/raw_rna_counts_female_4M.csv',select =  c('gene_name','gene_id',mtd[sex=='F'&age=='4M']$sample))


rna_mat<-Reduce(function(x,y)merge(x,y,by=c('gene_name','gene_id')),
                list(mat_male7m,
                     mat_male4m,
                     mat_female7m,
                     mat_female4m))

fwrite(rna_mat,'ref-data/raw_expr/raw_rna_counts_clean.csv.gz')

#filter for expressed genes
samples_rna<-colnames(rna_mat)[str_detect(colnames(rna_mat),'WC|WH')]
length(samples_rna)#34
isexpr <- rowSums(cpm(data.frame(rna_mat[,-'gene_name'],row.names = 'gene_id'))>1) >= 0.1 * length(samples_rna)
sum(isexpr) #14k

rna_matf <- rna_mat[isexpr,]
fwrite(rna_matf,'ref-data/raw_expr/raw_rna_counts_clean_filtered.csv.gz')
rna_matf<-fread('ref-data/raw_expr/raw_rna_counts_clean_filtered.csv.gz')

#Outliers detection

#proteo####
#get mtd proteo
mtd_prot<-Reduce(rbind,
                 lapply(unique(mtd$file_key),function(k)data.table(sample=colnames(fread(paste0('ref-data/raw_expr/log2_prot_expr_',k,'.csv'),)[,.SD,.SDcols=!c('gene_name','gene_id')]),
                                                                   file_key=k,
                                                                   sex=unique(mtd[file_key==k]$sex),
                                                                   age=unique(mtd[file_key==k]$age))))
mtd_prot#24 samples
mtd_prot[,condition:=str_extract(sample,'WC|WH')]
fwrite(mtd_prot,'ref-data/raw_expr/metadata_prot_mice.csv')
mtd_prot<-fread('ref-data/raw_expr/metadata_prot_mice.csv')

#get expr 
proteo_mat<-Reduce(function(x,y)merge(x,y,by=c('gene_name','gene_id')),
                   lapply(unique(mtd$file_key),function(k)fread(paste0('ref-data/raw_expr/log2_prot_expr_',k,'.csv'))))

fwrite(proteo_mat,'ref-data/raw_expr/log2_prot_expr_clean.csv.gz')
proteo_mat<-fread('ref-data/raw_expr/log2_prot_expr_clean.csv.gz')
proteo_mat # 2790 prot


#filter prot for too many missing values
?goodSamplesGenes
qcs<-goodSamplesGenes(t(data.frame(proteo_mat[,-'gene_name'],row.names = 'gene_id')),minNSamples = 10,minFraction = 0.5,minNGenes = 100)
qcs #all OK

#distr NA
protmat<-as.matrix(data.frame(proteo_mat[,-'gene_name'],row.names = 'gene_id'))
plot(hist(rowSums(is.na(protmat)),breaks = 8)) #all are <9/25 samples NA 


#Outliers detection

#using hierarchical clustering
sampleTree = hclust(dist(t(protmat)[str_replace_all(paste0('X',mtd_prot$sample),'-','.'),]),
                    method = "average")
plot(sampleTree)
traits<-c('sex','age','condition')
traitsnum<-paste(traits,'num',sep = '_')
mtd_prot[,(traits):=lapply(.SD,as.factor),.SDcols = traits]
mtd_prot[,(traitsnum):=lapply(.SD,as.numeric),.SDcols = traits]
traitColors = numbers2colors(data.frame(mtd_prot[,.SD,.SDcols = c('sample',traitsnum)],row.names = 'sample'), signed = FALSE);
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = traits,
                    main = "Sample dendrogram and trait heatmap")

#age bias
#log2 values bias?
mtd_prot[,tot.log2expr:=colSums(protmat[,str_replace_all(paste0('X',mtd_prot$sample),'-','.')],na.rm = T)]
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=tot.log2expr))
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=tot.log2expr))

mtd_prot[,avg.log2expr:=colMeans(protmat[,str_replace_all(paste0('X',mtd_prot$sample),'-','.')],na.rm = T)]
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=avg.log2expr))
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=avg.log2expr,col=sex))
#=> batch effect
#need normalized or QCs
#divide by med or mean ?
#dep of distribution
plot(density(na.omit(protmat[,str_replace_all(paste0('X',mtd_prot[age=='4M']$sample),'-','.')])))
plot(density(na.omit(protmat[,str_replace_all(paste0('X',mtd_prot[age=='7M']$sample),'-','.')])))
#just because of zero values?
#median
mtd_prot[,med.log2expr:=apply(protmat[,str_replace_all(paste0('X',mtd_prot$sample),'-','.')],2,median,na.rm = T)]
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=med.log2expr))
ggplot(mtd_prot)+geom_boxplot(aes(x=age,y=med.log2expr,col=sex))

#rm 7M zero values : need 75% non zero of NA values
protmatna<-protmat[(rowSums(is.na(protmat))/ncol(protmat))<0.25,]
nrow(protmatna)
#trans NA to 0
protmatna[is.na(protmatna)]<-0
protmatf<-protmatna[(rowSums(protmatna==0)/ncol(protmatna))<0.25,]
nrow(protmatf)#2680 /2790
#filter the main dataset
proteo_matf<-proteo_mat[gene_id%in%rownames(protmatf),]

fwrite(proteo_matf,'ref-data/raw_expr/log2_prot_expr_clean_filtered.csv.gz')



#annotate TFs prot
tfs_targ<-fread('ref-data/FT dérégulés en MSMS_WC vs WH 7M.csv',select = 1:2,skip = 1,col.names = c('TF','pred.target'))
tfs<-unique(tfs_targ$TF)
tfs
proteo_matf[,TF:=gene_name%in%tfs]
proteo_matf[(TF)]#52
nrow(proteo_matf[(TF)])

#merge prot and trscripto
proteo_matf[,mol:=paste0(gene_name,'_prot')]
rna_matf[,mol:=paste0(gene_name,'_trs')]

#rm duplicates
rna_matff<-unique(rna_matf[,avg.expr:=rowMeans(.SD),.SDcols=is.numeric][order(mol,-avg.expr)],by='mol')

mols<-rbind(proteo_mat[,.SD,.SDcols=intersect(colnames(rna_matff),colnames(proteo_mat))],rna_matff[,.SD,.SDcols=intersect(colnames(rna_matff),colnames(proteo_mat))])
dim(mols ) #20 common samples
samples_com<-colnames(mols)[str_detect(colnames(mols),'WC|WH')]

#which samples not common?
setdiff(colnames(rna_matff),colnames(proteo_mat))
#[1] "286-WC"    "288-WC"    "387-WH"    "845-WC-4M" "810-WC-4M" "891-WH-4M" "373-WC"
# [8] "241-WC"    "290-WC"    "363-WH"    "364-WH"    "391-WH"    "850-WC-4M" "870-WC-4M"
# [15] "895-WC-4M" "871-WH-4M" "884-WH-4M"

setdiff(colnames(proteo_mat),colnames(rna_matff)) #"849-WC-4M" "809-WC-4M" "282-WC"    "294-WH"


fwrite(mols,fp(out,'full_prot_trscripto_matrix.csv'))


# Trscriptomic change asso to a specific TF?####
#bulk refined SCENIC approach: TF-gene weight, > module identif > module filtering and gene filtering based on TF motif enrichment > merging all filtered gene module by TF

# 1) TF-gene expr weight using GENIE3
#rm NA in prot mat
proteo_matf[,n.na:=rowSums(is.na(.SD)),.SDcols=samples_com]
proteo_matff<-proteo_matf[n.na==0]
tfs<-proteo_matff[(TF)]$mol
length(tfs)#45

#perform GENIE coregulatory network analysis
weight_mat<-GENIE3(as.matrix(data.frame(mols[,-c('gene_name','gene_id')],row.names = 'mol')),
                   regulators =tfs,
                   targets = rna_matff$mol,
                   nCores = 16)
dim(weight_mat)
weight_mat[1:10,1:10]
## Get ranking of edges
linkList <- getLinkList(weight_mat)
head(linkList)
plot(density(linkList$weight))
linkList<-data.table(linkList)
linkList[regulatoryGene=='Nono_prot'][weight>0.04]

fwrite(linkList,fp(out,'tf_gene_weight.csv'))


#add sens of correlation using spearman correl
linkList<-fread(fp(out,'tf_gene_weight.csv'))
mols<-fread(fp(out,'full_prot_trscripto_matrix.csv'))
samples<-colnames(mols)[str_detect(colnames(mols),'WC|WH')]

#test 1
cor(unlist(mols[mol=='A2m_prot',.SD,.SDcols=samples]),
    unlist(mols[mol=='mt-Tl1_trs',.SD,.SDcols=samples]),
    method = 'spearman'
)

linkList[,r:=cor(unlist(mols[mol==regulatoryGene,.SD,.SDcols=samples]),
                 unlist(mols[mol==targetGene,.SD,.SDcols=samples]),
                 method = 'spearman'),by=.(regulatoryGene,targetGene)]

fwrite(linkList,fp(out,'tf_gene_weight_and_spearman.csv.gz'))
linkList<-fread(fp(out,'tf_gene_weight_and_spearman.csv.gz'))


#2) module identif: use differnt threshold
#In all the diff thresholds, only the links with IM > x were taken into account.
#determine x according to LinkScore distribution
#+ module by pos or neg correlation
linkList[,reg:=ifelse(r>0,'up','down')]
plot(density(linkList$weight))
abline(v=0.02)


unique(linkList$regulatoryGene)
plot(density(linkList[regulatoryGene=='Apex1_prot']$weight))
plot(density(linkList[regulatoryGene=='Arl6ip5_prot']$weight))
plot(density(linkList[regulatoryGene=='Bub3_prot']$weight))
plot(density(linkList[regulatoryGene=='Chd4_prot']$weight))
plot(density(linkList[regulatoryGene=='Cops5_prot']$weight))
plot(density(linkList[regulatoryGene=='Crtc1_prot']$weight))
plot(density(linkList[regulatoryGene=='Ctnnb1_prot']$weight))
plot(density(linkList[regulatoryGene=='Ddx5_prot']$weight))
abline(v=0.02) 
#x= 0.02
linkListf<-linkList[weight>0.02]
#different thresolds for TF modules: 
#(i) setting several IM thresholds (IM > 0.02 and IM > 0.05),
#(ii) taking the 50 targets with highest IM for each TF

tfs<-unique(linkListf$regulatoryGene)
mods_groups<-c('up','down')

mods<-lapply(mods_groups,function(s){
  mod_tfs<-lapply(tfs,function(tf){
    return(list(thr0.02=linkListf[reg==s&regulatoryGene==tf][weight>0.02]$targetGene,
                thr0.05=linkListf[reg==s&regulatoryGene==tf][weight>0.05]$targetGene,
                top50=linkListf[reg==s&regulatoryGene==tf][head(order(-weight),50)]$targetGene))
    
  })
  names(mod_tfs)<-tfs
  return(mod_tfs)
})
names(mods)<-mods_groups

mods

#(iii) keeping only the top 5, and 10  TFs for each target gene (then, split by TF).
mods2<-list(tartop5=5,
            tartop10=10)

mods2dt<-lapply(mods2,function(t){
  return(linkListf[,.SD[head(order(-weight),t)],by=.(targetGene,reg)])
})

for(s in names(mods)){
  
  for(tf in names(mods[[s]])){
    
    for(t in names(mods2dt)){
      mods[[s]][[tf]][[t]]<-mods2dt[[t]][reg==s&regulatoryGene==tf]$targetGene
      
    }
    
  }
}

mods$up
unlist(lapply(mods$up,function(tfmods)lapply(tfmods,function(x)length(x))))

saveRDS(mods,fp(out,'tf_modules_before_filtering.rds'))
mods<-readRDS(fp(out,'tf_modules_before_filtering.rds'))

#3) module filtering and gene filtering based on TF motif enrichment
#get mouse genes ranking by motif
system('wget -O ref-data/cistarget/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')

BiocManager::install("RcisTarget")
library(RcisTarget)

#rm '_trs' to have similar gene names
mods<-lapply(mods, function(lsens)lapply(lsens,function(lTF)lapply(lTF,function(lTHR)str_remove(lTHR,'_trs'))))
#rm '_prot' to have similar prot names
mods<-lapply(mods, function(lsens)setNames(lsens,str_remove(names(lsens),'_prot')))


#load ranking
motifRankings<-importRankings('ref-data/cistarget/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
head(motifRankings@rankings)


#load motif annotations
data(motifAnnotations_mgi)
#motif enrichment using fgsea
#devtools::install_github('ctlab/fgsea')

#annotate motifs for tested tf
tfs<-names(mods$up)
motifAnnof<-motifAnnotations[TF%in%tfs]
motifAnnof
motifs_ranks<-melt(data.table(t(data.frame(motifRankings@rankings,row.names = 'motifs')[motifAnnof$motif,]),keep.rownames = 'gene'),id.vars = 'gene',variable.name = 'motif',value.name = 'rank')
motifs_ranks<-merge(motifs_ranks,motifAnnof[,-'description'])

#filter for expressed genes
rna_matf<-fread('ref-data/raw_expr/raw_rna_counts_clean_filtered.csv.gz')

genes<-rna_matf$gene_name
length(genes) #14030
motifs_ranksf<-motifs_ranks[gene%in%genes]
rm(motifs_ranks)
motifs<-unique(motifs_ranksf$motif)
length(motifs)#298

#test motif enrichment for tested tfs 
#rm non unique motif-TF ranking
motifs_ranksf<-unique(motifs_ranksf[order(motif,TF,-directAnnotation,-inferred_Orthology,rank)],by=c('motif','TF','gene'))
res_up<-motifs_ranksf[,fgsea(pathways = mods[['up']][[TF[1]]],stats = setNames(rank,gene),minSize = 5,scoreType = "pos"),by=.(motif,TF)]
res_up[pval<0.01] #

unique(res_up[pval<0.01]$TF) #9: "Fhl2"  "Ube2k" "Mecp2" "Myef2" "Pura"  "Hdac2" "Ilf2"  "Apex1" "Nono" 
unique(res_up[pval<0.05]$TF) #11

fwrite(res_up,fp(out,'res_modules_up_tf_motif_enrichment.csv.gz'))

linkList[regulatoryGene=='Mecp2_prot'&targetGene=='Scand1_trs']
summary(sapply(mods[['up']],function(tfmod)sapply(tfmod,length)))
sum(duplicated(sort(setNames(motifs_ranksf[motif=='hdpi__FHL2'&TF=='Fhl2']$rank,motifs_ranksf[motif=='hdpi__FHL2'&TF=='Fhl2']$gene)))) #0


res_dn<-motifs_ranksf[,fgsea(pathways = mods[['down']][[TF[1]]],stats = setNames(rank,gene),minSize = 5,scoreType = "pos"),by=.(motif,TF)]
res_dn[pval<0.01] #
unique(res_dn[pval<0.01]$TF) #9: "Nf1"     "Pura"    "Ilf2"    "Smarca5" "Ctnnb1"  "Ilf3"    "Gtf2i"   "Trim28"  "Sfpq" 

fwrite(res_dn,fp(out,'res_modules_down_tf_motif_enrichment.csv.gz'))

#summary(sapply(mods[['down']],function(tfmod)sapply(tfmod,length)))


#(i) module filtering: 
#add NES threshold?
plot(density(res_up$NES))

#thr: pval<0.01
res_upsig<-res_up[pval<0.01]
res_dnsig<-res_dn[pval<0.01]

#(ii) gene filtering
#get all Leading Edges of signif modules by TFs
TFregsup<-res_upsig[,Reduce(union,leadingEdge),by='TF']
setnames(TFregsup,'V1','gene')
fwrite(TFregsup,fp(out,'TF_regulons_up.csv'))

TFregsdn<-res_dnsig[,Reduce(union,leadingEdge),by='TF']
setnames(TFregsdn,'V1','gene')
fwrite(TFregsdn,fp(out,'TF_regulons_down.csv'))

#stats
head(TFregsup[order(TF,gene)],100)
TFregsup[,size.regulon:=.N,by='TF']
table(TFregsup$TF)
# Apex1  Fhl2 Hdac2  Ilf2 Mecp2 Myef2  Nono  Pura Ube2k 
# 4883  4228   340   175  9702   608  2127  1270  1094 

TFregsdn[,size.regulon:=.N,by='TF']
table(TFregsdn$TF)
# Ctnnb1   Gtf2i    Ilf2    Ilf3     Nf1    Pura    Sfpq Smarca5  Trim28 
# 991    4149    5825     426    2670      93    2123    2094     678 


regsup<-fread('outputs/01-tf-downstream_gene-inference/TF_regulons_up.csv')
unique(regsup$TF)

regsdn<-fread('outputs/01-tf-downstream_gene-inference/TF_regulons_down.csv')
unique(regsdn$TF)



#ChipSeq Validation of regulon ####
#test MECP2 regulon
#MEPC2: MBD prot, binding specifically to methylated DNA.
# #MECP2 is dispensible in stem cells, but is essential for embryonic development. 
# #MECP2 gene mutations are the cause of most cases of Rett syndrome,
# #a progressive neurologic developmental disorder and one of the most common causes of cognitive disability in females
library(data.table)
out1<-fp(out,'validation_Mecp2')
dir.create(out1)
bedrep1<-fread('ref-data/validation/Mecp2/GSE71126_RAW/mecp2_rep1_c5.0_l200_g30_peaks.narrowPeak',
               select = c(1:5,7:9),
               col.names = c('chr','start','end','name','score','signalValue','pvalue','qvalue'))

bedrep1 #1.2M

bedrep2<-fread('ref-data/validation/Mecp2/GSE71126_RAW/mecp2_rep2_c5.0_l200_g30_peaks.narrowPeak',
               select = c(1:5,7:9),
               col.names = c('chr','start','end','name','score','signalValue','pvalue','qvalue'))
bedrep2 #1.4M
#merge common peak
#source('../utils/r_utils.R')
bedinter<-bed_inter(bedrep1,bedrep2)
bedinter#1209697

#annotate Mecp2 bed +\- 10kb TSS
gtfs<-fread('ref-data/validation/gencode.vM32.annotation.gtf.gz',select = c(1,2,3,4,5,7,9),col.names = c('chr','source','type','start','end','strand','gene_names'))
gtfs<-gtfs[type=='gene']
gtfs[,gene_id:=str_extract(gene_names,'ENSMUSG[0-9]+')]
gtfs[,gene_type:=str_remove_all(sapply(gene_names,function(x)str_split(x,';')[[1]][2]),'gene_type|"')]
gtfs[,gene_name:=str_remove_all(sapply(gene_names,function(x)str_split(x,';')[[1]][3]),'gene_name|"')]
#filter for genes expressed
expressed_genes<-fread('ref-data/raw_expr/raw_rna_counts_clean_filtered.csv.gz')$gene_id
length(expressed_genes)#14k
gtfsf<-gtfs[intersect(expressed_genes,gene_id), on='gene_id']
gtfsf[,chr:=factor(chr,levels = paste0('chr',c(1:22,'X','Y')))]
gtfsf[,startwin:=ifelse(strand=='+',start-10000,end-10000)]
gtfsf[startwin<0,startwin:=1]
gtfsf[,endwin:=ifelse(strand=='+',start+10000,end+10000)]

bedinter_anno<-bed_inter(bedinter[,.(V1,V2,V3,V4,V5)],
                         gtfsf[,.(chr,startwin,endwin,gene_id,gene_name)][order(chr,startwin)],
                         select =c(1:5,9:10),
                         col.names = c('chr','start','end','peak','score','gene_id','gene_name'))
bedinter_anno#185k peak in expressed genes
unique(bedinter_anno,by='gene_id') #in 13k genes / 14k genes

#overlap regulon signif ? => not too restrictive
Mecp2p_reg<-fread('outputs/01-tf-downstream_gene-inference/TF_regulons_up.csv')[TF=='Mecp2']
Mecp2p_reg #9.7kgenes 

#score enriched in regulon?
expressed_genes<-fread('ref-data/raw_expr/raw_rna_counts_clean_filtered.csv.gz')$gene_name

expressed_genes_mecp2pbind<-unique(bedinter_anno[expressed_genes,on='gene_name'][order(-score)],by = 'gene_name')
expressed_genes_mecp2pbind[is.na(score),score:=0]
expressed_genes_mecp2pbind[,in_regulon:=gene_name%in%Mecp2p_reg$gene]

expressed_genes_mecp2pbind[,`in Mecp2 regulon`:=gene_name%in%Mecp2p_reg$gene]

expressed_genes_mecp2pbind[,p:=wilcox.test(score[(in_regulon)],score[!(in_regulon)])$p.value]
fwrite(expressed_genes_mecp2pbind,fp(out1,'expressed_genes_ChIPscore_Mecp2.csv.gz'))

#ANNOT modules####
library(data.table)
#install.packages('gprofiler2')
library(gprofiler2)

source('../utils/r_utils.R')
source('../utils/visualisation.R')

#gos_pathways<-fread('/Users/alexandrepelletier/adpelle1@bu.edu - Google Drive/My Drive/P.Alexandre_Pelletier/Common_resources/MSigDB/all_CPandGOs_gene_and_genesets.csv.gz')

#Hypergeometric test on GO and Pathways
regulons<-rbind(fread('outputs/01-tf-downstream_gene-inference/validation_norm_data/TF_regulons_up.csv')[,regul_type:='up'],
                fread('outputs/01-tf-downstream_gene-inference/validation_norm_data/TF_regulons_down.csv')[,regul_type:='down'])
regulons[,regulon:=paste(TF,regul_type,sep='_')]

#on whole regulon
?gost
regs_to_annot<-unique(regulons$regulon)
res_enr<-gost(split(regulons[regs_to_annot,on='regulon']$gene,regulons[regs_to_annot,on='regulon']$regulon),
              organism = 'mmusculus',domain_scope = 'annotated')
res_enrdt<-data.table(res_enr$result)  

fwrite(res_enrdt,fp(out1,'res_annot_regulons.csv.gz'))




#REGULONS ENRICHMENT FOR HFD DEGSs#####
#Are the regulons enriched in DEs genes in a condition?
res_de<-Reduce(rbind,lapply(c('4M','7M'),function(m)Reduce(rbind,lapply(c('male','femelle'),function(s)fread(fp('ref-data',paste0('res_de_',m,'_',s,'.csv')),
                                                                                                             select = 1:5,
                                                                                                             col.names = c('gene_id','gene_name','log2FC','pvalue','padj'))[
                                                                                                               ,age:=m][
                                                                                                                 ,sex:=s][
                                                                                                                   ,compa:='HFDvsChow']))))
#filter for expressed genes and tested genes
res_def<-res_de[gene_name%in%genes&!is.na(padj)]
#rm duplicate
res_def<-unique(res_def[order(pvalue)],by = 'gene_name')
unique(res_def,by='gene_name') #14k
res_def[,gene_stat:=sign(log2FC)*-log10(pvalue),by=.(age,sex,compa)]


res_regup<-res_def[,fgsea(pathways = split(TFregsup$gene,TFregsup$TF),
                          stats = setNames(gene_stat,gene_name)),by=.(age,sex,compa)]

res_regup[padj<0.01] 
fwrite(res_regup,fp(out,'res_fgsea_degs_in_regulons_up.csv'))

res_regdn<-res_def[,fgsea(pathways = split(TFregsdn$gene,TFregsdn$TF),
                          stats = setNames(gene_stat,gene_name)),by=.(age,sex,compa)]

res_regdn[padj<0.01]
fwrite(res_regdn,fp(out,'res_fgsea_degs_in_regulons_down.csv'))

#regulon up in HFD have prot up ?
res_regup<-fread(fp(out,'res_fgsea_degs_in_regulons_up.csv'))
#up in HFD 7Mfemale
res_regup[padj<0.01&age=='7M'&sex=='femelle'&NES>0]
#9:  7M femelle HFDvsChow   Myef2 4.446688e-30 4.002019e-29 1.4245627  0.5498256  3.528174  197

myef2p<-melt(proteo_mat[gene_name=='Myef2'],id.vars = c('gene_name','gene_id'),variable.name = 'sample',value.name = 'log2expr')
myef2p<-merge(myef2p,mtd_prot)
ggplot(myef2p)+geom_boxplot(aes(x=sex,y=log2expr,fill=condition))+
  facet_grid(gene_name~age,scales = 'free')
#same for male
res_regup[padj<0.01&age=='7M'&sex=='male'&NES>0]
# 3:  7M male HFDvsChow   Mecp2 1.000000e-50 8.000000e-50


Mecp2p<-melt(proteo_mat[gene_name=='Mecp2'],id.vars = c('gene_name','gene_id'),variable.name = 'sample',value.name = 'log2expr')
Mecp2p<-merge(Mecp2p,mtd_prot)
ggplot(Mecp2p)+geom_boxplot(aes(x=sex,y=log2expr,fill=condition))+
  facet_grid(gene_name~age,scales = 'free')

#common in male and female ? no
res_regup[padj<0.01&age=='7M'&NES>0]
#reg down
res_regdn<-fread(fp(out,'res_fgsea_degs_in_regulons_down.csv'))
res_regdn[padj<0.01]
# 7:  7M    male HFDvsChow  Ctnnb1 1.000000e-50 2.500000e-50
# 8:  7M    male HFDvsChow    Sfpq 1.000000e-50 2.500000e-50
# 6:  7M femelle HFDvsChow Smarca5 1.227476e-22 1.104728e-21

downp<-melt(proteo_mat[gene_name%in%c('Ctnnb1','Sfpq','Smarca5')],id.vars = c('gene_name','gene_id'),variable.name = 'sample',value.name = 'log2expr')
downp<-merge(downp,mtd_prot)
ggplot(downp)+geom_boxplot(aes(x=sex,y=log2expr,fill=condition))+
  facet_grid(gene_name~age,scales = 'free')



