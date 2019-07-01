###########################################################################################
###########################################################################################
#aligning the transcripts from ARSUCD-1.2 and Bos indicus 1.0 to our primary brahman X chromosome to checking  protential X chormosome contigs we missed in primary assembly
###########################################################################################
###########################################################################################

load("~/Downloads/dam_agp_clean_assembly_to_salsa.RData")

damscaffoldinfo=read.delim("~/Downloads/dam_scaffold_salsa.agp",sep="\t",header=F)
damXscaffoldinfo=damscaffoldinfo[damscaffoldinfo[,1]=="X",]
damXscaffoldinfo=damscaffoldinfo[damscaffoldinfo[,5]=="A",]

load("~/Downloads/sire_agp_clean_assembly_to_salsa.RData")

sirescaffoldinfo=read.delim("~/Downloads/sire_scaffold_salsa.agp",sep="\t",header=F)
sireXscaffoldinfo=sirescaffoldinfo[sirescaffoldinfo[,1]=="X",]
sireXscaffoldinfo=sirescaffoldinfo[sirescaffoldinfo[,5]=="A",]

x2=read.delim("Downloads/ref_ARS-UCD1.2_top_level.gff3",sep="\t",header=F,comment.char="#")

refX=x2[x2[,1]=="NC_037357.1",]
refXgene=refX[refX[,3]=="mRNA",]
transciptid=apply(as.matrix(refXgene[,9]),1,function(x){strsplit(as.character(x),"transcript_id=")[[1]][2]})
genename=apply(as.matrix(refXgene[,9]),1,function(x){strsplit(as.character(x),"gene=")[[1]][2]})

genename=apply(as.matrix(genename),1,function(x){strsplit(as.character(x),";")[[1]][1]})

arsanno=cbind(transciptid,genename)
arsanno=cbind(refXgene[,1:8],arsanno)
###########################################################################################
###########################################################################################
damnotfit=read.delim("arsucdXvsDamFULLnotinfirst-88.gff",sep="\t",header=F,comment.char="#")
damnotfitgene=damnotfit[damnotfit[,3]=="gene",]

damtransciptid=apply(as.matrix(damnotfitgene[,9]),1,function(x){strsplit(as.character(x),"\\|")[[1]][2]})
damtransciptid=apply(as.matrix(damtransciptid),1,function(x){strsplit(as.character(x),"_cds")[[1]][1]})
damuniquetranscript=unique(damtransciptid)

dammatachedscaffold=match(damnotfitgene[,1],damXscaffoldinfo[,6])

damscaffoldinfo[dammatachedscaffold[!is.na(dammatachedscaffold)],]

dam_scaffolds_FINAL_agp_modi[dam_scaffolds_FINAL_agp_modi[,1]%in%unique(damnotfitgene[,1]),]

arsanno[match(damuniquetranscript,arsanno[,9]),]


sirenotfit=read.delim("arsucdXvsSireFULLnotinfirst-88.gff",sep="\t",header=F,comment.char="#")
sirenotfitgene=sirenotfit[sirenotfit[,3]=="gene",]

siretransciptid=apply(as.matrix(sirenotfitgene[,9]),1,function(x){strsplit(as.character(x),"\\|")[[1]][2]})
siretransciptid=apply(as.matrix(siretransciptid),1,function(x){strsplit(as.character(x),"_cds")[[1]][1]})
sireuniquetranscript=unique(siretransciptid)

sirematachedscaffold=match(sirenotfitgene[,1],sirescaffoldinfo[,6])
sirescaffoldinfo[sirematachedscaffold[!is.na(sirematachedscaffold)],]

match(sireuniquetranscript,damuniquetranscript)


###########################################################################################
###########################################################################################
bosindicus1damnotfit=read.delim("bosindicus1vsdamall-88-notlift.gff",sep="\t",header=F,comment.char="#")
bosindicus1damnotfitgene=bosindicus1damnotfit[bosindicus1damnotfit[,3]=="gene",]

bosindicus1damtransciptid=apply(as.matrix(bosindicus1damnotfitgene[,9]),1,function(x){strsplit(as.character(x),"\\|")[[1]][2]})
bosindicus1damtransciptid=apply(as.matrix(bosindicus1damtransciptid),1,function(x){strsplit(as.character(x),"_cds")[[1]][1]})
bosindicus1damuniquetranscript=unique(bosindicus1damtransciptid)

bosindicus1dammatachedscaffold=match(bosindicus1damnotfitgene[,1],damXscaffoldinfo[,6])

bosindicus1dammatachedscaffold[bosindicus1dammatachedscaffold[!is.na(bosindicus1dammatachedscaffold)],]

dam_scaffolds_FINAL_agp_modi[dam_scaffolds_FINAL_agp_modi[,1]%in%unique(bosindicus1damnotfitgene[,1]),]

arsanno[match(damuniquetranscript,arsanno[,9]),]
