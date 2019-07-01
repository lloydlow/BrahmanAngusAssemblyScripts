#goat PAR
x2=read.delim("Downloads/ref_ASM170441v1_top_level.gff3",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NW_017189516.1",]
xchr2=x2[x2[,1]=="NW_017189517.1",]
xchr1gene=xchr1[xchr1[,3]=="gene",]
xchr2gene=xchr2[xchr2[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9

xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]

xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[337:nrow(xchr1gene),]

xchr2genetype=apply(as.matrix(xchr2gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr2genetype=apply(as.matrix(xchr2genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr2genetype=="protein_coding"
xchr2gene=xchr2gene[sel,]

xchr2genename=apply(as.matrix(xchr2gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr2genename=apply(as.matrix(xchr2genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr2gene[337:nrow(xchr2gene),]
##################################
#ARsucd PAR

x2=read.delim("Downloads/ref_ARS-UCD1.2_top_level.gff3",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_037357.1",]
genes=which(xchr1[,3]=="gene")
mRNA=genes+1
xchr1gene=xchr1[xchr1[,3]=="gene",]

mRNA=xchr1[mRNA,]

mRNA= mRNA[mRNA[,3]=="mRNA",]



mRNAgenetype=apply(as.matrix(mRNA[,9]),1,function(x){strsplit(x,";")[[1]][5]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[815:nrow(xchr1gene),]

write.table(selectedgene,file="pargene-arsucd1.gff3",sep="\t",col.names=F,row.names=F,quote=F)

##################################
#waterbaffalo PAR

x2=read.delim("Downloads/GCF_003121395.1_UOA_WB_1_genomic.gff",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_037569.1",]
xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:50,]

write.table(selectedgene,file="pargene-waterbuffalo.gff3",sep="\t",col.names=F,row.names=F,quote=F)

##################################
#human PAR

x2=read.delim("Downloads/GCF_000001405.38_GRCh38.p12_genomic.gff",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_000023.11",]
#NC_000024.10
xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:45,]

write.table(selectedgene,file="pargene-human.gff3",sep="\t",col.names=F,row.names=F,quote=F)
##################################
#pig PAR

x2=read.delim("Downloads/GCF_000003025.6_Sscrofa11.1_genomic(1).gff",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_010461.5",]
#chr1=x2[x2[,1]=="NC_010462.3",] Y 

xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:45,]
pig=xchr1gene[1:18,]
write.table(pig,file="pargene_pig.gff",sep="\t",row.names=F,col.names=F,quote=F)
##################################
#dog PAR

x2=read.delim("Downloads/ref_CanFam3.1_top_level.gff3",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_006621.3",]
xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:50,]
isna=is.na(selectedgene[,1])
dog=selectedgene[!isna,]
write.table(dog,file="pargene_dog.gff",sep="\t",row.names=F,col.names=F,quote=F)
##################################
#sheep PAR

x2=read.delim("Downloads/ref_Oar_v4.0_top_level.gff3",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_019484.2",]
xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:33,]
isna=is.na(selectedgene[,1])
sheep=selectedgene[!isna,]
write.table(sheep,file="pargene_sheep.gff",sep="\t",row.names=F,col.names=F,quote=F)

##################################
#horse PAR

x2=read.delim("Downloads/ref_EquCab3.0_top_level.gff3",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_009175.3",]
xchr1gene=xchr1[xchr1[,3]=="gene",]


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][7]})
isna=is.na(xchr1genetype)

xchr1gene9=apply(as.matrix(xchr1gene[!isna,9]),1,function(x){paste(strsplit(as.character(x),";")[[1]][c(1,2,3,5,6,7)],collapse=";")})
xchr1gene[,9]=as.character(xchr1gene[,9])
xchr1gene[!isna,9]=xchr1gene9


xchr1genetype=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][6]})
xchr1genetype=apply(as.matrix(xchr1genetype),1,function(x){strsplit(x,"gene_biotype=")[[1]][2]})

sel=xchr1genetype=="protein_coding"
xchr1gene=xchr1gene[sel,]


xchr1genename=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(x,";")[[1]][3]})
xchr1genename=apply(as.matrix(xchr1genename),1,function(x){strsplit(x,"Name=")[[1]][2]})


selectedgene=xchr1gene[1:40,]
isna=is.na(selectedgene[,1])
horse=selectedgene[!isna,]
write.table(horse,file="pargene_horse.gff",sep="\t",row.names=F,col.names=F,quote=F)


##################################
#btau5.0.1 Y

x2=read.delim("Downloads/GCF_000003205.7_Btau_5.0.1_genomic.gff",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_016145.1",]
xchr1gene=xchr1[xchr1[,3]=="mRNA",]


transciptid=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(as.character(x),"transcript_id=")[[1]][2]})

write.table(transciptid,file="transcriptID_btau501_Y",sep="\t",quote=F,row.names=F,col.names=F)


#bos_indicus Y


x2=read.delim("Downloads/GCF_000247795.1_Bos_indicus_1.0_genomic.gff",sep="\t",header=F,comment.char="#")


xchr1=x2[x2[,1]=="NC_032680.1",]
xchr1gene=xchr1[xchr1[,3]=="mRNA",]


transciptid=apply(as.matrix(xchr1gene[,9]),1,function(x){strsplit(as.character(x),"transcript_id=")[[1]][2]})

write.table(transciptid,file="transcriptID_bosindicus1_Y",sep="\t",quote=F,row.names=F,col.names=F)

