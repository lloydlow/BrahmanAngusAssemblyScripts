#########################################################################################
# Hi-C and optical map scaffolding comparisons#
#########################################################################################
                                        #process bionano scaffolds
 angusmashmap=read.delim("bostaurus_angus_vs_sire_cleaned_assembly.mashmap",sep=" ",header=F)
 brahmanmashmap=read.delim("bostaurus_brahma_vs_dam_cleaned_assembly.mashmap",sep=" ",header=F)
 angusbionano=read.table("EXP_REFINEFINAL1_bppAdjust_cmap_bostaurus_angus_fasta_NGScontigs_HYBRID_SCAFFOLD.agp")
 brahmanbionano=read.table("EXP_REFINEFINAL1_bppAdjust_cmap_bostaurus_brahma_fasta_NGScontigs_HYBRID_SCAFFOLD.agp")
 
 load("Downloads/sire_agp_clean_assembly_to_salsa.RData")
 
load("Downloads/dam_agp_clean_assembly_to_salsa.RData")
#########################################################################################
                                        #X
#########################################################################################
x=read.table("contig_order_v4.txt",se="\t")
newx=apply(as.matrix(x[,1]),1,function(x){strsplit(x,".fa")[[1]][1]})
X=cbind(newx,x[,2])


brahmanmashmapX=brahmanmashmap[brahmanmashmap[,6]%in%X[,1],]
brahmanmashmapX=brahmanmashmapX[match(X[,1],brahmanmashmapX[,6]),]
brahmanmashmapXtig=apply(as.matrix(brahmanmashmapX[,1]),1,function(x){strsplit(x,"\\|")[[1]][1]})
brahmanmashmapXnew=cbind(brahmanmashmapX,brahmanmashmapXtig)

brahmanbionanotig=apply(as.matrix(brahmanbionano[,6]),1,function(x){strsplit(x,"\\|")[[1]][1]})
brahmanbionanotignew=cbind(brahmanbionano,brahmanbionanotig)

mergedbrahman=merge(brahmanmashmapXnew,brahmanbionanotignew,by.x="brahmanmashmapXtig",by.y="brahmanbionanotig")
supermergedbrahman=mergedbrahman[,c(1,3,4,5,6,7,12,13,14,15)]

m=match(X[,1],mergedbrahman[,7])
mergedbrahman=mergedbrahman[m,]


checking=brahmanbionano[!is.na(mm),]
checking1=as.matrix(checking[checking[,5]=="W",])
write.table(dam_scaffolds_FINAL_agp_modi[match(X[,1],dam_scaffolds_FINAL_agp_modi[,6]),][isna,] ,file="missingtable-Xmatchedbionano.txt",row.names=F,col.names=F,sep="\t")
write.table(dam_scaffolds_FINAL_agp_modi[match(X[,1],dam_scaffolds_FINAL_agp_modi[,6]),],file="finaltable-Xsalsa.txt",row.names=F,col.names=F,sep="\t")

write.table(supermergedbrahman,file="finaltable-Xmatchedbionano.txt",row.names=F,col.names=F,sep="\t")

write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00011795",],file="missedtig00011795.csv")
write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00003246",],file="missedtig00003246.csv")
write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00011886",],file="missedtig00011886.csv")
write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00011907",],file="missedtig00011907.csv")
write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00479033",],file="missedtig00479033.csv")
write.csv(brahmanbionanotignew[brahmanbionanotignew[,10]=="tig00011826",],file="missedtig00011826.csv")

#########################################################################################
                                        #dam full
#########################################################################################




colnames(brahmanmashmap)[6]="contig"
brahmanmashmapXnew2=merge(brahmanmashmap,dam_scaffolds_FINAL_agp_modi,by.x="contig",by.y="component_id")
colnames(brahmanmashmapXnew2)[2]="tig"

brahmanbionanotig=apply(as.matrix(brahmanbionano[,6]),1,function(x){strsplit(x,"_")[[1]][1]})
brahmanbionanotignew=cbind(brahmanbionano,brahmanbionanotig)


mergedbrahman=merge(brahmanmashmapXnew2,brahmanbionanotignew,by.x="tig",by.y="brahmanbionanotig")




#########################################################################################
                                        #Y
#########################################################################################
y=read.csv("sireY_scaffold_info.csv",header=F)
y=y[y[,5]=="A",]
newYcontig=sire_scaffolds_FINAL_agp_modi[sire_scaffolds_FINAL_agp_modi[,1]%in%y[,6],]


for(i in 1:nrow(Y)){
    if(i==1){
        newYcontigbase=sire_scaffolds_FINAL_agp_modi[sire_scaffolds_FINAL_agp_modi[,1]%in%y[1,6],]
        oriation=rep(as.character(y[i,9]),nrow(newYcontigbase))
    }else{
                newYcontigbase2=sire_scaffolds_FINAL_agp_modi[sire_scaffolds_FINAL_agp_modi[,1]%in%y[i,6],]
        oriation2=rep(as.character(y[i,9]),nrow(newYcontigbase2))

                      newYcontigbase=rbind(newYcontigbase,newYcontigbase2)
                                 oriation=c(oriation,oriation2)

    }
}



Y=cbind(newYcontigbase[,6],oriation)


angusmashmapY=angusmashmap[angusmashmap[,6]%in%Y[,1],]
angusmashmapY=angusmashmapY[match(Y[,1],angusmashmapY[,6]),]
angusmashmapYtig=apply(as.matrix(angusmashmapY[,1]),1,function(x){strsplit(x,"\\|")[[1]][1]})
angusmashmapYnew=cbind(angusmashmapY,angusmashmapYtig)

angusbionanotig=apply(as.matrix(angusbionano[,6]),1,function(x){strsplit(x,"\\|")[[1]][1]})
angusbionanotignew=cbind(angusbionano,angusbionanotig)

mergedangus=merge(angusmashmapYnew,angusbionanotignew,by.x="angusmashmapYtig",by.y="angusbionanotig")
supermergedangus=mergedangus[,c(1,3,4,5,6,7,12,13,14,15)]

m=match(Y[,1],supermergedangus[,6])
supermergedangus=supermergedangus[m,]

table(as.matrix(supermergedangus[,7]))

special=rownames(table(as.matrix(supermergedangus[,7])))[1:12]
> special 
 [1] "Super-Scaffold_100035" "Super-Scaffold_100072" "Super-Scaffold_100129" "Super-Scaffold_100133" "Super-Scaffold_100160" "Super-Scaffold_100217" "Super-Scaffold_100308" "Super-Scaffold_100367"
 [9] "Super-Scaffold_100404" "Super-Scaffold_100406" "Super-Scaffold_100477" "Super-Scaffold_100494"



checking=brahmanbionano[!is.na(mm),]
checking1=as.matrix(checking[checking[,5]=="W",])
write.table(sire_scaffolds_FINAL_agp_modi[match(X[,1],sire_scaffolds_FINAL_agp_modi[,6]),][isna,] ,file="missingtable-Xmatchedbionano.txt",row.names=F,col.names=F,sep="\t")
write.table(newYcontigbase,file="finaltable-Ysalsa.txt",row.names=F,col.names=F,sep="\t")

 write.table(supermergedangus,file="finaltable-Ymatchedbionano.txt",row.names=F,col.names=F,sep="\t")
superscaffold=angusbionano[angusbionano[,1]%in%rownames(table(as.character(as.matrix(supermergedangus[,7]))))[1:11],]
specialtig=superscaffold[superscaffold[,5]=="W",6]
specialtig=apply(as.matrix(specialtig),1,function(x){strsplit(x,"_")[[1]][1]})

specialtig=rownames(table(specialtig))
 sire_scaffolds_FINAL_agp_modi[match(angusmashmap[match(specialtig,angusmashmap[,1]),6],sire_scaffolds_FINAL_agp_modi[,6]),]


 rownames(table(as.character(as.matrix(supermergedangus[,7]))))[1:11]
angusbionano[angusbionano[,1]%in%rownames(table(as.character(as.matrix(supermergedangus[,7]))))[1:11],6]
