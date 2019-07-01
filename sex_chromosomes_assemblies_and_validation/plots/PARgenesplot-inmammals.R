
library(GenomicRanges)

library(Gviz)
library(GenomicAlignments)

library(rtracklayer)

#####################
gffRangedData<-import.gff("pargene-waterbuffalo.gff3")
myGranges<-as(gffRangedData, "GRanges")

trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="chrX",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=0,to=6800000)


#############
gffRangedData<-import.gff("pargene-arsucd.gff3")
myGranges<-as(gffRangedData, "GRanges")

trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="chrX",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=132000000,to=139000000)

pdf("arsucd-par.pdf",width=20,height=1)
plotTracks(trTrack,from=132000000,to=139000000)
 dev.off()

#############
gffRangedData<-import.gff("pargene-goat.gff3")
myGranges<-as(gffRangedData, "GRanges")

trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="chrX",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=58632000,to=65632000)

pdf("goat-par.pdf",width=20,height=1)
plotTracks(trTrack,from=58632000,to=65632000)
 dev.off()

#############
gffRangedData<-import.gff("pargene-human.gff3")
myGranges<-as(gffRangedData, "GRanges")

trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="chrX",genome="cow", fontsize=15, col="darkblue", fill="lightblue", name="",stacking="squish",col="lightblue")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=1,to=9949500)

pdf("human-par-v2.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=9949500)
 dev.off()

#############
gffRangedData<-import.gff("pargene_pig.gff")
myGranges<-as(gffRangedData, "GRanges")

trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="chrX",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=1,to=6500000)

pdf("pig-par.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=6500000)
 dev.off()
#############
gffRangedData<-import.gff("Ygene-withshroom2.gff")
myGranges<-as(gffRangedData, "GRanges")


trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",chromosome="sireX",genome="cow", fontsize=15, col="purple", fill="purple", name="",stacking="squish")

options(ucscChromosomeNames=FALSE) 

pdf("X-degenerate-region.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=8300000)
 dev.off()
#############
gffRangedData<-import.gff("pargene_sheep.gff")
myGranges<-as(gffRangedData, "GRanges")
options(ucscChromosomeNames=FALSE)
trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=1,to=7100000)

pdf("sheep-par.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=7100000)
 dev.off()
#############
gffRangedData<-import.gff("pargene_dog.gff")
myGranges<-as(gffRangedData, "GRanges")
options(ucscChromosomeNames=FALSE)
trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",genome="cow", fontsize=15, col="darkblue", fill="darkblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=1,to=7000000)

pdf("dog-par.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=7000000)
 dev.off()
#############
gffRangedData<-import.gff("pargene_horse.gff")
myGranges<-as(gffRangedData, "GRanges")
options(ucscChromosomeNames=FALSE)
trTrack <- GeneRegionTrack(myGranges, fontcolor.group="black",genome="cow", fontsize=15, col="lightblue", fill="lightblue", name="",stacking="squish")

symbol(trTrack) <-myGranges$gene
isna=is.na(symbol(trTrack))
trTrack=trTrack[!isna,]
plotTracks(trTrack,from=1,to=7000000)

pdf("horse-par-light.pdf",width=20,height=1)
plotTracks(trTrack,from=1,to=7000000)
 dev.off()
