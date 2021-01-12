# metagenome

## Cleanup Trimmomatic

```
$ export TRIM=/home/impact/biotools/local/assembly/Trimmomatic-0.39

#**** **** **** **** **** **** 
$ mkdir cleanup; \
 for i in ./fastq/*_L001_R1_001.fastq.gz; do a=$i; \
 j=${i##./*/}; k=${j%%_S*_L001_R1_001.fastq.gz}; echo $k; \
 ## trimmomatic
 java -jar $TRIM/trimmomatic-0.39.jar PE -threads 24 -phred33 $i ${i/_R1_/_R2_} ./cleanup/${k}_R1_PE.fq ./cleanup/${k}_R1_OP.fq ./cleanup/${k}_R2_PE.fq ./cleanup/${k}_R2_OP.fq ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100 HEADCROP:3 &> trimmomatic_log.txt; \
 done
#**** **** **** **** **** **** 

##==> ./cleanup/xxx_R1_PE.fq ./cleanup/xxx_R2_PE.fq

## COUNT READS

~/biotools/local/assembly/bin/seqkit stats ./cleanup/xxx_R1_PE.fq

```
## Assembly    SPades

```
$ export SPADE_HOME=/home/impact/biotools/local/assembly/SPAdes-3.14.0-Linux/bin

read1="../201211/cleanup/usoriyamakoNo5-1-M2_R1_PE.fq";
$SPADE_HOME/spades.py --pe-1 1 ${read1} --pe-2 1 ${read1/_R1/_R2} -m 200 --careful -k 21,33,55,77,99,127 -o usoriyamakoNo5-1-M2_spd_A1_assembly

```

## Statics of assemblies

```
##****** STATICS OF CONTIGS FROM FASTA ********
i="usoriyamakoNo5-1-M2_spd_A1_assembly.fas"; \
gc_contentSkew.pl -if $i -p gc; \
perl -F'\t' -anle 'BEGIN{$info={}; $info->{max}= 0; $info->{min}= 0} next if ($F[1] !~ /^\d+$/); ($con,$clen,$cgc,$ccov)=(@F[0,1,3],0); next if ($clen < 500 || $ccov < 0); $info->{count}++; push @{$info->{lens}},$clen; $info->{sumlen}+= $clen; $info->{sumgc}+= $clen*$cgc; $info->{sumcov}+= $clen*$ccov; END{ @sortlen= sort{$b<=>$a}@{$info->{lens}}; print "COUNT:\t$info->{count}"; print "SUM:\t$info->{sumlen}"; print "MIN:\t$sortlen[-1]"; print "MAX:\t@sortlen[0..9]"; print "AVE:\t".sprintf("%d",$info->{sumlen}/$info->{count}); foreach (@sortlen) {$info->{N50}+= $_; if ($info->{N50} > $info->{sumlen}/2) { $N50=$_; last; }  } print "N50:\t$N50"; print "GC:\t".sprintf("%.4f",$info->{sumgc}/$info->{sumlen}); print "Cov:\t".sprintf("%.4f",$info->{sumcov}/$info->{sumlen});}' outgc

##****** ****** ****** ****** ******
##****** STATICS OF CONTIGS FROM TABLE ********
# i="usoriyamakoNo5-1-M2_spd_A1_assembly_gc.txt"; \
perl -F'\t' -anle 'BEGIN{$info={}; $info->{max}= 0; $info->{min}= 0} next if ($F[1] !~ /^\d+$/); ($con,$clen,$ccov,$cgc)=@F[0,1,6,3]; next if ($clen < 500 || $ccov < 0); $info->{count}++; push @{$info->{lens}},$clen; $info->{sumlen}+= $clen; $info->{sumgc}+= $clen*$cgc; $info->{sumcov}+= $clen*$ccov; END{ @sortlen= sort{$b<=>$a}@{$info->{lens}}; print "COUNT:\t$info->{count}"; print "SUM:\t$info->{sumlen}"; print "MIN:\t$sortlen[-1]"; print "MAX:\t@sortlen[0..9]"; print "AVE:\t".sprintf("%d",$info->{sumlen}/$info->{count}); foreach (@sortlen) {$info->{N50}+= $_; if ($info->{N50} > $info->{sumlen}/2) { $N50=$_; last; } } print "N50:\t$N50"; print "GC:\t".sprintf("%.4f",$info->{sumgc}/$info->{sumlen}); print "Cov:\t".sprintf("%.4f",$info->{sumcov}/$info->{sumlen});}' $i
##****** ****** ****** ****** ******

```

## Distribution of contig

```
$ R
> source("/home/impact/biotools/local/statics/R/Library/frequency.R",encoding="utf-8")
# > library(KernSmooth)
 lrd<-log10(dat[,6]+1)
 clen<-dat[,1]
 summary(lrd)
# Min. 1st Qu.  Median  Mean 3rd Qu.   Max. 
#0.4575  2.3108  2.4507  2.3224  2.6382  3.5308 

hmin<- 0
hmax<- 3.6
hwide<- 0.1
stepnum<-(hmax-hmin)/hwide
freqLensum_1<-freqLensum_2P(lrd,clen,hmin,hmax,hwide) ## log10Redun CLen
xmin<-hmin
xmax<-hmax
xmin_stepnum<-(xmin-hmin)/hwide
xmax_stepnum<-(xmax-hmin)/hwide

#postscript("xxx.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 9, width = 9)

barplot(freqLensum_1[,"SubtotalRate"],xlab="LOG10(Coverage)",ylab="%rate of Subtotal length ",axes=F,axisnames=F,xlim=c(xmin_stepnum,xmax_stepnum),ylim=c(0,40),col="#ff00ff40",border="#ff00ff",space=0)
axis(2)
par(new = T)
plot(seq(0.5,stepnum-0.5,1),freqLensum_1[,"CumulSubtotalRate"],
 type="b", pch=1, lwd=0.5,cex=0.5, 
 xlim=c(xmin_stepnum,xmax_stepnum), ylim=c(0, 100),
 xlab="", ylab="", axes=FALSE)
axis(1, at=seq(0.5,stepnum-0.5,5), labels=seq(hmin,hmax-hwide,hwide*5),las=2)
axis(4)

```

## 
