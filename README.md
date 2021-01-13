# metagenome

### NAME: xxx

## Cleanup Trimmomatic

```
$ export TRIM=$HOME/biotools/local/assembly/Trimmomatic-0.39

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
$ cd ./assembly

export SPADE_HOME=$HOME/biotools/local/assembly/SPAdes-3.14.0-Linux/bin

read1="../201211/cleanup/xxx_R1_PE.fq";
$SPADE_HOME/spades.py --pe-1 1 ${read1} --pe-2 1 ${read1/_R1/_R2} -m 200 --careful -k 21,33,55,77,99,127 -o xxx_spd_A1_assembly

```

## Statics of assemblies

```
##****** STATICS OF CONTIGS FROM FASTA ********
i="xxx_spd_A1_assembly.fas"; \
gc_contentSkew.pl -if $i -p gc; \
perl -F'\t' -anle 'BEGIN{$info={}; $info->{max}= 0; $info->{min}= 0} next if ($F[1] !~ /^\d+$/); ($con,$clen,$cgc,$ccov)=(@F[0,1,3],0); next if ($clen < 500 || $ccov < 0); $info->{count}++; push @{$info->{lens}},$clen; $info->{sumlen}+= $clen; $info->{sumgc}+= $clen*$cgc; $info->{sumcov}+= $clen*$ccov; END{ @sortlen= sort{$b<=>$a}@{$info->{lens}}; print "COUNT:\t$info->{count}"; print "SUM:\t$info->{sumlen}"; print "MIN:\t$sortlen[-1]"; print "MAX:\t@sortlen[0..9]"; print "AVE:\t".sprintf("%d",$info->{sumlen}/$info->{count}); foreach (@sortlen) {$info->{N50}+= $_; if ($info->{N50} > $info->{sumlen}/2) { $N50=$_; last; }  } print "N50:\t$N50"; print "GC:\t".sprintf("%.4f",$info->{sumgc}/$info->{sumlen}); print "Cov:\t".sprintf("%.4f",$info->{sumcov}/$info->{sumlen});}' outgc
##****** ****** ****** ****** ******

##****** STATICS OF CONTIGS FROM TABLE ********
# i="xxx_spd_A1_assembly_gc.txt"; \
perl -F'\t' -anle 'BEGIN{$info={}; $info->{max}= 0; $info->{min}= 0} next if ($F[1] !~ /^\d+$/); ($con,$clen,$ccov,$cgc)=@F[0,1,6,3]; next if ($clen < 500 || $ccov < 0); $info->{count}++; push @{$info->{lens}},$clen; $info->{sumlen}+= $clen; $info->{sumgc}+= $clen*$cgc; $info->{sumcov}+= $clen*$ccov; END{ @sortlen= sort{$b<=>$a}@{$info->{lens}}; print "COUNT:\t$info->{count}"; print "SUM:\t$info->{sumlen}"; print "MIN:\t$sortlen[-1]"; print "MAX:\t@sortlen[0..9]"; print "AVE:\t".sprintf("%d",$info->{sumlen}/$info->{count}); foreach (@sortlen) {$info->{N50}+= $_; if ($info->{N50} > $info->{sumlen}/2) { $N50=$_; last; } } print "N50:\t$N50"; print "GC:\t".sprintf("%.4f",$info->{sumgc}/$info->{sumlen}); print "Cov:\t".sprintf("%.4f",$info->{sumcov}/$info->{sumlen});}' $i
##****** ****** ****** ****** ******

```

## Distribution of contig

```
$ R
> source("$HOME/biotools/local/statics/R/Library/frequency.R",encoding="utf-8")
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

## Gene (Protein) detection    Prodigal


```
$ cd ./genefind

ref="../assembly/xxx_spd_A1_assembly.fas"; \
name="xxx"; \
$HOME/biotools/local/genome/Prodigal/prodigal -i $ref -m -a ${name}.orfs.faa -d ${name}.orfs.fna -o ${name}.gbk -f gbk -p meta -q

```
##==> xxx.gbk
           xxx.orfs.faa
           xxx.orfs.fna
           
### CONVERT gbk TO GENELIST  ## FIX complete=1 partial=3

```
 name="xxx"; \
  echo -e "#Con\tlen\ttype\tID\tGENE\tdir\tST\tED\tseq_id\tcds_id\tCAT\tFIX\tEC\tNOTE\tAA" > ${name}_genelist; \
 perl -F'\t' -anle 'BEGIN{ $info={}; $no= 0; } if (/^DEFINITION  seqnum\=(\d+)\;seqlen\=(\d+)\;seqhdr\=\"(\S+?)\"\;/) { $ctg= $3; $cid= $1; $clen= $2; } elsif (/^\s+CDS\s+/) {  $no++;  if(/<?(\d+)\.\.>?(\d+)/) { ($st,$ed)=($1,$2); } $dir= /complement\(/ ? "-" : "+"; push @{$info->{$ctg}->{cds}},[$cid,$ctg,$clen,"CDS",$dir,$st,$ed]; } elsif (/note\=\"ID\=(.+?)\;/) { push @{$info->{$ctg}->{cds}->[-1]},$1; $comp= /partial=00/ ? 1 : 3; push @{$info->{$ctg}->{cds}->[-1]},$comp; }  END{ foreach $ctg (sort keys %{$info} ) { foreach (@{$info->{$ctg}->{cds}}) { $cds_no= $1 if ($$_[7] =~ /\d+_(\d+)/); $aalen=($$_[6]-$$_[5]-2)/3; print (join "\t",(@$_[1..3],"$$_[1]_$cds_no","",@$_[4,5,6],$$_[0],$$_[7],"A",$$_[8],"","",$aalen)); } }   }' ${name}.gbk >> ${name}_genelist

```

#==> xxx_genelist


##   Checkm

```
$ export PATH="(Your environment)/miniconda3/bin:$PATH"

#bioconda(link)python3 仮想環境に入れる
$ conda create -n checkm -c bioconda -y checkm-genome
  conda activate checkm

## GET DATABASE
checkm data setRoot checkm_db
cd checkm_db
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz

##  開始
conda activate checkm

## HMMER/prodigal/pplacer
$ export PATH=$PATH:$HOME/biotools/local/homology/hmmer-3.2.1/bin:$HOME/biotools/local/genome/Prodigal:$HOME/biotools/local/genome/pplacer-Linux-v1.1.alpha19

(checkm) $ checkm taxon_list   ## TAXON候補を表示

(checkm) $ checkm taxon_set genus TargetTax TargetTax.ms
(checkm) $ checkm analyze --genes -t 12 -x faa TargetTax.ms bins results
(checkm) $ checkm qa TargetTax.ms results

```


##   Comparem

仮想環境 comparem

```
$ export PATH="(Your environment)/miniconda3/bin:$PATH"

#bioconda(link)python3 仮想環境に入れる
$ conda create -n comparem

##  確認
$ conda info -e

##  開始
$ conda activate comparem
(comparem) $ pip install comparem

   OR
(comparem) $conda install -c bioconda comparem

##  comparem 環境開始
$ conda activate comparem
(comparem) $ export PATH="$HOME/biotools/local/genome/Prodigal:$HOME/biotools/local/homology/diamond-0.9.17:$PATH"

##  ゲノム毎のアミノ酸配列ファイルのリスト作成　（comp_genome_pep.lst）
###  comp_genome_pep.lst
./data/GCA_000016725.faa
./data/GCA_900156265.faa
./data/GCA_000202835.faa
./data/GCF_000701585.faa
./data/usoriyamakoNo5.faa

$ glist="comp_genome_pep.lst"

### QUERY_LIST TARGET_LIST OUTPUTS
  (comparem) $ output="aai_output3"; \
     comparem similarity -c 1 -x pep $glist $glist $output
     comparem aai -c 12 ./$output/query_genes.faa ./$output/hits_sorted.tsv $output

```

##==> ./aai_output1/aai_summary.tsv



## BLAST+


###   BLASTN
```
QUERY="RF00195_AY780887.1_24-141.fasta"; \
  DB="methR7_smrt_A1_pilon_IV.fa"; \


  $HOME/biotools/local/homology/ncbi-blast/bin/makeblastdb -in $DB -out $DB -dbtype nucl; \
  $HOME/biotools/local/homology/ncbi-blast/bin/blastn -query $QUERY -db $DB -out m_blastnout -evalue 1e-5 -dust no -num_threads 12; \
 exe_blastn_hit_plus_n.pl -r m_blastnout -s F -lh 70 -l 100 -n 1 -o hitlist_blastnout.txt; \
 awk '{FS="\t"}$1 ~ /^(A|b|s|[0-9]|#Query)/{print}' hitlist_blastnout.txt > hitlist_blastnoutS.txt

```

###  tBLASTn
```
  QUERY="lanmodulin.txt"; \
  DB="methR7IV.fasta"; \

  $HOME/biotools/local/homology/ncbi-blast/bin/makeblastdb -in $DB -out $DB -dbtype nucl; \
  $HOME/biotools/local/homology/ncbi-blast/bin/tblastn -query $QUERY -db $DB -out m_blastpout -evalue 1e-5 -seg no -num_threads 12; \
 exe_blastx_hit_plus.pl -r m_blastpout -s T -lh 70 -n 0 -o hitlist_blastxout.txt
 awk '{FS="\t"}$1 ~ /^(b|I|O|#Query)/{print}' hitlist_blastpout.txt > hitlist_blastpoutS.txt
 more hitlist_blastxout.txt

```

###  BLASTP
```
  QUERY="lanmodulin.txt";
 DB="methR7IV_cds.pep"; \

  $HOME/biotools/local/homology/ncbi-blast/bin/makeblastdb -in $DB -out $DB -dbtype prot; 
  $HOME/biotools/local/homology/ncbi-blast/bin/blastp -query $QUERY -db $DB -out m_blastpout -evalue 1e-5 -seg no -num_threads 12 -num_descriptions 100 -num_alignments 100; 
 exe_blastp_hit_plus.pl -r m_blastpout -s F -lh 20 -n 0 -o hitlist_blastpout.txt
 awk '{FS="\t"}$1 ~ /^(B|b|O|[0-9]|#Query)/{print}' hitlist_blastpout.txt > hitlist_blastpoutS.txt

```
