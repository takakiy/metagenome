# metagenome

## Trimmomatic

```
$ export TRIM=/home/impact/biotools/local/assembly/Trimmomatic-0.36

#**** **** **** **** **** **** 
$ mkdir cleanup; \
 for i in ./fastq/*_L001_R1_001.fastq.gz; do a=$i; \
 j=${i##./*/}; k=${j%%_S*_L001_R1_001.fastq.gz}; echo $k; \
 ## trimmomatic
 java -jar $TRIM/trimmomatic-0.36.jar PE -threads 24 -phred33 $i ${i/_R1_/_R2_} ./cleanup/${k}_R1_PE.fq ./cleanup/${k}_R1_OP.fq ./cleanup/${k}_R2_PE.fq ./cleanup/${k}_R2_OP.fq ILLUMINACLIP:$TRIM/adapters/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:100 HEADCROP:3 &> trimmomatic_log.txt; \
  done
#**** **** **** **** **** **** 


```


