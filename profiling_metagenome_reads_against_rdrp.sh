#!/bin/bash

INPATH="/PATH/TO/FQ/FILES"
OUTDIR="/PATH/TO/OUTPUTS"
REFORMATTER="/bbmap/reformat.sh"
DIAMOND="/Diamond/2.0.9/diamond"
DIAMOND_DB="/PATH/TO/rdrp_Orthornavirae.clustered.dmnd"
BEDTOOLS="/bedtools2/bin/bedtools"
PROTLIST="/PATH/TO/rdrp_Orthornavirae.clustered.length.tab"

ACCESSIONLIST="/PATH/TO/samples.list"

cat $ACCESSIONLIST | while read ACCESSION;
do

echo -e "start to process fastq reads from $ACCESSION at $(date)" 
mkdir $OUTDIR/$ACCESSION

OUTPATH="$OUTDIR/$ACCESSION"

echo -e "start to merge paired end reads of $ACCESSION at $(date)"

$REFORMATTER in1=${ACCESSION// /}_rmrRNA_fwd.fq.gz  in2=${ACCESSION// /}_rmrRNA_fwd.fq.gz  out=$OUTPATH/${ACCESSION// /}.merged.fq.gz

echo -e "fninishing merging paired end reads of $ACCESSION at $(date), results stored as $OUTPATH/$ACCESSION.merged.fq.gz"

echo -e "start to align merged reads of $ACCESSION at $(date)"

$DIAMOND blastx  -d  $DIAMOND_DB  -q $OUTPATH/$ACCESSION.merged.fq.gz -o $OUTPATH/$ACCESSION.merged.tab -f 6 --threads 8 

echo -e "finishing alignment of $OUTPATH/$ACCESSION.merged.fq, start to genrate coverage at $(date)"

cat $OUTPATH/$ACCESSION.merged.tab | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n  > $OUTPATH/$ACCESSION.merged.bed 

cat $OUTPATH/$ACCESSION.merged.tab | awk '$12 >= 50'  | awk '{if($9 < $10){print($2"\t"$9-1"\t"$10)}else{print($2"\t"$10-1"\t"$9)}}' |  sort -k1,1 -k2,2n > $OUTPATH/$ACCESSION.filtered.merged.bed

$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.merged.bed -g $PROTLIST > $OUTPATH/$ACCESSION.merged.cov.txt

$BEDTOOLS  genomecov -i $OUTPATH/$ACCESSION.filtered.merged.bed -g $PROTLIST > $OUTPATH/$ACCESSION.filtered.merged.cov.txt 

rm $OUTPATH/$ACCESSION.merged.fq

echo -e "finishing alignment on $ACCESSION at $(date)"

done

echo -e "job completed at $(date)"
