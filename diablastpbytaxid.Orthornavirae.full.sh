#!/bin/bash

DIAMOND="/hwfssz5/ST_INFECTION/GlobalDatabase/user/fengqikai/software/Diamond/WHICHDIAMOD"
DIAMOND_DB="/jdfssz1/ST_HEALTH/P20Z10200N0206/fengqikai/NCBI-BLAST/2021-05/NR/Diamond/WHICHNRDMND"
INDIR="/hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/new_taxdump"
OUTDIR="/hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/WHICHOUTPATH"
familylist="/hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/new_taxdump/WHICHFAMLIST"
Infile="fromSM.rdrp20201015.fasta"
OUTFAMDIR="WHICHDIRTOSTOREFAMS"

echo -e "Starting working on $Infile at $(date)"

$DIAMOND  blastp --db  $DIAMOND_DB -q $INDIR/$Infile -f 6 qseqid sseqid qlen slen length qcovhsp pident evalue bitscore mismatch staxids sscinames salltitles full_sseq  --threads 16 --very-sensitive -k 0  -e 0.00001 --taxonlist 2732396 -o $OUTDIR/$Infile.blastpout

echo -e "blastp completed at $(date)"

cat $OUTDIR/$Infile.blastpout | tr " " "+" | cut -f2,12,11,14 | sort -u -T /ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/liqian6 |while read f1 f2 f3 f4; do echo -e ">${f1// /}\t|\t${f2// /}\t|\t${f3// /}\n${f4// /}"; done | tr "+" " " > $OUTDIR/$Infile.sigblastp.faa

echo -e "extraction from blastp results completed at $(date)"

echo -e "Starting to compare CDDs against selected sequences with psi-blast at $(date)"

/hwfssz1/ST_META/AP/zhaohailong/ncbi-blast-2.10.1+/bin/psiblast -query  $OUTDIR/$Infile.sigblastp.faa -outfmt 6 -evalue 1E-5 -out $OUTDIR/$Infile.cdds.psiblast -db /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/CDD/cdd.smp.pn/all.html.accession  -num_threads 16

echo -e "Starting to compare CDDs against selected sequences with rps-blast at $(date)"
$
/hwfssz1/ST_META/AP/zhaohailong/ncbi-blast-2.10.1+/bin/rpsblast -query  $OUTDIR/$Infile.sigblastp.faa -outfmt 6 -evalue 1E-5 -out $OUTDIR/$Infile.cdds.rpsblast -db /hwfssz5/ST_INFECTION/GlobalDatabase/share/public_database/rdrp_AAseq/CDD/cdd.smp.pn/all.html.accession  -num_threads 16

cat $OUTDIR/$Infile.cdds.psiblast $OUTDIR/$Infile.cdds.rpsblast | cut -f1 | sort | uniq > $OUTDIR/$Infile.selected.list

echo -e "start to extract sequence from cdd blast results at $(date)"

cat $OUTDIR/$Infile.selected.list | while read acc; do grep -w "$acc" -A1  $OUTDIR/$Infile.sigblastp.faa >> $OUTDIR/$Infile.allsig.faa ; done

grep ">" $OUTDIR/$Infile.allsig.faa | sed -e 's/>//' > $OUTDIR/$Infile.allsig.list

echo -e "start to separate taxid at $(date)"

cat $OUTDIR/$Infile.allsig.list | grep ";" | head | while read f1 f2 f3 f4 f5; do echo -e "$f3" > $OUTDIR/tmpf1 ; echo -e "$f5" > $OUTDIR/tmpf2; cat $OUTDIR/tmpf1 | tr ";" "\n" > $OUTDIR/tmpf11; cat $OUTDIR/tmpf2 | tr ";" "\n" > $OUTDIR/tmpf22; paste -d "*" $OUTDIR/tmpf11  $OUTDIR/tmpf22 | while read line; do echo -e "$f1\t|\t$line"; done; done | sed -e "s/*/\t|\t/" > $OUTDIR/$Infile.allsig.list.alltaxa 

echo -e "finishing at $(date)"

cat $familylist | while read fam; do cp /hwfssz5/ST_INFECTION/GlobalDatabase/user/liqian6/new_taxdump/rmrdseqbyfam.sh $OUTDIR/rmrdseqbyfam.${fam// /}.tmp; sed -i "s/TESTFAM/$(echo ${fam}  | tr "/" "@")/" $OUTDIR/rmrdseqbyfam.${fam// /}.tmp ; sed -i "s/TESTWORKDIR/$(echo ${OUTFAMDIR} | tr "/" "@")/" $OUTDIR/rmrdseqbyfam.${fam// /}.tmp; sed -i "s/TESTINPATH/$(echo ${OUTDIR}  | tr "/" "@")/" $OUTDIR/rmrdseqbyfam.${fam// /}.tmp; sed -i "s/INPUTFASTA/$(echo ${Infile} | tr "/" "@")/" $OUTDIR/rmrdseqbyfam.${fam// /}.tmp ; cat $OUTDIR/rmrdseqbyfam.${fam// /}.tmp | tr "@" "/" > $OUTDIR/rmrdseqbyfam.${fam// /}.sh; rm $OUTDIR/rmrdseqbyfam.${fam// /}.tmp ; echo -e "script for $fam is prepared, please run it in $OUTDIR" ; done
 
