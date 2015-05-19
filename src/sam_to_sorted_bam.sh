#! /bin/bash

### try to catch a tonne of errors (h/t Michael Hoffman)
set -o nounset -o pipefail -o errexit

### grab genomic locations:
rep_prefix=$1
rep1=$2
rep2=$3

### convert sam to sorted bam files, filtering the poor quality reads

sam_to_bam(){
	outfile=`echo $1 | sed -e 's/.sam/.bam/g'`
	samtools view -bSq $2 $1 > $outfile 
}

export -f sam_to_bam

echo "convert sam to bam?"
select convert in "y" "n";
do
	if [ "$convert" == "y" ]; then

		echo "input scaled phred probability quality cutoff value"
		read qval

		cd "$rep_prefix$rep1"
		echo "currently in "`pwd`
		echo "using cutoff score $qval"
		parallel -j10 --dry-run --progress --xapply sam_to_bam ::: $qval ::: `ls -1 *.sam`
		

		cd "$rep_prefix$rep2"
		echo "currently in "`pwd`
		parallel -j10 --dry-run --progress --xapply sam_to_bam ::: $qval ::: `ls -1 *.sam`
	fi
	break;
done

echo "delete sam files?"
select rmsam in "y" "n";
do
	if [ "$rmsam" == "y" ]; then
		cd "$rep_prefix$rep1"
		rm *.sam

		cd "$rep_prefix$rep2"
		rm *.sam
	fi
	break;
done	

sort_index_bam(){
	outfile=`echo $1 | sed -e 's/.bam/_sorted.bam/g'`
	samtools sort $1 $outfile
	samtools index $outfile
}

echo "Sort and index bam files?"
select sortindex in "y" "n";
do
	if [ "$sortindex" == "y" ]; then
		cd "$rep_prefix$rep1"
		parallel -j10 --dry-run --progress --xapply sort_index_bam ::: `ls -1 *.bam`

		cd "$rep_prefix$rep2"
		parallel -j10 --dry-run --progress --xapply sort_index_bam ::: `ls -1 *.bam`

	fi
done

echo "Remove non-indexed bam files?"
select rmbam in "y" "n";
do
	if [ "$rmbam" == "y" ]; then
		cd "$rep_prefix$rep1"
		rm *[0-9].bam

		cd "$rep_prefix$rep2"
		rm *[0-9].bam
	fi
done	

