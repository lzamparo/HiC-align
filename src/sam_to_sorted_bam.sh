#! /bin/bash

### try to catch a tonne of errors (h/t Michael Hoffman)
set -o nounset -o pipefail -o errexit

### grab genomic locations:
rep_prefix=$1
rep1=$2
rep2=$3

### convert sam to sorted bam files, filtering the poor quality reads

sam_to_bam(){
	outfile=`echo $2 | sed -e 's/.sam/.bam/g'`
	samtools view -bSq $1 $2 > $outfile 
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
		parallel -j10 --progress --xapply sam_to_bam ::: $qval ::: `ls -1 *.sam`
		

		cd "$rep_prefix$rep2"
		echo "currently in "`pwd`
		parallel -j10 --progress --xapply sam_to_bam ::: $qval ::: `ls -1 *.sam`
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

sort_bam(){
	outfile=`echo $1 | sed -e 's/.bam/_sorted/g'`
	samtools sort $1 $outfile
}
export -f sort_bam

echo "Sort bam files?"
select dosort in "y" "n";
do
	if [ "$dosort" == "y" ]; then
		cd "$rep_prefix$rep1"
		ls -1 *.bam | parallel -j20 --progress sort_bam 

		cd "$rep_prefix$rep2"
		ls -1 *.bam | parallel -j20 --progress sort_bam

	fi
	break;
done

echo "Remove unsorted bam files?"
select rmbam in "y" "n";
do
	if [ "$rmbam" == "y" ]; then
		cd "$rep_prefix$rep1"
		rm *[0-9].bam

		cd "$rep_prefix$rep2"
		rm *[0-9].bam
	fi
	break;
done	


echo "Index sorted bam files?"
select doindex in "y" "n";
do
	if [ "$doindex" == "y" ]; then
		cd "$rep_prefix$rep1"
		ls -1 *_sorted.bam | parallel -j20 --progress samtools index 

		cd "$rep_prefix$rep2"
		ls -1 *_sorted.bam | parallel -j20 --progress samtools index 

	fi
	break;
done

