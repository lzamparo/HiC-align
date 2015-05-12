#! /bin/bash

### try to catch a tonne of errors (h/t Michael Hoffman)
set -o nounset -o pipefail -o errexit

### grab genomic locations:
hg=$1
rep_prefix=$2
rep1=$3
rep2=$4

echo "starting pipeline..."
echo "build index?"
select build_index in "y" "n";
do
if [ "$build_index" == 'n' ]; then
	break;
elif [ "$build_index" == 'y' ]
then
	### Build index of HG (in $1)
	echo $"trying to find hg in $hg"$
	if [ ! -d "$hg" ]; then
		echo "did not find $hg"
		exit;
	fi

	### remove *random.fa, *hap?.fa, chrUn*.fa, chrM.fa ?
	if [ -d "$hg" ]; then
  		cd "$hg"
  		echo "remove *random.fa, *hap?.fa, chrUn*.fa, chrM.fa ?"
  		select yn in "y", "n";
		do
			case $yn in
				y ) rm *random.fa *hap?.fa chrUn*.fa chrM.fa;;
				n ) echo "ok then" ;;
			esac
		done
	fi

	### create whole genome fasta file
	echo "creating whole genome fasta file"
	cat `ls -1v` > genome.fa

fi
done

### is there an index file already?
cd "$hg"
if [ ! -f "genome.fa.pac" ]; then
  	echo "did not detect an index, generating one..."
  	bwa index -a bwtsw genome.fa
else
	echo "detected an index in $hg..."
fi

### fastqc results will appear here
fastqc_results="/home/zamparol/data/fastqc_results/"
outdir_rep1="$fastqc_results$rep1"
outdir_rep2="$fastqc_results$rep2"

echo "run fastqc tests?"
select do_fastqc in "y" "n";
do
	if [ "$do_fastqc" == "y" ]; then
		# do fastqc, check for red flags (failed tests) in both reps
		cd "$rep_prefix$rep1"
		echo "processing $rep1..."
		ls -1 *.gz | parallel --progress fastqc -q --outdir=/home/zamparol/data/fastqc_results/Sample_HiC_LY1_1
		cd "$rep_prefix$rep2"
		echo "processing $rep2..."
		ls -1 *.gz | parallel --progress fastqc -q --outdir=/home/zamparol/data/fastqc_results/Sample_HiC_LY1_2

	fi
	break;
done

### check test output for failed tests
check_fail(){
	unzip -q $1
	filename=$(basename "$1")
	prefix="${filename%.*}"
	fails=`cat $prefix/summary.txt | grep -c FAIL`
	echo $fails
}

cleanup_dirs(){
	filename=$(basename "$1")
        prefix="${filename%.*}"
	if [ -d "$prefix" ]; then
		rm -rf "$prefix"
	fi
}

export -f check_fail
export -f cleanup_dirs

echo "analyze fastqc results?"
select analyze in "y" "n";
do
	if [ "$analyze" == "n" ]; then
		break;
	elif [ "$analyze" == "y" ]
	then
		cd "$outdir_rep1"
		num_files=`ls -1 *.zip | wc -l`
		echo "found $num_files files in rep1..."
		echo "unzipping, counting failures"
		fails=`ls -1 *.zip | parallel --progress check_fail | paste -sd+ | bc`
		echo "found $fails failures in rep 1"

		cd "$outdir_rep2"
		num_files=`ls -1 *.zip | wc -l`
		echo "found $num_files files in rep2..."
		echo "unzipping, counting failures"
		fails=`ls -1 *.zip | parallel --progress check_fail | paste -sd+ | bc`
		echo "found $fails failures in rep 2"
		break;
	fi
done

echo "clean up directories?"
select cleanup in "y" "n";
do
	if [ "$cleanup" == "y" ]; then
		cd $outdir_rep1
		ls -1 *.zip | parallel cleanup_dirs
		cd $outdir_rep2
		ls -1 *.zip | parallel cleanup_dirs
	fi
	break;
done
