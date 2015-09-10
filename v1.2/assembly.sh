#/bin/bash
exec 6>HenanBreast.log.txt
BASH_XTRACEFD=6
set -x
cd $1
echo "samples $2 processing..."
read_1="$2"_1.fq
read_2="$2"_2.fq
#QC_before
if [[ $recover -eq 1 && $(check $2 QC_before $Ehome) -gt 0 ]];then
	echo -e "\033[34m..Skip QC_before..\033[0m"
else
	if [ ! -d 01_QC_before_filter ];then mkdir 01_QC_before_filter;fi
	echo $(report $2 QC_before starting $Ehome)
	fastqc --outdir 01_QC_before_filter --threads $threads_num $read_1 $read_2 2>&1 | tee -a $Ehome/run.log 
	echo $(verdict $2 QC_before $Ehome)
fi
#filter#
if [[ "$recover" = "1" && $(check $2 filter $Ehome) > 0 ]];then
echo -e "\033[34m..Skip filter..\033[0m"
else
	if [ ! -d 02_Filter ];then mkdir 02_Filter;fi
	echo $(report $2 filter starting $Ehome)
	if [ -z $adapter ];then
		IlluQC_PRLL.pl -pe $read_1 $read_2 N A -cutOffReadLen4HQ $filter_percent -cpus $threads_num -outputFolder 02_Filter  | tee -a $Ehome/run.log
		echo $(verdict $2 filter $Ehome)
	else
		echo $(report $2 filter adapter_find $Ehome)
		echo $adapter | sed -nr 's/,/\n/gp' > 02_Filter/adapter.seq
		#0.2 增加echo $adapter | sed -nr 's/,/\n/gp' > 02_Filter/adapter.seq
		IlluQC_PRLL.pl -pe $read_1 $read_2 02_Filter/adapter.seq A -cutOffReadLen4HQ $filter_percent -cpus $threads_num -outputFolder 02_Filter  | tee -a $Ehome/run.log
		#0.2 02_Filter/adapter.seq 替换 $adapter 
		echo $(verdict $2 filter $Ehome)
	fi
	mv 02_Filter/`echo $read_1`_filtered 02_Filter/$read_1
	mv 02_Filter/`echo $read_2`_filtered 02_Filter/$read_2
	mv 02_Filter/`echo $read_1`_`echo $read_2`_unPaired_HQReads 02_Filter/unpaired.fq
fi
read_1=02_Filter/$read_1
read_2=02_Filter/$read_2
unpaired=02_Filter/unpaired.fq
#QC_after
if [[ "$recover" = "1" && $(check $2 QC_after $Ehome) > 0 ]];then
echo -e "\033[34m..Skip QC_after..\033[0m"
else
	echo $(report $2 QC_after starting $Ehome)
	if [ ! -d 03_QC_after_filter ];then mkdir 03_QC_after_filter;fi
	fastqc --outdir 03_QC_after_filter --threads $threads_num $read_1 $read_2  | tee -a $Ehome/run.log
	echo $(verdict $2 QC_after $Ehome)
fi
#insert-size
# if [ -z $insert_size ];then
	# if [[ "$recover" = "1" && $(check $2 insert-size $Ehome) > 0 ]];then
	# echo -e "\033[34m..Skip insert-size..\033[0m"
	# else
		# if [ ! -d 04_Insert_size ];then mkdir 04_Insert_size;fi
		# echo $(report $2 insert-size starting $Ehome)
		# echo $(report $2 insert-size starting_bowtie $Ehome)
		# bowtie2 -q -I 0 -X 500 -x $bt_genome -p $threads_num -1 $read_1 -2 $read_2 -S 04_Insert_size/insert.sam  2>&1 | tee -a $Ehome/run.log
		# echo $(report $2 insert-size sam_to_bam_ing $Ehome)
		# samtools view -bS 04_Insert_size/insert.sam > 04_Insert_size/insert.bam
		# echo $(report $2 insert-size sort_sam_ing $Ehome)
		# multi java^-jar^$PICARD_PATH/SortSam.jar^I=04_Insert_size/insert.bam^O=04_Insert_size/insert.sort.bam^SO=coordinate rm^04_Insert_size/insert.sam
		# echo $(report $2 insert-size Collect_ing $Ehome)
		# java -jar $PICARD_PATH/CollectInsertSizeMetrics.jar I=04_Insert_size/insert.sort.bam O=04_Insert_size/lib_length H=04_Insert_size/graph.pdf
		# echo $(verdict $2 insert-size $Ehome)
		# rm 04_Insert_size/insert.bam	
	# fi
	# lib_mean_length=`sed -nr '/^([0-9.]+\t){7}FR/p' 04_Insert_size/lib_length | awk '{print $5}' | sed -nr 's/([0-9]+)\.[0-9]+/\1/p'`
	# insert_size=$[$lib_mean_length-2*$reads_length]
# fi

#tophat
if [[ "$recover" = "1" && $(check $2 tophat $Ehome) > 0 ]];then
	echo -e "\033[34m..tophat..\033[0m"
else
	echo $(report $2 tophat starting $Ehome)
	if [ "$model" = "fast" ]
		then
		tophat --read-mismatches 1  --output-dir 05_TopHat --read-realign-edit-dist 0 --num-threads $threads_num --transcriptome-index $bt_trans $bt_genome $read_1 $read_2  | tee -a $Ehome/run.log
	elif [ "$model" = "medium" ]
		then
		tophat --read-mismatches 1  --output-dir 05_TopHat --read-realign-edit-dist 0 --num-threads $threads_num --prefilter-multihits --transcriptome-index $bt_trans $bt_genome $read_1 $read_2  | tee -a $Ehome/run.log
	else
		tophat --read-mismatches 1 --output-dir 05_TopHat --read-realign-edit-dist 0 --num-threads $threads_num --prefilter-multihits --coverage-search --transcriptome-index $bt_trans $bt_genome $read_1 $read_2  | tee -a $Ehome/run.log
	fi
	echo $(verdict $2 tophat $Ehome)
fi
#mask duplicates
if [[ "$recover" = "1" && $(check $2 mask_duplicates $Ehome) > 0 ]];then
echo -e "\033[34m..mask_duplicates..\033[0m"
else
	echo $(report $2 mask_duplicates starting $Ehome)
	java -jar $PICARD_PATH/MarkDuplicates.jar INPUT=05_TopHat/accepted_hits.bam OUTPUT=05_TopHat/accepted_hits_dupmaked.bam METRICS_FILE=05_TopHat/dup_metrics_file REMOVE_DUPLICATES=false
	echo $(verdict $2 mask_duplicates $Ehome)
fi
#assembly by cufflinks
if [[ "$recover" = "1" && $(check $2 cufflinks $Ehome) > 0 ]];then
echo -e "\033[34m..cufflinks..\033[0m"
else
	echo $(report $2 cufflinks starting $Ehome)
	cufflinks --output-dir 06_Cufflinks --num-threads $threads_num --GTF-guide $annotation 05_TopHat/accepted_hits_dupmaked.bam
	echo $(verdict $2 cufflinks $Ehome)
fi
echo > single_sample_completed
echo $(report $2 assembly completed $Ehome)