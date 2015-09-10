#/bin/bash
exec 6>HenanBreast.log.txt
BASH_XTRACEFD=6
set -x
get()
{
        out=$(awk -F "=" '$1~"'$1'"{print $2}' $2)
		#0.2 ~替换==
		#0.2 允许对配置文件进行注释
        echo $out
}

report()
{
#threads step info home
#1.0[main] report Wed_Jul_15_04:37:21_EDT_2015 [类] [信息]
	echo -e "\033[36m[$1]\treport\t`date | sed -n 's/\s/_/gp'`\t$2\t$3\033[0m" | tee -a $4/run.log
}
export -f report
check()
{
	awk '$1~"'$1'" && $2=="passed" && $4=="'$2'" && $5~"finished"{print}' $3/run.log | wc -l
}
export -f check
verdict()
# threads step home_dir
{
	#1.0[main] passed/error Wed_Jul_15_04:37:21_EDT_2015 [类] finished/failed
	if [ $? -eq 0 ]
	then
		echo -e "\033[36m[$1]\tpassed\t`date | sed -n 's/\s/_/gp'`\t$2\tfinished\033[0m" | tee -a $3/run.log
	else
		echo -e "\033[31m[$1]\terror\t`date | sed -n 's/\s/_/gp'`\t$2\tfailed\033[0m" | tee -a $3/run.log
		exit
	fi
}
export -f verdict
#check arg.
recover=0 #1.0 normal model default
config=$1
case $1 in
        *"config.ini")
		Ehome=`echo $(pwd)/$1 | sed -n 's/\/config\.ini//p'`
		echo > $Ehome/run.log;;
        "-R")
                if [ `echo $2 | grep "config.ini" |wc -l` = 0 ];then
					echo -e "\033[36mNo configure file\nUsage:$0 <config.ini> or $0 -R <config.ini>\033[0m"
					exit 0
					else Ehome=`echo $(pwd)/$2 | sed -n 's/\/config\.ini//p'`
					if [ -f $Ehome/run.log ];then
						$(report main log_file log_file_finded,recovering_previous_work $Ehome)
						recover=1
						config=$2
					else 
						echo -e "\033[31m no log file \033[0m\n"
					fi
                fi;;
        *)echo -e "\033[36mNo configure file\nUsage:$0 <config.ini> or $0 -R <config.ini>\033[0m";exit 0;;
esac

export Ehome
export recover


JAVA_PATH=$(get Java $config)
NGSQCToolkit_PATH=$(get NGSQCToolkit $config)
BOWTIE_PATH=$(get Botiew2 $config)
SAMTOOLS_PATH=$(get Samtools $config)
TOPHAT_PATH=$(get Tophat $config)
CUFF_PATH=$(get Cufflinks $config)
FastaQC_PATH=$(get FastaQC $config)
export PATH=$JAVA_PATH:$NGSQCToolkit_PATH:$BOWTIE_PATH:$SAMTOOLS_PATH:$TOPHAT_PATH:$CUFF_PATH:$FastaQC_PATH:$PATH
PICARD_PATH=$(get Picard $config)
export PICARD_PATH
model=$(get Model $config)
export model
filter_percent=$(get filter_percent $config)
export filter_percent
annotation=$(get annotation $config)
export annotation
genome=$(get genome $config)
export genome
threads_num=$(get usable_threads_number $config)
export threads_num
declare -a Samples
samples_dir=$(get Samples_Dir $config)
Samples=(`ls $samples_dir`)
#optional
adapter=$(get Adapter_seqs $config)
export adapter
insert_size=$(get insert_szie $config)
export insert_size
#1.0增加reads长度配置
reads_length=$(get reads_length $config)
export reads_length

#required check
if [[ $PICARD_PATH && $model && $filter_percent && $annotation && $reads_length ]];then
echo -e "\033[34m..OK..\nStarting Pipeline...\033[0m"
else
	echo -e "\033[33mconfigure file error:
	detail: 
		PICARD_PATH $PICARD_PATH
		model $model
		filter_percent $filter_percent
		annotation $annotation
		reads_length $reads_length
	\033[0m"
	exit 0
fi
#0.3 check index dir
# 自动识别建立index
bt_genome=`echo -n $genome | sed -nr 's/\.\w+$//p'`
if [ `ls -l $bt_genome.*.bt2 | wc -l ` != 6 ]
then
	echo $(report main building_genome_bt starting $Ehome) 
	bowtie2-build $genome $bt_genome 2>&1 | tee -a $Ehome/run.log
	echo $(vardict main building_genome_bt $Ehome)
fi
export bt_genome
#0.3 check transcriptome file in annotation
#0.3 命名规则：annotation_trans.fa
#0.3 自动识别建立index
bt_trans=`echo -n $annotation|sed -nr 's/\.\w+$/.trans/p'`
if [ `ls -l $bt_trans.*.bt2 | wc -l` != 6 ]
then
#0.3 建立转录组，构建index
	echo $(report main building_trans_bt starting $Ehome)
	tophat -G $annotation --transcriptome-index=$bt_trans $bt_genome 2>&1 | tee -a $Ehome/run.log
	echo $(vardict main building_trans_bt $Ehome)
fi
export bt_trans
if [ ! -f `echo $0 | sed -nr 's/main.sh/assembly.sh/p'` ]
then
	echo -e "\033[36mno assembly.sh error
	main path:$0\033[0m"
	exit 0
fi
#创建命令
for i in ${Samples[@]}; do
	if [ "$recover" = "0" ];then
		if [ ! -f $samples_dir/$i/"$i"_1.fq ];then mv $samples_dir/$i/*_1.fq $samples_dir/$i/"$i"_1.fq;fi;
		if [ ! -f $samples_dir/$i/"$i"_2.fq ];then mv $samples_dir/$i/*_2.fq $samples_dir/$i/"$i"_2.fq;fi;
	fi
	cmdlines=bash^`echo $0 | sed -n 's/main.sh/assembly.sh/p'`"^$samples_dir/$i^$i $cmdlines"
	#0.2 `echo $0 | sed -n 's/main.sh/assembly.sh/p'`" 替换 `pwd`"/assembly.sh
	#0.2 从主程序中获取子程序位置，实现程序与数据分离
done
#并行
multi $cmdlines
if [ ! -d exp ];then mkdir exp;fi
cd exp
if [ -f merge.list ]
then
	rm merge.list
fi

for i in ${Samples[@]}; do
	while :
	do
		if [ -f $samples_dir/$i/single_sample_completed ];then break;fi;
		sleep 300s;
		echo $(report main check waiting $Ehome)
	done
	awk '$3=="CDS" || $3=="exon"{print}' $samples_dir/$i/06_Cufflinks/transcripts.gtf > $i.gtf
	echo $i.gtf >> merge.list
	if [ -f $samples_dir/$i/EG ];then
		eglab="$eglab,$i"
		egpath="$egpath,$samples_dir/$i/05_TopHat/accepted_hits_dupmaked.bam"
	elif [ -f $samples_dir/$i/CG ];then
		cglab="$cglab,$i"
		cgpath="$cgpath,$samples_dir/$i/05_TopHat/accepted_hits_dupmaked.bam"
	else
		echo "no group info."
		exit 1
	fi
done
eglab=`echo $eglab |sed -nr 's/^,(.+)/\1/p'`
egpath=`echo $egpath |sed -nr 's/^,(.+)/\1/p'`
cgpath=`echo $cgpath |sed -nr 's/^,(.+)/\1/p'`
echo $(report main exp_diff_info $cglab,$eglab,$egpath,$cgpath $Ehome)
if [[ "$recover" = "1" && $(check main cuffmerge $Ehome) > 0 ]];then
echo -e "\033[34m..Skip cuffmerge..\033[0m"
else
	echo $(report main cuffmerge starting $Ehome)
	cuffmerge -o merged -p $threads_num -g $annotation -s $genome merge.list  2>&1 | tee -a $Ehome/run.log
	echo $(verdict main cuffmerge $Ehome)
fi
echo $(report main cuffdiff starting $Ehome)
if [ "$model" = "fast" ];then
	cuffdiff -labels $eglab$cglab -p $threads_num -output-dir cuffdiff merged/merged.gtf $egpath $cgpath  2>&1 | tee -a $Ehome/run.log
else
	cuffdiff -labels $eglab$cglab -multi-read-correct -p $threads_num -output-dir cuffdiff -frag-bias-correct $genome --total-hits-norm -quiet merged/merged.gtf $egpath $cgpath  2>&1 | tee -a $Ehome/run.log
fi
awk '{a[$14]++}END{for(i in a){print i,a[i]}}' $Ehome/exp/cuffdiff/gene_exp.diff > $Ehome/exp/exp_check
echo $(verdict main cuffdiff $Ehome)
