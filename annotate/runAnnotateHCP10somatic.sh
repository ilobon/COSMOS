#!/bin/bash

pathIn=/scratch/devel/irenel/parkinson/haplotypeCallerP10
pathOut=/scratch/devel/irenel/parkinson/haplotypeCallerP10annotated
inBed=/home/devel/irenel/scratch/parkinson/somatic
qu=${pathOut}/qu
out=${pathOut}/out
mkdir -p $qu $out

xhmm=/home/devel/irenel/scratch/parkinson/XHMM/DV.HCorder.vcf.gz

beds="1000Gstrict=${inBed}/hs37d5.1000G.NonPassStrictMask.taggedN.bed.gz;mappability=${inBed}/hg19.100mers.non1.taggedN.ordered.bed.gz;segdups=${inBed}/hs37d5.segdups.taggedN.bed.gz;dbsnp=${inBed}/dbsnp_common.taggedN.bed.gz;homopolymer=${inBed}/strs/reference_set_minScore5_all.homopol.nuc.bed.gz;onTarget=/home/devel/irenel/scratch/parkinson/SureSelect/Target.hs37d5.bed.gz"

bams=""
for bam in $(ls /home/devel/irenel/scratch/parkinson/bwaMapping/exomeBams/*bam)
do
	if [ "$bams" == "" ]
	then
		bams=${bam}
	else
		bams=${bams}";${bam}"
	fi
done

#for chr in $(seq 1 1 21)
for chr in 3.missed
do 
	jobName=${qu}/annotateHCP10somatic.${chr}.sh
	vcfIn=${pathIn}/All_HP_P10.chr${chr}.vcf.gz
	vcfOut=${pathOut}/All_HP_P10.chr${chr}.annotated.vcf
	echo "#!/bin/bash" > ${jobName}
	echo "./annotateHCP10somatic.sh $vcfIn $vcfOut \"$beds\" \"$bams\" $xhmm" >> ${jobName}
	chmod +x $jobName
	~/bin/submit.py -c ${jobName} -n ${chr}annotate -o ${out}/annotateHCP10somatic.${chr}.out -e ${out}/annotateHCP10somatic.${chr}.err -p main,genB -w "35:00:00" -r lowprio
done

