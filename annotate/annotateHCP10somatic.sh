#Script to add fields relevant for somatic calling filtering 
#It adds fields to HaplotypeCaller --ploidy 10 vcfs

module unload gcc;module load gcc/6.3.0; module load BCFTOOLS/1.6

vcfIn=$1 
vcfOut=$2
beds=$3
bams=$4
xhmm=$5
pathOut=$(dirname $vcfOut)
name=$(echo $vcfOut | sed "s/annotated.vcf/annotated/g")
qu=$pathOut/qu
out=$pathOut/out
mkdir -p $qu $out

echo $vcfOut $name

#: <<'END'
#0. Previously: Run XHMM, get homopolymers bed

#1. Remove multiallelic positions and indels + indels to bed
#Remove multiallelic positions and indels from VCF
zcat $vcfIn | grep "^#" > ${name}.snps.0.vcf
zcat $vcfIn | grep -v "^#" | awk '{split($5,a,",");if(length($4)==1 && length(a[1])==1 && length(a[2])==0 && a[2]!="*")print $0}' >> ${name}.snps.0.vcf; bgzip ${name}.snps.0.vcf
#Create a bed with indels +-5bp
zcat $vcfIn | grep -v "^#" | awk '{OFS="\t";split($5,a,",");if(length($4)>1 || length(a[1])>1 || length(a[2])>1 || a[2]=="*")print $1,$2-6,$2+4,$1":"$2}' > ${name}.indels.bed; bgzip ${name}.indels.bed; tabix -p bed ${name}.indels.bed.gz

#2. Annotate position BEDs: 1000G strict, Segdups, Mappability, dbsnp, on target, indels, homopolymer
nbeds=$(echo $beds | awk '{split($1,a,";");print length(a)}')
for i in $(seq 1 1 $nbeds);do tag=$(echo $beds | cut -f$i -d";" | cut -f1 -d"=");bed=$(echo $beds | cut -f$i -d";" | cut -f2 -d"=");header=$(echo "##INFO=<ID=$tag,Number=1,Type=String,Description="Overlaps $tag bed">");bcftools annotate -a $bed -c CHROM,FROM,TO,${tag} -h <(echo "$header") ${name}.snps.$((i-1)).vcf.gz | bgzip > ${name}.snps.${i}.vcf.gz;done
bcftools annotate -a ${name}.indels.bed.gz -c CHROM,FROM,TO,indel -h <(echo "##INFO=<ID=indel,Number=1,Type=String,Description="Overlaps indel bed">") ${name}.snps.${nbeds}.vcf.gz | bgzip > ${name}.snps.beds.vcf.gz
rm ${name}.snps.[0-9]*

#3. Python features and pvalues + Annotate CNVs from XHMM
python3.6 ./annotateHCP10somatic.py ${name}.snps.beds.vcf.gz $bams $xhmm $vcfOut
bgzip $vcfOut
rm ${name}.indels.bed.gz* ${name}.snps.beds.vcf.gz 


#END
