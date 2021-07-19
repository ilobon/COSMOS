![cosmos_logo_rounded](https://user-images.githubusercontent.com/17475181/125962073-12640e62-ffeb-499b-9a95-e86d60278681.png)

# COSMOS

COSMOS is an annotation and filtering strategy that accurately indentifies somatic SNVs from sequencing data.

## Inputs:

* VCF file with SNV calls
* VCF file with indels
* VCF file with copy number calls
* BAM files
* sample_ind file (tab sep file with sample and individual in VCF order)

## Recommended annotation BED files
Test files can be found in resources_test

* 1000G strict mask 
* Mappability <1 for read length
* dbsnp common variants
* Homopolymer track
* Segmental duplications track
* Capture BED (if not WGS)

## Dependencies

* HTSlib
* BCFtools
* Python >=3.6
#### Python modules
* Pysam
* natsort
* scipy


## Usage

COSMOS has been designed to detect somatic mutations in the absence of a "normal/germline" reference.
Sensitivity and specificity improve when multiple samples from the same individual are used.

We recommend calling variants with HaplotypeCaller --ploidy 10 
This allows for the detection of somatic mutations without a paired normal, but retrieves mostly noise.



#### The first step is to annotate with the BED file info and read features
`./annotateHCP10somatic.sh sample.vcf.gz sample.bedanno.vcf "anno1=anno1.bed.gz;anno2=anno2.bed.gz" "bam1.bam;bam2.bam;bam3.bam"`

#### Filtering by annotated flags and number of individuals
`python3.6 ./excludeRequireTags_otherInds.py --exclude 1000G,Mappability,SegDups,indel,Homopolymer,dbsnp --require onTarget --samplesInds sample_ind.txt --nOtherInds 1 -a VCF --includesSelf YES --vcfPath vcfPath/ --overwriteSamples YES --out sample.bedanno.flags.ninds.vcf.gz`

#### Filtering by read features
`python3.6 ./COSMOS.py -i sample.bedanno.flags.ninds.vcf.gz -o sample.bedanno.flags.ninds.cosmos.vcf.gz -is sample_ind.txt -c BOTH -ns 4 -ad 2 -adss 3 -vaf 0.5 -dp1 20 -dp2 100 -sr 2 -pr 4 -sp 0.05 -b 0.05 -nrl 4 -hap 4 -cnv NO -pir 4,4 -vafq 0.4 -clip 0.9 -mq 0.05 -mm 0.05`

