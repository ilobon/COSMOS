![cosmos_logo_rounded](https://user-images.githubusercontent.com/17475181/125962073-12640e62-ffeb-499b-9a95-e86d60278681.png)

# COSMOS

COSMOS is an annotation and filtering strategy that accurately indentifies somatic SNVs from sequencing data.

## Inputs:

* VCF file with SNV calls
* VCF file with indels
* VCF file with copy number calls
* BAM files

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


./annotateHCP10somatic.sh sample.vcf.gz sample.bedanno.vcf "1000Gstrict=resources_fragment/hs37d5.1000GNonPassStrictMask.tagged.fragment.bed.gz;mappability=resources_fragment/hs37d5.100mersMappability.tagged.fragment.bed.gz;segdups=resources_fragment/hs37d5.segdups.tagged.fragment.bed.gz;dbsnp=resources_fragment/hs37d5.dbsnpCommon.tagged.fragment.bed.gz;homopolymer=resources_fragment/hs37d5.HomopolminScore5.fragment.bed.gz;onTarget=resources_fragment/exomeSureSelect.hs37d5.fragment.bed.gz" "data/bams/DV10B.fragment.bam;data/bams/DV10C.fragment.bam;data/bams/DV10E.fragment.bam;data/bams/DV10N.fragment.bam;data/bams/DV10S.fragment.bam;data/bams/DV1B.fragment.bam;data/bams/DV1C.fragment.bam;data/bams/DV1E.fragment.bam;data/bams/DV1N.fragment.bam;data/bams/DV1S.fragment.bam;data/bams/DV2B.fragment.bam;data/bams/DV2C.fragment.bam;data/bams/DV2E.fragment.bam;data/bams/DV2N.fragment.bam;data/bams/DV2S.fragment.bam;data/bams/DV3B.fragment.bam;data/bams/DV3C.fragment.bam;data/bams/DV3E.fragment.bam;data/bams/DV3N.fragment.bam;data/bams/DV3S.fragment.bam;data/bams/DV4B.fragment.bam;data/bams/DV4C.fragment.bam;data/bams/DV4E.fragment.bam;data/bams/DV4N.fragment.bam;data/bams/DV4S.fragment.bam;data/bams/DV5B.fragment.bam;data/bams/DV5C.fragment.bam;data/bams/DV5E.fragment.bam;data/bams/DV5N.fragment.bam;data/bams/DV5S.fragment.bam;data/bams/DV6B.fragment.bam;data/bams/DV6C.fragment.bam;data/bams/DV6E.fragment.bam;data/bams/DV6N.fragment.bam;data/bams/DV6S.fragment.bam;data/bams/DV7B.fragment.bam;data/bams/DV7C.fragment.bam;data/bams/DV7E.fragment.bam;data/bams/DV7N.fragment.bam;data/bams/DV7S.fragment.bam;data/bams/DV8B.fragment.bam;data/bams/DV8C.fragment.bam;data/bams/DV8E.fragment.bam;data/bams/DV8N.fragment.bam;data/bams/DV8S.fragment.bam;data/bams/DV9B.fragment.bam;data/bams/DV9C.fragment.bam;data/bams/DV9E.fragment.bam;data/bams/DV9N.fragment.bam;data/bams/DV9S.fragment.bam" data/XHMM.vcf.gz
