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



