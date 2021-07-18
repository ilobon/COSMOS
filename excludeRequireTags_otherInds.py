import sys
import argparse
import gzip
import os
import io
from natsort import natsorted

##This script filters a VCF to retrieve higher confidence somatic SNVs
##It works on with the annotations you have previously added to the VCF
##With --e you can pass the tags to filter out calls, such as segdups etc, that indicate low confidence or non-callability
##The second feature removes calls present in too many individuals (germline heterozygous or recurrent errors)
##As finding somatic mutations in different tissues/replicates from the same individual is expected, 
##you need to input a file with the sample ID in the first column and the individual ID in the second (-s file)
##If you trust your caller is sensitive enough, you can use -a VCF which will use the calls in the VCF to determine how many inds have been called
##Otherwise, you can retrive the information differently and annotate it in the INFO field under any tag such as "samples"
##It has to have the following format: samples=G:sample1,sample4,sample8;T=sample6,sample7
##For both uses, with -n you can determine the number of other individuals allowed to carry support from the same variant

#Parse arguments
parser=argparse.ArgumentParser()

parser.add_argument('-e','--exclude', type=str, required=False, default="EXCLUDE", dest="EXCLUDE",help="Keep calls that do NOT include these tags")
parser.add_argument('-r','--require', type=str, required=False, default="", dest="REQUIRE", help="Keep calls IF they include these tags ")
parser.add_argument('-s','--samplesInds', type=str, required=True, dest="SAMPLEIND", help="File with sample IDs and the corresponding individual")
parser.add_argument('-n','--nOtherInds', type=int, required=False, default=1, dest="NOTHERINDS", help="Keep calls with support in maximum this number of other individuals")
parser.add_argument('-f','--overwriteSamples', type=str, required=False, default="NO", dest="OW", help="Do you want to overwrite sample names in VCF and use those in the samplesInds file? Do this if annotation samples match the sampleInds file")
parser.add_argument('-a','--sampleAnnotation', type=str, required=True, default="VCF", dest="SAMPLEANN", help="Use VCF to calculate support from other samples directly from VCF. Use annotation name to get information from the INFO column. Annotation has to be Allele:sample,sample,sample;Allele:sample")
parser.add_argument('-i','--includesSelf',type=str,required=False, default="NO", dest="SELF",help="Indicate if the VCF samples are included in the annotation")
parser.add_argument('-m','--maxMultipleSamples',type=int, required=False, default=1, dest="MAXMULTI",help="Keep calls with support from multiple samples coming form maximum this number of individuals")
parser.add_argument('-p','--vcfPath', type=str, required=True, dest="VCFPATH", help="Filter all vcfs from this path")
parser.add_argument('-o','--out',type=str, required=True, dest="OUTF", help="Name for out gzipped vcf")

args=parser.parse_args()
if args.SAMPLEANN!="VCF" and args.SELF is None:
	sys.exit("ERROR! --sampleAnnotation not 'VCF' requires --includeSelf")

if args.SELF!="YES" and args.SELF!="NO":
	sys.exit("ERROR! --includeSelf can only be YES or NO")
if args.SELF=="YES":
	args.NOTHERINDS+=1

#Get and sort vcf names
vcfs=[]
for r,d,f in os.walk(args.VCFPATH):
	for file in f:
		if file.split(".")[-2:]==['vcf', 'gz']:
			vcfs.append(args.VCFPATH+"/"+file)

vcfs=natsorted(vcfs)

#Get sample names
samples=[]
for line in io.TextIOWrapper(gzip.open(vcfs[0],"rb")):
	if line.startswith("#CHROM"):
		for s in line.split("\t")[9:]:
			samples.append(s.strip("\n"))

#Check samples match file and get ind info
sampleToInd={}
for line in open(args.SAMPLEIND,'r'):
	sampleToInd[line.split()[0]]=line.split()[1]

if not list(sampleToInd.keys())==samples:
	if args.OW!="YES": 
		sys.exit("ERROR! Sample names in VCF do not match sample names in samplesInd file. If you want to overwrite the VCF names run the program with -f YES")
		

#Open output file
out=gzip.GzipFile(args.OUTF,"wb")

#Write header
for line in io.TextIOWrapper(gzip.open(vcfs[0],"rb")):
	if line.startswith("#"):
		out.write(line.encode())


#Filter lines by exclude tags and number of samples called from a diff ind (assuming sample=indX)
for vcf in vcfs:
	for line in io.TextIOWrapper(gzip.open(vcf,"rb")):
		if not line.startswith("#") and not any(x in line for x in args.EXCLUDE.split(",")) and all(x in line for x in args.REQUIRE.split(",")):
			if args.SAMPLEANN not in line:
				out.write(line.encode())
			else:
				inds=[]
				if args.SAMPLEANN=="VCF":
					coln=0
					for col in line.split()[9:]:
						if "1" in col.split(":")[0]:
							inds.append(sampleToInd[samples[coln]])
						coln+=1
				else:
					for allele in line.split("\t")[7].split(args.SAMPLEANN+"=")[1].split(";"):
						if line.split("\t")[4]==allele.split(":")[0]:
							for x in allele.split(":")[1].split(","):
								try:
									inds.append(sampleToInd[x])
								except KeyError as err:
									print("Samples in sampleinds file do not match annotations in VCF:", err, "is missing from sampleinds")
									sys.exit(1)
				#If annotation in line but with a different allele
				if len(inds)==0:
					out.write(line.encode())
				#Filter according to n of inds
				elif len(set(inds))<=args.NOTHERINDS:
					if args.SELF=="NO":
						out.write(line.encode())
					else:
					#And only 1 ind with multiple samples (self)
						y=0
						for ind in set(inds):
							if inds.count(ind)>1:
								y+=1
						if y<=args.MAXMULTI:
							out.write(line.encode())

out.close()




