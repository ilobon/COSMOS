import sys
import gzip
import io
import argparse
import SOMEXfunctions

## Parse arguments ##

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--mode', type=str, required=True, dest='MODE',help='')
parser.add_argument('-p', '--pos', type=str,required=True, dest='p',help='')
parser.add_argument('-i', '--vcf', type=str, required=True, dest='VCF', help='Input VCF file')
##parser.add_argument('-o', '--outfile', type=str, required=True, dest='OUTF', help='Output file name')
#Maybe not require then each sample one ind, combined NO?
parser.add_argument('-is', '--samplesInds', type=str, required=True, dest='SAMPLEIND', help='File with sample IDs and the corresponding individual')
parser.add_argument('-c', '--combined',type=str, required=True, default='NO', dest='COMBINED', help='Different tissues can pass each filter: YES/NO/BOTH')
#Only id combined BOTH or YES
parser.add_argument('-ns', '--nSamplesPerInd', type=int, required=False, default=99999, dest='NPERIND', help='Keep calls passing each filter in at least int samples per individual')
#Filters:
parser.add_argument('-ad', '--minAD', type=int, required=False, default=99999, dest='MIN_AD', help='')
parser.add_argument('-adss', '--minADsingleSample', type=int, required=False, default=99999, dest='MIN_AD_SS', help='')
parser.add_argument('-vaf', '--maxVAF', type=float, required=False, default=999.99, dest='MAX_VAF', help='')
parser.add_argument('-dp1', '--minDP', type=int, required=False, default=99999, dest='MIN_DP', help='')
parser.add_argument('-dp2', '--maxDP', type=int, required=False, default=99999, dest='MAX_DP', help='')
parser.add_argument('-sr', '--strandRatio', type=int, required=False, default=99999, dest='S_RATIO', help='')
parser.add_argument('-pr', '--readpairRatio', type=int, required=False, default=99999, dest='P_RATIO', help='')
parser.add_argument('-ser1', '--seqErrR1', type=float, required=False, default=999.99, dest='SER1', help='')
parser.add_argument('-ser2', '--seqErrR2', type=float, required=False, default=999.99, dest='SER2', help='')
parser.add_argument('-sp', '--stradPOISSON', type=float, required=False, default=999.99, dest='SPOIS', help='')
parser.add_argument('-sf', '--strandFET', type=float, required=False, default=999.99, dest='SFET', help='')
parser.add_argument('-b', '--ghBinomial', type=float, required=False, default=999.99, dest='BINOM', help='')
parser.add_argument('-nrl', '--maxNreadlength', type=int, required=False, default=99999, dest='MAX_NRL', help='')
parser.add_argument('-hap', '--maxHaplotypes', type=int, required=False, default=99999, dest='MAX_NHAP', help='')
parser.add_argument('-cnv', '--copyNumberVariant', type=str, required=False, default="empty", dest='CNV', help='')
parser.add_argument('-pir', '--positionInReads', type=str, required=False, default="empty", dest='PIR', help='')
parser.add_argument('-vafq', '--maxVAFq20', type=float, required=False, default=999.99, dest='MAX_VAFQ20', help='')
parser.add_argument('-clip', '--maxDifClipping', type=float, required=False, default=999.99, dest='MAX_CLIP', help='')
parser.add_argument('-bq', '--baseQualMW', type=float, required=False, default=999.99, dest='MWBQ', help='')
parser.add_argument('-mq', '--mappingQualMW', type=float, required=False, default=999.99, dest='MWMQ', help='')
parser.add_argument('-mm', '--mismatchesMW', type=float, required=False, default=999.99, dest='MWMM', help='')
parser.add_argument('-pon', '--PONbbPV', type=float, required=False, default=999.99, dest='PON', help='')

args=parser.parse_args()




##Main

samples=SOMEXfunctions.getSamplesVCF(args.VCF)
sampleToInd=SOMEXfunctions.getSampleInd(args.SAMPLEIND)

i=0
for line in io.TextIOWrapper(gzip.open(args.VCF,'rb')):
	if not line.startswith("#"):
		while i==0:
			fieldNumbers=SOMEXfunctions.getFieldNumbers(line)
			i+=1

for line in io.TextIOWrapper(gzip.open(args.VCF,"rb")):
	if not line.startswith("#") and line.split("\t")[1]==args.p:
		myline=line

if args.MODE=="failed":
	result=SOMEXfunctions.getFailed(args.COMBINED,args.NPERIND,myline,sampleToInd,fieldNumbers,samples,args.MIN_AD,args.MIN_AD_SS,args.MAX_VAF,args.MIN_DP,args.MAX_DP,args.S_RATIO,args.P_RATIO,args.SER1,args.SER2,args.SPOIS,args.SFET,args.BINOM,args.MAX_NRL,args.MAX_NHAP,args.CNV,args.PIR,args.MAX_VAFQ20,args.MAX_CLIP,args.MWBQ,args.MWMQ,args.MWMM,args.PON)
	for x in result:
		print(x)
elif args.MODE=="all":
	result=SOMEXfunctions.getFiltResultsAllSamples(args.COMBINED,args.NPERIND,myline,sampleToInd,fieldNumbers,samples,args.MIN_AD,args.MIN_AD_SS,args.MAX_VAF,args.MIN_DP,args.MAX_DP,args.S_RATIO,args.P_RATIO,args.SER1,args.SER2,args.SPOIS,args.SFET,args.BINOM,args.MAX_NRL,args.MAX_NHAP,args.CNV,args.PIR,args.MAX_VAFQ20,args.MAX_CLIP,args.MWBQ,args.MWMQ,args.MWMM,args.PON)
	for x in result:
		print(x)


