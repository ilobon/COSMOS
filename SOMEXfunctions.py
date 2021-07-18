import sys
import gzip
import io


filtNames=["0/1","MIN_AD","MIN_AD_SS","MAX_VAF","DPrange","S_RATIO","P_RATIO","SER1","SER2","SPOIS","SFET","BINOM","MAX_NRL","MAX_NHAP","CNV","PIR","MAX_VAFQ20","MAX_CLIP","MWBQ","MWMQ","MWMM","PON"]


## Functions ##



##Filtering

#1 0/1
def isvariant():
	if "0" in CV['GT'] and "1" in CV['GT']:
		return True
	else:
		return False

#2 AD
def minAD(MIN_AD):
	if MIN_AD!=99999:
		if int(CV['AD'].split(",")[1])>=MIN_AD:
			return True
		else:
			return False

def minADss(MIN_AD_SS):
	if MIN_AD_SS!=99999:
		if int(CV['AD'].split(",")[1])>=MIN_AD_SS:
			return True
		else:
			return False



#3 VAF
def maxVAF(MAX_VAF):
	if MAX_VAF!=999.99:
		if CV['DP']!="." and int(CV['AD'].split(",")[1])/int(CV['DP'])<MAX_VAF:
			return True
		else:
			return False

#4 DP
def DPrange(MIN_DP,MAX_DP):
	if MIN_DP!=99999 and MAX_DP!=99999:
		if CV['DP']!="." and int(CV['DP'])>MIN_DP and int(CV['DP'])<MAX_DP:
			return True
		else:
			return False

#5 SAC ratio
def sratio(S_RATIO):
	if S_RATIO!=99999:
		if CV['SAC']!=".":
			s=[int(i) for i in CV['SAC'].split(',')]
			ssum=[s[0]+s[2],s[1]+s[3]]
			ssum.sort()
			if ssum[0]/ssum[1] >= 1/S_RATIO:
				return True
			else:
				return False

#6 and 7 RAC ratios
def pratio(P_RATIO):
	if P_RATIO!=99999:
		r=[int(i) for i in CV['RAC'].split(',')]
		rR=r[0:2]; rA=r[2:4]
		rR.sort(); rA.sort()
		if (sum(rR)==0 or rR[0]/rR[1] >= 1/P_RATIO) and (sum(rA)==0 or rA[0]/rA[1] >= 1/P_RATIO or (rA[0]==0 and rA[1]<=P_RATIO)):
			return True
		else:
			return False

#8 SER1
def ser1(SER1):
	if SER1!=999.99:
		if CV['SER1']!="." and float(CV['SER1'])<SER1:
			return True
		else:
			return False

#9 SER2
def ser2(SER2):
	if SER2!=999.99:
		if CV['SER2']!="." and float(CV['SER1'])<SER2:
			return True
		else:
			return False

#10 strand Poisson
def spois(SPOIS):
	if SPOIS!=999.99:
		if CV['SPOIS']!="." and float(CV['SPOIS'])>SPOIS:
			return True
		else:
			return False

#11 strand FET
def sfet(SFET):
	if SFET!=999.99:
		if CV['SFET']!="." and float(CV['SFET'])>SFET:
			return True
		else:
			return False

#12 Binomial g het
def binom(BINOM):
	if BINOM!=999.99:
		if CV['BINOM']!="." and float(CV['BINOM'])<BINOM:
			return True
		else:
			return False

#13 n within read length
def maxnrl(MAX_NRL):
	if MAX_NRL!=99999:
		if int(CV['nRL'])<=MAX_NRL:
			return True
		else:
			return False

#14 n haplotypes
def nhap(MAX_NHAP):
	if MAX_NHAP!=99999:
		if int(CV['PHA'].split(",")[0])<=MAX_NHAP:
			return True
		else:
			return False

#15 CNV call
def cnv(CNV):
	if CNV!="empty":
		if CNV=="NO":
			if CV['SQ']==".":
				return True
		else:
			return False

#16 PIR
def pir(PIR):
	if PIR!="empty":
		if CV['PIR90'] in [i for i in PIR.split(";")]:
			return True
		else:
			return False

#17 VAF Q20
def maxVAFq20(MAX_VAFQ20):
	if MAX_VAFQ20!=999.99:
		if CV['Q20VAF']!="." and float(CV['Q20VAF'])<MAX_VAFQ20:
			return True
		else:
			return False

#18 Clipping 
def clip(MAX_CLIP):
	if MAX_CLIP!=999.99:
		c=[i for i in CV['CLIP'].split(",")]
		if c[0]=="." or c[1]=="." or float(c[0])-float(c[1])>(-MAX_CLIP) and float(c[0])-float(c[1])<MAX_CLIP:
			return True
		else:
			return False

#19 MWBQ
def mwbq(MWBQ):
	if MWBQ!=999.99:
		if CV['MWBQ']=="." or float(CV['MWBQ'])>MWBQ:
			return True
		else:
			return False

#20 MWMQ
def mwmq(MWMQ):
	if MWMQ!=999.99:
		if CV['MWMQ']=="." or float(CV['MWMQ'])>MWMQ:
			return True
		else:
			return False

#21 MWMM
def mwmm(MWMM):
	if MWMM!=999.99:
		if CV['MWMM']=="." or float(CV['MWMM'])>MWMM:
			return True
		else:
			return False

##22 Sanger PON
def spon(PON):
	if PON!=999.99:
		if CV['PON']=="." or float(CV['PON'])<PON:
			return True
		else:
			return False




#Get field numbers
def getFieldNumbers(vcfline):
	fieldNumbers={}
	i=0
	for x in vcfline.split()[8].split(":"):
		fieldNumbers[x]=i
		i+=1
	return fieldNumbers




def filterColumn(samplecol,fieldNumbers,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON):
	global CV
	CV={}
	for name, field in fieldNumbers.items():
		CV[name]=samplecol.split(":")[field]
	return [isvariant(), minAD(MIN_AD),minADss(MIN_AD_SS), maxVAF(MAX_VAF), DPrange(MIN_DP,MAX_DP), sratio(S_RATIO), pratio(P_RATIO), ser1(SER1), ser2(SER2), spois(SPOIS), sfet(SFET), binom(BINOM), maxnrl(MAX_NRL), nhap(MAX_NHAP), cnv(CNV), pir(PIR), maxVAFq20(MAX_VAFQ20), clip(MAX_CLIP), mwbq(MWBQ), mwmq(MWMQ), mwmm(MWMM), spon(PON)]



#Get sample names from VCF
def getSamplesVCF(vcfgz):
	samples=[]
	for line in io.TextIOWrapper(gzip.open(vcfgz,'rb')):
		if line.startswith('#CHROM'):
			for s in line.split('\t')[9:]:
				samples.append(s.strip('\n'))
	return samples





#Get sample and ind names from file with 2 cols
def getSampleInd(file):
	sampleToInd={}
	for line in open(file,'r'):
		sampleToInd[line.split()[0]]=line.split()[1]
	return sampleToInd


def filterLine(COMBINED,NPERIND,line,sampleToInd,fieldNumbers,samples,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON):
	index=0
	result=[]
	for eachind in set(sampleToInd.values()):
		cols=[]
		filtersInd=[]
		#For a certain ind record filter result per sample + cols of that ind
		for sample, ind in sampleToInd.items():
			if ind == eachind:
				col=index+9
				filtersSample=filterColumn(line.split()[col],fieldNumbers,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON)
				filtersInd.append(filtersSample)
				cols.append(col)
#				print(col, filtersSample)
				index+=1
		#count samples passing BINOM and Q20VAF
		nbinom=0; nq20vaf=0
		for fs in filtersInd:
			if fs[11]!=False:
				nbinom+=1
			if fs[16]!=False:
				nq20vaf+=1
		#samplesPassing already controls for MIN_AD_SS - check after
		if COMBINED=='NO':
#			print("inside NO")
			samplesPassing=[]
			for i in range(0,len(filtersInd)):
				if  filtersInd[i].count(False)==0:
					samplesPassing.append(samples[cols[i]-9])
			if len(samplesPassing)>0 and nbinom==len(filtersInd) and nq20vaf==len(filtersInd):
				for s in samplesPassing:
					result.append("\t".join([s]+line.split()[0:9]+[line.split()[samples.index(s)+9]]))
		#Now we need to calculate n of samples passing each tissue
		else:
#			print("calculate filtPassed")
			filtPassed=[]
			for f in range(0, len(filtersInd[0])):
				passing=0
				for s in filtersInd:
					if s[f]!=False:
						passing+=1
				filtPassed.append(passing)
#			print(filtPassed)
#			print(all(n>=NPERIND for n in filtPassed[0:2]+filtPassed[3:]))
		#Here we do not care about MIN_AD_SS so always remove it
			if COMBINED=='YES' and all(n>=NPERIND for n in filtPassed[0:2]+filtPassed[3:]):
#				print("inside YES")
				result.append("\t".join([sampleToInd[samples[index-1]]]+line.split()[0:9]+line.split()[cols[0]:(cols[-1]+1)]))
		#Now for BOTH
			elif COMBINED=='BOTH':
#				print("inside BOTH")
				samplesPassing=[]
				for i in range(0,len(filtersInd)):
					if  filtersInd[i].count(False)==0:
						samplesPassing.append(samples[cols[i]-9])
				if MIN_AD_SS==99999:
					if all(n>=NPERIND for n in filtPassed):
						result.append("\t".join([",".join([sampleToInd[samples[index-1]]]+samplesPassing)]+line.split()[0:9]+line.split()[cols[0]:(cols[-1]+1)]))
					elif len(samplesPassing)>0 and nbinom==len(filtersInd) and nq20vaf==len(filtersInd):
						result.append("\t".join([",".join(samplesPassing)]+line.split()[0:9]+line.split()[cols[0]:(cols[-1]+1)]))
				elif  MIN_AD_SS!=99999:
					if all(n>=NPERIND for n in filtPassed[0:2]+filtPassed[3:]):
						result.append("\t".join([",".join([sampleToInd[samples[index-1]]]+samplesPassing)]+line.split()[0:9]+line.split()[cols[0]:(cols[-1]+1)]))
					elif len(samplesPassing)>0 and nbinom==len(filtersInd) and nq20vaf==len(filtersInd):
						result.append("\t".join([",".join(samplesPassing)]+line.split()[0:9]+line.split()[cols[0]:(cols[-1]+1)]))
	return("\n".join(result))



def getFailed(COMBINED,NPERIND,line,sampleToInd,fieldNumbers,samples,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON):
	index=0
	result=[]
	cols=[]
	filtersInd=[]
	#For a certain ind record filter result per sample + cols of that ind
	for sample, ind in sampleToInd.items():
		col=index+9
		filtersSample=filterColumn(line.split()[col],fieldNumbers,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON)
		filtersSampleDict={}
		i=0
		for x in filtNames:
			filtersSampleDict[x]=filtersSample[i]
			i+=1
		filtersInd.append(filtersSampleDict)
		cols.append(col)
		index+=1
	j=0
	falseFilts=[]
	for x in filtersInd:
		sampleFalseFilts=[samples[j]]
		for filt,value in x.items():
			if value==False:
				sampleFalseFilts.append(filt)
		falseFilts.append(sampleFalseFilts)
		j+=1
	return(falseFilts)


def getFiltResultsAllSamples(COMBINED,NPERIND,line,sampleToInd,fieldNumbers,samples,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON):
	index=0
	result=[]
	cols=[]
	filtersInd=[]
	#For a certain ind record filter result per sample + cols of that ind
	for sample, ind in sampleToInd.items():
		col=index+9
		filtersSample=filterColumn(line.split()[col],fieldNumbers,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM,PON)
		filtersSampleDict={}
		i=0
		for x in filtNames:
			filtersSampleDict[x]=filtersSample[i]
			i+=1
		filtersInd.append(filtersSampleDict)
		cols.append(col)
		index+=1
	j=0
	for x in filtersInd:
		result.append([samples[j],x])
		j+=1
	return(result)

















































### OLD just in case

#def filterLine(COMBINED,NPERIND,line,sampleToInd,fieldNumbers,samples,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM):
#	index=0
#	for eachind in set(sampleToInd.values()):
#	filtersInd=[]
#	samplesPassing=[]
#	for sample, ind in sampleToInd.items():
#		if ind == eachind:
#		col=index+9
#		filtersSample=filterColumn(line.split()[col],fieldNumbers,MIN_AD,MIN_AD_SS,MAX_VAF,MIN_DP,MAX_DP,S_RATIO,P_RATIO,SER1,SER2,SPOIS,SFET,BINOM,MAX_NRL,MAX_NHAP,CNV,PIR,MAX_VAFQ20,MAX_CLIP,MWBQ,MWMQ,MWMM)
#		filtersInd.append(filtersSample)
#		if filtersSample.count(False)==0:
#			samplesPassing.append(samples[index])
#		index+=1
#	filtPassed=[]
#	for f in range(0, len(filtersInd[0])):
#		passing=0
#		for s in filtersInd:
#		if s[f]!=False:
#			passing+=1
#		filtPassed.append(passing)
#	if all(n>=NPERIND for n in filtPassed) and COMBINED!='NO':
#		indPassing=sampleToInd[samples[index-1]]
#		if COMBINED=="YES":
#		return("\t".join([indPassing]+line.split()[0:9]+line.split()[col-4:col+1]))
#		elif COMBINED=="BOTH":
#		printPassing=",".join([indPassing]+samplesPassing)
#		return("\t".join([printPassing]+line.split()[0:9]+line.split()[col-4:col+1]))
#	elif len(samplesPassing)>0 and filtPassed[filtNames.index("BINOM")]>=NPERIND:
#		if COMBINED=="NO":
#		multipleLines=[]
#		for k in samplesPassing:
#			multipleLines.append("\t".join([k]+line.split()[0:9]+[line.split()[list(sampleToInd.keys()).index(k)+9]]))
#		return("\n".join(multipleLines))
#		elif COMBINED=="BOTH":
#		return("\t".join([",".join(samplesPassing)]+line.split()[0:9]+line.split()[col-4:col+1]))
#
#



