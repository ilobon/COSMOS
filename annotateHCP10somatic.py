import sys
import pysam
import gzip
import io
import scipy.stats


vcf=sys.argv[1]
bams=sys.argv[2].split(";")
vcfxhmm=sys.argv[3]
outf=sys.argv[4]

nucleotides=["A","C","G","T"]
readlength=150


#### Functions #####

## RAC ##
def read1(read):
	if read.is_read1:
		return True
	else:
		return False

def read2(read):
	if read.is_read2:
		return True
	else:
		return False

## PIR ## 
def n1t(read):
        if p-read.pos<=(readlength*0.33):
                return True
        else:
                return False

def n2t(read):
        if p-read.pos>(readlength*0.33) and p-read.pos<(readlength*0.66):
                return True
        else:
                return False

def n3t(read):
        if p-read.pos>=(readlength*0.66):
                return True
        else:
                return False

## nRL and PHA ##

#Get non homozygous lines per column 
def getLines(col):
        lines=[]
        for line in io.TextIOWrapper(gzip.open(vcf,"rb")):
                if not line.startswith("#") and "0" in line.split("\t")[col].split(":")[0] and "1" in line.split("\t")[col].split(":")[0]:
                        lines.append(line)
        return lines

#Make groups of lines within read length
def makeGroups(lines):
        groups=[]
        subgroup=[]
        pchr=lines[0].split("\t")[0]
        pp=int(lines[0].split("\t")[1])-50
        for line in lines:
                cchr=line.split("\t")[0]
                cp=line.split("\t")[1]
                if cchr==pchr and int(pp)+(readlength-1)>=int(cp):
                        subgroup.append(line)
                else:
                        groups.append(subgroup)
                        subgroup=[]
                        subgroup.append(line)
                pchr=cchr
                pp=cp
        return groups

#Get nRL
def countVars(groups):
        result=[]
        for subgroup in groups:
                n=len(subgroup)
                for line in subgroup:
                        result.append([line.split("\t")[0],int(line.split("\t")[1]),n])
        return result

#Get haplotype info
def getHaplotypes(groups,bamfile):
        result=[]
        for g in range(0,len(groups)-1):
                if len(groups[g])>1:
                        for i in range(1,len(groups[g])):
                                chr=groups[g][0].split("\t")[0]
                                p=int(groups[g][i-1].split("\t")[1])-1
                                p2=int(groups[g][i].split("\t")[1])-1
                                bases=[]
                                for read in bamfile.fetch(chr,p,p+1):
                                        if p-read.pos<read.query_alignment_end-read.query_alignment_start:
                                                bases.append(read.query_alignment_sequence[p-read.pos])
                                ubases=[]
                                for x in set(bases): ubases.append(x)
                                poscounts=[]
                                for allele in ubases:
                                        x=[]
                                        for read in bamfile.fetch(chr,p,p+1):
                                                if p-read.pos<read.query_alignment_end-read.query_alignment_start and read.query_alignment_sequence[p-read.pos]==allele and (read.pos+read.query_alignment_end-read.query_alignment_start)>p2 and (p2-read.pos)>145:
                                                        x.append(read.query_alignment_sequence[p2-read.pos])
                                        if len(x)>0:
                                                for a in set(x):
                                                        if x.count(a)>1:
                                                                poscounts.append([allele,a,x.count(a)])
                                if len(poscounts)>1:
                                        result.append([chr,p+1,p2+1,len(poscounts)])
        return result



#### 1. Make header with new fields #####


#I want to add the lines after the last of ##FORMAT lines
out=open(outf,"w")
done=0
nformat=0
for line in io.TextIOWrapper(gzip.open(vcf,"rb")):
    if line.startswith("#"):
            #Before getting to ##FORMAT lines, just write them
            if not line.startswith("##FORMAT") and nformat==0:
                    out.write(line)
            else:
                    #Then as long as lines start with ##FORMAT keep writting them
                    if line.startswith("##FORMAT"):
                            nformat+=1
                            out.write(line)
                    #We've gotten to ##FORMAT lines and this line doesn't start with ##FORMAT so write new lines
                    else:
                            if done==0:
                                    out.write('##FORMAT=<ID=RAC,Number=.,Type=Integer,Description="Number of R1 and R2 reads supporting each allele (R1,R2,A1,A2)">'+'\n')
                                    out.write('##FORMAT=<ID=SER1,Number=.,Type=Float,Description="Fisher exact test p-value for R1 Ref and Alt counts vs error for the main allele nucleotide in the R1">'+'\n')
                                    out.write('##FORMAT=<ID=SER2,Number=.,Type=Float,Description="Fisher exact test p-value for R2 Ref and Alt counts vs error for the main allele nucleotide in the R2">'+'\n')
                                    out.write('##FORMAT=<ID=SPOIS,Number=.,Type=Float,Description="Poisson test for read counts for each strand">'+'\n')
                                    out.write('##FORMAT=<ID=SFET,Number=.,Type=Float,Description="Fisher exact test p-value for Ref+ and Ref- vs Alt+ and Alt- read counts">'+'\n')
                                    out.write('##FORMAT=<ID=BINOM,Number=.,Type=Float,Description="Binomial p-value for Ref and total read counts">'+'\n')
                                    out.write('##FORMAT=<ID=PIR90,Number=.,Type=Integer,Description="Position in read per allele, where 0 means no coverage, 1,2 or 3 >0.9 of the reads carrying the allele in that third of the read and 4 means no bias">'+'\n')
                                    out.write('##FORMAT=<ID=MC,Number=.,Type=Integer,Description="Maximum coverage in +- read length region">'+'\n')
                                    out.write('##FORMAT=<ID=nRL,Number=.,Type=Integer,Description="Number of variants called in the individual +- read length around, including the variant">'+'\n')
                                    out.write('##FORMAT=<ID=PHA,Number=.,Type=String,Description="Number of haplotypes this variant has in fase with, chr, position">'+'\n')
                                    out.write('##FORMAT=<ID=SQ,Number=2,Type=Float,Description="XHMM phred-scaled Qualities of Some CNV event of allele types: DEL,DUP">'+"\n"+'##FORMAT=<ID=RDX,Number=1,Type=Float,Description="XHMM mean Read Depth over region">'+"\n")
                                    out.write('##FORMAT=<ID=Q20VAF,Number=.,Type=Float,Description="Variant allele frequency of reads with BQ>20 with pysam">'+'\n')
                                    out.write('##FORMAT=<ID=CLIP,Number=.,Type=Float,Description="Ratio of soft clipped reads per allele (ref,alt)">'+'\n')
                                    out.write('##FORMAT=<ID=MWBQ,Number=.,Type=Float,Description="Mann-Withney U test p-value for base quality scores (ref greater than alt)">'+'\n')
                                    out.write('##FORMAT=<ID=MWMQ,Number=.,Type=Float,Description="Mann-Withney U test p-value for mapping quality scores (ref greater than alt)">'+'\n')
                                    out.write('##FORMAT=<ID=MWMM,Number=.,Type=Float,Description="Mann-Withney U test p-value for number of mismatches per read (ref less than alt)">'+'\n')
                                    #We record we've done it so that it isn't written for every line after ##FORMAT that doesn't start with ##FORMAT 
                                    done+=1
                                    out.write(line)
                            #Now output the lines after
                            else:
                                    out.write(line)

out.close()


#### 2. Run cluster and phasing info per sample ####

lines={}
groups={}
resultn={}
resulth={}
for col in range(9,len(bams)+9):
        lines[col]=[]
        groups[col]=[]
        resultn[col]=[]
        resulth[col]=[]
        #Get lines for that column
        for i in getLines(col):
                lines[col].append(i)
        if lines[col]==[]:
                break
        #Make groups with the lines
        for i in makeGroups(lines[col]):
                groups[col].append(i)
        #Get nRL
        for i in countVars(groups[col]):
                resultn[col].append(i)
        #Open bam
        bam=bams[col-9]
        bamfile=pysam.AlignmentFile(bam,"rb")
        #Get haplotypes 
        for i in getHaplotypes(groups[col],bamfile):
                resulth[col].append(i)

#### 3. Read XHMM vcf ####

## XHMM ##

#Get position SQ and RD from reading xhmm vcf
resultxhmm=[]
for line in io.TextIOWrapper(gzip.open(vcfxhmm,"rb")):
	if not line.startswith("#"):
		lan=line.split("\t")[0:2]
		for col in range(9,len(bams)+9):
			lan.append(line.split("\t")[col].split(":")[4]+':'+line.split("\t")[col].split(":")[9])
		resultxhmm.append(lan)

#### 4. Add new fields ####


out=open(outf,"a")

#Errors from Schirmer et al 2016 for HiSeq per nucleotide and R1/R2
errors1=dict(zip(["A","C","G","T"],[0.0004,0.0004,0.0005,0.0008]))
errors2=dict(zip(["A","C","G","T"],[0.0008,0.0008,0.0010,0.0015]))

for line in io.TextIOWrapper(gzip.open(vcf,"rb")):
	if not line.startswith("#"):
		chr=line.split("\t")[0]
		p=int(line.split("\t")[1])-1 #Python is 0 based
		lbroken=line.split("\t")
		alleles=[lbroken[3],lbroken[4]]
		lbroken[8]=lbroken[8]+":RAC:SER1:SER2:SPOIS:SFET:BINOM:PIR90:MC:nRL:PHA:SQ:RDX:Q20VAF:CLIP:MWBQ:MWMQ:MWMM" #Add new fields names
		col=9
		ns=[lbroken[3],lbroken[4]]	
		for bam in bams:
			print(chr,p,bam,col)
			bamfile=pysam.AlignmentFile(bam,"rb") #Read bam with pysam
			#RAC
			countsR1=[]
			for i in bamfile.count_coverage(chr,p,p+1,quality_threshold=10,read_callback=read1):
				countsR1.append("".join(str(x) for x in i)) #Counts on R1 per nucleotide in list
			countsR2=[]
			for i in bamfile.count_coverage(chr,p,p+1,quality_threshold=10,read_callback=read2):
				countsR2.append("".join(str(x) for x in i))
			dictR1=dict(zip(nucleotides,countsR1))
			dictR2=dict(zip(nucleotides,countsR2))
			RACounts=[]
			for i in alleles:
				RACounts.append(",".join(str(x) for x in [dictR1[i],dictR2[i]])) #Get R1 and R2 counts for alleles on VCF
			RAC=[]
			for i in RACounts:
				for x in i.split(","):
					RAC.append(int(x))
			print("RAC annotated",p,col)
			#Tests
			if "0" not in lbroken[col].split(":")[0] and "1" not in lbroken[col].split(":")[0]:
				SER1=".";SER2=".";SPOIS=".";SFET=".";BINOM="."
			else:
				#Alleles depth
				AD=list(map(int,lbroken[col].split(":")[1].split(",")))
				#HaplotypeCaller's SAC can be "." or have 0,0 counts for a non reported allele
				SAC="." if lbroken[col].split(":")[5].replace("\n","")=="." else list(map(int,lbroken[col].replace("\n","").split(":")[5].split(",")))
				SACsums=[]
				if SAC!=".":
					for i in range(0,len(SAC),2): SACsums.append(SAC[i]+SAC[i+1])
				while len(SACsums)>2 and len(AD)<len(SACsums):
					empty=SACsums.index(min(SACsums))
					del SAC[empty*2+1]
					del SAC[empty*2]
					SACsums=[]
					for i in range(0,len(SAC),2): SACsums.append(SAC[i]+SAC[i+1])
					#Since we've corrected SAC, add it to the column corrected
				lbroken[col]=":".join([item for sublist in [lbroken[col].split(":")[:5],[",".join(str(x) for x in SAC)],lbroken[col].split(":")[6:]]for item in sublist])				
				RACsums=[]
				for i in range(0,len(RAC),2): RACsums.append(RAC[i]+RAC[i+1])
				#Index of most frequent and second allele for tests orientation 
				n=AD.index(max(AD))
				n2=AD.index(sorted(AD,reverse=True)[1])
				if n2==n: n2=n+1	
				#Sequencing error test for R1 and R2
				SER1=scipy.stats.fisher_exact([[RAC[n*2],RAC[n2*2]],[(RAC[n*2]+RAC[n2*2])*(1-errors1[alleles[n]]),(RAC[n*2]+RAC[n2*2])*errors1[alleles[n]]]])[1]
				SER2=scipy.stats.fisher_exact([[RAC[n*2+1],RAC[n2*2+1]],[(RAC[n*2+1]+RAC[n2*2+1])*(1-errors2[alleles[n]]),(RAC[n*2+1]+RAC[n2*2+1])*errors2[alleles[n]]]])[1]
				BINOM=scipy.stats.binom_test(AD[n],(AD[n]+AD[n2]))
				if SAC=="." or sum(SAC)==0: SPOIS=".";SFET="."
				else:
					SACp=SAC[0]+SAC[2]
					SACm=SAC[1]+SAC[3]
					SPOIS=scipy.stats.poisson.cdf(min([SACp,SACm]),sum(SAC)/2)
					SFET=scipy.stats.fisher_exact([[SAC[n*2],SAC[n*2+1]],[SAC[n2*2],SAC[n2*2+1]]])[1]
			print("Tests annotated",p,col)
			#PIR
			cov1=dict(zip(nucleotides,[str(x) for y in bamfile.count_coverage(chr,p,p+1,quality_threshold=10,read_callback=n1t) for x in y]))
			cov2=dict(zip(nucleotides,[str(x) for y in bamfile.count_coverage(chr,p,p+1,quality_threshold=10,read_callback=n2t) for x in y]))
			cov3=dict(zip(nucleotides,[str(x) for y in bamfile.count_coverage(chr,p,p+1,quality_threshold=10,read_callback=n3t) for x in y]))
			PIR=[]
			for allele in alleles:
				tot=int(cov1[allele])+int(cov2[allele])+int(cov3[allele])
				if tot==0:
					pir=0
				elif int(cov1[allele])/tot>0.9:
					pir=1
				elif int(cov2[allele])/tot>0.9:
					pir=2
				elif int(cov3[allele])/tot>0.9:
					pir=3
				else:
					pir=4
				PIR.append(pir)
			print("PIR annotated",p,col)
			#Max cov
			covs=[]
			for pileupcolumn in bamfile.pileup(chr,p-readlength,p+readlength):
					if pileupcolumn.pos>=p-readlength and pileupcolumn.pos<=p+readlength:
						covs.append(pileupcolumn.n)
			MC=0 if len(covs)==0 else max(covs)
			print("MC annotated",p,col)
			#nRL and PHA
			if col in resultn:
				inlistn=[i for i, e in enumerate(resultn[col]) if e[0]==chr and e[1]==(p+1)]
				if len(inlistn)==1:
					annotationn=resultn[col][inlistn[0]][2]
				else:
					annotationn=0
			else:
				annotationn=0
			if col in resulth:
				inlisth=[i for i, e in enumerate(resulth[col]) if e[0]==chr and e[1]==(p+1)]
				if len(inlisth)==1:
					annotationh=",".join(str(x) for x in [resulth[col][inlisth[0]][3],resulth[col][inlisth[0]][0],resulth[col][inlisth[0]][2]])	
				else:
					annotationh=0
			else:
				annotationh=0
			print("nRL and PHA annotated",p,col)
			#XHMM fields
			inresult=[i for i, e in enumerate(resultxhmm) if e[0]==chr and int(e[1])==(p+1)]
			if len(inresult)==1:
				XHMM=resultxhmm[inresult[0]][col-7]
			else:
				XHMM="."+":"+"."
			print("XHMM annotated",p,col)
			#Q20VAF,CLIP,MWBQ,MWMQ,MWMM
			covq20=dict(zip(nucleotides,[x for y in bamfile.count_coverage(chr,p,p+1,quality_threshold=20) for x in y]))
			q20VAF=0 if covq20[ns[1]]==0 and covq20[ns[0]]==0 else covq20[ns[1]]/(covq20[ns[0]]+covq20[ns[1]])
			refN=0; refclip=0; altN=0; altclip=0; mm=0
			rBQ=[]; aBQ=[]; rMQ=[]; aMQ=[]; mmR=[]; mmA=[]
			for read in bamfile.fetch(chr,p,p+1):
				for read_pos, ref_pos, base in read.get_aligned_pairs(with_seq=True):
					if base is not None and base in "acgt":
						mm+=1
				if p-read.pos<read.query_alignment_end-read.query_alignment_start and read.seq[p-read.pos]==ns[0]:
					refN+=1
					if "S" in read.cigarstring:
						refclip+=1
					rBQ.append(read.query_alignment_qualities[p-read.pos])
					rMQ.append(read.mapping_quality)
					mmR.append(mm)
				elif p-read.pos<read.query_alignment_end-read.query_alignment_start and read.seq[p-read.pos]==ns[1]:
					altN+=1
					if "S" in read.cigarstring:
						altclip+=1
					aBQ.append(read.query_alignment_qualities[p-read.pos])
					aMQ.append(read.mapping_quality)
					mmA.append(mm-1)
				mm=0
			RCLIP="." if refN==0 else str("%.2f" % (refclip/refN))
			ACLIP="." if altN==0 else str("%.2f" % (altclip/altN))
			MWBQ="." if (len(set(rBQ+aBQ))==1 or (len(rBQ)==0 or len(aBQ)==0)) else str("%.4f" % scipy.stats.mannwhitneyu(rBQ,aBQ,alternative='greater')[1])
			MWMQ="." if (len(set(rMQ+aMQ))==1 or (len(rMQ)==0 or len(aMQ)==0)) else str("%.4f" % scipy.stats.mannwhitneyu(rMQ,aMQ,alternative='greater')[1])
			MWMM="." if (len(set(mmR+mmA))==1 or (len(mmR)==0 or len(mmA)==0)) else str("%.4f" % scipy.stats.mannwhitneyu(mmR,mmA,alternative='less')[1])
			print("Q20VAF,CLIP,MWBQ,MWMQ,MWMM annotated",p,col)
			#Add all to line	
			if "0" not in lbroken[col].split(":")[0] and "1" not in lbroken[col].split(":")[0]:
				lbroken[col]=":".join(str(x) for x in [lbroken[col].replace("\n",""),",".join(str(x) for x in RAC),SER1,SER2,SPOIS,SFET,BINOM,",".join(str(x) for x in PIR),MC,annotationn,annotationh,XHMM,str("%.3f" % q20VAF),",".join([RCLIP,ACLIP]),":".join([MWBQ,MWMQ,MWMM])])
			else:
				if SAC=="." or sum(SAC)==0:
					lbroken[col]=":".join(str(x) for x in [lbroken[col].replace("\n",""),",".join(str(x) for x in RAC),float(("%0.4f"%SER1)),float(("%0.4f"%SER2)),SPOIS,SFET,float(("%0.4f"%BINOM)),",".join(str(x) for x in PIR),MC,annotationn,annotationh,XHMM,str("%.3f" % q20VAF),",".join([RCLIP,ACLIP]),":".join([MWBQ,MWMQ,MWMM])])
				else:
					lbroken[col]=":".join(str(x) for x in [lbroken[col].replace("\n",""),",".join(str(x) for x in RAC),float(("%0.4f"%SER1)),float(("%0.4f"%SER2)),float(("%0.4f"%SPOIS)),float(("%0.4f"%SFET)),float(("%0.4f"%BINOM)),",".join(str(x) for x in PIR),MC,annotationn,annotationh,XHMM,str("%.3f" % q20VAF),",".join([RCLIP,ACLIP]),":".join([MWBQ,MWMQ,MWMM])])
			print(col,"annotated")
			col+=1
		lbroken[col-1]=lbroken[col-1]+"\n" #Add new line to last field
		out.write("\t".join(str(x) for x in lbroken))



out.close()






