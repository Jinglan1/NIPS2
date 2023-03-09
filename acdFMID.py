'''    
    <acdFMID.py>
    <ACD & FMID Filter On SNV Only VCF>
    Copyright (C) <2021>  <Jing Ruilin>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any
    later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
'''


# !/usr/bin/python3
import sys
import os
import time
import math
import re
import gzip
import argparse

# Third-party imports
import pysam
import numpy
import scipy.stats

usage = '''

=============================================================================================================

Version: 1.1.20201106   Date: Dec 11 CST 2020
Author: Jing Ruilin     Contact: ruilin\.jing\@gmail.com

Description

   ACD & FMID Filter On SNV Only VCF
           
                    
=============================================================================================================

'''

epilog = r'''
                  ___       ___           ___           ___ 
      ___        /\__\     /\  \         |\__\         /\__\
     /\  \      /:/  /    /::\  \        |:|  |       /:/  /
     \:\  \    /:/  /    /:/\:\  \       |:|  |      /:/  / 
     /::\__\  /:/  /    /::\~\:\  \      |:|__|__   /:/  /  
  __/:/\/__/ /:/__/    /:/\:\ \:\__\ ____/::::\__\ /:/__/   
 /\/:/  /    \:\  \    \:\~\:\ \/__/ \::::/~~/~    \:\  \   
 \::/__/      \:\  \    \:\ \:\__\    ~~|:|~~|      \:\  \  
  \:\__\       \:\  \    \:\ \/__/      |:|  |       \:\  \ 
   \/__/        \:\__\    \:\__\        |:|  |        \:\__\
                 \/__/     \/__/         \|__|         \/__/
                 
'''

parser = argparse.ArgumentParser(description=usage, epilog=epilog, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-v', '--invcf', help='Input VCF File', action='store')
parser.add_argument('-b', '--inbam', help='Input BAM File', action='store')
parser.add_argument('-l', '--libid', help='Sample Library ID (In accord with VCF FORMAT)', action='store')
parser.add_argument('-q', '--lib2qc', help='LibID to QC Info File', action='store')
parser.add_argument('-i', '--inssm', help='InsertSize Smoothing Instead of Filtering', default='on', action='store')
parser.add_argument('-c', '--rinsc', help='Ref InsertSize Correction By Alt', default='on', action='store')
parser.add_argument('-d', '--outdir', help='General Output Directory', default='./', action='store')
args = parser.parse_args()

if not (os.path.exists(args.outdir)):
	os.makedirs(args.outdir)

print("\nFeature Activation Status:\n")
print("# InsertSize Smoothing: %s" %(args.inssm))
print("# Ref InsertSize Correction By Alt: %s" %(args.rinsc))
print("\n")


# Globle Varibles

minIns = 20
maxIns = 600

dnaMolecule_10ng = 3333
avgInputdna = 23


# Function Define

def binarySearch(inlist, x):
	# Empty List
	if(not inlist or len(inlist) == 0):
		return -1
	# 1 element list
	# x lower than minimum
	if(len(inlist) == 1 or x < inlist[0]):
		return 0
	# x greater than maximum
	if(x > inlist[-1]):
		return len(inlist) - 1
	
	low = 0
	high = len(inlist) - 1
	
	while True:
		# if low <= high:
		if(high-low > 1):
			mid = int((low + high) / 2)
			if inlist[mid] == x:
				return mid
			elif x < inlist[mid]:
				# high = mid - 1
				high = mid
			else:
				# low = mid + 1
				low = mid
		else:
			return low



# Extract File Prefix and Dir

basename = os.path.basename(args.invcf)
basedir = os.path.abspath(args.invcf)
searchObj = re.search('(.*?)\.(.*)', basename)
prefix = searchObj.group(1)


# LibID to Fetal Fraction and Average Depth
lib2ff = {}
lib2avgdp = {}
with open(args.lib2qc, 'r') as lib2qcObj:
	for line in lib2qcObj:
		line = line.rstrip(os.linesep)
		lines = line.split("\t")
		lib2ff.setdefault(lines[0], lines[1])
		lib2avgdp.setdefault(lines[0], lines[2])
lib2qcObj.close()


insInfoFile = args.outdir + "/" + prefix + ".acd_fmid.xls";
outbam = args.outdir + "/" + prefix + ".addVariantInsInfo.bam";

inVcfObj = pysam.VariantFile(args.invcf, "r");
inBamObj = pysam.AlignmentFile(args.inbam, "rb")

insInfoFileObj = open(insInfoFile, 'w')
insInfoFileObj.write("VarString\tBetabinomLogCDF\t")
insInfoFileObj.write("RefDepth\tAltDepth\tIns_tTest\tIns_ksTest\tIns_hTest\tIns_uTest\tMinTestPvalue\t")
insInfoFileObj.write("RefInsMid\tAltInsMid\tInsMidDiff\n")
outBamObj = pysam.AlignmentFile(outbam, "wb", template=inBamObj)

# Variant Info Dictionary
varinfo = {}

for variant in inVcfObj.fetch():
	
	varSTR = "%s:%d:%s>%s" %(variant.chrom, variant.pos, variant.ref, variant.alts[0])
	varinfo[varSTR] = {}
	varinfo[varSTR].setdefault('refD', 0)
	varinfo[varSTR].setdefault('altD', 0)
	
	var_reads = inBamObj.fetch(contig=variant.chrom, start=variant.pos-1, stop=variant.pos)
	for varRead in var_reads:	
		for (readPos, refPos) in varRead.get_aligned_pairs():
			
			if(readPos == None):
				# print("\nNo readPos for Read: ", varRead.query_name, "\n")
				continue
			currentBase = varRead.query_sequence[readPos]
			
			if(refPos == variant.pos-1):
				
				insertsize = abs(varRead.template_length)
				# Abnormal InsertSize Correction
				if(args.inssm == 'on'):
					if(insertsize > maxIns):
						# print("\nSoomthing: %d => %d\n" %(insertsize, maxIns))
						insertsize = maxIns
					if(insertsize < minIns):
						# print("\nSoomthing: %d => %d\n" %(insertsize, minIns))
						insertsize = minIns
				
				# Ref Read
				if(currentBase.upper() == variant.ref.upper()):
					varinfoTag = ("vr", varSTR, "Z")
					if(insertsize >= minIns and insertsize <= maxIns):
						varinfo[varSTR].setdefault('refIns', []).append(insertsize)
					varinfo[varSTR]['refD'] = varinfo[varSTR]['refD'] + 1
			  # Alt Read
				elif(currentBase.upper() == variant.alts[0].upper()):
					varinfoTag = ("va", varSTR, "Z")
					if(insertsize >= minIns and insertsize <= maxIns):
						varinfo[varSTR].setdefault('altIns', []).append(insertsize)
					varinfo[varSTR]['altD'] = varinfo[varSTR]['altD'] + 1
			  # Other Alt Read
				else:
					continue
				
				# Update Tags and Print
				updatedTags = varRead.get_tags(with_value_type=True)
				updatedTags.append(varinfoTag)
				varRead.set_tags(updatedTags)
				outBamObj.write(varRead)
	
	if 'refIns' not in varinfo[varSTR]:
		varinfo[varSTR].setdefault('refIns', []).append(0)
	if 'altIns' not in varinfo[varSTR]:
		varinfo[varSTR].setdefault('altIns', []).append(0)
		
	altInsList = varinfo[varSTR]['altIns']
	refInsList = varinfo[varSTR]['refIns']
	
	# Drop Fetal Reads Insertsize for Ref Reads List by Fetal Insertsize Distribution
	if(args.rinsc == 'on' and len(refInsList)-len(altInsList) > 0):
		altInsList.sort()
		refInsList.sort()
		for altIns in altInsList:
			dropIdx = binarySearch(refInsList, altIns)
			refMaxIdx = len(refInsList)-1
			# print("\nDropIndex: %d/%d\n" %(dropIdx,refMaxIdx))
			if(dropIdx == refMaxIdx):
				# print("\nAltIns vs RefIns: %d vs %d, RefIdx:%d/%d\n" %(altIns,refInsList[dropIdx],dropIdx,refMaxIdx))
				del refInsList[dropIdx]
			elif(abs(refInsList[dropIdx]-altIns) > abs(refInsList[dropIdx+1]-altIns)):
				# print("\nAltIns vs RefIns: %d vs %d, RefIdx:%d/%d\n" %(altIns,refInsList[dropIdx+1],dropIdx,refMaxIdx))
				del refInsList[dropIdx+1]
			else:
				# print("\nAltIns vs RefIns: %d vs %d, RefIdx:%d/%d\n" %(altIns,refInsList[dropIdx],dropIdx,refMaxIdx))
				del refInsList[dropIdx]
	
	# Odd Number Median: (n1+n2)/2
	refInsMid = numpy.median(refInsList)
	altInsMid = numpy.median(altInsList)
	insMidDiff = refInsMid-altInsMid
	
	# Low Ref or Alt Support, Set Pvalue = -1
	if(len(altInsList) < 2 or len(refInsList) < 2):
		tTest_p = ksTest_p = hTest_p = uTest_p = minPvalue = -1
	else:
		# Welch's T-Test
		tTest_p = scipy.stats.ttest_ind(altInsList, refInsList, equal_var=False)[1]
		# Kolmogorov-Smirnov Test
		ksTest_p = scipy.stats.ks_2samp(altInsList, refInsList)[1]
		# Kruskal-Wallis H Test
		hTest_p = scipy.stats.kruskal(altInsList, refInsList)[1]
		# Mann-Whitney U Test
		uTest_p = scipy.stats.mannwhitneyu(altInsList, refInsList, alternative='two-sided')[1]
		
		minPvalue = min(tTest_p, ksTest_p, hTest_p, uTest_p)
	
	
	# Beta-Binomial, LOG CDF Calculation
	dnaMolecule = int(dnaMolecule_10ng*avgInputdna/10)
	
	varDepth = variant.info["DP"]
	fetalFraction = float(lib2ff[args.libid])
	
	nfactor = varDepth/float(lib2avgdp[args.libid])
	alpha = float(nfactor*dnaMolecule*fetalFraction/2)
	beta = float(nfactor*dnaMolecule-alpha)
	aaf = float(variant.samples[args.libid]["AF"][0])
	
	altDepth = int(varDepth*aaf)
	
	betabinom_logcdf = scipy.stats.betabinom.logcdf(altDepth, varDepth, alpha, beta)
	
	# Final Output
	insInfoFileObj.write("%s\t%.4f\t" %(varSTR, betabinom_logcdf))
	insInfoFileObj.write("%d\t%d\t" %(varinfo[varSTR]['refD'], varinfo[varSTR]['altD']))
	insInfoFileObj.write("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t" %(tTest_p, ksTest_p, hTest_p, uTest_p, minPvalue))
	insInfoFileObj.write("%.2f\t%.2f\t%.2f\n" %(refInsMid, altInsMid, insMidDiff))

inVcfObj.close()
inBamObj.close()
outBamObj.close()
insInfoFileObj.close()