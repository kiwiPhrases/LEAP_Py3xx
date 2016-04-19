#!/usr/bin/python

import random
import math
import os
import subprocess
random.seed()
from optparse import OptionParser
from argparse.ArgumentParser


import numpy as np
np.set_printoptions(precision=3, linewidth=200)

parser = ArgumentParser()
parser.add_argument()
parser.add_argument('--bfile', metavar='bfile', type=str, default = "bfile_test")
parser.add_argument('--csnps', metavar='csnps', type=int, default=500)
parser.add_argument('--fst', metavar='fst',  type=float, default=0)
parser.add_argument('--her', metavar='her', type=float, default=0.5)
parser.add_argument('--prev', metavar='prev',  type=float, default = 0.01)
parser.add_argument('--caseFrac', metavar='caseFrac', type=float, default=0.5)
parser.add_argument('--hidCoeff', metavar='hidCoeff', type=float, default=5.0)
parser.add_argument('--numSnps', metavar='numSnps', type=int, default=5000)
parser.add_argument('--numIndividuals', metavar='numIndividuals', type=int, default=2000)
parser.add_argument('--numDiff', metavar='numDiff', type=int, default=0)
parser.add_argument('--numSiblingsRatio', metavar='numSiblingsRatio', type=float, default=0)
parser.add_argument('--numHiddenSnps', metavar='numHiddenSnps', type=int, default=0)
parser.add_argument('--numLD', metavar='numLD', type=int, default=0)
parser.add_argument('--minFreq', metavar='minFreq', type=float, default=0.05)
parser.add_argument('--covarStd', metavar='covarStd', type=float, default=0)
parser.add_argument('--diffSize', metavar='diffSize', type=float, default=-1)
parser.add_argument('--numRare', metavar='numRare', type=int, default=0)
parser.add_argument('--rareFreq', metavar='rareFreq', type=float, default=0.01)
parser.add_argument('--geneSize', metavar='geneSize', type=int, default=10)
parser.add_argument('--effStd', metavar='effStd', type=float, default=1.0)
parser.add_argument('--numSibs', metavar='numSibs', type=int, default=0)
parser.add_argument('--diffCause', metavar='diffCause', type=int, default=1)
parser.add_argument('--mateSibs', metavar='mateSibs', type=int, default=1)
parser.add_argument('--numEpi', metavar='numEpi', type=int, default=0)
parser.add_argument('--epiEff', metavar='epiEff', type=float, default=5.0)
(options, args) = parser.parse_args()

"""
#parser = OptionParser()
parser.add_option('--bfile', metavar='bfile', type=str, default = "bfile_test")
parser.add_option('--csnps', metavar='csnps', type=int, default=500)
parser.add_option('--fst', metavar='fst',  type=float, default=0)
parser.add_option('--her', metavar='her', type=float, default=0.5)
parser.add_option('--prev', metavar='prev',  type=float, default = 0.01)
parser.add_option('--caseFrac', metavar='caseFrac', type=float, default=0.5)
parser.add_option('--hidCoeff', metavar='hidCoeff', type=float, default=5.0)
parser.add_option('--numSnps', metavar='numSnps', type=int, default=5000)
parser.add_option('--numIndividuals', metavar='numIndividuals', type=int, default=2000)
parser.add_option('--numDiff', metavar='numDiff', type=int, default=0)
parser.add_option('--numSiblingsRatio', metavar='numSiblingsRatio', type=float, default=0)
parser.add_option('--numHiddenSnps', metavar='numHiddenSnps', type=int, default=0)
parser.add_option('--numLD', metavar='numLD', type=int, default=0)
parser.add_option('--minFreq', metavar='minFreq', type=float, default=0.05)
parser.add_option('--covarStd', metavar='covarStd', type=float, default=0)
parser.add_option('--diffSize', metavar='diffSize', type=float, default=-1)
parser.add_option('--numRare', metavar='numRare', type=int, default=0)
parser.add_option('--rareFreq', metavar='rareFreq', type=float, default=0.01)
parser.add_option('--geneSize', metavar='geneSize', type=int, default=10)
parser.add_option('--effStd', metavar='effStd', type=float, default=1.0)
parser.add_option('--numSibs', metavar='numSibs', type=int, default=0)
parser.add_option('--diffCause', metavar='diffCause', type=int, default=1)
parser.add_option('--mateSibs', metavar='mateSibs', type=int, default=1)
parser.add_option('--numEpi', metavar='numEpi', type=int, default=0)
parser.add_option('--epiEff', metavar='epiEff', type=float, default=5.0)
(options, args) = parser.parse_args()
"""

sqrt12 = math.sqrt(12.0)
globalCovarCoeff = np.random.randn(1)[0] * options.covarStd
print( 'Covar coeff:', globalCovarCoeff)

bfile = options.bfile
numSnps = options.numSnps
numCausativeSNPs = options.csnps
numDiffSnps = options.numDiff
N = options.numIndividuals
casesFrac = options.caseFrac
hidCoeff = options.hidCoeff
numCases = int(N*casesFrac)
numControls = N - numCases
if (options.numSibs == 0): options.numSiblingsRatio = 0

fst = [options.fst, options.fst]
epsilon = [math.sqrt(1.0-options.her), math.sqrt(1.0-options.her)]
prevalence = options.prev


def runCommand(command, arguments, showOutput=True):
	# if (getOutput): p.StartInfo.RedirectStandardOutput = True
	#print 'Running:', command, arguments
	# p.StartInfo.FileName = command
	# p.StartInfo.Arguments  = arguments
	# p.Start()
	# if (getOutput): return p.StandardOutput.ReadToEnd()
	# p.WaitForExit()

	#print '\n'+command + " " + arguments+'\n'
	if showOutput: subprocess.call([command]+arguments.split())
	#else:
	#	with open(os.devnull, "w") as fnull:
	#		subprocess.call([command]+arguments.split(), stdout = fnull, stderr = fnull)


	
def createFrqFile(freqs, fileName, numCausativeSNPs, numHiddenSnps):
	f = open(fileName, 'w')	
	numCausativeSnpsCreated=0
	numSnpsCreated=0
	for i in range(len(freqs[0])):
		snpName = "snp"+str(i+1)
		if (numCausativeSnpsCreated < numCausativeSNPs-numHiddenSnps):
			snpName = 'c'+snpName
			numCausativeSnpsCreated += 1
			isCausative = True
		if (i >= numCausativeSNPs-numHiddenSnps and i < numCausativeSNPs):continue
		if (numSnpsCreated >= numSnps+numCausativeSNPs-numHiddenSnps): snpName = 'd'+snpName	
		f.write(' '.join([str(x) for x in [snpName, freqs[0][i], freqs[1][i], abs(freqs[0][i]-freqs[1][i])]])+'\n')
		numSnpsCreated+=1
	f.close()

def generateIndividualPhenotype(individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon, freqs, covars):	

	normalizedSNPs = (individualMafCounts[:numCausativeSNPs]-2*freqs[:numCausativeSNPs]) / np.sqrt(2*freqs[:numCausativeSNPs]*(1-freqs[:numCausativeSNPs]))
	if (options.numEpi == 0):
		sumBeta2_np = snpEffectSizes[:numCausativeSNPs].dot(snpEffectSizes[:numCausativeSNPs])
		
	#epistasis handling
	else:	
		sumBeta2_np = snpEffectSizes[(2*options.numEpi):numCausativeSNPs].dot(snpEffectSizes[(2*options.numEpi):numCausativeSNPs])
		for i in range(0, 2*options.numEpi, 2):
			#dominant model			
			f = 2*freqs[i:i+2]*(1-freqs[i:i+2]) + freqs[i:i+2]**2	#the probability for either one or two risk alleles in the two SNPs
			epiSNPFreq = f.prod()			
			epiSNPStd = np.sqrt(epiSNPFreq*(1-epiSNPFreq))
			epiSNPMean = epiSNPFreq #the mean is the frequency				
			normalizedSNPs[i] = ((1 if (individualMafCounts[i]>0 and individualMafCounts[i+1]>0) else 0) - epiSNPMean) / epiSNPStd
			normalizedSNPs[i+1] = 0
			sumBeta2_np += (snpEffectSizes[i]**2)		
			
			# #recessive model
			# f = freqs[i:i+2]**2	#the probability for two risk alleles in the two SNPs
			# epiSNPFreq = f.prod()
			# epiSNPStd = np.sqrt(epiSNPFreq*(1-epiSNPFreq))
			# epiSNPMean = epiSNPFreq #the mean is the frequency				
			# normalizedSNPs[i] = ((1 if (individualMafCounts[i]>1 and individualMafCounts[i+1]>1) else 0) - epiSNPMean) / epiSNPStd
			# normalizedSNPs[i+1] = 0
			# sumBeta2_np += (snpEffectSizes[i]**2)		
			
	
	liability_np = snpEffectSizes[:numCausativeSNPs].dot(normalizedSNPs)	
	if (options.covarStd > 0):
		liability_np += np.sum(covars * globalCovarCoeff)
		sumBeta2_np += globalCovarCoeff**2 * covars.shape[0]
		
	if (sumBeta2_np>0.0): liability_np *= np.sqrt((1.0-epsilon*epsilon) / sumBeta2_np)
	residual = random.gauss(0, epsilon)
	liability_np += residual
	return liability_np
	
def generateIndividualSnps(numSnps, mafs):
	#randomNumbers = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in mafs]
	#return randomNumbers
	return np.random.binomial(2, mafs)
	
	
def generateCovars():	
	return np.random.randn(1)
	
def createCovarsFile(covarsList, fileName):
	f = open(fileName, 'w')
	# f.write('FID IID ')
	# for i in range(len(covarsList[0])): f.write('covar'+str(i+1)+ ' ')
	# f.write('\n')
	for personCounter in range(len(covarsList)):
		line = 'FAM1 person'+str(personCounter+1) + ' ' + ' '.join([str(x) for x in covarsList[personCounter]])
		f.write(line+'\n')
	f.close()
	
def createConfoundersFile(snpsList, liabilitiesList, fileName, numCausativeSNPs, numHiddenSnps, popList):
	f = open(fileName, 'w')
	f.write('FID IID ')
	for i in range(numHiddenSnps): f.write('confounder'+str(i+1)+ ' ')
	f.write('\n')
	
	for personCounter in range(len(liabilitiesList)):
		line = 'FAM1 person'+str(personCounter+1) + ' ' + str(popList[personCounter])
		i = -1
		for maf in snpsList[personCounter]:
			i+=1
			#if (i>= (numCausativeSNPs-numHiddenSnps) and i<numCausativeSNPs): line += ' ' + str(maf)
				
		f.write(line+'\n')
	f.close()
	

	
def createPedFile(liabilitiesList, snpsList, fileName, numCausativeSNPs, numHiddenSnps, numLD):
	f = open(fileName, "w")	
	for personCounter in range(len(liabilitiesList)):
		line = 'FAM1 person'+str(personCounter+1)+' 0 0 1 '+str(liabilitiesList[personCounter])
		i = -1
		randomDraws = [random.random()<0.02 for r in range(2*numLD*len(snpsList[personCounter]))]
		drawIndex=0
		for maf in snpsList[personCounter]:
			i+=1
			if (i>= (numCausativeSNPs-numHiddenSnps) and i<numCausativeSNPs): continue
			if (maf==0): line += ' 1 1'
			elif (maf==1): line += ' 1 2'
			elif (maf==2): line += ' 2 2'

			#create LD markers
			if (i >= numCausativeSNPs or True):				
				for LDMArkerInd in range(numLD):
					if (maf==0): snp1, snp2 = 0, 0
					elif (maf==2): snp1, snp2 = 1, 1
					else: snp1, snp2 = 1, 0
					p1 = randomDraws[drawIndex]
					drawIndex+=1
					p2 = randomDraws[drawIndex]
					drawIndex+=1
					if p1: snp1=1-snp1
					if p2: snp2=1-snp2
					line += ' ' + str(snp1+1) + ' ' +str(snp2+1)
		f.write(line+'\n')
		
	f.close()
		

def createMapFile(numSnps, numDiffSnps, numCausativeSNPs, fileName, numHiddenSnps, numLD):
	f = open(fileName, 'w')
	
	numCausativeSnpsCreated = 0
	numSnpsCreated = 0
	
	chromosome = '1'
	pos = 0.00
	base = 0
	for i in range(numSnps+numCausativeSNPs+numDiffSnps):
		isCausative = False
		snpName = "snp"+str(i+1)
		if (numCausativeSnpsCreated < numCausativeSNPs):
			snpName = 'c'+snpName
			numCausativeSnpsCreated += 1
			isCausative = True
		if (numSnpsCreated >= numSnps+numCausativeSNPs): snpName = 'd'+snpName		
		lineArr = [chromosome, snpName, str(pos), str(base)]
		if (i< (numCausativeSNPs-numHiddenSnps) or i>=numCausativeSNPs):
			f.write(' '.join(lineArr)+'\n')		

			#create LD markers
			if (not isCausative or True):
				for LDMArkerInd in range(numLD):
					ldsnpName = snpName + '_'+str(LDMArkerInd+1)
					lineArr = [chromosome, ldsnpName, str(pos+0.01*(LDMArkerInd+1)), str(base+1+LDMArkerInd)]
					f.write(' '.join(lineArr)+'\n')
		
		
		pos += 5
		base += (numLD+10) 		
		numSnpsCreated += 1
	f.close()
	
def writeGenesFile(fileName, numSnps, numDiffSnps, numCausativeSNPs, chrom, posDiff):
	f = open(fileName, 'w')
	totalNumSnps = numSnps+numDiffSnps+numCausativeSNPs
	maxDist = ((totalNumSnps-1)//options.geneSize)*posDiff
	pos = posDiff-10
	while True:
		f.write(str(chrom) + ' ' + str(pos) + '\n')
		pos += posDiff
		if (pos >= maxDist):
			if (pos-posDiff < maxDist): f.write(str(chrom) + ' ' + str(pos) + '\n')			
			break
	f.close()
	
def createPhenotypeFile(liabilitiesList, filename, binary=False, threshold=0.0, addColumns=False):
	f = open(filename, 'w')
	for personCounter in range(len(liabilitiesList)):	
		phenotype = liabilitiesList[personCounter]
		if (binary): phenotype = [1,2][phenotype>=threshold]
		if (addColumns): lineArr = ["FAM1", "person"+str(personCounter+1), '0', '0', '1', str(phenotype)]
		else: 			 lineArr = ["FAM1", "person"+str(personCounter+1), str(phenotype)]
		f.write(' '.join(lineArr)+'\n')
	f.close()
	
	
	
def findThreshold(diseasePrevalence, sampleSize, numCausativeSNPs, snpEffectSizes, epsilon, freqs, realFreqs):
	liabilitiesArr = []	
	for i in range(sampleSize):
		if (i % 100000 == 0 and i>0): print('generating individual', i, 'out of', sampleSize	)
		individualMafCounts = generateIndividualSnps(numCausativeSNPs, freqs)
		covars = generateCovars()
		individualLiability = generateIndividualPhenotype(individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon, realFreqs, covars)
		liabilitiesArr.append(individualLiability)	
	liabilitiesArr.sort()
	prevIndex = int(round(sampleSize * (1-diseasePrevalence)))
	return liabilitiesArr[prevIndex]
	
	
def generateBrothers(numBrothers, freqs, popIndex, numCausativeSNPs):

	#fatherSnps = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
	#motherSnps = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]	
	fatherSnps = np.random.binomial(2, freqs[popIndex][numCausativeSNPs:])
	motherSnps = np.random.binomial(2, freqs[popIndex][numCausativeSNPs:])
	
	brotherSnps = []
	for i in range(numBrothers):
		# childFatherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in fatherSnps]
		# childMotherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in motherSnps]
		# brotherSnps.append(map(int.__add__, childFatherSnps, childMotherSnps))
		# brotherSnps.append(list(np.array(childFatherSnps) + np.array(childMotherSnps)))
		
		fatherOne = (fatherSnps==1)
		fatherSnps[fatherOne] = np.random.randint(0,2,size=np.sum(fatherOne))
		fatherSnps[fatherSnps==2] = 1
		motherOne = (motherSnps==1)
		motherSnps[motherOne] = np.random.randint(0,2,size=np.sum(motherOne))
		motherSnps[motherSnps==2] = 1
		brotherSnps.append(list(fatherSnps + motherSnps))
		
	return (fatherSnps, motherSnps, brotherSnps)

	

	
#Generate SNP effect sizes
if (options.effStd == 0.0): snpEffectSizes = [1.0 for x in range(numCausativeSNPs)]
else: snpEffectSizes = [random.gauss(0.0, options.effStd) for x in range(numCausativeSNPs)]

#Effect size for hidden SNPs
try:
	for i in range(numCausativeSNPs-options.numHiddenSnps, numCausativeSNPs):		
		if (options.effStd > 0.0): snpEffectSizes[i] = random.gauss(0.0, options.effStd*hidCoeff)
		else: snpEffectSizes[i] = hidCoeff
except: pass

#Effect size for epistatic SNPs
for i in range(2*options.numEpi):			
	snpEffectSizes[i] = options.epiEff

snpEffectSizes = np.array(snpEffectSizes)
print('first SNP effect sizes:', snpEffectSizes[:10])

#Generate allele frequencies for the old and new populations
ancestralFreqs = [random.random()*(0.5-options.minFreq) + options.minFreq for x in range(numSnps+numCausativeSNPs)]
for i in range(1,min(options.numRare+1, options.csnps-1)): ancestralFreqs[i] = options.rareFreq
freqs = [[] for x in range(len(fst))]
fstCoeffs = [(1.0-fst[popIndex]) / fst[popIndex] for popIndex in range(len(fst)) if fst[popIndex]>0]
for snpNum in range(numSnps+numCausativeSNPs):	
	ancientFreq = ancestralFreqs[snpNum]	
	for popIndex in range(len(fst)):	
		snpFreq = 0.0
		if (fst[popIndex] < 1e-20 or (options.diffCause == 0 and snpNum < numCausativeSNPs-options.numHiddenSnps)): snpFreq = ancientFreq
		else:	#if SNP is hidden and causal, or non-causal
			while (snpFreq < options.minFreq or snpFreq > 1-options.minFreq):
				#if (snpNum >= 1 and snpNum < min(options.numRare+1, options.csnps-1)): break #rare variants should remain as is
				snpFreq = random.betavariate(ancientFreq*fstCoeffs[popIndex], (1.0-ancientFreq)*fstCoeffs[popIndex])				
		freqs[popIndex].append(snpFreq)
		#if (snpNum < 5): print snpNum, snpFreq

#Generate allele frequencies for unusually differentiated SNPs
for snpNum in range(numDiffSnps):	
	freq1 = random.random()*0.4
	freq2 = freq1 + 0.6
	if (options.diffSize >= 0):
		freq1 = 0.5 - options.diffSize
		freq2 = 0.5 + options.diffSize
	diffFreqs = [freq1, freq2]
	for popIndex in range(len(fst)): freqs[popIndex].append(diffFreqs[popIndex])	

for popIndex in range(len(fst)): freqs[popIndex] = np.array(freqs[popIndex])
ancestralFreqs = np.array(ancestralFreqs)
	
	
#Find the affection threshold
threshold = 0
for popIndex in range(len(fst)):	
	popThreshold = findThreshold(prevalence, int(3000.0/prevalence), numCausativeSNPs, snpEffectSizes, epsilon[popIndex], freqs[popIndex][:numCausativeSNPs], ancestralFreqs[:numCausativeSNPs])
	threshold = max(threshold, popThreshold)
	
#Generate individuals
casesCounter = 0
controlsCounter = 0
popCasesSnps = []
popControlsSnps = []
popCasesLiabilities = []
popControlsLiabilities = []
popControlsCovars = []
popCasesCovars = []


for popIndex in range(len(fst)):
	popCasesSnps.append([])
	popControlsSnps.append([])
	popCasesLiabilities.append([])
	popControlsLiabilities.append([])
	popControlsCovars.append([])
	popCasesCovars.append([])
	
while (casesCounter < numCases or controlsCounter < numControls):
	for popIndex in range(len(fst)):
		individualMafCounts = generateIndividualSnps(numCausativeSNPs, freqs[popIndex][:numCausativeSNPs])
		
		#Add covariates		
		covars = generateCovars()
		
		individualLiability = generateIndividualPhenotype(individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon[popIndex], ancestralFreqs[:numCausativeSNPs], covars)
		if (controlsCounter < numControls and individualLiability < threshold):
			#print 'control', individualMafCounts[:2], individualLiability
			controlsCounter += 1
			popControlsSnps[popIndex].append(list(individualMafCounts))
			popControlsLiabilities[popIndex].append(individualLiability)
			popControlsCovars[popIndex].append(covars)
		elif (casesCounter < numCases and individualLiability >= threshold):
			#print 'case', individualMafCounts[:2], individualLiability
			casesCounter += 1
			popCasesSnps[popIndex].append(list(individualMafCounts))
			popCasesLiabilities[popIndex].append(individualLiability)
			popCasesCovars[popIndex].append(covars)
			
print( [len(popControlsLiabilities[0]), len(popControlsLiabilities[1])], [len(popCasesLiabilities[0]), len(popCasesLiabilities[1])])

#Find the subpopulation with the required number of cases
numCaseCaseSiblings = int(options.numSiblingsRatio * numCases)
numCaseControlSiblings = numCaseCaseSiblings
numControlControlSiblings = numCaseCaseSiblings
siblingPop = -1
maxPopCases = -1
for popIndex in range(len(fst)):
	if (len(popCasesLiabilities[popIndex]) > maxPopCases):
		siblingPop = popIndex
		maxPopCases = len(popCasesLiabilities[popIndex])
		
#Generate random alleles
print('Generating alleles...')
for popIndex in range(len(fst)):
	caseIndex = 0
	controlIndex = 0
	#Generate siblings
	if (options.mateSibs==0 and (popIndex == siblingPop or True)):
		while (numCaseCaseSiblings > 0):
			(fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
			for i in range(options.numSibs):
				popCasesSnps[popIndex][caseIndex] += sibSnps[i]
				caseIndex+=1
				numCaseCaseSiblings -= 1
		while (numControlControlSiblings > 0):
			(fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
			for i in range(options.numSibs):
				popControlsSnps[popIndex][controlIndex] += sibSnps[i]
				controlIndex+=1
				numControlControlSiblings -= 1
		while (numCaseControlSiblings > 0 and controlIndex < len(popControlsSnps[popIndex])):
			(fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
			for sibIndex in range(0, options.numSibs, 2):
				popCasesSnps[popIndex][caseIndex] += sibSnps[sibIndex]
				caseIndex+=1
				popControlsSnps[popIndex][controlIndex] += sibSnps[sibIndex+1]
				controlIndex+=1
				numCaseControlSiblings -= 2
				
	
	#Generate non-siblings controls
	while (controlIndex < len(popControlsSnps[popIndex])):
		#popControlsSnps[popIndex][controlIndex] += [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
		popControlsSnps[popIndex][controlIndex] += list(np.random.binomial(2, freqs[popIndex][numCausativeSNPs:]))
		controlIndex+=1
	#Generate non-siblings cases
	while (caseIndex < len(popCasesSnps[popIndex])):
		#popCasesSnps[popIndex][caseIndex] += [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
		popCasesSnps[popIndex][caseIndex] += list(np.random.binomial(2, freqs[popIndex][numCausativeSNPs:]))
		caseIndex+=1
		
		
#Generate mated sibs
if (options.mateSibs==1):
	numRemainingSibs = int(options.numSiblingsRatio * numCases * 3)
	while (numRemainingSibs > 0):
		parent1 = random.randint(0, len(popCasesSnps[siblingPop])+len(popControlsSnps[siblingPop])-1)
		parent2 = parent1
		while (parent1 == parent2): parent2 = random.randint(0, len(popCasesSnps[siblingPop])+len(popControlsSnps[siblingPop])-1)
		if (parent1 < len(popCasesSnps[siblingPop])): parent1SNPs = popCasesSnps[siblingPop][parent1]
		else: parent1SNPs = popControlsSnps[siblingPop][parent1-len(popCasesSnps[siblingPop])]
		if (parent2 < len(popCasesSnps[siblingPop])): parent2SNPs = popCasesSnps[siblingPop][parent2]
		else: parent2SNPs = popControlsSnps[siblingPop][parent2-len(popCasesSnps[siblingPop])]
		#parent1SNPs = list(parent1SNPs)
		#parent2SNPs = list(parent2SNPs)
		
		np.set_printoptions(precision=3, linewidth=200)
		parent1SNPs = np.array(parent1SNPs)
		parent2SNPs = np.array(parent2SNPs)
		
		
		for childNum in range(options.numSibs):
			# childFatherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in parent1SNPs]
			# childMotherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in parent2SNPs]
			# childSNPs = map(int.__add__, childFatherSnps, childMotherSnps)
			
			childFatherSnps = parent1SNPs.copy()
			childMotherSnps = parent2SNPs.copy()
			fatherOne = (childFatherSnps==1)
			childFatherSnps[fatherOne] = np.random.randint(0,2,size=np.sum(fatherOne))
			childFatherSnps[childFatherSnps==2] = 1
			motherOne = (childMotherSnps==1)
			childMotherSnps[motherOne] = np.random.randint(0,2,size=np.sum(motherOne))
			childMotherSnps[childMotherSnps==2] = 1
			childSNPs = childFatherSnps + childMotherSnps
			
			covars = generateCovars()
			childlLiability = generateIndividualPhenotype(np.array(childSNPs[:numCausativeSNPs]), numCausativeSNPs, snpEffectSizes, epsilon[siblingPop], ancestralFreqs[:numCausativeSNPs], covars)
			childSNPs = list(childSNPs)
			if (childlLiability < threshold):
				indToRemove = parent1
				while (indToRemove == parent1 or indToRemove == parent2):
					indToRemove = random.randint(0, len(popControlsSnps[siblingPop])-1)
				popControlsSnps[siblingPop][indToRemove] = childSNPs
				popControlsLiabilities[siblingPop][indToRemove] = childlLiability
				popControlsCovars[siblingPop][indToRemove] = covars
		
			else:
				indToRemove = parent1
				while (indToRemove == parent1 or indToRemove == parent2):
					indToRemove = random.randint(0, len(popCasesSnps[siblingPop])-1)
				popCasesSnps[siblingPop][indToRemove] = childSNPs
				popCasesLiabilities[siblingPop][indToRemove] = childlLiability
				popCasesCovars[siblingPop][indToRemove] = covars
				
			numRemainingSibs-=1
			
		#Regenerate parents SNPs
		#parent1NewSNPs = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[siblingPop][numCausativeSNPs:]]
		#parent2NewSNPs = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[siblingPop][numCausativeSNPs:]]		
		parent1NewSNPs = list(np.random.binomial(2, freqs[siblingPop][numCausativeSNPs:]))
		parent2NewSNPs = list(np.random.binomial(2, freqs[siblingPop][numCausativeSNPs:]))
		
		if (parent1 < len(popCasesSnps[siblingPop])): popCasesSnps[siblingPop][parent1][numCausativeSNPs:] = parent1NewSNPs
		else: popControlsSnps[siblingPop][parent1-len(popCasesSnps[siblingPop])][numCausativeSNPs:] = parent1NewSNPs
		if (parent2 < len(popCasesSnps[siblingPop])): popCasesSnps[siblingPop][parent2][numCausativeSNPs:] = parent2NewSNPs
		else: popControlsSnps[siblingPop][parent2-len(popCasesSnps[siblingPop])][numCausativeSNPs:] = parent2NewSNPs
		
		


#Create combined lists for all individuals
allLiabilities = []
allSnps = []
allCovars = []
popList = []
for popIndex in range(len(fst)):
	allLiabilities += popControlsLiabilities[popIndex]
	allLiabilities += popCasesLiabilities[popIndex]
	popList += [popIndex+1 for i in range(len(popControlsLiabilities[popIndex]) + len(popCasesLiabilities[popIndex]))]
	allSnps += popControlsSnps[popIndex]
	allSnps += popCasesSnps[popIndex]	
	allCovars += popControlsCovars[popIndex]
	allCovars += popCasesCovars[popIndex]
	

#Create Plink Files
print('Creating Plink files...')
createMapFile(numSnps, numDiffSnps, numCausativeSNPs, bfile+'.map', options.numHiddenSnps, options.numLD)
writeGenesFile(bfile+'.genes', numSnps, numDiffSnps, numCausativeSNPs, 1, 10*options.geneSize)
createPedFile(allLiabilities, allSnps, bfile+'.ped', numCausativeSNPs, options.numHiddenSnps, options.numLD)
createPhenotypeFile(allLiabilities, bfile+'.phe', binary=True, threshold=threshold)
createPhenotypeFile(allLiabilities, bfile+'.phe.liab', binary=False)
createPhenotypeFile(allLiabilities, bfile+'_pca.ped', binary=True, threshold=threshold, addColumns = True)
createCovarsFile(allCovars, bfile+'.cov')
createFrqFile(freqs, bfile+'.frq', options.csnps, options.numHiddenSnps)
createConfoundersFile(allSnps, allLiabilities, bfile+'.confounders', numCausativeSNPs, options.numHiddenSnps, popList)

print(bfile)
#args = ['C:/Users/SpiffyApple/Downloads/plink_win64/plink1_09.exe', '--noweb', '--silent',' --file ' + bfile ,' --make-bed',' --out ' + bfile]
#subprocess.call(args)

