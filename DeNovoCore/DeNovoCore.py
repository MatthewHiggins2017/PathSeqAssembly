import subprocess
import regex as re
import pandas as pd
import numpy as np
import json
import os
import sys
import random
from scipy.optimize import minimize
import urllib.request
import glob
import shutil


def CreateLoggingFiles(args):

	global CommandLogFile
	global ErrorLogFile
	global StdOutLogFile

	CommandLogFile = open(f'{args.Prefix}_Command_Log.txt','a')
	ErrorLogFile = open(f'{args.Prefix}_Error_Log.txt','a')
	StdOutLogFile = open(f'{args.Prefix}_StdOut_Log.txt','a')



# VALIDATE CHANGES TO CAPTURING OUTPUTS
def Run(Command):
	"""
	Run Command using subprocess.
	"""

	CommandLogFile.write(Command+'\n')
	CommandLogFile.flush()


	RO = subprocess.run(Command,
						shell=True,
						capture_output=True)


	ErrorLogFile.write(RO.stderr.decode("utf-8")+'\n')
	ErrorLogFile.flush()

	StdOutLogFile.write(RO.stdout.decode("utf-8")+'\n')
	StdOutLogFile.flush()

	return RO


def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False


def FastaToDict(InputFile):
	# Converts Fasta File to Dict
	fastadict = {}
	with open(InputFile) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				active_sequence_name = line[1:].split()[0]
				if active_sequence_name not in fastadict:
					fastadict[active_sequence_name] = ''
				continue
			sequence = line
			# If Uracil Present, Substitute for T (I.e. stick with 4 DNA bases)
			fastadict[active_sequence_name] += re.sub('U','T',sequence.upper())
	return fastadict


def BWAIndex(Ref):
	"""
	Index reference file
	"""
	if nofile(f"{Ref}.bwt"):
		Run(f"bwa index {Ref}")



def TrimShort(args):
	"""
	Denovo Trimming
	"""
	if nofile(f"{args.Prefix}_P1"):
		TrimC = f'trimmomatic PE -phred33 {args.ShortRead1} {args.ShortRead2} \
				 -threads {args.Threads} -baseout {args.Prefix} \
				 LEADING:3 TRAILING:3 \
				 SLIDINGWINDOW:4:20 MINLEN:36'

		Run(TrimC)




def ContamRemoval(args,ReadType='Short'):
	"""
	Decontamination Reads via Mapping
	"""

	ThreadSubset = int(args.Threads/2)


	# If short read i.e. Illumina.
	if ReadType == 'Short':

		# If file present already end function
		if nofile(f'{args.Prefix}_IlluminaNoContam_1.fq') == False:
			return True



		ContamCommand = f'bwa mem -t {ThreadSubset} {args.Contam} {args.ShortRead1} {args.ShortRead2} | samtools view -@ {ThreadSubset} -b -o  {args.Prefix}.ShortContam.bam ; \
		samtools sort -@ {ThreadSubset} -n  {args.Prefix}.ShortContam.bam | samtools view -@ {ThreadSubset} -F 2 -o {args.Prefix}.ShortContam.unmapped.bam ; \
		samtools fastq -N {args.Prefix}.ShortContam.unmapped.bam -1 {args.Prefix}_IlluminaNoContam_1.fq -2 {args.Prefix}_IlluminaNoContam_2.fq ;'

	# If long read i.e. Minion.
	else:

		# If file present already end function
		if nofile(f'{args.Prefix}_Minion_NoContam_1.fq') == False:
			return True

		ContamCommand = f'minimap2 -t {ThreadSubset} -ax map-ont {args.Contam} {args.LongReads} | samtools view -@ {ThreadSubset} -b -o {args.Prefix}.LongContam.bam ; \
		 samtools sort -@ {ThreadSubset} -n {args.Prefix}.LongContam.bam | samtools view -@ {ThreadSubset} -f 4 -o {args.Prefix}.LongContam.unmapped.bam ; \
		 samtools fastq -N {args.Prefix}.LongContam.unmapped.bam -0 {args.Prefix}_Minion_NoContam_1.fq ;'

	# Run command to remove contamination
	Run(ContamCommand)


'''
def FlyeAssembly(args):
	"""
	Long-Read Flye Denovo Assembly
	"""
	# Ensure file
	if nofile(f'{args.Prefix}_Flye/assembly.fasta'):

		FlyeC = f'flye  --nano-raw {args.LongReads} \
				--out-dir {args.Prefix}_Flye \
				--min-overlap {args.FlyeOverlap} \
				--scaffold \
				--threads {args.Threads}'

		# See if meta parameter is present
		if args.FlyeMeta:
			FlyeC+= ' --meta'

		print(FlyeC)
		Run(FlyeC)

	return True
'''



def PilonPolish(args,Iteration):

	"""
	Pilon Short-Read based Polishing
	"""

	# Check if Pilon Polish has already been done.
	if nofile(f'{args.Prefix}_Polish_{Iteration}/'+
			  f'{args.Prefix}_FlyePolish{Iteration}.fasta')== False:
		return True

	# Index Working Assembly
	BWAIndex(args.WorkingAssembly)

	ThreadSubset = int(args.Threads / 2)

	# Create Polishing Bam
	PilonPreAlign = f'bwa mem -t {ThreadSubset} {args.WorkingAssembly} {args.ShortRead1} {args.ShortRead2} | samtools view -@ {ThreadSubset} -b -o  {args.Prefix}_Polishing{Iteration}_Mapped.bam ; \
	samtools sort  -@ {ThreadSubset} {args.Prefix}_Polishing{Iteration}_Mapped.bam -o {args.Prefix}_Polishing{Iteration}_Mapped.sorted.bam ; \
	samtools index {args.Prefix}_Polishing{Iteration}_Mapped.sorted.bam ;'


	Run(PilonPreAlign)

	# Run Pilon
	PilonCommand = f'pilon --genome {args.WorkingAssembly} \
					--frags {args.Prefix}_Polishing{Iteration}_Mapped.sorted.bam \
					--output {args.Prefix}_Polish{Iteration} \
					--outdir ./{args.Prefix}_Polish_{Iteration}/ \
					--tracks \
					--verbose > {args.Prefix}_Pilon_Polish_{Iteration}.log 2>&1'

	Run(PilonCommand)



def SSPACEBasic(args,Iteration):

	"""
	SSPACE Basic Contig Extension
	"""


	# Run Picard Collect Insert Size Metrics
	# to get mean insert size for SSPACE-Basic.
	PicardCommandISM = f'picard CollectInsertSizeMetrics -I {args.WorkingBam} \
						-O {args.Prefix}_{Iteration}_Pilon1_ISM.txt \
						-H {args.Prefix}_{Iteration}_Pilon1_ISM.pdf'

	Run(PicardCommandISM)

	ISMDataframe = pd.read_csv(f'{args.Prefix}_{Iteration}_Pilon1_ISM.txt',
								comment='#',
								sep='\t')

	# Extract mean insert size and error
	MeanInsertSize = int(ISMDataframe.loc[0,'MEAN_INSERT_SIZE'])
	InsertSizeError = round(((ISMDataframe.loc[0,'STANDARD_DEVIATION'] * 2) \
					  /MeanInsertSize),2)

	# Create Test Library File change name
	LibFile = open(f'{args.Prefix}_{Iteration}_SSPACE_Library.txt','w')
	LibFile.write(f'Lib1\t{args.ShortRead1}\t{args.ShortRead2}\t' +
				  f'{MeanInsertSize}\t{InsertSizeError}\tFR')
	LibFile.close()


	ContigFile = args.WorkingAssembly

	## SSPACE basic iterations.
	for depth in list(range(args.MinDepth,
							args.MaxDepth+args.DepthStep,
							args.DepthStep))[::-1]:

		# Check if file exist already
		if os.path.isfile(f'{depth}_{args.Prefix}_{Iteration}_SSPACE' + \
						  f'_Depth_File.final.scaffolds.fasta'):
			continue


		# Run SSPACE Basic command
		SSPaceBasicRun =   f'''perl {args.SSPACEBasicPath} \
							-l {args.Prefix}_{Iteration}_SSPACE_Library.txt \
							-s {ContigFile} \
							-x 1 \
							-o {depth} \
							-T {args.Threads} \
							-b {depth}_{args.Prefix}_{Iteration}_SSPACE_Depth_File
							'''

		Run(SSPaceBasicRun)
		# Update the contig file
		ContigFile = f'{depth}_{args.Prefix}_{Iteration}_SSPACE' + \
					 f'_Depth_File.final.scaffolds.fasta'


	return (MeanInsertSize, InsertSizeError)



def UnmapReadExtraction(args):
	"""
	Extract Unmapped Reads
	"""

	BWAIndex(args.WorkingAssembly)

	ThreadsSubset = int(args.Threads/2)

	# Extract unmapped reads
	SpadesPriorCommand = f'bwa mem -t {ThreadsSubset} {args.WorkingAssembly} {args.ShortRead1} {args.ShortRead2} | samtools view -@ {ThreadsSubset} -b -o {args.Prefix}.SpadesPrior.bam ; \
	samtools sort -@ {ThreadsSubset} -n {args.Prefix}.SpadesPrior.bam | samtools view -@ {ThreadsSubset} -F 2 -o {args.Prefix}.SpadesPrior.unmapped.bam ; \
	samtools fastq -N {args.Prefix}.SpadesPrior.unmapped.bam -1 {args.Prefix}_IlluminaSpadesPrior_1.fq -2 {args.Prefix}_IlluminaSpadesPrior_2.fq ;'

	Run(SpadesPriorCommand)



def Spades(args):
	"""
	Short Read Denovo Assembly
	"""

	# Run Spades to generate additional contigs
	Spades = f'spades.py -1 {args.UnmappedReads1} -2 {args.UnmappedReads2} \
			 -t {args.Threads} -o {args.Prefix}_unmapped_illumina_spades_output \
			 > {args.Prefix}_Spades.log 2>&1'

	Run(Spades)



def RagTag(args):
	"""
	Reference Guided Scaffolding
	"""
	# Create directory to store ragtag output
	os.mkdir(f'{args.Prefix}_RagTag_Optimisation')

	# Define F and Q initial guesses
	x = np.array((args.RTMinUnique,
				  args.RTMinQual))

	# Define fixed parameters.

	p0=(args.GuideRef,
		args.WorkingAssembly,
		args.Prefix,
		args.Threads,
		False,
		args)

	# Run Parameter Optimisation
	RagTagOptimise = minimize(RagTagOptimisation,
							  x,
							  args=(p0),
							  method='Nelder-Mead',
							  bounds=[(0,np.Inf),(0,40)],
							  options={'maxiter':int(args.RTIterations)})



	# Run Final RagTag optimise with optimal F and Q
	# parameters
	RagTagOptimisation(RagTagOptimise.x,
					   args.GuideRef,
				   		args.WorkingAssembly,
				   		args.Prefix,
				   		args.Threads,
				   		True,
						args)






def RagTagOptimisation(Param,
					   RefScaffold,
					   ContigsFiles,
					   OutputPrefix,
					   Threads,
					   Final,
					   args):

	"""
	Ragtag function used in Optimisation
	"""


	F,Q = Param
	F = int(F)
	Q = int(Q)
	if Final == False:
		RC = random.randint(0, 10**10)
	else:
		RC = 'Final'

	RTComm = f'ragtag.py scaffold -u --debug -t {Threads} -f {F} -q {Q} \
			 {RefScaffold} {ContigsFiles} \
			 -o {OutputPrefix}_RagTag_Optimisation/{RC}_RagTag_Fixed_Gap'

	Run(RTComm)



	# Run contig renaming, this will add Unplaced
	ContigRenaming(args,f'{OutputPrefix}_RagTag_Optimisation/{RC}_RagTag_Fixed_Gap/ragtag.scaffold.fasta')



	# Run BUSCO Analysis
	BUSCOScore = BUSCOAssessment(f'{OutputPrefix}_RagTag_Optimisation/{RC}_RagTag_Fixed_Gap/ragtag.scaffold.fasta',
								 f'{OutputPrefix}_RagTag_Optimisation/{RC}_RagTag_Fixed_Gap/{OutputPrefix}',
								 args,
								 True)


	# Write BUSCO TO File
	BuscoScoreFile = open(f'{OutputPrefix}_RagTag_Optimisation/{RC}_RagTag_Fixed_Gap/Busco_Score.txt','w')
	BuscoScoreFile.write(str(BuscoScore))
	BuscoScoreFile.close()


	'''
	# Extract Unplaced BP count
	RagStat = pd.read_csv(f'{OutputPrefix}_RagTag_Optimisation/'+
						  f'{RC}_RagTag_Fixed_Gap/ragtag.scaffold.stats',
						  sep='\t')

	UnplacedBp = RagStat.loc[0,'unplaced_bp']
	'''

	return BUSCOScore



def ContigRenaming(args,FastaPath=False):
	"""
	Replace Fasta IDs
	"""

	# Load in Working Assembly to Dict and create file.
	if FastaPath==False:
		WD = FastaToDict(args.WorkingAssembly)
		OF = open(f'{args.Prefix}_Renamed.fasta','w')
	else:
		WD = FastaToDict(FastaPath)
		OF = open(FastaPath,'w')


	# Load in Json File, creating dictionary with replacement IDs
	with open(args.RenameJson) as json_file:
		RD = json.load(json_file)
		ReplaceD = RD['Replacements']


	# Clean working assembly (removing _RagTags suffix)
	# and then replace with corrected IDs and write to
	# a cleaned fasta file
	for I,S in WD.items():
		III = I.replace('_RagTag','')
		if III in list(ReplaceD.keys()):
			HO = ReplaceD[III]
		else:
			HO = 'Unplaced_'+III

		OF.write(f'>{HO}\n{S}\n')

	OF.close()


def Circularisation(args):
	"""
	Run Circularisation
	"""
	# Load in working scaffold to dictionary
	WD = FastaToDict(args.WorkingAssembly)

	# Extract contigs based on provided IDs and save
	# to tempoary fasta file
	CTF = open(f'{args.Prefix}_Circ_Temp.fasta','w')

	for ID in args.CircIDs.split(':'):
		ID = ID.strip()
		Seq = WD[ID]
		CTF.write(f'>{ID}\n{Seq}\n')
	CTF.close()

	# Run Circularisation with temp fasta and long reads
	CirC = f'circlator all  {args.Prefix}_Circ_Temp.fasta  \
			{args.LongReads}  {args.Prefix}_CirculatorOutput \
			--threads {args.Threads}'

	Run(CirC)


	# Read in Circularisation output and update
	# the fasta dictionary
	CS = FastaToDict(f'{args.Prefix}_CirculatorOutput/06.fixstart.fasta')
	for I,S in CS.items():
		WD[I]=S


	# Write updated dictionary to fasta file.
	FOW = open(f'{args.Prefix}_Renamed_Circ.fasta','w')
	for I,S in WD.items():
		FOW.write(f'>{I}\n{S}\n')
	FOW.close()



def GapFiller(args):
	"""
	Gap Filling using Short-reads
	"""

	# Create Library Guide File
	GFGuide = open(f'{args.Prefix}_GF_Library.txt','w')
	GFGuide.write(f'Lib1\tbowtie\t{args.ShortRead1}\t{args.ShortRead2}' + \
				  f'\t{args.MeanInsertSize}\t{args.InsertSizeError}\tFR')
	GFGuide.close()

	InputFile = args.WorkingAssembly


	# Download geoplts file required for gap filler
	# later find more elegant solution
	urllib.request.urlretrieve("https://raw.githubusercontent.com/AGCCAL/perl5/master/getopts.pl","getopts.pl")


	# Loop through viable depth range to close gaps where necessary.
	for depth in list(range(args.MinDepth,
							args.MaxDepth+args.DepthStep,
							args.DepthStep))[::-1]:

		# Check if file exist already
		if nofile(f'{args.Prefix}_{depth}_Extended_GapFiller/'+ \
				  f'{args.Prefix}_{depth}_' + \
				  'Extended_GapFiller.gapfilled.final.fa')==False:
			continue


		GFCommand = f'perl {args.GapFillerPath} \
					-l {args.Prefix}_GF_Library.txt -s {InputFile} \
					-m 30 -o {depth} -r 0.8 -n 10 -d 100000 \
					-t 10 -T {args.Threads} -i 10 \
					-b {args.Prefix}_{depth}_Extended_GapFiller'

		Run(GFCommand)

		# Update Input file
		InputFile = f'{args.Prefix}_{depth}_Extended_GapFiller/'+\
					f'{args.Prefix}_{depth}_' +\
					'Extended_GapFiller.gapfilled.final.fa'


	return True


def AbyssSealer(args):
	"""
	Sealer - Gap Closing
	"""
	if nofile(f'{args.Prefix}_Sealer_scaffold.fa')==False:
		return True

	ASCommand = f'abyss-sealer -b20G -k64 -k80 -k96 -k112 -k128 \
				-o {args.Prefix}_Sealer \
				-S {args.WorkingAssembly} {args.ShortRead1} {args.ShortRead2}'

	subprocess.run(ASCommand,shell=True)

	return True


def CombineFasta(FastaFileList,
				 OutPrefix):
	"""
	Concatenate Fasta Files
	"""

	AF = ' '.join(FastaFileList)
	CFC = f'cat {AF} > {OutPrefix}.fasta'
	Run(CFC)



def MissassemblyCovCorrection(args):


	# Index Wokring Assembly
	BWAIndex(args.WorkingAssembly)

	# Determine number of threads
	ThreadsSubset = int(args.Threads/2)

	# Generate coverage files (with samtools filter bwa sets multiple mapped
	# reads as having q = 0 so can use q=1 threshold as filter.)

	MCC = f'bwa mem -t {ThreadsSubset} {args.WorkingAssembly} {args.MisAssembRead1} {args.MisAssembRead2} | samtools view -@ {ThreadsSubset} -b -o {args.Prefix}.MissassemblyCheck.bam ; \
	samtools sort -@ {ThreadsSubset} {args.Prefix}.MissassemblyCheck.bam -o {args.Prefix}.MissassemblyCheck.sorted.bam ; \
	samtools depth -a {args.Prefix}.MissassemblyCheck.sorted.bam > {args.Prefix}.MissassemblyCheck.coverage ; \
	samtools view -F 2048 -b -q 1 -@ 18 {args.Prefix}.MissassemblyCheck.sorted.bam  -o {args.Prefix}.MissassemblyCheck_UniqMap.sorted.bam  ; \
	samtools depth -a {args.Prefix}.MissassemblyCheck_UniqMap.sorted.bam  > {args.Prefix}.MissassemblyCheck_UniqMap.coverage ; \
	'

	Run(MCC)


	# Loop through coverage files
	BadRegionsList = []

	# Chromosomes to exclude based on total coverage to skip.
	FullChromoToExclude = []


	# Boolean to define what coverage file to use
	# i.e. for contigs only to switch too

	ChromoOnlyBoolean = True


	for CF in [f'{args.Prefix}.MissassemblyCheck.coverage',
			   f'{args.Prefix}.MissassemblyCheck_UniqMap.coverage']:


		ComboDF = pd.read_csv(CF,
						  	  sep='\t',
						  	  header=None)

		ComboDF.columns = ['Chromo','Pos','Cov']


		# Filter based on chromo boolean = Later tidy this up not to
		# rely on scaffold name tag.
		if ChromoOnlyBoolean == True:
			ComboDF = ComboDF[ComboDF['Chromo'].str.contains('Unplaced')==False]
			ChromoOnlyBoolean = False
		else:
			ComboDF = ComboDF[ComboDF['Chromo'].str.contains('Unplaced')]


		# Check if full chromosome is below the coverage limit.
		# if coverage threshold is set as 0 then adapt and set it
		# too 1.
		if args.MisAssemCovThresh == 0:
			FCO = 1
		else:
			FCO = args.MisAssemCovThresh


		# Loop through Chromosome.
		for C in ComboDF['Chromo'].unique():

			# subset and extract all mean coverage metrics per 50bp
			SubCombo = ComboDF[ComboDF['Chromo']==C]
			SubMean = SubCombo['Cov'].tolist()
			MaxChromoSize = len(SubCombo)


			# Check if full chromosome is below the limit and if so add to list
			if SubCombo['Cov'].mean() <= FCO:
				FullChromoToExclude.append(C)


			# Complete assessment along chromosome according to user defined
			# windows
			else:

				# If no coverage in this place
				if os.path.isfile(f'{args.Prefix}_Bad_Region_Coverage.csv')==False:

					# Reset to 0 for every new chromosome.
					Tstart = 0
					Tend = 0

					# Loop through all positions in the chromosomes up to window
					# defined limit.
					for i in range(MaxChromoSize-args.MisAssemWindow):

						ROI = SubMean[i:i+args.MisAssemWindow]
						SMC = np.mean(ROI)

						# Check if region is below the limit.
						if SMC <= args.MisAssemCovThresh:

							# Calculate Start and End Site
							PS = i
							PE = i+args.MisAssemWindow

							# If first region in loop
							if Tstart == 0:
								Tstart = PS
								Tend = PE

							# If not first in loop but overlapping with existing
							# region then update TEND the end site.
							elif (PS > Tstart and PS < Tend):
								Tend = PE

							# If not overlapping then add existing TStart TEnd
							# to dataframe and reset Tstart and TEnd
							elif (Tstart != 0 and Tend != 0):

								BadRegionsList.append({'Chromo':C,
														'Start':Tstart,
														'End':Tend})

								Tstart = PS
								Tend = PE

					# Once I have looped through all windows add existing TStart and TEnd to
					# dataframe this may create duplicate entry but that can be correct for.
					if (Tstart != 0 and Tend != 0):

						BadRegionsList.append({'Chromo':C,
												'Start':Tstart,
												'End':Tend})

					# Reset Tstart and Tend as 0 for next iteration
					Tstart=0
					Tend=0



	# Convert into a DataFrame
	# or read in existing dataframe.
	if os.path.isfile(f'{args.Prefix}_Bad_Region_Coverage.csv')==False:
		BadRegionDF = pd.DataFrame(BadRegionsList)
		BadRegionDF = BadRegionDF.drop_duplicates()
		BadRegionDF.to_csv(f'{args.Prefix}_Bad_Region_Coverage.csv',index=None)
	else:
		BadRegionDF = pd.read_csv(f'{args.Prefix}_Bad_Region_Coverage.csv')



	# Load in reference file as dictionary
	RefDict = FastaToDict(args.WorkingAssembly)

	# Determine length of each chromosome or contig
	LenDict={}
	for a,b in RefDict.items():
		LenDict[a]=len(b)

	CutDF = BadRegionDF.sort_values(by=['Chromo','Start'])


	SplitUnplacedContigsDict = {}
	TOBEREMOVED = []

	# Loop through regions which need to be removed.
	for SubChromo in CutDF['Chromo'].unique().tolist():

		Sub = CutDF[CutDF['Chromo']==SubChromo]
		FinalEnd = len(RefDict[SubChromo])

		# Check if in nuclear genome as if so will simply look to replace
		# gaps with NNNN values.
		if 'Unplaced' not in SubChromo:

			# Generating new chromsome
			NewChromo = ''

			POI = list(zip(Sub['Start'].tolist(),
							Sub['End'].tolist()))

			TempStart=0
			for I in POI:
				if TempStart == 0:
					NewChromo+=RefDict[SubChromo][:I[0]-1]
				else:
					NewChromo+=RefDict[SubChromo][TempStart:I[0]-1]

				NewChromo+='N'*100
				TempStart = I[1]

			NewChromo+=RefDict[SubChromo][TempStart:FinalEnd]
			NewChromo = NewChromo.rstrip('N').lstrip('N')
			RefDict[SubChromo]=NewChromo

		# If Bad regions in unplaced contigs then need to split and add to
		# the SplitUnplacedContigsDict dictionary.
		else:

			TOBEREMOVED.append(SubChromo)


			# Loop over positions in Sub DF and split iplaced contig accordingly.
			SplitsContigsFasta = []
			FirstSub = True
			for i in range(len(Sub)):
				if FirstSub == True:
					FirstSub = False
					SplitsContigsFasta.append(RefDict[SubChromo][:Sub.iloc[i,1]])
				else:
					SplitsContigsFasta.append(RefDict[SubChromo][Sub.iloc[i,1]:Sub.iloc[i-1,2]])

			# Add in end.
			SplitsContigsFasta.append(WorkingAssemblyFasta[UC][SubDrop.iloc[-1,3]:])


			# Filter SplitsContigsFasta if length is greater than filter and add
			# to SplitContigsDict with unique ID.
			TCCID = 1
			SplitContigMinLen = 200
			for f in SplitsContigsFasta:
				if len(f)>= SplitContigMinLen:
					SplitUnplacedContigsDict[SubChromo+f'_Split_{str(TCCID)}'] = f
					TCCID+=1


	# Update dictionary to remove
	TCD = {a,b for a,b in RefDict.items() if a not in TOBEREMOVED}
	UpdatedRefDict = {**TCD,**SplitUnplacedContigsDict}

	# Write to file but exclude the full chromosomes necessary which do not
	# pass the full-chromosome coverage assessment.
	MisAssemOutID = f'{args.Prefix}_Misassembly_Cov_Filtered.fa'
	MisAssemOut = open(MisAssemOutID,'w')
	for a,b in UpdatedRefDict.items():
		if a not in FullChromoToExclude:
			MisAssemOut.write(f'>{a}\n{b}\n')
	MisAssemOut.close()


	# Write what chromosomes can be excluded
	DroppedChromoFile = open(f'{args.Prefix}_Chromo_Dropped.txt','w')
	for DC in FullChromoToExclude:
		DroppedChromoFile.write(f'{DC}\n')
	DroppedChromoFile.close()

	return True



def ExtractInsertMetrics(args,Iteration):

	"""
	Extract Pilon derived insert metrics
	"""

	ISMDataframe = pd.read_csv(f'{args.Prefix}_{Iteration}_Pilon1_ISM.txt',
								comment='#',
								sep='\t')

	# Extract mean insert size and error
	MeanInsertSize = int(ISMDataframe.loc[0,'MEAN_INSERT_SIZE'])
	InsertSizeError = round(((ISMDataframe.loc[0,'STANDARD_DEVIATION'] * 2) \
					  /MeanInsertSize),2)

	return (MeanInsertSize,InsertSizeError)



def FlyeOptimisation(args):
	"""
	Flye Optimisation - Long read assembly
	"""

	if 'FlyeMeta' in vars(args):
		FOM = 1
	else:
		FOM = 0


	# Initial guesses
	x = np.array((args.FlyeOverlap,
				  FOM))

	# Define fixed parameters.
	p0=(args.Threads,
		args.Prefix,
		args.LongReads,
		False,
		args)


	# Run Parameter Optimisation
	FlyeOptAss = minimize(FlyeAssembly,
						  x,
						  args=(p0),
						  method='Nelder-Mead',
						  bounds=[(1000,np.Inf),(0,2)],
						  options={'maxiter':int(args.FlyeIterations)})



	# Identify best Flye scaffolding
	MetricDf = pd.DataFrame()
	for FlyeFile in glob.glob(f'*_{args.Prefix}_Flye'):

		# Dictionary to store all metric values
		FFMD = {}

		# Extract Random ID
		RI = FlyeFile.split('_')[0]

		FFMD['ID']=RI

		# Identify and extract flye assembly metrics
		TFF = open(f'{FlyeFile}/flye.log','r')
		for TL in TFF.readlines():
			for KW in ['Total length:',
					   'Fragments:',
					   'Fragments N50:',
					   'Largest frg:',
					   'Scaffolds:',
					   'Mean coverage:']:

				if KW in TL:
					FFMD[KW.replace(':','')] = int(TL.split(':')[-1].replace(' ',''))



		# Extract BUSCO Score and  add to dictionary.
		BSF = open(f'{FlyeFile}/Busco_Score.txt','r')
		BSFOutStr = ''
		for B in BSF.readlines():
			BSFOutStr+=B
		FFMD['Busco Score'] = float(BSFOutStr.replace('\n',''))

		MetricDf = MetricDf.append(FFMD,ignore_index=True)


	# Sort Metrics DF (based on Busco score and make so smallest is at the
	# start)
	MetricDf = MetricDf.sort_values(by='Busco Score',ascending=True)
	MetricDf = MetricDf.reset_index(drop=True)
	MetricDf.to_csv(f'{args.Prefix}_Flye_Optimisation_Metrics.csv',index=None)


	# Extract top Random ID and Rename directory to be Final
	RIOI = MetricDf.loc[0,'ID']
	shutil.move(f'{RIOI}_{args.Prefix}_Flye',
				f'Final_{args.Prefix}_Flye')




def FlyeAssembly(OptimisationCombo,
				 Threads,
				 Prefix,
				 LongReads,
				 Final,
				 args):
	"""
	Long-Read Flye Denovo Assembly
	"""

	MinOverlap, MetaReview = OptimisationCombo

	MinOverlap = int(MinOverlap)

	if Final == False:
		RC = random.randint(0, 10**10)
	else:
		RC = 'Final'


	FlyeC = f'flye  --nano-raw {LongReads} \
			--out-dir {RC}_{Prefix}_Flye \
			--min-overlap {MinOverlap} \
			--scaffold \
			--threads {Threads}'

	# See if meta parameter is present
	if MetaReview >= 1:
		FlyeC+= ' --meta'

	Run(FlyeC)

	# Definte Flye Fasta Path
	FlyeOutputFasta = f'{RC}_{Prefix}_Flye/assembly.fasta'

	# Run BUSCO Command
	BuscoScore = BUSCOAssessment(FlyeOutputFasta,
								f'{RC}_{Prefix}_Flye/{Prefix}',
								 args)

	# Write BUSCO TO File
	FlyeBuscoScoreFile = open(f'{RC}_{Prefix}_Flye/Busco_Score.txt','w')
	FlyeBuscoScoreFile.write(str(BuscoScore))
	FlyeBuscoScoreFile.close()


	'''
	#####################################
	OLD Optimisation Delete When Ready
	#####################################
	FlyeLog = open(f'{RC}_{Prefix}_Flye/flye.log','r')
	for l in FlyeLog.readlines():
		if 'Total length' in l:
			TL = int(l.split(':')[-1].strip().replace(' ',''))

	FlyeScore = (10**15)-TL
	#####################################
	'''

	return BuscoScore




def RemoveSmallContigs(args):

	# Load in existing working assembly
	WAD = FastaToDict(args.WorkingAssembly)

	# Complete assessment of working assembly
	NWA = {}
	for i,s in WAD.items():
		if len(s) >= args.MinContigLen:
			NWA[i]=s

	# Write new assembly
	Out = open(f'{args.Prefix}_Post_Contig_Filt.fasta','w')
	for a,b in NWA.items():
		Out.write(f'>{a}\n{b}\n')
	Out.close()




def BUSCOAssessment(InputFasta,
					OutputPrefix,
					args,
					CoreGenomeOnly=False):

	BUSCO_Command = f'busco -i {InputFasta} -o {OutputPrefix}_Busco \
	 				  -m genome -l {args.BUSCOLineage} -c {args.Threads}'

	Run(BUSCO_Command)

	BuscoTable = pd.read_csv(f'{OutputPrefix}_Busco/run_{args.BUSCOLineage}/full_table.tsv',
							 sep='\t',
							 comment='#',
							 names=range(10))

	# If core genome is true remove any tags with scaffold in it. (Future
	# may change this to Unplaced but would need to modify post RagTag naming
	# etc do this later.)
	if CoreGenomeOnly == True:
		BuscoTable = BuscoTable[BuscoTable[2].str.contains('Unplaced')==False]


	# Extract BUSCO table Metrics
	BUSCOOutput = BuscoTable[1].value_counts().to_dict()

	# Determine total number
	with open(f'{OutputPrefix}_Busco/run_{args.BUSCOLineage}/short_summary.json') as json_file:
		SSJson = json.load(json_file)
		TotalBuscos = int(SSJson['dataset_total_buscos'])


	# Determine Score
	BUSCOOptimisationScore = ((TotalBuscos-BUSCOOutput['Complete'])*args.BUSCOWeight1) + \
							 (BUSCOOutput['Duplicated'] * args.BUSCOWeight2) + \
							 (BUSCOOutput['Missing'] * args.BUSCOWeight3) + \
							 (BUSCOOutput['Fragmented'] * args.BUSCOWeight4)


	return BUSCOOptimisationScore



def MisassemblyBuscoCorrection(args):

	# Run Standard BUSCO assessment command on WorkingAssembly.
	BUSCOAssessment(args.WorkingAssembly,
					f'{arg.Prefix}_BUSCO_Misassembly/{arg.Prefix}',
					args)

	# Analyses BuscoTable and subset so only contains Scaffold
	BuscoTable = pd.read_csv(f'{arg.Prefix}_BUSCO_Misassembly/{arg.Prefix}_Busco/run_{args.BUSCOLineage}/full_table.tsv',
							 sep='\t',
							 comment='#',
							 names=range(10))



	# Subset table to only contain duplicates
	DuplicateTable = BuscoTable[BuscoTable[1]=='Duplicated']

	# Store all sites which need to be droped
	DropTable = pd.DataFrame()

	# Loop through unique BIDs
	for BID in DuplicateTable[0].unique():

		# Subset based on BID
		TempSub = DuplicateTable[DuplicateTable[0]==BID]

		# 1) If duplicates are both in core genome then should keep both.
		# 2) If duplicates are split between core genome and unscaffolded, then
		# delete the scaffolded ones.
		# 3) If duplicates both in scaffolded keep the first occurance and remove
		# the rest.

		# Add column identifying if scaffold genome or not
		TempSub[10] = (TempSub[2].str.contains('Unplaced')==False)

		# Check if core genome duplicate
		DuplicateIDs = [i for i in TempSub.index.tolist() if i not in TempSub[[0,10]].drop_duplicates().index.tolist()]

		# Subset regions to drop and make sure that the none are in the core genome.
		# as dont want to remove these duplicates.
		PreDrop = TempSub.loc[DuplicateIDs,:]
		PreDrop = TempSub[TempSub[10]!=True]

		DropTable = pd.concat([DropTable,PreDrop])



	WorkingAssemblyFasta = FastaToDict(args.WorkingAssembly)


	# Loop through Drop table column 2 and use this to modify the fasta file.
	# 2 = chromosome, 3 = start, 4 = end (I assume starts from 1 but need to
	# check this)

	SplitContigsDict = {}


	for UC in DropTable[2].unique():

		# Subset and sort according to start so in order.
		SubDrop = DropTable[DropTable[2]==UC]
		SubDrop = SubDrop.sort_values(by=3)

		# Minus 1 from start and end site to 0-indexed for python instead
		# of 1-indexed.
		df[3] = df[3]-1
		df[4] = df[4]-1


		# Split contigs around BUSCO duplicate sites.
		SplitsContigsFasta = []
		FirstSub = True
		for i in range(len(SubDrop)):
			if FirstSub == True:
				FirstSub = False
				SplitsContigsFasta.append(WorkingAssemblyFasta[UC][:SubDrop.iloc[i,2]])
			else:
				SplitsContigsFasta.append(WorkingAssemblyFasta[UC][SubDrop.iloc[i,2]:SubDrop.iloc[i-1,3]])

		# Add in end.
		SplitsContigsFasta.append(WorkingAssemblyFasta[UC][SubDrop.iloc[-1,3]:])


		# Filter SplitsContigsFasta if length is greater than filter and add
		# to SplitContigsDict with unique ID.
		TCCID = 1
		SplitContigMinLen = 200
		for f in SplitsContigsFasta:
			if len(f)>= SplitContigMinLen:
				SplitContigsDict[UC+f'_Split_{str(TCCID)}'] = f
				TCCID+=1


	# Filting existing working assemby for unplaced contigs which have undergone
	# spliting and then add in new split contigs.
	TempFullDict = {a:b for a,b in WorkingAssemblyFasta.items() if a not in DropTable[2].unique()}
	NewWorkingAssembly = {**TempFullDict,**SplitContigsDict}

	# Write new Working Assembly
	OutWorking = open(f'{args.Prefix}_BUSCO_Corrected_Assembly.fa','w')
	for a,b in NewWorkingAssembly.items():
		OutWorking.write(f'>{a}\n{b}\n')
	OutWorking.close()

	return True
