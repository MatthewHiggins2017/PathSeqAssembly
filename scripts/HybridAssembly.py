import os
import sys
import argparse
import DeNovoCore as dnc


# Pipeline Stage Gauge
def get_step_num(args):

	files = {f"{args.Prefix}_1P":1,
			 f"{args.Prefix}_IlluminaNoContam_1.fq":2,
			 f"{args.Prefix}_Minion_NoContam_1.fq":2.5,
			 f'Final_{args.Prefix}_Flye/assembly.fasta':3,
			 f'{args.Prefix}_Polish_1/' + \
			 f'{args.Prefix}_Polish1.fasta':4,
			 f'{args.MinDepth}_{args.Prefix}_1' + \
			 '_SSPACE_Depth_File.final.scaffolds.fasta':5,
			 f'{args.Prefix}_IlluminaSpadesPrior_1.fq':6,
			 f'{args.Prefix}_unmapped_illumina_spades' + \
			 '_output/scaffolds.fasta':7,
			 f'{args.Prefix}_Combo_contigs.fasta':8,
			 f'{args.Prefix}_Polish_2/' +
			 f'{args.Prefix}_Polish2.fasta':9,
			 f'{args.MinDepth}_{args.Prefix}_2' + \
			 '_SSPACE_Depth_File.final.scaffolds.fasta':10,
			 f'{args.Prefix}_RagTag_Optimisation/' +
			 'Final_RagTag_Fixed_Gap/ragtag.scaffold.fasta':11,
			 f'{args.Prefix}_Post_Contig_Filt.fasta':11.5,
			 f'{args.Prefix}_Renamed.fasta':12,
			 f'{args.Prefix}_Renamed_Circ.fasta':13,
			 f'{args.Prefix}_{args.MinDepth}_Extended_'+ \
			 f'GapFiller/{args.Prefix}_{args.MinDepth}_' + \
			 'Extended_GapFiller.gapfilled.final.fa':14,
			 f'{args.Prefix}_Sealer_scaffold.fa':15,
			 f'{args.Prefix}_Polish_3/' + \
			 f'{args.Prefix}_Polish3.fasta':16,
			 f'{args.Prefix}_BUSCO_Corrected_Assembly.fa':17,
			 f'{args.Prefix}_Misassembly_'+ \
	 		 'Cov_Filtered.fa':18}

	step = 0
	WorkingFile = False
	for f in files:
		if os.path.isfile(f):
			step = files[f]
			WorkingFile = f
			sys.stderr.write(f"Found {f}\n")

	return (step,WorkingFile)




# Controller Function
def DenovoControl(args):

	# Create Log files
	dnc.CreateLoggingFiles(args)

	# If read decontamination necessary, ensure
	# reference is indexed and class updated
	# to ensure optional sections run.
	if args.Contam != 'False':
		dnc.BWAIndex(args.Contam)
		args.decontaminate = True
	else:
		args.decontaminate = False


	#  Keep track of pipeline stage and working
	# assembly.
	args.step, args.WorkingAssembly = get_step_num(args)

	# If step prior to 2, then reset working assembly
	# as nothing.
	if args.step < 2:
		args.WorkingAssembly = ''

	# Trimming
	if args.redo or args.step<1:
		dnc.TrimShort(args)

	args.ShortRead1 = f"{args.Prefix}_1P"
	args.ShortRead2 = f"{args.Prefix}_2P"


	# Decontamination [OPTIONAL]
	if args.decontaminate == True:
		if args.redo or args.step<2:
			dnc.ContamRemoval(args)
		args.ShortRead1 = f"{args.Prefix}_IlluminaNoContam_1.fq"
		args.ShortRead2 = f"{args.Prefix}_IlluminaNoContam_2.fq"

		if args.redo or args.step<2.5:
			dnc.ContamRemoval(args,"Long")
		args.LongReads = f"{args.Prefix}_Minion_NoContam_1.fq"


	# Long Read Assembly
	if args.redo or args.step<3:
		dnc.FlyeOptimisation(args)

	args.WorkingAssembly = f'Final_{args.Prefix}_Flye/assembly.fasta'


	# Pilon Polishing
	Iteration=1
	if args.redo or args.step<4:

		# Mark first iteration of pilon polishing
		dnc.PilonPolish(args,Iteration)

	args.WorkingAssembly =  f'{args.Prefix}_Polish_{Iteration}/' + \
							f'{args.Prefix}_Polish{Iteration}.fasta'

	args.WorkingBam = f'{args.Prefix}_Polishing{Iteration}_Mapped.sorted.bam'


	# SSPACE Extension
	Iteration=1
	if args.redo or args.step<5:

		args.MeanInsertSize, args.InsertSizeError = dnc.SSPACEBasic(args,
																	Iteration)

	args.WorkingAssembly = f'{args.MinDepth}_{args.Prefix}_{Iteration}' + \
							'_SSPACE_Depth_File.final.scaffolds.fasta'


	if 'MeanInsertSize' not in vars(args):
		args.MeanInsertSize, args.InsertSizeError = dnc.ExtractInsertMetrics(args,
																		 Iteration)




	# Unmapped Read Extraction
	if args.redo or args.step<6:
		dnc.UnmapReadExtraction(args)

	args.UnmappedReads1 = f'{args.Prefix}_IlluminaSpadesPrior_1.fq'
	args.UnmappedReads2 = f'{args.Prefix}_IlluminaSpadesPrior_2.fq'

	# Spades contig assembly
	if args.redo or args.step<7:
		dnc.Spades(args)

	args.AddOnAssembly = f'{args.Prefix}_unmapped_illumina_spades' + \
						  '_output/scaffolds.fasta'


	# COMBINE CONTIGS
	if args.redo or args.step<8:
		dnc.CombineFasta([args.WorkingAssembly,
						  args.AddOnAssembly],
						 f'{args.Prefix}_Combo_contigs')

	args.WorkingAssembly = f'{args.Prefix}_Combo_contigs.fasta'

	# Pilon Polishing (Round 2)
	Iteration=2
	if args.redo or args.step<9:
		dnc.PilonPolish(args,Iteration)

	args.WorkingAssembly =  f'{args.Prefix}_Polish_{Iteration}/' + \
								f'{args.Prefix}_Polish{Iteration}.fasta'

	args.WorkingBam = f'{args.Prefix}_Polishing{Iteration}_Mapped.sorted.bam'


	# SSPACE Extension (Round 2)
	Iteration=2
	if args.redo or args.step<10:
		args.MeanInsertSize, args.InsertSizeError = dnc.SSPACEBasic(args,
																	Iteration)

	args.WorkingAssembly = f'{args.MinDepth}_{args.Prefix}_{Iteration}' + \
							'_SSPACE_Depth_File.final.scaffolds.fasta'


	if 'MeanInsertSize' not in vars(args):
		args.MeanInsertSize, args.InsertSizeError = dnc.ExtractInsertMetrics(args,
																		 Iteration)

	# RagTag Scaffolding [OPTIONAL]
	if args.GuideRef != False:
		if args.redo or args.step<11:
			dnc.RagTag(args)

	args.WorkingAssembly = f'{args.Prefix}_RagTag_Optimisation/' + \
							'Final_RagTag_Fixed_Gap/ragtag.scaffold.fasta'


	# Remove small contigs
	if args.redo or args.step<11.5:
		dnc.RemoveSmallContigs(args)


	args.WorkingAssembly = f'{args.Prefix}_Post_Contig_Filt.fasta'


	# Contig Name Conversion [OPTIONAL]
	if args.RenameJson != False:
		if args.redo or args.step<12:
			dnc.ContigRenaming(args)

	args.WorkingAssembly = f'{args.Prefix}_Renamed.fasta'


	# Circularisation [OPTIONAL]
	if args.CircIDs != False:
		if args.redo or args.step<13:
			dnc.Circularisation(args)

	args.WorkingAssembly = f'{args.Prefix}_Renamed_Circ.fasta'


	# GapFiller
	if args.redo or args.step<14:
		dnc.GapFiller(args)


	args.WorkingAssembly = f'{args.Prefix}_{args.MinDepth}_Extended_'+ \
						   f'GapFiller/{args.Prefix}_{args.MinDepth}_' + \
						   'Extended_GapFiller.gapfilled.final.fa'


	# Abyss Sealer
	if args.redo or args.step<15:
		dnc.AbyssSealer(args)

	args.WorkingAssembly = f'{args.Prefix}_Sealer_scaffold.fa'


	# Pilon Polishing
	Iteration=3
	if args.redo or args.step<16:
		dnc.PilonPolish(args,Iteration)

	args.WorkingAssembly =  f'{args.Prefix}_Polish_{Iteration}/' + \
							f'{args.Prefix}_Polish{Iteration}.fasta'


	################ NEW NEW NEW ################
	# Busco Based Missassembly Correction.
	if args.redo or args.step<17:
		dnc.MisassemblyBuscoCorrection(args)

	args.WorkingAssembly = f'{args.Prefix}_BUSCO_Corrected_Assembly.fa'
	################ NEW NEW NEW ################


	# Missassembly Correction [OPTIONAL]
	if args.MisAssembRead1 != 'False':
		# Coverage based correction
		if args.redo or args.step<17:
			dnc.MissassemblyCovCorrection(args)

		args.WorkingAssembly = f'{args.Prefix}_Misassembly_'+ \
								'Cov_Filtered.fa'


	# Score Assembly [Optional]
	if args.Score:
		AssemblyScore = dnc.Assembly(args.Prefix,
									args.WorkingAssembly,
									log=True)


	return True

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

parser = argparse.ArgumentParser(description='Hybrid Denovo Assembly Pipeline')
subparsers = parser.add_subparsers(help="Options")

parser_sub = subparsers.add_parser('all',
									help='Full Denovo Assembly')


parser_sub.add_argument('--Prefix','-P',
						help='Prefix for all files generated',
						required=True)


parser_sub.add_argument('--redo',
						action="store_true",
						help='Redo everything')

parser_sub.add_argument('--ShortRead1',
						'-1',
						help='First illumina read file',
						type=str)

parser_sub.add_argument('--ShortRead2',
						'-2',
						help='Second illumina read file',
						type=str)

parser_sub.add_argument('--LongReads',
						'-3',
						help='Long read file (minion)',
						type=str)

parser_sub.add_argument('--Threads',
						'-T',
						default=4,
						type=int,
						help='Number of threads')

parser_sub.add_argument('--FlyeOverlap',
						default=1000,
						type=int,
						help='Minimum overlap parameter for Fyle assembly.')

parser_sub.add_argument('--FlyeMeta',
						action="store_true",
						help='Flye Meta parameter for long-read assembly. May \
							  be required for samples which have undergone \
							  SWGA due to intrinsic uneven coverage of the \
							  target genome.')


parser_sub.add_argument('--FlyeIterations',
						default=6,
						type=int,
						help='Number of iterations to optimise Flye.')


parser_sub.add_argument('--RTIterations',
						default=20,
						type=int,
						help='Number of iterations to optimise RagTag \
							  scaffolding over.')


parser_sub.add_argument('--RTMinUnique',
						default=1000,
						help='RagTag minimum unique bp strech (Q)')


parser_sub.add_argument('--RTMinQual',
						default=10,
						help='RagTag minimum alignment quality')


parser_sub.add_argument('--SpadesMinCov',
					help='Spades Minimum Coverage',
					default=2,
					type=int)

parser_sub.add_argument('--MinDepth',
					help='Min Depth for SSPACE and GapFiller',
					default=4,
					type=int)

parser_sub.add_argument('--MaxDepth',
					help='Max Depth for SSPACE and GapFiller',
					default=30,
					type=int)

parser_sub.add_argument('--DepthStep',
					help='Depth Step for SSPACE & GapFiller iterations',
					default=2,
					type=int)

parser_sub.add_argument('--Contam',
						'-C',
						help='Path to contamination reference fasta',
						default='False',
						type=str)


parser_sub.add_argument('--GuideRef','-R',
						help='Guide Reference',
						default=False,
						type=str)

parser_sub.add_argument('--RenameJson',
						help='Path to renaming json file',
						default=False,
						type=str)


parser_sub.add_argument('--CircIDs',
						help='Fasta IDs, seperated by (:)',
						default=False,
						type=str)

parser_sub.add_argument('--SSPACEBasicPath',
						help='Path to SSPACE Basic Script',
						required=True,
						type=str)


parser_sub.add_argument('--GapFillerPath',
						help='Path to GapFiller Script',
						required=True,
						type=str)

parser_sub.add_argument('--Score',
						action="store_true",
						help='Score assembly generated and write output file')


parser_sub.add_argument('--MisAssembRead1',
						help='First read file for misassembly correction',
						default='False',
						type=str)

parser_sub.add_argument('--MisAssembRead2',
						help='Second read file for misassembly correction',
						type=str)

parser_sub.add_argument('--MisAssemCovThresh',
						help='Misassembly Coverage Threshold',
						type=float,
						default=0)

parser_sub.add_argument('--MisAssemWindow',
						help='Misassembly Window e.g. 100bp',
						type=int,
						default=100)

parser_sub.add_argument('--MinContigLen',
						help='Minimum contig length to keep post scaffolding',
						type=int,
						default=1000)



parser_sub.add_argument('--BUSCOLineage',
						help='Minimum contig length to keep post scaffolding',
						type=str,
						default='auto_lineage')


parser_sub.add_argument('--BUSCOWeight1',
						help='BUSCO Weight1 for optimisation algorithm',
						type=int,
						default=1)


parser_sub.add_argument('--BUSCOWeight2',
						help='BUSCO Weight1 for optimisation algorithm',
						type=int,
						default=1)

parser_sub.add_argument('--BUSCOWeight3',
						help='BUSCO Weight1 for optimisation algorithm',
						type=int,
						default=1)

parser_sub.add_argument('--BUSCOWeight4',
						help='BUSCO Weight1 for optimisation algorithm',
						type=int,
						default=1)


parser_sub.set_defaults(func=DenovoControl)

################################################################################
################################################################################
################################################################################

# If no options specified add help
if len(sys.argv) <= 1:
	sys.argv.append('--help')

args = parser.parse_args()
args.func(args)
