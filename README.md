# Pathogen Hybrid Assembly Pipeline

**Description**: Combining short read ( Illumina ) and long-read ( Minion | Pacbio ) sequencing data for de-novo hybrid assembly of novel pathogens.


## Steps
* Short read QC & trimming
* Decontamination **[Optional]**
* Long read assembly via Flye
* Error correction with Pilon
* Contig Extension with SSPACE
* Unmapped Short read extraction
* Short Read Assembly via SPADES
* Contig Extension with SSPACE
* Error correction with Pilon
* Scaffolding with RagTag **[Optional]**
* Gap closure with Gapfiller
* Gap closure with Sealer
* Error correction via Pilon
* Missassmbly check and correction **[Optional]**
* Assembly scoring **[Optional]**


Some steps are run iteratively such as contig extension and gap filling, decreasing stringency with each round following best-practise established by Sanger sequencing. Use the --help option to identify purpose of each parameter 
