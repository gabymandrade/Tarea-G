# TAREA OPCIONAL G

### SUMMARY OF THE PIPELINE AND SOFTWARE USED IN MY THESIS PROJECT


*Martínez Andrade Elsa Gabriela. Junio, 2020*


The pipeline I'm going to explain is about data cleansing to analyze samples of *Capsicum annuum* var. *glabriusculum* or "chiltepín", and how to download and index the *Capsicum* reference genome.

The samples were sequenced by **Illumina** High Throughput Sequencing, by **GBS** technique these samples are paired-end, and the enzime for the cutter site was ApeKI.


* **TRIMMOMATIC**

The fist step is use Trimmomatic to trim adapters.
Trimmomatic is bioinformatic program is fast, and have multithreaded command line tool that can be used to trim and crop Illumina(FASTQ) data as well as to remove adapters (Bolger *et al*., 2014).

First I decided to use the Single End Mode, because i am processing GBS data and not complete genomes. The parameters that i used for the program, i defined them with to the graphs of MultiQC, for the analysis i decided to occupy 8 threads (Bolger *et al*., 2014).

ILLUMINACLIP: TruSeq3-SE.fa:2:30:10 that used to find and remove Illumina adapter. This is done by Iientifying adapter or other contaminant sequences within a dataset is inherently a trade off between sensitivity and specificity (Bolger *et al*., 2014).

Finally MINLEN:70 , that removes reads that fall below the specified minimal length, i select this number when i see the graphs of my data. These were the only parameters because my reads had good quality (Bolger *et al*., 2014).

In this part, the parameters used in Trimmomatic must be reported in the results, which are in agreement, as I mentioned in the MultoQC charts. This is relevant since the parameters used depend on the quality of extraction and sequencing of the samples.


* **PROCESS RADTAGS**

It is a program implemented in Stacks. Process radtags examines raw reads from an Illumina sequencing, checks that the barcode and the RAD cutsite are intact, and demultiplexes the data (Catchen *et al*., 2011) and (Catchen *et al*., 2013).

I run process radtags whit Generic FASTQ Data, in paired-end mode. In this analysis use my fastq.gz files forward (-1) and reverse (-2) or the output of trimmomatic, --inline_null that means barcode is inline with sequence (Catchen *et al*., 2011) and (Catchen *et al*., 2013).

The next -o or the path to the output what is a directory, -b path to a file containing barcodes for this run, -e or provide the restriction enzyme used in my case apeKI (Catchen *et al*., 2011) and (Catchen *et al*., 2013).

By last the flags -r rescue barcodes and RAD-Tags, -c clean data, remove any read with an uncalled base and -q discard reads with low quality scores (Catchen *et al*., 2011) and (Catchen *et al*., 2013).

In the case of this analysis, the .log file is very important, since we can graph the information provided in R. In this case, the number of readings per sample.

The graph is one of the most important results, since it will show us the quality of the readings per sample. We can also decide which sequences can be discarded by having few reads.

In case you have a blank sample, we would rather perceive few or no reads. This is also important since it tells us if there was contamination in any part of the process, up to sequencing.


* **FASTQC**

FastQC check output of process radtags, therefore is a quality control tool for high throughput sequence data. FastQC also give a quick impression of whether your data has any problems in further analysis (Andrews, 2010).

The HTML report contains important information like Per Base Sequence Quality that shows an overview of the range of quality values across all  bases at each position in the FastQ file (Andrews, 2010).

Per Base Sequence Content of proportion of each base position in a file for which each of the four normal DNA bases has been called. Overrepresented Kmers, this willspot an increase in any exactly duplicated sequences, but there are a different subset of problems where it will not work (Andrews, 2010).

These graphics with some of the most important since they tell us about the quality of sequencing. This can be seen with the graph of Per Base Sequence Content.

Per Base Sequence Quality, provides us with a visual form of quality in Phred Score for samples. This value is of great importance in determining the quality of sequencing.

Therefore it is important to mention these results in our work to demonstrate the quality of the samples.





* **DOWNLOAD THE REFERENCE GENOME**

First we have to choose our reference genome, in my case is: Reference quality assembly of the 3.5GB genome of *Capsicum annuum* from a single linked-read library. 2018. Amanda M. Hulse-Kemp *et al*. University
of California, Davis.

Since the genome is in NCBI, the code is PRJNA376668. Then follow the next steps: go to access, Assembly details, Assembly, access the data, FTP directory for GenBank assembly, .fna.gz (for complete genome) right click, copy link address. 

Next in server `wget` and the link obtained from NCBI.

In the methods of my work, I must emphasize the importance of mentioning where I obtained the reference genome. This will determine the quality of the alignment against the reference genome and will have a great impact on the subsequent analyzes.


* **BWA-INDEX GENOME**

BWA is a software package for mapping low-divergent sequences against a large reference genome (Li and Durbin, 2009).

This part it is necessary to construct the FM-index for the reference genome, what is necessary to run BWA.

Also we need the path from bin to the reference genome, index database sequences in the FASTA format and and the path to the output.

We have to run  BWA with bwtsw algoritm, because *Capsicum* has a large genome. Therefore apply the flag -a bwtsw, and the flag -p, or prefix of the output database.

The indexing part is of great importance since they cannot run the BWA program. In addition to maximizing the effectiveness of alignment with BWA.

The parameters that we put in the scrip must be properly selected to avoid errors when alienating against the reference genome.



**REFERENCES**

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

J. Catchen, P. Hohenlohe, S. Bassham, A. Amores, and W. Cresko. Stacks: an analysis tool set for population genomics. Molecular Ecology. 2013. [reprint] 

J. Catchen, A. Amores, P. Hohenlohe, W. Cresko, and J. Postlethwait. Stacks: building and genotyping loci de novo from short-read sequences. G3: Genes, Genomes, Genetics, 1:171-182, 2011. [reprint]

Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] 