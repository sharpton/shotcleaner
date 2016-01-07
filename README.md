shotcleaner
=======

A high-throughput and modular workflow to quality control shotgun metagenomic DNA sequence libraries

Version: 0.1

Overview
--------

Shotcleaner is an automated workflow that takes host-associated shotgun metagenomes and prepares them for 
analysis by conducting the following quality control measures:

1. Trim and filter low quality sequences
2. Remove host DNA
3. Collapse identical duplicate sequences

Shotcleaner can work with both single end and paired end data and provides access to several methodological
options at most quality control steps listed above, though the default methods are recommended. This software
was designed to be run on a multi-core, UNIX OS machine.

Quickstart
----------

**(1) Build a host genome index:**

    perl shotcleaner/index_genome -i <host_genome.fasta> -d <index_database> -n <index_basename> -m <filter_method>

For example:

    perl shotcleaner/index_genome -i /data/genomes/hg19.fa -d ~/projects/index_database -n h_sapien -m bowtie2

You only need to do this once per host genome per filter method

**(2) Execute shotcleaner:**

For paired end data:

     perl shotcleaner/shotcleaner.pl -1 <forward_reads.fq.gz> -2 <reverse_reads.fq.gz> -o <output_directory> -d <index_database> -n <index_basename> --nprocs <number_of_processors>

For singled end data:

     perl shotcleaner/shotcleaner.pl -1 <forward_reads.fq.gz> -o <output_directory> -d <index_database> -n <index_basename> --nprocs <number_of_processors>

shotcleaner will produce a directory as specified by -o that contains the results of each step of the workflow, with the final output
being a cleaned fasta file in <output_directory/fasta_cleaned>. 

INSTALLATION
------------

You can use the auto-installer, install.pl, which attempts to install all dependencies for shotcleaner (see Dependencies). Note that this may not work on all systems (e.g., it assumes an x86 architecture), which means dependencies will need to be installed by hand.

 To implement this installer:

cd shotcleaner/
perl install.pl &> install.log

This may take some time. 

Options
-------

### index_genome.pl ###

* **-i, /PATH/TO/HOST/GENOME/FASTA** (REQUIRED argument) NO DEFAULT VALUE

Location of the host genome fasta file that will be used to filter host reads from the metagenome. The script index_genome.pl will
build an index using this file.

* **-d, /PATH/TO/INDEX_DATABASE/** (REQUIRED argument) NO DEFAULT VALUE

Location of the directory that will contain the index of the host genome that is produced by this script

* **-n, HOST GENOME NAME** (REQUIRED argument) NO DEFAULT VALUE

User specified base name of the host genome index files. Will take on the form of <BASENAME>.bt2, for example, with the extension
varying depending on the method used.

* **-m, INDEX METHOD** (OPTIONAL argument) DEFAULT = bowtie2

The type of index that should be built. You must build an index that is appropriate for the method you will use to filter host DNA.
For example, if you will use bowtie2 to filter host DNA, then you must run index_genome.pl with -m bowtie2 beforehand.

### shotcleaner.pl ###

* **-1, /PATH/TO/FORWARD_READS.fq.gz** (REQUIRED argument) NO DEFAULT VALUE

Location of the forward reads to be processed by shotcleaner.pl, in fastq format. Gzipped compressed files can be specified.
Files must end with either .fastq, .fq, .fastq.gz, .fq.gz. Use this option for either paired end or single end sequences

* **-2, /PATH/TO/REVERSE_READS.fq.gz** (REQUIRED argument) NO DEFAULT VALUE

Location of the reverse, paired reads to be processed by shotcleaner.pl. File format and naming conventions are consistent with
-1. Not used for single end sequences.

* **-d, /PATH/TO/INDEX_DATABASE/** (REQUIRED argument) NO DEFAULT VALUE

Location of the directory that contains the index of the host genome that is produced by this script. You must build an index
prior to running shotcleaner.pl by using index_genome.pl.

* **-n, HOST GENOME NAME** (REQUIRED argument) NO DEFAULT VALUE

User specified base name of the host genome index files. Ttake on the form of <BASENAME>.bt2, for example, with the extension
varying depending on the method used.

* **--trim, LIST OF TRIM METHODS** (OPTIONAL) DEFAULT = "trimmomatic"

A comma separated (no white space) list of methods to use for quality trimming. Must be one or more of the following:

trimmomatic prinseq

* **--filter, LIST OF FILTER METHODS** (OPTIONAL) DEFAULT = "bowtie2"

A comma separated (no white space) list of methods to use for host filtering. Must be one or more of the following:

bowtie2 bmtagger

* **--derep, LIST OF DEREPLICATION METHODS** (OPTIONAL) DEFAULT = "fastuniq" for paired end, "prinseq" for single end

A comma separated (no white space) list of methods to use for collapsing duplicates. Must be one or more of the following:

fastuniq prinseq

Note that fastuniq can only work with paired end data.

* **--compress (OPTIONAL)** DEFAULT ACTIVE

Should the results be compressed? This removes data produced during intermediary steps and temporary files. Only the cleaned
fastq, fasta, and fastqc results will be retained, and all files will be gzipped. Silence by using --nocompress.

Dependencies
------------

###Perl Modules
*Carp
*File::Basename
*File::Copy
*File::Path
*File::Spec
*File::Spec::Functions
*Parallel::ForkManager
*Getopt::Long

###External Software
*bowtie2
*bmtagger
*FastQC
*FastUnique
*prinseq
*seqret (EMBOSS)
*trimmomatic