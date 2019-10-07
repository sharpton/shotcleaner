#!/usr/bin/env perl

use warnings;
use strict;
use File::Spec;
use Cwd 'abs_path';
use File::Basename;
use File::Path;

my $nprocs = $ARGV[0];

if( !defined( $nprocs ) ){
    $nprocs = 1;
}

my $host_name = "e_coli";

my $master = dirname( abs_path($0) );
my $data   = File::Spec->catdir( $master, "test" );
my $idx_db = File::Spec->catdir( $data, "index_db" );
my $host   = File::Spec->catdir( $data, $host_name );
my $meta   = File::Spec->catdir( $data, "metagenome" );
my $output = File::Spec->catdir( $data, "test_out" );

my $host_genome = File::Spec->catfile( $host, "CP001163.fna" );

if( -e $output ){
    rmtree( $output );
}

if( ! -d $output ) {
    `mkdir -p $output`;
}

my $ap = $master . "/pkg/Trimmomatic-0.39/adapters/NexteraPE-PE.fa";

my $nocompress = 1;

####
# BOWTIE2
####
#test building a database index
my $cmd1 = "perl ${master}/index_genome.pl -i $host_genome -d $idx_db -n $host_name -m bowtie2";
print $cmd1 . "\n";
if( system($cmd1) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}
#test on a single end data
##fastq
my $out2 = File::Spec->catdir( $output, "single_fastq" ); 
my $cmd2 = "perl ${master}/shotcleaner.pl -1 ${meta}/single/single.fq.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out2 --adapt-path $ap";
if( $nocompress ){
    $cmd2 .= " --nocompress ";
}
print $cmd2 . "\n";
if( system($cmd2) ){ #system returns 0 on success
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
} 
##fasta
my $out3 = File::Spec->catdir( $output, "single_fasta" ); 
my $cmd3 = "perl ${master}/shotcleaner.pl -1 ${meta}/single/single.fa.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out3";
if( $nocompress ){
    $cmd3 .= " --nocompress ";
}
print $cmd3 . "\n";
if( system($cmd3) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}
#test on paired end data
##fastq
my $out4 = File::Spec->catdir( $output, "paired_fastq" ); 
my $cmd4 = "perl ${master}/shotcleaner.pl -1 ${meta}/paired/fastq/mate_1.fq.gz -2 ${meta}/paired/fastq/mate_2.fq.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out4 --adapt-path $ap";
if( $nocompress ){
    $cmd4 .= " --nocompress ";
}
print $cmd4 . "\n";
if( system($cmd4) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}
##fasta
my $out5 = File::Spec->catdir( $output, "paired_fasta" ); 
my $cmd5 = "perl ${master}/shotcleaner.pl -1 ${meta}/paired/fasta/mate_1.fa.gz -2 ${meta}/paired/fasta/mate_2.fa.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out5";
if( $nocompress ){
    $cmd5 .= " --nocompress ";
}
print $cmd5 . "\n";
if( system($cmd5) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}

##################
# ADD TESTS HERE
#test on paired end data
##fastq
my $out6 = File::Spec->catdir( $output, "paired_fastq-trimmomatic" );
my $cmd6 = "perl ${master}/shotcleaner.pl -1 ${meta}/paired/fastq/mate_1.fq.gz -2 ${meta}/paired/fastq/mate_2.fq.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out6 --adapt-path $ap --trim trimmomatic";
if( $nocompress ){
    $cmd6 .= " --nocompress ";
}
print $cmd6 . "\n";
if( system($cmd6) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}

##prinseq
my $out7 = File::Spec->catdir( $output, "paired_fastq-prinseq" );
my $cmd7 = "perl ${master}/shotcleaner.pl -1 ${meta}/paired/fastq/mate_1.fq.gz -2 ${meta}/paired/fastq/mate_2.fq.gz -d $idx_db -n $host_name -m bowtie2 --nprocs $nprocs -o $out7 --adapt-path $ap --trim prinseq";
if( $nocompress ){
    $cmd7 .= " --nocompress ";
}
print $cmd7 . "\n";
if( system($cmd7) ){
    die( "GOT AN ERROR RUNNING  SHOTCLEANER!\n" );
}

print "TESTS COMPLETE!\n";
