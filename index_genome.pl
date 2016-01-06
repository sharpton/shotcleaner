#!/usr/bin/perl -w

use strict;

use File::Basename;
use File::Copy;
use File::Path qw( make_path remove_tree );
use File::Spec;
use File::Spec::Functions;
use Parallel::ForkManager;
use Carp;
use Data::Dumper;
use Getopt::Long;

my $h_genome;           #host genome fasta file
my $index_db;           #directory containing indexed genomes
my $method = "bowtie2"; #search method, used to determine how to index
my $index_name;         #basename of the index

GetOptions(
    "i=s" => \$h_genome,
    "d=s" => \$index_db,
    "m:s" => \$method,
    "n=s" => \$index_name,
    );

if( $method eq "bowtie2" ){
    my $cmd = "bowtie2-build ${h_genome} ${index_db}/${index_name}";
    print "$cmd\n";
    system( $cmd );	    
}

if( $method eq "deconseq" ){
    my $cmd = "bwa index -p ${index_db}/${index_name} ${h_genome}";
    print "$cmd\n";
    system( $cmd );	    
}

if( $method eq "bmtagger" ){
    #bmfilter index
    my $cmd = "bmtool -d ${h_genome} -o ${index_db}/${index_name}.bitmask -w 18";
    print "$cmd\n";
    system( $cmd );
    #srprism index
    $cmd = "srprism mkindex -i ${h_genome} -o ${index_db}/${index_name}.srprism -M 7168";
    print "$cmd\n";
    system( $cmd );
    #blast index
    $cmd = "makeblastdb -in ${h_genome} -dbtype nucl -out ${index_db}/${index_name}.blast";
    print "$cmd\n";
    system( $cmd );	
}
