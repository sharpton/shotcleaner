#!/usr/bin/perl -w

use strict;
use File::Basename;
use File::Path qw( make_path );
use File::Spec;
use Getopt::Long;
use POSIX;

my $in_file;
my $n_splits;
my $out_dir;
my $limit_seqs = 0;
my $seq_limit = 0; #set to another value to downsample sequences
GetOptions(
    "i=s"  => \$in_file,
    "n=i"  => \$n_splits,
    "o=s"  => \$out_dir,
    "s:i"  => \$seq_limit,
    "sub!" => \$limit_seqs,
    );

if( !(defined( $in_file ) ) ){
    die( "You must specify a fastq file to process with -i\n" );
}
if( !(defined( $n_splits ) ) ){
    die( "You must specify the number of splits to create with -n\n" );
}
if( !(defined( $out_dir ) ) ){
    die( "You must specify an output directory with -o\n" );
}

make_path( $out_dir );

my $is_compressed   = _compression_check( $in_file );
my $nseqs_per_split = _calculate_split_size( $in_file, $n_splits, $is_compressed );
_split_sequence_file( $in_file, $nseqs_per_split, $out_dir, $is_compressed, $limit_seqs, $seq_limit );

###############
# SUBROUTINES #
###############

sub _calculate_split_size{
    my $in_file  = shift;
    my $n_splits = shift;
    my $is_compressed = shift;
    my $n_lines;   
    if( $is_compressed ){
        $n_lines = `zcat $in_file | wc -l`;
    } else {
	$n_lines = `wc -l $in_file`;
    }
    my $nseqs = $n_lines / 4; #for standard fastq
    $nseqs_per_split = ceil( $nseqs / $n_splits  ); #round up to nearest integer to be sure we get all reads
    return $nseqs_per_split;
}

#this compression check is a bit crude - 
#assumes file is gzipped and ends in .gz extension
sub _compression_check{
    my $in_file = shift;
    my $compressed = 0;
    if( $in_file =~ m/\.gz$/ ){
	$compressed = 1;
    } else {
	$compressed = 0;
    }
    return $compressed;
}

sub _split_sequence_file{
    my $in_file          = shift;
    my $nseqs_per_split  = shift;
    my $split_dir        = shift;
    my $is_compressed    = shift;
    my $limit_seqs       = shift;
    my $seq_limit        = shift;
    my @suffix           = ( ".fastq", ".fq", ".fastq.gz", ".fq.gz" );
    my $basename         = basename( $in_file, @suffix );
    my @output_names     = (); #a list of filenames
    if( $is_compressed ){
	open( SEQS, "zcat $in_file|" ) || die "Can't open $in_file for read\n"; 
    } else {
	open( SEQS, $in_file ) || die "Can't open $in_file for read\n";
    }
    my $counter  = 1;
    my $outname  = $basename . "_" . $counter . ".fq";
    my $splitout = File::Spec->catfile( $split_dir,  $outname );
    open( OUT, ">$splitout" ) || die "Can't open $splitout for write in Shotmap::DB::split_sequence_file_no_bp\n";
    push( @output_names, $outname );
    my $seq_ct   = 0;
    my $seq_count_across_splits = 0;
    my $line_counter = 0;
    while( <SEQS> ){
	$line_counter++;
	print OUT $_;
	if( $line_counter == 4 ){
	    $seq_ct++;
	    $seq_count_across_splits++;
	    $line_counter = 0;
	}
	#have we reached the prerarefy sequence count, if that is set?
	if( $limit_seqs &&
	    $seq_count_across_splits == $seq_limit ){
	    close OUT;
	    last;	    
	}	
	if( eof ){
	    close OUT;
	}
	#do we need to process a split?
	if( $seq_ct == $nseqs_per_split ){	
	    close OUT;
	    $counter++;
	    $outname  = $basename . "_" . $counter . ".fq";
	    $splitout = File::Spec->catfile( $split_dir,  $outname );	    
	    unless( eof ){
		open( OUT, ">$splitout" ) || die "Can't open $splitout for write\n";
		push( @output_names, $outname );
		$seq_ct = 0;		
	    }
	}
    }
    close SEQS;
    return \@output_names;
}
