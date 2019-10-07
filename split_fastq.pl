#!/usr/bin/env perl

use warnings;
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
my $in_format = "fastq"; #fasta | fastq
GetOptions(
    "i=s"  => \$in_file,
    "n=i"  => \$n_splits,
    "o=s"  => \$out_dir,
    "s:i"  => \$seq_limit,
    "sub!" => \$limit_seqs,
    "f:s"  => \$in_format
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
if( $in_format ne "fasta" &&
    $in_format ne "fastq" ){
    die( "You must use fasta or fastq with -f\n" );
}

make_path( $out_dir );

my $is_compressed   = _compression_check( $in_file );
my $nseqs_per_split = _calculate_split_size( $in_file, $n_splits, $is_compressed, $in_format );
_split_sequence_file( $in_file, $nseqs_per_split, $out_dir, $is_compressed, $limit_seqs, $seq_limit, $in_format );

###############
# SUBROUTINES #
###############

sub _calculate_split_size{
    my $in_file  = shift;
    my $n_splits = shift;
    my $is_compressed = shift;
    my $in_format     = shift;
    my $nseqs;
    if( $in_format eq "fastq" ){
	my $n_lines;   
	if( $is_compressed ){
	    $n_lines = `zcat $in_file | wc -l`;
	} else {
	    $n_lines = `wc -l $in_file`;
	}
	$nseqs = $n_lines / 4; #for standard fastq
    } else { #is fasta
	if( $is_compressed ){
	    $nseqs = `zgrep -c ">" $in_file`;
	} else {
	    $nseqs = `grep -c ">" $in_file`;
	}
    }
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
    my $in_format        = shift;
    my @suffix           = ( ".fastq", ".fq", ".fastq.gz", ".fq.gz",
			     ".fasta", ".fa", ".fasta.gz", ".fa.gz" );
    my $basename         = basename( $in_file, @suffix );
    my @output_names     = (); #a list of filenames
    if( $is_compressed ){
	open( SEQS, "zcat $in_file|" ) || die "Can't open $in_file for read\n"; 
    } else {
	open( SEQS, $in_file ) || die "Can't open $in_file for read\n";
    }
    my $counter  = 1;
    my $suf;
    if( $in_format eq "fastq" ){
	$suf = ".fq";
    } else {
	$suf = ".fa";
    }   
    my $outname  = $basename . "_" . $counter . $suf;
    my $splitout = File::Spec->catfile( $split_dir,  $outname );
    open( OUT, ">$splitout" ) || die "Can't open $splitout for write in Shotmap::DB::split_sequence_file_no_bp\n";
    push( @output_names, $outname );
    my $seq_ct                  = 0;
    my $seq_count_across_splits = 0;
    my $line_counter            = 0;
    while( <SEQS> ){
	if( $in_format eq "fasta" ){
	    #do we need to process a split?
	    if( $seq_ct == $nseqs_per_split &&
		$_ =~ m/\>/ ){	
		close OUT;
		$counter++;
		my $suf;
		if( $in_format eq "fasta" ){
		    $suf = ".fa";
		} else {
		    $suf = ".fq"
		}
		$outname  = $basename . "_" . $counter . $suf;
		$splitout = File::Spec->catfile( $split_dir,  $outname );	    
		unless( eof ){
		    open( OUT, ">$splitout" ) || die "Can't open $splitout for write\n";
		    push( @output_names, $outname );
		    $seq_ct = 0;		
		}	    
	    }
	    if( $limit_seqs &&
		$seq_count_across_splits == $seq_limit ){
		close OUT;
		last;	    
	    }	
	    if( $_ =~ m/\>(.*?)\s/ ){
		$seq_ct++;
		$seq_count_across_splits++;		
		#bowtie doesn't play nice with whitespace in headers, at least in
		#version we tested with (2.2.3)
		print OUT ">${1}\n";
	    } else {
		print OUT uc( $_ );
	    }
	    if( eof ){
		close OUT;
	    }
	}
  	elsif( $in_format eq "fastq" ){
	    print OUT uc( $_ );
	    $line_counter++;
	    #count sequences 
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
		my $suf;
		if( $in_format eq "fasta" ){
		    $suf = ".fa";
		} else {
		    $suf = ".fq"
		}
		$outname  = $basename . "_" . $counter . $suf;
		$splitout = File::Spec->catfile( $split_dir,  $outname );	    
		unless( eof ){
		    open( OUT, ">$splitout" ) || die "Can't open $splitout for write\n";
		    push( @output_names, $outname );
		    $seq_ct = 0;		
		}	    
	    }
	}
    }
    close SEQS;
    return \@output_names;
}
