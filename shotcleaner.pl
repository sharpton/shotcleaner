#!/usr/bin/perl -w

# Inspiration: http://www.hmpdacc.org/doc/ReadProcessing_SOP.pdf
# And: http://bioinformatics.oxfordjournals.org/content/suppl/2011/12/12/btr669.DC1/SupplementaryFile1.pdf

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

my $root = _get_root();
my $in_fastq;
my $pair_fastq; #the mate paired file. Optional.
my $out_dir;
my $index_dir;
my $index_basename;


#defaults
my $nprocs             = 1;
my $paired_end         = 0;
my $compress           = 1;
my $qc_method_list     = ( "fastqc" );
#can use comma separated (no whitespace) string to do multiple methods per variable
my $trim_method_list   = "trimmomatic"; #trimmomatic | prinseq
my $filter_method_list = "bowtie2"; #bowtie2 | bmtagger (deconseq obsolete)
#note: we might add fastx_collapser or some other
#tool for single end data
my $derep_method_list  =  "fastuniq"; #fastuniq | prinseq-derep
my $overwrite          = 1; #for now, assume the user is aware
my @in_suffixes        = ( ".fq", ".fastq", ".fq.gz", ".fastq.gz" ); #useful for basename

GetOptions(
    "1=s"         => \$in_fastq,
    "2:s"         => \$pair_fastq,
    "o=s"         => \$out_dir,
    "d=s"         => \$index_dir,
    "n=s"         => \$index_basename,
    "nprocs=i"    => \$nprocs,
    "compress!"   => \$compress,
    "qc:s"        => \$qc_method_list,
    "trim:s"      => \$trim_method_list,
    "filter:s"    => \$filter_method_list,
    "derep:s"     => \$derep_method_list,
    );

### INITIALIZATION
_check_vars( $in_fastq,  $pair_fastq,    $out_dir, 
	     $index_dir, $index_basename );
my @qc_methods     = @{ _parse_method_list(  $qc_method_list,    "qc"     ) };
my @trim_methods   = @{ _parse_method_list( $trim_method_list,   "trim"   ) };
my @filter_methods = @{ _parse_method_list( $filter_method_list, "filter" ) };
my @derep_methods  = @{ _parse_method_list( $derep_method_list,  "derep"  ) };

my $mate_basename;
if( defined( $pair_fastq ) ){
    $paired_end = 1;
    $mate_basename = basename( $pair_fastq, @in_suffixes );
}

#fastuniq only works on paired end data
if( $derep_method_list =~ m/fastuniq/ &&
    !$paired_end ){
    die( "You cannot invoke fastuniq on non-paired end data. Please instead only invoke " .
	 "--derep prinseq-derep\n" );
}

my $logdir    = File::Spec->catdir( $out_dir, "logs" );
if( ! defined( $logdir ) ){
    $logdir    = $out_dir . "/logs/";
}
my $tmp_dir   = File::Spec->catdir( $out_dir, "_tmp");
make_path( $out_dir );
make_path( $logdir );
make_path( $tmp_dir );

### GET PIPELINE SETTINGS

my $settings = _set_settings( \@qc_methods, \@trim_methods, 
			      \@filter_methods, \@derep_methods );
my ( $run_raw_qc, $split_reads, $run_trim, $run_filter,
  $cat_reads, $derep, $check_qc, $fasta_cleaned ) = @{ $settings->{"parameters"} };
 

### START PIPELINE

print scalar(localtime()) . "\n";

if( $run_raw_qc ){
    foreach my $method( @qc_methods ){
	if( $method eq "fastqc" ){
	    my $fastqc_result_dir = File::Spec->catdir( $out_dir, "/fastqc_raw/" );
	    my $fastqc_log_dir    = File::Spec->catdir( $logdir, "/fastqc_raw/" );
	    my $fastqc_bin        = File::Spec->catfile( $root, "bin", "fastqc" );
	    my @files; 
	    make_path( $fastqc_result_dir );
	    make_path( $fastqc_log_dir );
	    if( $paired_end ){
		@files = ( $in_fastq, $pair_fastq );
	    } else {
		@files = ( $in_fastq );
	    }
	    print "FASTQC: " . scalar(localtime()) . "\n";
	    foreach my $file( @files ){
		_run_fastqc( {
		    in_file     => $file, 
		    result_dir  => $fastqc_result_dir, 
		    log_dir     => $fastqc_log_dir, 
		    fastqc      => $fastqc_bin,
			     }
		    );
	    }
	}
    }
}

my @reads = ();

if( $split_reads ){
    my $split_reads_input = $in_fastq;
    my $split_reads_dir   = File::Spec->catdir( $out_dir, $settings->{"split_reads"}->{"output"} );
    my $split_log_dir     = File::Spec->catdir( $logdir, "split_reads" );
    my $split_reads_bin   = File::Spec->catfile( $root, "split_fastq.pl" );
    make_path( $split_reads_dir );
    make_path( $split_log_dir );
    print "SPLIT READS: " . scalar(localtime()) . "\n";
    _run_split_reads( {
	in_seqs     => $split_reads_input,
	result_dir  => $split_reads_dir,
	n_splits    => $nprocs,	
	log_dir     => $split_log_dir,
	split_reads => $split_reads_bin,
		      }
	);
    if( $paired_end ){
	my $split_pair_input = $pair_fastq;
	my $split_pair_dir = File::Spec->catdir( $out_dir, $settings->{"split_reads"}->{"output"} );
	_run_split_reads( {
	    in_seqs     => $split_pair_input,
	    result_dir  => $split_pair_dir,
	    n_splits    => $nprocs,	
	    log_dir     => $split_log_dir,
	    split_reads => $split_reads_bin,
			  }
	    );
    }
    if( defined( $split_reads_dir ) ){
	@reads = @{ _get_fastq_file_names( $split_reads_dir, 0 ) };
    } else {
	die "Couldn't get split read file names!\n";
    }
}

if( $run_trim ){
    foreach my $method( @trim_methods ){
	if( $method eq "prinseq" ){
	    my $prinseq_trim_in_dir      = File::Spec->catdir( $out_dir, $settings->{"prinseq"}->{"input"} );
	    my $prinseq_trim_results_dir = File::Spec->catdir( $out_dir, $settings->{"prinseq"}->{"output"} );
	    my $prinseq_trim_log_dir     = File::Spec->catdir( $logdir, $settings->{"prinseq"}->{"log"} );
	    my $prinseq_bin              = File::Spec->catfile( $root, "bin", "prinseq-lite.pl" );
	    make_path( $prinseq_trim_results_dir );
	    make_path( $prinseq_trim_log_dir );
	    print "PRINSEQ, TRIM&FILTER: " . scalar(localtime()) . "\n";
	    _run_prinseq( {
		prinseq       => $prinseq_bin,
		in_dir        => $prinseq_trim_in_dir,
		file_names    => \@reads,
		log_dir       => $prinseq_trim_log_dir,
		result_dir    => $prinseq_trim_results_dir,
		nprocs        => $nprocs,
		overwrite     => $overwrite,
		derep         => 0,
		paired_end    => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
		  });	
	}
	if( $method eq "trimmomatic" ){
	    my $trim_in_dir      = File::Spec->catdir( $out_dir, $settings->{"trimmomatic"}->{"input"} );
	    my $trim_results_dir = File::Spec->catdir( $out_dir, $settings->{"trimmomatic"}->{"output"} );
	    my $trim_log_dir     = File::Spec->catdir( $logdir, $settings->{"trimmomatic"}->{"log"} );
	    my $trim_bin              = File::Spec->catfile( $root, "bin", "trimmomatic-0.35.jar" );
	    make_path( $trim_results_dir );
	    make_path( $trim_log_dir );
	    print "TRIMMOMATIC: " . scalar(localtime()) . "\n";
	    _run_trimmomatic( {
		trimmomatic   => $trim_bin,
		in_dir        => $trim_in_dir,
		file_names    => \@reads,
		log_dir       => $trim_log_dir,
		result_dir    => $trim_results_dir,
		nprocs        => $nprocs,
		overwrite     => $overwrite,
		paired_end    => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
		  });		
	}
    }
}

if( $run_filter ){
    foreach my $method( @filter_methods ){
	if( $method eq "bowtie2" ){
	    my $bowtie_in_dir      = File::Spec->catdir( $out_dir, $settings->{"bowtie2"}->{"input"});
	    my $bowtie_results_dir = File::Spec->catdir( $out_dir, $settings->{"bowtie2"}->{"output"});
	    my $bowtie_log_dir     = File::Spec->catdir( $logdir, $settings->{"bowtie2"}->{"log"} );
	    my $bowtie_bin         = File::Spec->catfile( $root, "bin", "bowtie2" );
	    my @db_names             = ( $index_basename );
	    make_path( $bowtie_results_dir );
	    make_path( $bowtie_log_dir );
	    print "BOWTIE2: " . scalar(localtime()) . "\n";
	    _run_bowtie2( {
		bowtie2     => $bowtie_bin,
		in_dir      => $bowtie_in_dir, 
		file_names  => \@reads, 
		result_dir  => $bowtie_results_dir, 
		log_dir     => $bowtie_log_dir, 
		nprocs      => $nprocs,
		index_dir   => $index_dir,
		tmp_dir     => $tmp_dir,
		db_names    => \@db_names,
		overwrite   => $overwrite,
		paired_end  => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
			   });       	    
	}
	if( $method eq "bmtagger" ){
	    my $bmtagger_in_dir      = File::Spec->catdir( $out_dir, $settings->{"bmtagger"}->{"input"});
	    my $bmtagger_results_dir = File::Spec->catdir( $out_dir, $settings->{"bmtagger"}->{"output"});
	    my $bmtagger_log_dir     = File::Spec->catdir( $logdir, $settings->{"bmtagger"}->{"log"} );
	    my $bmtagger_bin         = File::Spec->catfile( $root, "bin", "bmtagger.sh" );
	    my $extract              = 1; #YOU PROBABLY WANT THIS, prints only non-host sequences in output file
	    my @db_names             = ( $index_basename );
	    make_path( $bmtagger_results_dir );
	    make_path( $bmtagger_log_dir );
	    print "BMTAGGER: " . scalar(localtime()) . "\n";
	    _run_bmtagger( {
		bmtagger    => $bmtagger_bin,
		in_dir      => $bmtagger_in_dir, 
		file_names  => \@reads, 
		result_dir  => $bmtagger_results_dir, 
		log_dir     => $bmtagger_log_dir, 
		nprocs      => $nprocs,
		bitmask_dir => $index_dir,
		srprism_dir => $index_dir,
		tmp_dir     => $tmp_dir,
		db_names    => \@db_names,
		extract     => $extract,
		overwrite   => $overwrite,
		paired_end  => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,

			   });       
	}
	#note that deconseq doesn't currently support mate_paired data
	#also, it is challenging to automate the deconseq config file
	#finally, an update to bwa now breaks deconseq in our hands (change in the output produced
	#by bwa index). 
	#given the above, we do not currently support deconseq.
	if( $method eq "deconseq" ){
	    my $deconseq_in_dir      = File::Spec->catdir( $out_dir, $settings->{"deconseq"}->{"input"} );
	    my $deconseq_results_dir = File::Spec->catdir( $out_dir, $settings->{"deconseq"}->{"output"} );
	    my $deconseq_log_dir     = File::Spec->catdir( $logdir, $settings->{"deconseq"}->{"log"} );
	    my $deconseq_bin         = File::Spec->catfile( $root, "bin", "deconseq.pl" );
	    my @db_names             = ( $index_basename );
            #note that decon_db_names MUST point to a key in the deconseq perl module
            #where database locations are specified. A change the array below
            #may require revising said perl module. You can probably find the module
            #in some place like /src/deconseq-standalone-0.4.3/DeconSeqConfig.pm
	    make_path( $deconseq_results_dir );
	    make_path( $deconseq_log_dir );
	    
	    print "DECONSEQ: " . scalar(localtime()) . "\n";
	    _run_deconseq( {
		deconseq    => $deconseq_bin,
		in_dir      => $deconseq_in_dir, 
		file_names  => \@reads, 
		result_dir  => $deconseq_results_dir,
		log_dir     => $deconseq_log_dir, 
		nprocs      => $nprocs,
		tmp_dir     => $tmp_dir,
		db_names    => \@db_names,
		overwrite   => $overwrite,
		paired_end  => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
			   });
	    
	}
    }
}

# Merge results
if( $cat_reads ){
    my $cat_reads_in_dir      = File::Spec->catdir( $out_dir, $settings->{"cat_reads"}->{"input"} );
    my $cat_reads_results_dir = File::Spec->catdir( $out_dir, $settings->{"cat_reads"}->{"output"});
    make_path( $cat_reads_results_dir );
    my $cat_reads_log = File::Spec->catdir( $logdir, $settings->{"cat_reads"}->{"log"} );
    make_path( $cat_reads_log );
        
    print "CATTING READS: " . scalar(localtime()) . "\n";
    _cat_reads( {
	in_dir     => $cat_reads_in_dir,
	file_names => \@reads,
	log_dir    => $cat_reads_log,
	result_dir => $cat_reads_results_dir,
	nprocs     => $nprocs,	
	paired_end => $paired_end,
	mate_basename => $mate_basename,
	in_suffixes   => \@in_suffixes,
		});
}

if( @derep_methods ){
    foreach my $method( @derep_methods ){
	if( $method eq "fastuniq" ){
	    my $fu_in_dir      = File::Spec->catdir( $out_dir, $settings->{"fastuniq"}->{"input"} );
	    my $fu_results_dir = File::Spec->catdir( $out_dir, $settings->{"fastuniq"}->{"output"} );
	    my $fu_log_dir     = File::Spec->catdir( $logdir, $settings->{"fastuniq"}->{"log"});
	    my $fu_bin         = File::Spec->catfile( $root, "bin", "fastuniq" );
	    make_path( $fu_results_dir );
	    make_path( $fu_log_dir );
	    
	    print "FASTUNIQ: " . scalar(localtime()) . "\n";
	    _run_fastuniq( {
		fastuniq   => $fu_bin,
		in_dir     => $fu_in_dir,
		file_names => \@reads,
		log_dir    => $fu_log_dir,
		result_dir => $fu_results_dir,
		nprocs     => $nprocs,
		overwrite  => $overwrite,
		tmp_dir    => $tmp_dir,
		paired_end => $paired_end,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
			  });	
	}
	if( $method eq "prinseq-derep" ){
	    my $prinseq_derep_in_dir     = File::Spec->catdir( $out_dir, $settings->{"prinseq-derep"}->{"input"} );
	    my $prinseq_derep_results_dir = File::Spec->catdir( $out_dir, $settings->{"prinseq-derep"}->{"output"} );
	    my $prinseq_derep_log_dir        = File::Spec->catdir( $logdir, $settings->{"prinseq-derep"}->{"log"});
	    my $prinseq_bin               = File::Spec->catfile( $root, "bin", "prinseq-lite.pl" );
	    make_path( $prinseq_derep_results_dir );
	    make_path( $prinseq_derep_log_dir );
        
	    print "PRINSEQ, DEREP: " . scalar(localtime()) . "\n";
	    _run_prinseq( {
		prinseq    => $prinseq_bin,
		in_dir     => $prinseq_derep_in_dir,
		file_names => \@reads,
		log_dir    => $prinseq_derep_log_dir,
		result_dir => $prinseq_derep_results_dir,
		nprocs     => $nprocs,
		overwrite  => $overwrite,
		derep      => 1,
		mate_basename => $mate_basename,
		in_suffixes   => \@in_suffixes,
			  });
	}
    }
}

if( $check_qc ){
    my $fastqc_clean_in_dir      = File::Spec->catdir( $out_dir, $settings->{"check_qc"}->{"input"}  );
    my $fastqc_clean_results_dir = File::Spec->catdir( $out_dir, $settings->{"check_qc"}->{"output"});
    my $fastqc_clean_log_dir     = File::Spec->catdir( $logdir, $settings->{"check_qc"}->{"log"} );
    my $fastqc_bin               = File::Spec->catfile( $root, "bin", "fastqc" );
    my @files; 
    make_path( $fastqc_clean_results_dir ); 
    make_path( $fastqc_clean_log_dir );

    print "FASTQC, CLEAN: " . scalar(localtime()) . "\n";
    my $in_base = basename( $in_fastq, @in_suffixes ); 
    my $in_fastq = File::Spec->catfile( $fastqc_clean_in_dir, $in_base . ".fq" );
    if( $paired_end ){
	my $in_pair = $in_fastq;
	$in_pair    =~ s/${in_base}/${mate_basename}/;	   
	@files = ( $in_fastq, $in_pair );
    } else {
	@files = ( $in_fastq );
    }
    print "FASTQC: " . scalar(localtime()) . "\n";
    foreach my $file( @files ){
	_run_fastqc( {
	    in_file     => $file, 
	    result_dir  => $fastqc_clean_results_dir, 
	    log_dir     => $fastqc_clean_log_dir, 
	    fastqc      => $fastqc_bin,
		     }
	    );
    }
}

if( $fasta_cleaned ){    
    my $fasta_in_dir      = File::Spec->catdir( $out_dir, $settings->{"fasta_cleaned"}->{"input"}  );
    my $fasta_results_dir = File::Spec->catdir( $out_dir, $settings->{"fasta_cleaned"}->{"output"} );
    my $fasta_log_dir     = File::Spec->catdir( $logdir, $settings->{"fasta_cleaned"}->{"log"} );
    my $seqret_bin        = File::Spec->catfile( $root, "bin", "seqret" );
    make_path( $fasta_results_dir );
    make_path( $fasta_log_dir );
    
    print "SEQRET: " . scalar(localtime()) . "\n";
    _run_seqret( {
	in_dir      => $fasta_in_dir,
	file_names  => \@reads, 
	result_dir  => $fasta_results_dir, 
	log_dir     => $fasta_log_dir, 
	nprocs      => 1,
	seqret      => $seqret_bin,
	mate_basename => $mate_basename,
	in_suffixes   => \@in_suffixes,

		 });
}

die;

if( $compress ){
    print "COMPRESSING DATA: " . scalar(localtime()) . "\n";
    _compress_results( { 
	in_dir => File::Spec->catdir( $out_dir, $settings->{"fasta_cleaned"}->{"output"}  )
		       });
    foreach my $method( @derep_methods ){
	_compress_results( { 
	    in_dir => File::Spec->catdir( $out_dir, $settings->{$method}->{"output"}  )
			   });
    }
    _compress_results( { 
	in_dir => File::Spec->catdir( $out_dir, $settings->{"cat_reads"}->{"output"}  ),
	delete => 1,
		       });
    foreach my $method( @trim_methods ){
	_compress_results( { 
	    in_dir => File::Spec->catdir( $out_dir, $settings->{$method}->{"output"}  ),
	    delete => 1,
		       });
    }
    foreach my $method( @filter_methods ){
	_compress_results( { 
	    in_dir => File::Spec->catdir( $out_dir, $settings->{$method}->{"output"}  ),
	    delete => 1,
			   })
    };
}

print scalar(localtime()) . "\n";


########
#
# SUBROUTINES
#
#######

sub _run_seqret{
    my $args = shift;
    
    my $in_dir     = $args->{"in_dir"};
    my @reads      = @{ $args->{"file_names"} };
    my $result_dir = $args->{"result_dir"};
    my $log_dir    = $args->{"log_dir"};
    my $nprocs     = $args->{"nprocs"};
    my $seqret     = $args->{"seqret"};
    my $mate_base  = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };

    for( my $i=0; $i<$nprocs; $i++ ){
	my $cmd;
	#do some housekeeping
	my $in_file   = $reads[0]; #includes split value!
	$in_file      =~ s/\_\d+\.fq$/\.fq/; #no longer does
	my $f_in  = File::Spec->catfile( $in_dir, $in_file );
	my $f_base = basename( $in_file, @in_suffixes ); #this contains the split value!
	my $f_log = File::Spec->catfile( $log_dir, $f_base . ".log" );
	my $outseq = File::Spec->catfile( $result_dir, $f_base . ".fa" );	
	$cmd = "$seqret -sequence $f_in -outseq $outseq -sformat1 fastq -osformat2 fasta &> $f_log";
	print "$cmd\n";
	system( $cmd );
	if( $paired_end ){
	    my $r_in;
	    $r_in  = $f_in;
	    $r_in  =~ s/$f_base/${mate_base}/;
	    my $r_base  = basename( $r_in, @in_suffixes );
	    my $r_log   = File::Spec->catfile( $log_dir, $r_base . ".log" );
	    my $routseq = File::Spec->catfile( $result_dir, $r_base . ".fa" );	
	    my $rcmd = "$seqret -sequence $r_in -outseq $routseq -sformat1 fastq -osformat2 fasta &> $r_log";
	    print "$rcmd\n";
	    system( $rcmd );
	}
    }
}


sub _compress_results{
    my $args = shift;
    
    my $in_dir = $args->{"in_dir"};
    my $delete = $args->{"delete"};

    opendir( DIR, $in_dir ) || die "Can't opendir $in_dir for read: $!\n";
    my @files = readdir( DIR );
    closedir DIR;
    foreach my $file( @files ){
	next if ( $file =~ m/^\./ );
	my $path = $in_dir . "/" . $file;
	if( defined( $delete ) && $delete ){
	    print "Deleting $path\n";
	    unlink( $path )
	}
	else{
	    next if( -d $path );
	    print "Compressing $path\n";
	    system( "gzip $path" );
	}
    }
}

sub _get_fastq_file_names{
    my $masterdir = shift;
    my $is_clean_check = shift; #are we running this on dereped data rather than raw reads?
    opendir( MASTER, $masterdir) || die "Can't open $masterdir for read: $!\n";
    my @files = readdir( MASTER );
    closedir MASTER;
    my @reads = (); 
    foreach my $read( @files ){
	if( !$is_clean_check ){
	    next if( $read !~ m/\.fastq/ &&
		     $read !~ m/\.fq/    && 
		     $read !~ m/\_R1\_/ );
	} else {
	    next if( $read !~ m/\.fastq/ &&
		     $read !~ m/\.fq/    );
	} 
	my @suffixlist = ( ".gz" );
	my( $name, $path, $suffix ) = fileparse( $read, @suffixlist );
	push( @reads, $name );
    }
    return \@reads;
}

sub _run_fastqc{
    my( $args ) = @_;
 
    #args->{in_file} is the fastq file to process
    #args->{result_dir} is output directory
    #args->{log_dir} is directory containing log files for each task
    #args->{nprocs} is number of tasks to run in parallel
    #args->{fastqc} is path to fastqc binary
    
    my $in_file    = $args->{"in_file"};
    my $result_dir = $args->{"result_dir"};
    my $log_dir    = $args->{"log_dir"};
    my $fastqc     = $args->{"fastqc"};
    my $is_clean   = $args->{"is_clean_check"};

    my $in_base  = basename( $in_file );
    my $log_file = File::Spec->catfile( $log_dir, $in_base . ".log" );
    my $cmd = "${fastqc} -o=$result_dir $in_file > $log_file 2>&1";
    print "$cmd\n";
    system( $cmd );
}

sub _run_split_reads{
    my( $args ) = @_;
    my $in_file     = $args->{"in_seqs"};
    my $result_dir  = $args->{"result_dir"};
    my $n_splits    = $args->{"n_splits"};
    my $log_dir     = $args->{"log_dir"};
    my $split_reads = $args->{"split_reads"};
    my $file_basename = basename( $in_file );
    my $log_file = File::Spec->catfile( $log_dir, $file_basename . ".log" );
    my $cmd = "${split_reads} -i ${in_file} -n ${n_splits} -o ${result_dir} > ${log_file} 2>&1";
    print "$cmd\n";
    system( $cmd );
}


sub _run_bowtie2{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $index_dir   = $args->{"index_dir"};
    my @db_names    = @{ $args->{"db_names"} };
    my $overwrite   = $args->{"overwrite"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my $bowtie2     = $args->{"bowtie2"};
    my $paired_end  = $args->{"paired_end"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };


    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];          
	my $f_mate = $read;
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	my $f_base = basename( $read, @in_suffixes ); #this contains the split value!
	my $r_in;
	if( $paired_end ){
	    $r_in  = $f_in;
	    $r_in     =~ s/$f_base/${mate_base}_${i}/;
	}
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	#prepare the output
        #might need to consider if we want a single output file or not.  
	my $out_stem = $read;
	my $out_path = File::Spec->catfile( $result_dir, $out_stem );
	#loop over dbs and run bmtagger
	foreach my $db( @db_names ){
	    my $database = File::Spec->catfile( $index_dir,  $db );
	    my $input_string = "";
	    if( $paired_end ){
	        $input_string .= "-1 $f_in ";
		$input_string .= "-2 $r_in ";
	    } else {
		$input_string .= "-U $f_in ";
	    }
	    my $out_string = "";
	    #with this notation, we keep all reads that don't concordantly align. 
	    #not necessairly conservative, but maybe best for metagenomes? Need to
	    #test accuracy.
	    if( $paired_end ){
		$out_string .= "--un-conc ${out_path} ";
	    } else { 
		$out_string .= "--un ${out_path} ";
	    }	    
	    $cmd =  "$bowtie2 ";
	    $cmd .= "${out_string} -x ${database} ${input_string} -S /dev/null ";
	    $cmd .= " &> $f_log";
	    print "$cmd\n";
	    system( $cmd );	    
	    #trim the name of the output file to standard to ease next steps. Use a move
	    if( !$paired_end ){
		#nothing needs to be done - names are preserved
	    } 
	    if( $paired_end ){ 	    #a little awkward with how bowtie renames files
	    #read 1.
		my $out_base = basename( $out_path, ".fq" );
		print $out_base . "\n";
		my $to_move = File::Spec->catfile( $result_dir, $out_base . ".1.fq" );
		print $to_move . "\n";
		my $new     = $out_path;
		print $new . "\n";
		move( $to_move, $new );
		#read 2
		my $to_move2 = File::Spec->catfile( $result_dir, $out_base . ".2.fq" );
		my $new2     = $out_path;
		$new2        =~ s/$f_base/${mate_base}_${i}/;
		move( $to_move2, $new2 );
	    }
	}
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
}

sub _run_bmtagger{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $bitmask_dir = $args->{"bitmask_dir"};
    my $sprism_dir  = $args->{"srprism_dir"};
    my @db_names    = @{ $args->{"db_names"} };
    my $extract     = $args->{"extract"};
    my $overwrite   = $args->{"overwrite"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my $bmtagger    = $args->{"bmtagger"};
    my $paired_end  = $args->{"paired_end"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };


    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];          
	my $f_mate = $read;
	my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	my $f_base = basename( $read, @in_suffixes ); #this contains the split value!
	my $r_in;
	if( $paired_end ){
	    $r_in  = $f_in;
	    $r_in     =~ s/$f_base/${mate_base}_${i}/;
	}
	my $f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	#prepare the output
        #might need to consider if we want a single output file or not.  
	my $out_stem = $read;
	my $out_path = File::Spec->catfile( $result_dir, $out_stem );
	#loop over dbs and run bmtagger
	foreach my $db( @db_names ){
	    my $bitmask  = File::Spec->catfile( $bitmask_dir, $db . ".bitmask" );
	    my $srprism  = File::Spec->catfile( $sprism_dir,  $db . ".srprism"  );
	    my $database = File::Spec->catfile( $sprism_dir,  $db . ".fa" );
	    my $f_string = "-1 " . $f_in;
	    my $r_string = "-2 " . $r_in;
	    if( $extract ){
		#bmtagger.sh -X -b $BITMASK -x $SPRISM -T $TMP -q1 $FREAD $RREAD -o $OUTPUT  >> $LOGS/bmtagger/${JOB_ID}.all 2>&1
#		$cmd = "run_bmtagger.sh -X -b $bitmask -x $sprism -T $tmp_dir -q1 $f_string $r_string -o $out_path > $f_log 2>&1";
		$cmd =  "$bmtagger ";
		$cmd .= "-X ";
		#$cmd .= "-b $bitmask -x $srprism -T $tmp_dir -d $database ";
		$cmd .= "-b $bitmask -x $srprism -d $database ";
		if( $paired_end ){
		    $cmd .= "-q 1 $f_string $r_string -o $out_path"; 
		} else {
		    $cmd .= "-q 1 $f_string -o $out_path" ;
		}
	    }
	    else{
                #bmtagger.sh -b $BITMASK -x $SPRISM -T $TMP -q1 $FREAD $RREAD -o $OUTPUT  >> $LOGS/bmtagger/${JOB_ID}.all 2>&1
#		$cmd = "run_bmtagger.sh -b $bitmask -x $sprism -T $tmp_dir -q1 $f_string $r_string -o $out_path > $f_log 2>&1";
		$cmd =  "$bmtagger ";
		$cmd .= "-b $bitmask -x $srprism -T $tmp_dir -d $database ";
		if( $paired_end ){
		    $cmd .= "-q 1 $f_string $r_string -o $out_path"; 
		} else {
		    $cmd .= "-q 1 $f_string -o $out_path"; 
		}
	    }
	    $cmd .= " &> $f_log";
	    print "$cmd\n";
	    system( $cmd );	    
	    #trim the name of the output file to standard to ease next steps. Use a move
	    #read 1
	    my $to_move = $out_path . "_1.fastq";
	    my $new     = $out_path;
	    move( $to_move, $new );
	    #read 2
	    if( $paired_end ){
		my $to_move = $out_path . "_2.fastq";
		my $new     = $to_move;
		$new        =~ s/$f_base/${mate_base}_${i}/;
		$new        =~ s/_2\.fastq$//;
		move( $to_move, $new );
	    }
	}
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
}

sub _run_deconseq{
    my( $args )  = @_;

    my $deconseq    = $args->{ "deconseq" };
    my $in_dir      = $args->{ "in_dir" };
    my @reads       = @{ $args->{ "file_names" } }; 
    my $result_dir  = $args->{ "result_dir" };
    my $log_dir     = $args->{"log_dir"}; 
    my $nprocs      = $args->{"nprocs"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my @db_names    = @{ $args->{"db_names"} };
    my $overwrite   = $args->{"overwrite"};
    my $paired_end  = $args->{"paired_end"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };

    my $db_string   = join(",", @db_names );

    #forward reads
    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	my $cmd;
	#do some housekeeping
	my $read = $reads[$i-1];
	my $f_mate     = $read;
	my  $f_out_name = $read;
	my( $f_in, $f_log );
	$f_in  = File::Spec->catfile( $in_dir, $f_mate );
	$f_log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	#Start threads
	$cmd = "$deconseq -id $f_out_name -out_dir $result_dir -f $f_in -dbs $db_string";
	print "$cmd\n";
	system( $cmd );
	#trim the name of the output file to standard to ease next steps. Use a move
	#read 1
	my $result_stem = File::Spec->catfile( $result_dir, $f_out_name );
	my $to_move = $result_stem . "_clean.fastq";
	my $new     = $result_stem;
	move( $to_move, $new );
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
    #reverse reads
    if( $paired_end ){
	my $pm2 = Parallel::ForkManager->new($nprocs);
	for( my $i=1; $i<=$nprocs; $i++ ){
	    my $pid = $pm2->start and next;
	    my $cmd;
	    #do some housekeeping
	    my $read = $reads[$i-1];
	    my $f_mate = $read;
	    my $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    my $f_base = basename( $read, @in_suffixes ); #this contains the split value!
	    my $r_in;
	    if( $paired_end ){
		$r_in  = $f_in;
		$r_in     =~ s/$f_base/${mate_base}_${i}/;
	    }
	    my $r_out_name = "${mate_base}_${i}";
	    my $r_log = File::Spec->catfile( $log_dir, "${mate_base}_${i}" . ".log" );
	    #Start threads
	    $cmd = "$deconseq -id $r_out_name -out_dir $result_dir -f $r_in -dbs $db_string";
	    print "$cmd\n";
	    system( $cmd );
	    #read 2
	    if( $paired_end ){
		my $result_stem = File::Spec->catfile( $result_dir, $r_out_name );
		my $to_move = $result_stem . "_clean.fastq";
		my $new     = $to_move;
		$new        =~ s/$f_base/${mate_base}_${i}/;
		$new        =~ s/_clean\.fastq$//;
		move( $to_move, $new );
	    }
	    $pm2->finish; # Terminates the child process
	}
	$pm2->wait_all_children;    
    }
}


sub _run_prinseq{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $overwrite   = $args->{"overwrite"};
    my $prinseq     = $args->{"prinseq"};
    my $derep       = $args->{"derep"};
    my $paired_end  = $args->{"paired_end"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };

    if( ! $derep ){
	my $pm = Parallel::ForkManager->new($nprocs);
	for( my $i=1; $i<=$nprocs; $i++ ){
	    my $pid = $pm->start and next;
	    #do some housekeeping
	    my $read = $reads[$i-1];          
	    my ( $f_in, $r_in );
	    my $f_mate = $read;
	    $f_in  = File::Spec->catfile( $in_dir, $f_mate );
	    my $f_base = basename( $read, @in_suffixes ); #this contains the split value!
	    if( $paired_end ){
		$r_in  = $f_in;
		$r_in     =~ s/$f_base/${mate_base}_${i}/;
	    }
	    my $log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	    my $out_path = File::Spec->catfile( $result_dir, $read );
	    my $compressed = 0;	   
	    my $gz_file = File::Spec->catfile( $f_in . ".gz" );
	    if( -e $gz_file ){
		$compressed = 1;
	    }
	    my $cmd = "";
	    if( $compressed ){
		$cmd .= "zcat ${f_in}.gz | ";
	    }
	    #you might want to turn derep back on here for downstream efficiency, but we turned off for
	    #courtney's analysis
	    #$cmd =  "$prinseq -verbose -derep 14 -derep_min 2 -no_qual_header "; #do we want -exact_only?
	    $cmd .=  "$prinseq -verbose -no_qual_header "; #do we want -exact_only?
	    $cmd .= "-min_len 60 -max_len 200 -min_qual_mean 25 -ns_max_n 0 ";
	    $cmd .= "-lc_method entropy -lc_threshold 60 -trim_qual_left 20 -trim_qual_right 20 ";
	    if( $paired_end ){
		$cmd .= "-out_good $out_path -fastq $f_in -fastq2 $r_in -log $log ";
	    } else {
		if( $compressed ){
		    $cmd .= "-fastq stdin ";
		} else {
		    $cmd .= "-fastq $f_in ";
		}
		$cmd .= "-out_good $out_path -log $log ";
	    }
	    $cmd .= "-out_bad null ";
	    $cmd    .= "&> $log ";
	    print "$cmd\n";
	    system( $cmd );
	    #trim the name of the output file to standard to ease next steps. Use a move
	    #read 1
	    my $to_move = $out_path . "_1.fastq";
	    my $new     = $out_path;
	    move( $to_move, $new );
	    #read 2
	    if( $paired_end ){
		my $to_move = $out_path . "_2.fastq";
		my $new     = $to_move;
		$new        =~ s/$f_base/${mate_base}_${i}/;
		$new        =~ s/_2\.fastq$//;
		move( $to_move, $new );
	    }
	    $pm->finish; # Terminates the child process
	}
	$pm->wait_all_children;    
    } else {
	my $in_file   = $reads[0]; #includes split value!
	$in_file      =~ s/\_\d+\.fq$/\.fq/; #no longer does
	my $f_in  = File::Spec->catfile( $in_dir, $in_file );
	my $f_base = basename( $in_file, @in_suffixes ); #this contains the split value!
	my $r_in;
	if( $paired_end ){
	    $r_in  = $f_in;
	    $r_in  =~ s/$f_base/${mate_base}/;
	}
	my $log      = File::Spec->catfile( $log_dir, $in_file . ".log" );
	my $out_file = $in_file; #just the file name 
	my $out_path = File::Spec->catfile( $result_dir, $out_file );
	my $cmd =  "$prinseq -verbose -derep 14 -derep_min 2 -no_qual_header "; #do we want -exact_only?
	#$cmd    .= "-out_good $out_path -fastq $in_path -log $log_path ";
	if( $paired_end ){
	    $cmd .= "-out_good $out_path -fastq $f_in -fastq2 $r_in -log $log ";
	} else {
	    $cmd .= "-fastq $f_in ";	   
	    $cmd .= "-out_good $out_path -log $log ";
	}
	$cmd    .= "-out_bad null ";
	$cmd    .= "&> $log ";
	print "$cmd\n";
	system( $cmd );	    
	#trim the name of the output file to standard to ease next steps. Use a move
	#read 1
	my $to_move = $out_path . ".fastq";
	my $new     = $out_path;
	$new        =~ s/\.fastq$//;
	move( $to_move, $new );
	#read 2
	if( $paired_end ){
	    my $to_move = $out_path . ".fastq";
	    my $new     = $to_move;
	    $new        =~ s/$f_base/${mate_base}/;
	    $new        =~ s/\.fastq$//;
	    move( $to_move, $new );
	}
    }
}

#Currently ignores sequences in the prinseq singletons file! This is the conservative move, it would seem.
#Also, F and R reads are pushed into the same file.
sub _cat_reads{
    my( $args )    = shift;
    my $in_dir     = $args->{"in_dir"};
    my $paired_end = $args->{"paired_end"};
    my @reads      = @{ $args->{"file_names"} };
    my $log_dir    = $args->{"log_dir"};
    my $result_dir = $args->{"result_dir"};
    my $mate_base  = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };
    #prep the array of files for cat
    my @sorted     = sort( @reads );    
    @sorted        = map { $in_dir . "/" . $_ } @sorted;
    my @f_sorted   = ();
    my @r_sorted   = ();
    #do some work for easy living below
    my $f_base  = basename( $reads[0], @in_suffixes ); #this contains the split value!
    $f_base =~ s/\_\d+$//; #no longer has split value
    #do the work
    if( $paired_end ) {
	@f_sorted   = @sorted;
	@r_sorted   = @f_sorted;
	for( @r_sorted ){
	    s/${f_base}/${mate_base}/;
	}
	my $f_out_stem = $f_base;
	my $f_out_path = File::Spec->catfile( $result_dir, $f_out_stem . ".fq" );
	print("cat @f_sorted > $f_out_path\n");
	system("cat @f_sorted > $f_out_path");
	my $r_out_stem = $mate_base;
	my $r_out_path = File::Spec->catfile( $result_dir, $r_out_stem . ".fq" );
	print("cat @r_sorted > $r_out_path\n");
	system("cat @r_sorted > $r_out_path");
    } else {
	@f_sorted    = @sorted;
	my $out_stem = $f_base;
	my $out_path = File::Spec->catfile( $result_dir, $out_stem );
	print("cat @f_sorted > $out_path\n");
	system("cat @f_sorted > $out_path");
    }
}

# settings is a hashref with the following pointers:
# setting->method->data_type
# where method is a function name (e.g., bmtagger)
# data_type is input, output, log, or something custom to function

sub _set_settings{
    my ($ra_qc_methods, $ra_trim_methods,
	$ra_filter_methods, $ra_derep_methods ) = @_;
    my $run_raw_qc   = 1;
    my $split_reads  = 1;
    my $run_trim     = 1;
    my $run_filter   = 1;
    my $cat_reads    = 1;
    my $derep        = 1;
    my $check_qc     = 1;
    my $fasta_cleaned= 1;
    my $settings = _build_settings($ra_qc_methods, $ra_trim_methods,
				   $ra_filter_methods, $ra_derep_methods ) ;
    my @params   = ( $run_raw_qc, $split_reads, $run_trim, $run_filter,
		     $cat_reads, $derep, $check_qc, $fasta_cleaned );
    $settings->{"parameters"} = \@params;
    return $settings;
}

sub _build_settings{
    my ($ra_qc_methods, $ra_trim_methods,
	$ra_filter_methods, $ra_derep_methods ) = @_;
    my @qc_methods     = @{ $ra_qc_methods };
    my @trim_methods   = @{ $ra_trim_methods };
    my @filter_methods = @{ $ra_filter_methods };
    my @derep_methods  = @{ $ra_derep_methods };
    my $settings     = ();
    #the order matters here!
    my @keys = ( @qc_methods, "split_reads", @trim_methods,
		      @filter_methods, "cat_reads", @derep_methods,
		      "check_qc", "fasta_cleaned" );
    for( my $i=0; $i<scalar(@keys); $i++ ){
	my $key = $keys[$i];
	if( $key eq "fastqc" ){
	    $settings->{$key}->{"input"}  = "";
	    $settings->{$key}->{"output"} = "fastqc_raw";
	    $settings->{$key}->{"log"}    = "fastqc_raw_log";
	    next;
	} elsif( $key eq "fasta_cleaned" ){
	    my $prior = $keys[$i-2]; #check_qc is always just before fasta_cleaned
	    $settings->{$key}->{"input"}  = $prior;
	    $settings->{$key}->{"output"} = $key;
	    $settings->{$key}->{"log"}    = $key . "_log";
	    next;
	} else {
	    my $next = $keys[$i+1];
	    my $prior = $keys[$i-1]; 
	    if( $prior eq "fastqc" ){
		$settings->{$key}->{"input"}  = "";       
	    } else {
		$settings->{$key}->{"input"}  = $prior;
	    }
	    $settings->{$key}->{"output"} = $key;
	    $settings->{$key}->{"log"}    = $key . "_log";
	    next;
	}
    }
    return $settings;
}

sub _get_root{
    my $exec_path = $0;
    my ($name, $root) = fileparse( $exec_path );
    return $root;
}

sub _check_vars{
    my ( $in_fastq,  $pair_fastq, $out_dir, 
	 $index_dir, $index_basename ) = @_;
    if( !defined $in_fastq ){
	die( "You must point me to a fastq file using -1\n" );
    }
    if(	! -e $in_fastq ) {
	die( "You must point me to a fastq file using -1. You specified:\n" .
	     $in_fastq . "\n which I cannot locate\n" );
    }
    if( defined( $pair_fastq ) &&
	! -e $pair_fastq ){
	die( "I can't located the mate pair file you specified with -2. You specified:\n" .
	     $pair_fastq . "\n" );
    }
    if( !defined( $out_dir ) ){
	die( "You must specify an output directory with -o\n" );
    }
    if( !defined( $index_dir ) ){
	die( "You must specify an index_directory with -d\n" );
    }
    if( ! -d $index_dir ){
	die( "I cannot access the index_directory. You specified:\n" .
	     $index_dir . "\n" );
    }
    if( !defined( $index_basename) ){
	die( "You must provide an index_basename -n so that I know which genome in index_db to use." );
    }
    my @files = glob( $index_dir . "/" . $index_basename . "*" );
    if( !@files ){
	die( "I cannot located an indexed genome database in the index_directory " .
	     "using the index_basename you specified, which was\n" .
	     $index_basename . "\n" );	     
    }
}

sub _parse_method_list{
    my $method_list = shift;
    my $type        = shift;
    my @methods = ( $method_list );
    if( $method_list =~ m/\,/ ){
	@methods = split( "\,", $method_list );
    }
    #check the naming convetion
    foreach my $method( @methods ){
	if( $type eq "qc" ){
	    if( $method ne "fastqc" ){
		die( "I don't know how to process $method\n" );
	    }
	} elsif( $type eq "trim" ){
	    if( $method ne "trimmomatic" && 
		$method ne "prinseq"     ){
		die( "I don't know how to process $method\n" );
	    }
	} elsif( $type eq "filter" ){
	    if( $method ne "bowtie2"  &&
		$method ne "bmtagger" &&
 		$method ne "deconseq" ){
		die( "I don't know how to process $method\n" );
	    }
	} elsif( $type eq "derep" ){
	    if( $method ne "fastuniq" &&
		$method ne "prinseq-derep" ){
		die( "I don't know how to process $method\n" );
	    }
	}
    }
    return \@methods;
}

sub _run_fastuniq{

    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $overwrite   = $args->{"overwrite"};
    my $fastuniq    = $args->{"fastuniq"};
    my $paired_end  = $args->{"paired_end"};
    my $tmp_dir     = $args->{"tmp_dir"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };

    my $in_file   = $reads[0]; #includes split value!
    $in_file      =~ s/\_\d+\.fq$/\.fq/; #no longer does
    my $f_in  = File::Spec->catfile( $in_dir, $in_file );
    my $f_base = basename( $in_file, @in_suffixes ); #this contains the split value!
    my $r_in = $f_in;
    #we have to create an input_list file for fastuniq
    my $input_list = File::Spec->catfile( $tmp_dir, "fu.tmp" );
    open( TMP, ">$input_list" ) || die "Can't write to $input_list: $!\n";
    print TMP $f_in . "\n" . $r_in;
    close TMP;
    my $log      = File::Spec->catfile( $log_dir, $in_file . ".log" );
    my $out_file = $in_file; #just the file name 
    my $f_out    = File::Spec->catfile( $result_dir, $out_file );
    my $r_out    = $f_out;
    $r_out       =~ s/$f_base/${mate_base}/;
    my $cmd = "${fastuniq} -i ${input_list} -t q -o ${f_out} -p ${r_out} ";
    $cmd    .= "&> $log ";
    print "$cmd\n";
    system( $cmd );	    
}

sub _run_trimmomatic{
    my( $args ) = @_;

    my $in_dir      = $args->{"in_dir"};
    my @reads       = @{ $args->{"file_names"} };
    my $result_dir  = $args->{"result_dir"};
    my $log_dir     = $args->{"log_dir"};
    my $nprocs      = $args->{"nprocs"};
    my $overwrite   = $args->{"overwrite"};
    my $trimmomatic = $args->{"trimmomatic"};
    my $paired_end  = $args->{"paired_end"};
    my $mate_base   = $args->{"mate_basename"};
    my @in_suffixes = @{ $args->{"in_suffixes"} };

    my $pm = Parallel::ForkManager->new($nprocs);
    for( my $i=1; $i<=$nprocs; $i++ ){
	my $pid = $pm->start and next;
	#do some housekeeping
	my $read = $reads[$i-1];          
	#prepare input files
	my ( $f_in, $r_in );
	my $f_mate = $read;
	$f_in  = File::Spec->catfile( $in_dir, $f_mate );
	my $f_base = basename( $read, @in_suffixes ); #this contains the split value!
	if( $paired_end ){
	    $r_in  = $f_in;
	    $r_in     =~ s/$f_base/${mate_base}_${i}/;
	}
	#prepare log file
	my $log = File::Spec->catfile( $log_dir, $f_mate . ".log" );
	#prepare output files
	my $f_out = File::Spec->catfile( $result_dir, $read );
	my $f_single = $f_out . ".singletons";
	my( $r_out, $r_single );
	if( $paired_end ){
	    $r_out = $f_out;
	    $r_out     =~ s/$f_base/${mate_base}_${i}/;
	    $r_single = $r_out . ".singletons";
	}
	my $compressed = 0;	   
	my $gz_file = File::Spec->catfile( $f_in . ".gz" );
	if( -e $gz_file ){
	    $compressed = 1;
	}
	my $cmd = "java -jar ${trimmomatic} ";
	if( $paired_end ){
	    $cmd .= "PE -phred33 ${f_in} ${r_in} ${f_out} ${f_single} ${r_out} ${r_single} ";
	} else {
	    $cmd .= "SE -phred33 ${f_in} ${f_out} ";
	}
	#now all of the thesholds
	$cmd .= "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
	$cmd    .= "&> $log ";
	print "$cmd\n";
	system( $cmd );
	$pm->finish; # Terminates the child process
    }
    $pm->wait_all_children;    
}
