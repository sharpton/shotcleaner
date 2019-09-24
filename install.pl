#!/usr/bin/perl -w 

use strict;
use File::Path;
use File::Spec;
use File::Basename;
use Cwd;
use Getopt::Long qw(GetOptionsFromString GetOptionsFromArray);
use Carp;
use Data::Dumper;

$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

my $testing = 0;

my $options = get_options( \@ARGV );

my $perlmods  = $options->{"perlmods"};
my $algs      = $options->{"algs"};
my $clean     = $options->{"clean"};
my $get       = $options->{"get"};
my $build     = $options->{"build"};
my $test      = $options->{"test"};
my $db        = $options->{"db"};
my $source    = $options->{"source"};
my $all       = $options->{"all"};
my $iso_alg;
if( defined( $options->{"iso-alg"} ) ){
    $iso_alg = $options->{"iso-alg"};
}

my $dir = _get_root();

our $ROOT = $dir;
my $root = $ROOT;

my $inc = "${root}/inc/"; #location of included source code (cpanm)
my $pkg = "${root}/pkg/"; #location that we'll place algorithms
my $ext = "${root}/ext/"; #location of installed external perl modules (separate so we can wipe)
my $bin = "${root}/bin/"; #location of installed executables (symlinks)

system( "mkdir -p $pkg" );

######################
# Install Perl Module Dependencies using cpanm
if( $perlmods ){
 my @mods = (
     "Carp",
     "File::Basename",
     "File::Copy",
     "File::Path",
     "File::Spec",
     "File::Spec::Functions",
     "Parallel::ForkManager",
     "Getopt::Long",
     );
 
 foreach my $mod( @mods ){
     my $cmd = "perl ${bin}/cpanm -L ${ext} ${mod}";
     system( $cmd );
 }
}

######################
######################
# Install 3rd Party Applications
if( $algs ){
    my $alg_data = parse_alg_data( $root . "/installer_alg_data.txt", $source, $test);
    foreach my $alg( keys( %$alg_data ) ){
	if( defined( $iso_alg ) ){
	    next if( $alg ne $iso_alg );
	}
	my $src   = $alg_data->{$alg}->{"src"};
	my $stem  = $alg_data->{$alg}->{"stem"};
	my $links = $alg_data->{$alg}->{"links"};
	my $build = $alg_data->{$alg}->{"build"};
	my $bins  = $alg_data->{$alg}->{"bins"};	    
	my $dload = $alg_data->{$alg}->{"download"};	    

	print "Installing ${alg}...\n";
	#need a better way to automate this step...
	my $alg_pkg = "${pkg}/${stem}/";
	if( $alg eq "blast" && $source ){
	    $alg_pkg = "${pkg}/${stem}/c++/";
	}
	if( $clean ){
	    clean_src( $links,
		       $alg_pkg, 
		       "${pkg}" . basename( $src ) 
		);
	}    
	if( $get ){
	    get_src( $src, $pkg, $dload );
	    decompress_src( $pkg, basename( $src ) );
	}	
	if( $build ){
	    unless( $build eq "NA" ){
		build_src( $alg_pkg, $build );	    
	    }
	}	
	foreach my $target( split( "\;", $bins ) ){
	    link_src( $alg_pkg . $target, 
		      $bin );
	}
	print "\t...done with ${alg}\n";
    }
    print( "The requested items have been installed. If you haven't already, please be sure " .
	   "to make the following additions to your ~/.bash_profile:\n" .
	   "###################################\n"               . 
	   "export PATH=\$PATH:${root}/bin/\n"         .
	   "export PERL5LIB=\$PERL5LIB:${root}/lib/\n" .
	   "###################################\n" 
	);
    print( "If you like, you can copy and paste the following command to automate the above entries:\n" . 
	   "###################################\n"                                         . 
	   "echo 'export PATH=\$PATH:${root}/bin/' >> ~/.bash_profile\n"         .
	   "echo 'export PERL5LIB=\$PERL5LIB:${root}/lib/' >> ~/.bash_profile\n" .
	   "###################################\n" 
	);
}



###############
###############
### SUBROUTINES

    
sub build_src{
    my $loc = shift; #path to the downloaded, decompressed source, relative to shotmap root
    my $cmds =  shift; #semicolon delimited list of commands

    print "building source...\n";
    print "$loc\n";
    chdir( $loc );
    my @commands = split( "\;", $cmds );
    print $loc . "\n";
    foreach my $command( @commands ){
	print "$command\n";
	system( "$command" );
    }
    chdir( $ROOT );
    return;
}

sub check_name{
    my $name = shift;
    if( ! defined( $name ) ){
	die "Could not find a name when parsing alg_data!";
    }
}

sub decompress_src{
    my $loc = shift; #path to directory containing download, relative to shotmap root
    my $stem = shift; #name of the downloaded file

    print "decompressing source...\n";
    chdir( $loc );
    #if( ! -e $stem ){
	#die "For some reason, I didn't download $stem into $loc. " . 
	#    "Please try to reproduce this error before contacting " . 
	#    "the author for assistance\n";
    #}
    if( $stem =~ m/\.tar\.gz/ ){
	system( "tar xzf $stem" );
    }
    elsif( $stem =~ m/\.zip/ ){
	system( "unzip $stem" );
    }
    else{
	warn "No decompression needed for $stem\n";
    }
    chdir( $ROOT );
    return;
}

sub dereference_options{
    my $options = shift; #hashref
    foreach my $key( keys( %$options )){
	my $value = $options->{$key};
	if( defined( $value ) ){
	    $options->{$key} = ${ $value };
	}
    }
    return $options;
}


sub get_src{
    my $url = shift;
    my $loc = shift;
    my $method = shift; #git, wget

    print "downloading source...\n";
    chdir( $loc );
    if( $method eq "git" ){
	system( "git clone $url" );
    }
    elsif( $method =~ m/(wget.*)/ ){ #bmtagger requires extra wget options
	system( "$1 $url" );
    }
    else{
	die( "I don't know how to deal with the download method $method\n");
    }
    if( $url =~ m/bmtagger/ ){ #need to move it to keep things clean...
	my $current = $loc . "ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/";
	my $move_to = $loc . "bmtagger/";
	system( "mv $current $move_to" );
    }
    chdir( $ROOT );
    return;
}

sub clean_src{
    my $link_names = shift;
    my $src_dir    = shift;
    my $dl_stem    = shift; #the downloaded compressed file, if it exists

    print "\tcleaning old install...\n";
    foreach my $link( split( "\;", $link_names ) ){
	print( "removing $link\n" );
	unlink( $link );
    }
    if( -d "${src_dir}" ){
	rmtree( "${src_dir}" );
    }
    if( defined( $dl_stem )) { 
	if( -e "${dl_stem}"){
	    unlink( "${dl_stem}" );
	}
    }

    chdir( $ROOT );
    return;
}

sub get_options{
    my $ra_args  = shift;
    my @args     = @{ $ra_args };

    #DEFAULT VALUES
    my $perlmods  = 0; #should we install perl modules
    my $rpackages = 0; #should we install R modules
    my $algs      = 0; #should we install 3rd party gene prediction/search algorithms?
    my $clean     = 0; #wipe old installations of algs?
    my $get       = 0; #download alg source code?
    my $build     = 0; #build alg source code?
    my $test      = 0; #should we run make checks during build?
    my $db        = 0; #should we install myql libraries? 
    my $all       = 1; 
    my $source    = 0; #should we build from source instead of x86 libraries
    my $iso_alg;       #let's only build one specific algoritm.

    my %ops = (
	"use-db"   => \$db, #try to build the libraries needed for mysql communication
	"source"   => \$source, #not currently implemented
	"clean"    => \$clean,
	"test"     => \$test,
	"build"    => \$build,
	"get"      => \$get,
	"algs"     => \$algs,
	"rpacks"   => \$rpackages,
	"perlmods" => \$perlmods,
	"all"      => \$all,
	"iso-alg"  => \$iso_alg,
	);
    
    my @opt_type_array = (
	"use-db!",
	"source!",
	"clean!",
	"test!",
	"build!",
	"get!",
	"algs!",
	"rpacks!",
	"perlmods!",
	"iso-alg:s"
	);
    
    #add --source options later..
    #print "Note that this installer attempts to install precompiled x86 binaries when possible. If ".
    #	"this doesn't work for your architecture, please rerun the installer and invoke --source\n";
    
    GetOptionsFromArray( \@args, \%ops, @opt_type_array );
    %ops = %{ dereference_options( \%ops ) };

    if( $ops{"perlmods"} || 
	$ops{"rpacks"}   ||
	$ops{"algs"}     ||
	$ops{"clean"}    ||
	$ops{"get"}      ||
	$ops{"build"}    ||
	$ops{"test"}     ||
	$ops{"source"}   ){
	print "You specified a modular aspect of installation, so I will not run the entire " .
	    "installation pipeline\n";
	$ops{"all"} = 0;
    }

    #If we clean and want to build, we must also get.
    if( $ops{"clean"} && 
	$ops{"build"} ){
	$ops{"get"} = 1;
    }

    if( $ops{"all"} ){
	$ops{"perlmods"} = 1;
	$ops{"rpacks"}   = 1;
	$ops{"algs"}     = 1;
	$ops{"clean"}    = 1;
	$ops{"get"}      = 1;
	$ops{"build"}    = 1;
	$ops{"test"}     = 1;
    }
    return \%ops;
}

sub link_src{
    my $targets = shift; #must be path to target relative to shotmap root, semicolon sep list
    my $linkdir = shift;
    
    print "creating symlinks...\n";
    chdir( $linkdir );
    foreach my $target( split( "\;", $targets ) ){
	print( $linkdir . "\n");
	print( $target  . "\n" );
	system( "chmod a+x ${target}" );
	system( "ln -s ${target}" );
    }
    chdir( $ROOT );
    return;
}

sub parse_alg_data{
    my $file     = shift;
    my $source   = shift;
    my $test     = shift;
    my $alg_data = ();
    open( IN, $file ) || die "Can't open $file for read: $!\n";
    my $name;
    while(<IN>){
	chomp $_;
	if( $_ =~ /^name/ ){
	    (my $string, $name ) = split( "\:", $_ );

	}  elsif( $_ =~ /^download/ ){
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"download"} = $value;
	    	   
	} elsif( $_ =~ /^x86\:/ ){
	    next if( $source );
	    check_name( $name );
	    #srcs may have : in the http/ftp string
	    my ($string, $value, @rest ) = split( "\:", $_ );
	    $value = join( ":", $value, @rest );
	    $alg_data->{$name}->{"src"} = $value;
	} elsif( $_ =~ /^src\:/ ){
	    next unless( $source );
	    check_name( $name );
	    my ($string, $value, @rest ) = split( "\:", $_ );
	    $value = join( ":", $value, @rest );
	    $alg_data->{$name}->{"src"} = $value;

	} elsif( $_ =~ /^x86stem/ ){
	    next if( $source );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"stem"} = $value;
	} elsif( $_ =~ /^srcstem/ ){
	    next unless( $source );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"stem"} = $value;

	} elsif( $_ =~ /^x86bins/ ){
	    next if( $source );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"bins"} = $value;
	} elsif( $_ =~ /^srcbins/ ){
	    next unless( $source );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"bins"} = $value;

	} elsif( $_ =~ /^testcmds/ ){
	    next unless( $test );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"build"} = $value;
	} elsif( $_ =~ /^srccmds/ ){
	    next if( $test );
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"build"} = $value;

	} elsif( $_ =~ /^installed/ ){
	    check_name( $name );
	    my ($string, $value ) = split( "\:", $_ );
	    $alg_data->{$name}->{"links"} = $value;
	} elsif( $_ =~ m/^$/ ){
	    $name    = undef;
	}
    }
    close IN;
    return $alg_data;
}

sub _get_root{
    my $exec_path = $0;
    #exec_path is relative - make absolute
    $exec_path = File::Spec->rel2abs( $exec_path );
    my ($name, $root) = fileparse( $exec_path );
    return $root;
}
