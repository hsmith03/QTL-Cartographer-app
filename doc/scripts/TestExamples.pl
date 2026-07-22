#
#  TestExamples.pl
#
#   Copyright (C) 2002 Christopher J. Basten 
#  Usage governed by the terms of the GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
#
local (  $opt_s, $opt_D, $opt_t,$opt_e,$opt_b,$opt_h) = (  "/bin", 0,"/Users/basten/test/examples","/Users/basten/software/qtlcart/QTLCartUnix/example", $TheBinDir,0);
#
my $stem;
my $bin;
my $model;
my $reps;
my $ihypo;
my $i;
my $j;
my $u;
my $ext;
my %parameters;
my $traits;
my $xtraits;
my $resource = "qtlcart.rc";
my $verbiage = "Verbiage.txt";
my $gobalopt;
#
$usage = $usage . "$0: [-D debug] [-t test dir ] [-e example dir ][-b bin dir ]   [-h] 
	The default values are:
		test dir      = $opt_t
		example dir   = $opt_e
		bin dir       = $TheBinDir\n";
getopts "s:D:t:e:b:h" ;
die $usage if ( $opt_h == 1 );
#
#  These are the real data files
#
my %datafiles = (
lauriem4c => lauriem,
lauriem6c => lauriem,
lauries4c => lauriem,
lauries6c => lauriem,
nuzhdinc => nuzhdinm,
rauhall => rauhmap,
rauham => rauhmap,
rauhratio => rauhmap,
rauhrl => rauhmap,
rauhrm => rauhmap,
realdatc => realdatm,
vieirac   => vieiram,
weber519c => weber519m,
weber701c => weber701m,
);
#
#   Create the test directory and cd to it
#
`$sysmkdir   $opt_t ` unless -d $opt_t ;
chdir $opt_t || die "Can't cd to $opt_t\n" ;
unlink($verbiage) if -f $verbiage ;
`$systouch $verbiage`; 
$globalopt = "-A -D$opt_D  >> ../$verbiage";
#
#  Copy, transform and anaylze the data files
#
foreach $u ( sort keys %datafiles ) {
  $ihypo = 10;
  $stem = $u;  
  &prep_space; 
  `$syscp $opt_e/$u.inp .`;
  `$syscp $opt_e/$datafiles{$u}.inp .`;
  `$opt_b/Rmap  -X $u -i $datafiles{$u}.inp $globalopt`;
  `$opt_b/Rcross -i $u.inp $globalopt`;
  &get_resource;
  $traits = $traits + 1;
  &do_basics;
  chdir ".." ;
}
#
#  mletest stuff
#
$ihypo = 10;
$stem="mletest"; $traits = 1;
&prep_space;
`$syscp $opt_e/$stem.map .`;
`$syscp $opt_e/$stem.cro .`;
&do_basics;
chdir "..";
#
#  sample stuff
#
$ihypo = 30;
$stem = "sample";  $traits = 1;
&prep_space;
`$syscp $opt_e/sample.mps .`;
`$syscp $opt_e/sample.raw .`;
`$opt_b/Rmap -X $stem -i sample.mps $globalopt`;
`$opt_b/Rcross -i sample.raw $globalopt`;
&do_basics;
unlink("qtlcart.rc");  #  check out emap.  
`$opt_b/Emap -X emap -i sample.raw -S 0.0 -M 13 $globalopt`;
chdir "..";
#
#  1 trait simulation
#
$ihypo = 30;
$stem = "sim1";  $traits = 1;
&prep_space;
`$opt_b/Rmap -X $stem -vd 2.0 -vm 2.0 -t 10.0 $globalopt`;
`$opt_b/Rqtl -q 3 $globalopt`;
`$opt_b/Rcross -c SF2 -H 0.9 $globalopt`;
&do_basics;
chdir "..";
#
#  2 trait simulation
#
$ihypo = 30;
$stem = "sim2";  $traits = 3;
&prep_space;
`$opt_b/Rmap -X $stem -vd 2.0 -vm 2.0 -t 10.0 $globalopt`;
`$opt_b/Rqtl -t 2 -q 3 $globalopt`;
`$opt_b/Rcross -c RF4 -H 0.9 $globalopt`;
&do_basics;
chdir ".."; 

sub prep_space {

`$sysrmdir $stem` if -d $stem ;
`$sysmkdir $stem` ;
chdir "$stem";
print "\n\tNow on $stem...";


}
#
#     Subroutine to get the parameters from qtlcart.rc
#
sub get_resource {
	#
	#  Need to get the number of traits...might as well get all the parameter values
	#  
	if ( -f $resource ) {
	  open RESOURCE, $resource or die "Can't find resource file $resource : $!\n";
	  while ( <RESOURCE> ) {  
		if ( /^-/ ) {
		  chomp;
		  @fields = split;
		  $parameters{$fields[0]} = $fields[1] ;     
		}  
	  }
	  close RESOURCE;
	  $traits = $parameters{"-traits"} ;
	}
	else {
	  $traits = 0;
	}
}

sub do_basics {
local $outputfile;
local $inputfile;
`$opt_b/Qstats -X $stem $globalopt`;
`$opt_b/LRmapqtl  -t $traits $globalopt`;
`$opt_b/SRmapqtl -t $traits   $globalopt`;
`$opt_b/Zmapqtl  -t $traits  $globalopt`;
`$opt_b/JZmapqtl -I $ihypo -t $traits $globalopt` if $traits > 1; 
`$opt_b/Eqtl   -H $ihypo $globalopt`;
`$opt_b/Preplot   $globalopt`;
`$opt_b/JZmapqtl -t $traits  -M 9 $globalopt`; 
`$opt_b/MultiRegress -t $traits   $globalopt`; 
$outputfile = $stem . "Phase0.mqt";
$inputfile = $stem . ".mr";
`$opt_b/Rqtl -i $inputfile -o $outputfile   $globalopt`; 
`$opt_b/MImapqtl -t $traits -p 1 -I sMPrtseC $globalopt`; 
`$opt_b/Prune  -b 1 $globalopt`;
}


#
print "\n";


exit(0);

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

TestExamples.pl - Run a basic analysis on the example files

=head1 SYNOPSIS

  TestExamples.pl [-D debug] [-t test dir ] [-e example dir ][-b bin dir ] [-h]   


=head1 DESCRIPTION

B<TestExamples.pl> is a script that will use the B<QTL Cartographer> commands to 
run basic analyses on all of the example data sets.  

=head1 OPTIONS



=over 4

=item  B<-b> 

This option requires the path to the B<QTL Cartographer> binaries and perl scripts.  

=item B<-t>

This option requires the path to the base directory for doing the analyses.
This directory will be created if it doesn't exist.  Each data set will have
a subdirectory in the base directory for analysis.

=item  B<-e> 

Use this to specify the location of the example subdirectory.

=item  B<-D> 

Allows you to set the debug level.   The default is zero.  You can set it to
any of 0, 1, 2, or 3.


=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 EXAMPLE

Suppose that the example files are in F</home/user/QTLCartUnix/example> and you want to 
run all the tests in a base directory F</home/user/test>.   Suppose further that
you have installed the B<QTL Cartographer> binaries in F</usr/local/bin>.  Use

    % TestExamples.pl -b /usr/local/bin -t /home/user/test \
        -e /home/user/QTLCartUnix/example

This will create the directory F</home/user/test> if it doesn't exist, and then
create subdirectories for each data set.  It will copy the data file and genetic linkage
map from the example subdirectory into the data set subdirectory of F</home/user/test>.
Then, it will convert the map and data files into the proper formats and run
B<Qstats>, B<LRmapqtl>, B<SRmapqtl>, B<Zmapqtl> on all traits.  If there are multiple traits,
B<JZmapqtl> will be run followed by B<Eqtl> and B<Preplot>.  Then, B<JZmapqtl> will convert
the data so that B<MultiRegress> can generate an initial model that B<Rqtl> translates and
B<MImapqtl> refines.   This will be done for all the data sets in the example directory that
have filename extension F<.inp>.   

The script also runs analyses for the F<mletest> and F<sample> data sets, and 
does two simulations:  One with a single trait and one with two traits.

At the end of the script, you will find a file F</home/user/test/Verbiage.txt> that 
contains all of the verbal output from the programs.  This file is deleted by the
script  if you run it again.  

=head1 NOTE

You will need at least 36 megabytes of disk space to run this program.   

=head1 SEE ALSO 

B<QTLcart(1)> and all programs referenced therein.  

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut

