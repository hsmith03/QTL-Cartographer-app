# 
#  Run Model 8 iteration
#
local ($opt_b,$opt_X,$opt_S,$opt_i,$opt_m,$opt_H,$opt_h) = ($TheBinDir,$TheStem,10.0,10,25,10,0);

my $stem;
my $bindir;
my $siglevel;
my $iterations;
my $maxnbp;
my $ihypo;
my $i;
my $j;
my $output;

$usage = $usage . "$0: [-b bindir ] [-X stem ] [-S sigthresh ] [-i iterations ] 
  [-H hypothesis ] [-m maxnbp ] [-h]  
	The default values are:
		bin dir    = $TheBinDir
		stem       = $TheStem
		sigthresh  = $opt_S
		iterations = $opt_i
		hypothesis = $opt_H
		maxnbp     = $opt_m  \n";

getopts "b:X:S:i:m:H:h" ;

die $usage if ( $opt_h == 1 );

$bindir = $opt_b;
$stem = $opt_X;
$siglevel = $opt_S;
$iterations = $opt_i;
$maxnbp = $opt_m;
$ihypo = $opt_H;

`$bindir/Qstats -X $stem -A -V` ;
`$bindir/Zmapqtl -A -V -M 3`;
`$bindir/Eqtl -A -V -S $siglevel -I Z -H $ihypo`;
#
#  Save the original files
#
rename "$stem.eqt", "$stem.eqt.0" ;
rename "$stem.z", "$stem.z.0" ;
#
#  Use model 8 iteratively with cofactors from previous run.
#
for ( $i=1; $i<=$iterations; $i++ ) {
  print "\nAt iteration $i  ";
  `$bindir/Zmapqtl -A -V   -M 8 -n $maxnbp`;
  $j = $i-1 ;
  rename "$stem.sr", "$stem.sr.$j";
  `$bindir/Eqtl -A -V -S $siglevel -I Z -H $ihypo`;
  rename "$stem.eqt", "$stem.eqt.$i";
  rename "$stem.z", "$stem.z.$i";
  $output = `$bindir/SRcompare.pl -f $stem.sr.$j < $stem.sr `;
  print $output ; 
}
$i -=1;
rename "$stem.sr", "$stem.sr.$i";
print "\n\nFinished\n";

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Model8.pl - Iterate Zmapqtl model 8 to determine a stable set of cofactors

=head1 SYNOPSIS

  Model8.pl  [-b bindir] [-X stem] [-S sig. threshold] [-i iterations] 
    [-m max nbp] [-H hypothesis] [-h]   


=head1 DESCRIPTION

B<Model8.pl> iterates using B<Zmapqtl> and  B<Eqtl> to determine a stable set of
cofactors for composite interval mapping.   First, interval mapping is run and the
nearest markers to significant peaks are identified.  These markers are used as 
cofactors in the first iteration of composite interval mapping.  A new set of 
cofactors are identified by proximity to the likelihood peaks and the process is
repeated.     

=head1 OPTIONS

=over 4

=item  B<-b> 

This option requires the path to the B<QTL Cartographer> binaries and perl scripts.  

=item B<-i>

requires an integer to control how many iterations you want to do.

=item  B<-X> 

This option allows you to specify the filename stem.

=item  B<-S> 

This option requires a real number to indicate the significance threshold for the likelihood ratio.

=item  B<-m> 

This option allows you to specify the maximum number of background parameters in composite interval 
mapping.  

=item  B<-H> 

Use this option to specify which hypothesis test you want to use.  The usual values are 10 or 30.  

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.



=back

=head1 EXAMPLE

Suppose the files F<mletest.map> and F<mletest.cro> are in the current working
directory.   

    % Model8.pl  -b /home/basten/bin  -X mletest -S 13.0 -m 25 -H 10 -i 15

Will assume that the B<QTL Cartographer> programs are in F</home/basten/bin>. 
It will use a significance threshold of 13.0 and allow for up to 25 markers in 
composite interval mapping.   It will iterate 15 times in an attempt to find a stable
set of cofactors.  This script uses another script called B<SRcompare.pl> which compares
the set of cofactors (in F<SRmapqtl.out> format) between consecutive runs and reports
how many cofactors have been added or deleted.  

=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Qstats(1)>, B<Eqtl(1)>, B<SRcompare.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
