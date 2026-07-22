#
#  Do the basic analyses for a QTL mapping data set
#
local ($opt_b,$opt_X,$opt_M, $opt_m,$opt_d,$opt_t,$opt_S,$opt_H,$opt_V,$opt_h) = ($TheBinDir,$TheStem,6,"","",1,10.0,10,0,0);

my $stem;
my $bindir;
my $siglevel;
my $mapfile;
my $datafile;
my $traits;
my $ihypo;
my $i;
my $j;
my $output;
my $model;

$usage = $usage . "$0: [-b bindir ] [-X stem ] [-m mapinput ] [-d datainput ] [-t traits ] 
  [-S sigthresh ]  [-H hypothesis ]  [-V] [-h]   
	The default values are:
		bin dir    = $TheBinDir
		stem       = $TheStem
		mapinput   = $opt_m
		datainput  = $opt_d
		hypothesis = $opt_H
		traits     = $opt_t
		sigthresh  = $opt_S  \n";

getopts "b:X:m:d:t:S:H:h" ;

die $usage if ( $opt_h == 1 );

$opt_V = 1 - $opt_V ;
$bindir = $opt_b;
$stem = $opt_X;
$siglevel = $opt_S;
$mapfile = $opt_m;
$datafile = $opt_d;
$traits = $opt_t;
$ihypo = $opt_H;
$model = $opt_M;

print "\n   Dobasics.pl running with the following parameters:  
		bin dir    = $TheBinDir
		stem       = $stem
		mapinput   = $opt_m
		datainput  = $opt_d
		hypothesis = $opt_H
		traits     = $opt_t
		sigthresh  = $opt_S\n   It is $date_time\n"   if $opt_V == 1 ;



print "\n   Converting the map in $mapfile" if $opt_V == 1 && length($mapfile) > 0 ;
`$bindir/Rmap   -X $stem -i $mapfile -A -V` if length($mapfile) > 0 ;
print "\n   Converting the data in $datafile" if $opt_V == 1 && length($datafile) > 0 ;
`$bindir/Rcross -X $stem -i $datafile -A -V` if length($datafile) > 0 ;
print "\n   Running Qstats " if $opt_V == 1 ;
`$bindir/Qstats -X $stem  -A -V` ;
print "\n   Running LRmapqtl " if $opt_V == 1 ;
`$bindir/LRmapqtl -t $traits  -A -V` ;
print "\n   Running SRmapqtl " if $opt_V == 1 ;
`$bindir/SRmapqtl -t $traits  -A -V` ;
print "\n   Running Zmapqtl with model $model" if $opt_V == 1 ;
`$bindir/Zmapqtl  -M $model -t $traits -A -V `;
print "\n   Running JZmapqtl with model $model" if $opt_V == 1 && $traits > 1;
`$bindir/JZmapqtl -M $model -t $traits -I $ihypo  -A -V ` if   $traits > 1   ;
print "\n   Running Eqtl for hypothesis $ihypo at significance level $siglevel " if $opt_V == 1 ;
`$bindir/Eqtl -I ZM  -M $model -S $siglevel -H $ihypo -A -V`;
print "\n   Running Preplot" if $opt_V == 1 ;
`$bindir/Preplot -H $ihypo -A -V `;
print "\n   Converting data with JZmapqtl for input into MultiRegress " if $opt_V == 1  ;
`$bindir/JZmapqtl -M 9 -t $traits -I $ihypo  -A -V `    ;
print "\n   Running MultiRegress" if $opt_V == 1  ;
`$bindir/MultiRegress  -t $traits -I $ihypo  -A -V `    ;
print "\n   Converting MultiRegress output for input into MImapqtl " if $opt_V == 1  ;
`$bindir/Rqtl  -i $stem.mr -o $stem.imq   -A -V `    ;
print "\n   Running MImapqtl" if $opt_V == 1  ;
`$bindir/MImapqtl  -E $stem.imq  -I sMPrtSEC -t $traits -A -V `    ;
print "\n   Dobasics.pl is finished\n" if $opt_V == 1 ;

exit 0 ;

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Dobasics.pl - Do a basic set of analyses for a QTL mapping data set

=head1 SYNOPSIS

  Dobasics.pl  [-b bindir] [-X stem] [-S sig. threshold] [-m mapinput] 
    [-d datainput]  [-t trait] [-H hypothesis] [-V] [-h]   


=head1 DESCRIPTION

B<Dobasics.pl> does a basic analysis of a QTL mapping data set.   
Map and data files are translated first.  Then, basic single marker
analyses are followed by interval mapping and composite interval mapping.
B<Eqtl> produces a summary of results and B<Preplot> prepares the results for
plotting with B<Gnuplot>.   

=head1 OPTIONS

=over 4

=item  B<-b> 

This option requires the path to the B<QTL Cartographer> binaries and perl scripts.  

=item B<-t>

requires an integer telling which trait to analyze.  If you want to do them all, set this to
one more than the number of traits.   

=item  B<-X> 

This option allows you to specify the filename stem.

=item  B<-S> 

This option requires a real number to indicate the significance threshold for the likelihood ratio.

=item  B<-m> 

This allows you to specify an input file for the map data.  It should be in the format of F<map.inp>
or B<MAPMAKER/EXP> output.  If not used, the script assumes that the map file already exists in the
correct format.  

=item  B<-d> 

This allows you to specify an input file for the mapping data.  It should be in the format of F<cross.inp>
or B<MAPMAKER/EXP> raw format.  If not used, the script assumes that the date file already exists in the
correct format.  

=item  B<-H> 

Use this option to specify which hypothesis test you want to use.  The usual values are 10 or 30.  

=item B<-V>

requires no operand.  If used, it turns off messages indicating the progress of the script.

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.



=back

=head1 ANALYSES

Assume that we will use the filename stem I<qtlcart> for the following.   The script does the tasks in
this order:

=over

=item B<Rmap>

Convert the map found in the file specified by the B<-m> option into B<QTL Cartographer> format with B<Rmap>.  
If the option is not used, then the script assumes that the map is already in the correct format and resides in
the file F<qtlcart.map>.

=item B<Rcross>

Convert the data found in the file specified by the B<-d> option into B<QTL Cartographer> format with B<Rcross>.  
If the option is not used, then the script assumes that the data is already in the correct format and resides in the
file F<qtlcart.cro>.

=item B<Qstats>

Run B<Qstats> to generate basic statistics.

=item B<LRmapqtl>

Do simple linear regression with B<LRmapqtl>.

=item B<SRmapqtl>

Do stepwise linear regression with B<SRmapqtl>.

=item B<Zmapqtl>

Do composite interval mapping with B<Zmapqtl> and the model specified on the command line.   The default is
to use model 6.  

=item B<JZmapqtl>

Do multiple trait composite interval mapping   with B<JZmapqtl>.  This will only be done if you speciefied more than
one trait with the B<-t> option.  


=item B<Eqtl>

Estimate QTL positions from the analyses above.   


=item B<Preplot>

Prepare the results for display with B<GNUPLOT>.   You will need to run B<GNUPLOT> on your own.

=item B<JZmapqtl>

This step converts the data into a format that can be analyzed by B<MultiRegress>.

=item B<MultiRegress>

Use least-squares analysis to estimate QTL postions.

=item B<Rqtl>

Convert the output of B<MultiRegress> into a format that can be read as an initial model for
multiple interval mapping.

=item B<MImapqtl>

Use the reformatted results from B<MultiRegress> to estimate new QTL and epistatic effects. 

=back

=head1 EXAMPLE

Suppose the files F<mletest.map> and F<mletest.cro> are in the current working
directory.   

    % Dobasics.pl  -b /home/basten/bin  -X mletest -S 13.0  -H 10  

Will assume that the B<QTL Cartographer> programs are in F</home/basten/bin>. 
It will use a significance threshold of 13.0.   It will do basic statistics, 
single marker analyses, interval and composite interval mapping, and prepare
the data for B<Gnuplot>.

Suppose the files F<realdatm.inp> and F<realdatc.inp> are in the current working
directory.  These files are formatted like F<map.inp> and F<cross.inp> and need to 
be translated.  Use

	% Dobasics.pl -b /home/basten/bin  -X realdat -m realdatm.inp -d realdatc.inp  -S 13.0  -H 30

for this analysis.

=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Qstats(1)>, B<Eqtl(1)>, B<LRmapqtl(1)>, B<SRmapqtl(1)>, B<Rmap(1)>, 
B<Rcross(1)>, B<Rqtl(1)>, B<MultiRegress(1)>, B<MImapqtl(1)>, B<Preplot(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
