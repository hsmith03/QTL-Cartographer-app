# 
#               SRcompare
#
local ( $opt_f,$opt_h) = ( "SRmapqtl.old",0);

my $indata;
my $lost;
my $same;
my $gained;
my $whichone  ;
my $chrom;
my $marker;
my %cofactors;
my @rest;


$usage = $usage . "$0:   [-f SRmapqtl.old ] [-h]   < SRmapqtl.new > output  \n";

getopts "f:h" ;

die $usage if ( $opt_h == 1 );

$lost = $gained = $same = 0;

$SRMAPQTL = $opt_f;
#
#   There must be an old file to process.  Die if not.
#
if ( length($opt_f) > 0 ) {   
    $indata = 0;
	open SRMAPQTL  or die "There is no old file $SRMAPQTL to compare with...\n";
    while ( <SRMAPQTL> ) {
      chomp;
      if ( /-start/  ) {
        $indata = 1;
      }
      elsif ( $indata == 1 ) {
        ($chrom, $marker, @rest ) = split;
        $whichone = "c" . $chrom . "m" . $marker ;

        $lost +=1  unless $cofactors{$whichone} == 1 ;
        $cofactors{$whichone} = 1 ;
        $indata = 0 if ( /-end/ );
      }
    }
 	close(SRMAPQTL);
}
else {
  die "There is no old file to compare with...\n";
}


$indata = 0;
while (<>) {

      chomp;
      if ( /-start/ && $indata == 0) {
        $indata = 1;
      }
      elsif ( $indata == 1 ) {
        ($chrom, $marker, @rest ) = split;
        $whichone = "c" . $chrom . "m" . $marker ;

        if  ( $cofactors{$whichone} == 1 ) {
          $cofactors{$whichone} = 2;
          $same +=1;
          $lost -=1;                
        }
        elsif ( !(  $cofactors{$whichone} == 2 ||  $cofactors{$whichone} == 3 )  ) {
          $cofactors{$whichone} = 3;   
          $gained +=1;     
        }
        $indata = 0 if ( /-end/ );
      }

}

print " $same are the same, $lost were lost and $gained were gained. ";


#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

SRcompare.pl - Compare a pair of SRmapqtl outut files to see how many cofactors were gained or lost

=head1 SYNOPSIS

  SRcompare.pl [-f SRmapqtl.old] [-h] < SRmapqtl.new > output


=head1 DESCRIPTION

B<SRcompare.pl> reads from the standard input and writes to the standard
output. It is meant to compare the set of cofactors in two B<SRmapqtl>
output files.   It will print out the number of cofactors that are
the same, lost or gained in the two files.  

=head1 OPTIONS

=over 4

=item  B<-f> 

This option requires an input filename that must exist.  It allows the user to specify 
the old F<SRmapqtl.out> file for processing.   If it is not given, then the script dies.

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 INPUT FILES

The input files should be of the same format as the output of B<SRmapqtl>.   
You should have only one set of such results in a file to be processed by
B<SRcompare.pl>.   These results are created by B<Eqtl> as well as B<SRmapqtl>.  

=head1 OUTPUT

At this time the output will be a single line telling how many cofactors 
were lost, gained or remained, between the two files.  For example:

  10 are the same, 2 were lost and 0 were gained.

Ideally, you would like to see 0 lost and 0 gained.  

=head1 EXAMPLE

Suppose we have a map in F<qtlcart.map> and a data file in F<qtlcart.cro>, and
that there is no resource file in the current directory.   For this example, 
assume that there are only two marker genotypes (that is a backcross or recombinant
inbred line) and one trait.   The following series of commands can be used to 
compare the cofactors chosen from interval mapping with those from composite interval
mapping:

    % Qstats -X qtlcart
	% Zmapqtl -M 3 -A
	% Eqtl -S 12.0 -H 10 -I Z -A
	% mv qtlcart.z qtlcart.z.1
	% Zmapqtl -M 8 -A
	% mv qtlcart.eqt qtlcart.eqt.1
	% mv qtlcart.sr qtlcart.sr.1
	% Eqtl -S 12.0 -H 10 -I Z -A 
	% SRcompare.pl -f qtlcart.sr.1 < qtlcart.sr

In the above,   B<Qstats> is run to set the filename stem.  Next, B<Zmapqtl> does an
interval mapping analysis and B<Eqtl> picks a set of markers closest to the peaks
from  interval mapping:  These peaks are only used if they have likelihood ratios greater
than 12.0.  Next, the interval mapping results are moved to an new file and 
composite interval mapping is run using the cofactors identified via interval mapping.
The previous output files of B<Eqtl> (F<qtlcart.eqt> and F<qtlcart.sr>) are renamed and 
new files are generated from the composite interval mapping results.   Then, B<SRcompare.pl> 
is used to compare the two sets of results, those from interval mapping (F<qtlcart.sr.1>)
and those from composite interval mapping (F<qtlcart.sr>).   

This type of analysis can be iterated until a stable set of cofactors are identified.


=head1 SEE ALSO 

B<Eqtl(1)>, B<Qstats(1)>, B<SRmapqtl(1)>, B<Zmapqtl(1)> 

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
