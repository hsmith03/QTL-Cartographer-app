#
#               EWThreshold
#
local ($opt_s,$opt_r,$opt_m,$opt_h) = ("0.05",0,0,0);

my $indata;
my @fields;
my @lratio;
my @sortedlratio;
my $col;
my $i;
my $header;

$usage = $usage . "$0: [-s size] [-r] [-m] [-h]  < input \n";

getopts "s:rmh" ;


die $usage if ( $opt_h == 1 );

$indata = 0;
$i = 0;
while (<>) {
  $indata = 0 if ( /^-e/ ) ;  # end of data
#
  if ( $indata == 1 ) {       # if we are in the data, do something.
    chomp;
    @fields = split;
    $lratio[$i] = $fields[1];
    $i++;
  }
#
  $indata = 1  if ( /-start/ ) ; # line prior to beginning of data
}

@sortedlratio = sort by_number @lratio;
$sigind = int ($i * (1.0 - $opt_s) ) ;
#$i++;

if ( $opt_m == 0 ) {
	if ( $opt_r == 1 ) {
	  print "\n-siglevel $sortedlratio[$sigind]    \#  Experimentwise threshold for $i permutations at size $opt_s\n";
	}
	else {
	  print "\n Experimentwise threshold at size $opt_s for $i permutations is $sortedlratio[$sigind]\n";
	}
}
else {

  $sig01 = int ($i * 0.99);
  $sig05 = int ($i * 0.95);
  $sig1  = int ($i * 0.9 );
  print "\n-s01 $sortedlratio[$sig01]   -s05 $sortedlratio[$sig05]   -s1 $sortedlratio[$sig1]";

}
#
sub by_number {
        $a <=> $b;
}


#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

EWThreshold.pl - Calculate the Experimentwise Threshold from a permutation test

=head1 SYNOPSIS

  EWThreshold [-s size] [-h]   < input  > output

=head1 DESCRIPTION

B<EWThreshold.pl> reads from the standard input and writes to the standard output.
It reads the results of permutation test that are in the format
of   B<ZipermE.out>.  It sorts the likelihood maxima from that file   
and prints out the 100(1 - size)th percentile.   

=head1 OPTIONS

=over 4

=item  B<-s> 

This option requires an a real value.  It should be the size of the test,
that is the Type I error probability.   The default is 0.05.   

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 EXAMPLE

After running the example outlined in the B<GetMaxLR> manpage, you can run 

	% EWThreshold.pl -s 0.01 < corn.z3.e

to get the experimentwise threshold for a 1% test.   The threshold will be printed to the
standard output.  


=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>, B<GetMaxLR.pl(1)>, B<Permute.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
