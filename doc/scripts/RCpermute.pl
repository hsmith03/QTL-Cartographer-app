#
#               CSSupdate
#
local ($opt_R,$opt_C,$opt_h) = ("","",0);

my $indata;
my @fields;
my %csum;
my %cssum;
my %ccounts;
my %rsum;
my %rssum;
my %rcounts;
my $i;
my $lr;
my $header;
my $window = 10.0;
my $nbp = 5;
my $Model = 3;
my $trait = 1;
my $cross = "BC1";
my $permutation = 0;

$date_time = `date`;

$usage = $usage . "$0: [-C col file] [-R row file ] [-h]   < input > output  \n";

getopts "C:R:h" ;
die $usage if ( $opt_h == 1 );
die $usage if ( length($opt_R) == 0 ); #  need a row file
die $usage if ( length($opt_C) == 0 ); #  need a column file


$ROWSSS = $opt_R;  #  row sum and sum of squares
$COLSSS = $opt_C;  #  col sum and sum of squares
#
#   If a CSSupdate.out file is given, then process it.
#
if ( length($opt_R) > 0 ) {   
    $indata = 0;
	open ROWSSS  or die "Can't open Row sum file $ROWSSS $!\n";
    while ( <ROWSSS> ) {
      chomp;
      ($first, $second, @rest) = split; 
       if ( $indata == 0 ) {
        if ( /-start/ ) {
          $indata = 1 ;
        }
      } 
      else {
        $row = $first;
        $rsum{$row} = $second;
        $rssum{$row} = $rest[0];
        $rcounts{$row} = $rest[1];
        $indata = -1 if (/-end/ );
      }    
    }
	close(ROWSSS);
}
$permutation = $row;


if ( length($opt_C) > 0 ) {   
    $indata = 0;
	open COLSSS  or die "Can't open Column results file $COLSSS $!\n";
    while ( <COLSSS> ) {
      chomp;
      ($first, $second, @rest) = split; 
       if ( $indata == 0 ) {
        if ( /-start/ ) {
          $indata = 1 ;
        }
      } 
      else {
        $row = $first;
        $csum{$row} = $second;
        $cssum{$row} = $rest[0];
        $ccounts{$row} = $rest[1];
        $indata = -1 if (/-end/ );
      }    
    }
	close(COLSSS);
}



#
#   The next line starts the data.   At this point, print the header
#   
$header = "\#      123456789   -filetype RCpermute.out
\#
\#       QTL Cartographer  
\#       This output file  was created by 
\#         '$0'
\#       It is $date_time\#
\#
\#
\#The position is from the left telomere on the chromosome
-window         $window     Window size for models 5 and 6
-background     $nbp     Background parameters in model 6
-Model          $Model     Model number
-trait          $trait    Analyzed trait  $traitname
-cross          $cross     Cross
-permutations   $permutation    Number of permutations in test";

print $header;



print "\n\#\n\#This section is for Row Statistics\n\#Perm     Mean       S.Var.    S.Size    -start  rows";
$sum = $ssum = $ss = 0;

foreach $u ( sort by_number keys %rsum ) {

  $sum = $sum + $rsum{$u};
  $ssum = $ssum + $rssum{$u};
  $ss = $ss + $rcounts{$u};
  $mean = $rsum{$u} / $rcounts{$u};
  
  $variance = ($rssum{$u} - $rsum{$u} * $rsum{$u} / $rcounts{$u} ) / ( $rcounts{$u} - 1 ) ;
  
  printf("\n%5d  %8.4f  %10.4f  %5d ", $u,   $mean,   $variance,   $rcounts{$u});


}

  $mean = $sum / $ss;  
  $variance = ($ssum - $sum * $sum / $ss ) / ( $ss - 1 ) ;
  print "    -end \n-------------------------------------------";
  printf("\n Total %8.4f  %10.4f  %5d ",     $mean,   $variance,   $ss);


$sum = $ssum = $ss = 0;
print "\n\#\n\#This section is for Column Statistics\n\#               Site     Mean       S.Var.    S.Size    -start columns";

foreach $u ( sort by_number keys %csum ) {

  $sum = $sum + $csum{$u};
  $ssum = $ssum + $cssum{$u};
  $ss = $ss + $ccounts{$u};
  $mean = $csum{$u} / $ccounts{$u};
  
  $variance = ($cssum{$u} - $csum{$u} * $csum{$u} / $ccounts{$u} ) / ( $ccounts{$u} - 1 ) ;
  
  printf("\n%20s  %8.4f  %10.4f  %5d ", $u,   $mean,   $variance,   $ccounts{$u});


}

  $mean = $sum / $ss;  
  $variance = ($ssum - $sum * $sum / $ss ) / ( $ss - 1 ) ;
  print "    -end \n----------------------------------------------------------";
  printf("\n Total                %8.4f  %10.4f  %5d ",     $mean,   $variance,   $ss);



  print "\n-quit\n";

#
sub by_number {
        $a <=> $b;
}



#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

RCpermute.pl - Compute statistics from a permutation test

=head1 SYNOPSIS

   RCpermute.pl [-R row file] [-R col file] [-h]   > output


=head1 DESCRIPTION

B<RCpermute.pl> reads the results from two special files produced
by B<Permute.pl>.   These files contain the sum and sum of squares 
for the rows or the columns of the permutation analysis.   We define
rows as the permutation, and columns as the sites.   

=head1 OPTIONS

=over 4

=item  B<-C> 

This option requires an input filename that must exist.  It allows the user to 
specify the column results file for processing.   If it is not given, then the
script will exit.

=item  B<-R> 

This option requires an input filename that must exist.  It allows the user to 
specify the row results file for processing.   If it is not given, then the
script will exit.

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.


=back

=head1 ANALYSIS

B<Zmapqtl> calculates likelihood ratios at numerous sites during each
permutation.   Think of a matrix of these likelihood ratios where 
columns are for sites in the genome and rows are for permutations.
If you run B<Permute.pl> with the B<-u> flag, then the row and column
sums will be stored in files.   Sums of squares are also saved. 
The row results file will store a sum, sum of squares and sample size
over the genome for each of the permutations.  The column results file
will store the sum, sum of squares and sample size over permutations for
each test site.   For the stem F<qtlcart>, model I<3> and I<r> permutations,
these results will be in F<qtlcart.z3.rss> and F<qtlcart.z3.css.r>, respectively.


=head1 EXAMPLE

For a row results file F<qtlcart.z3.rss> and a column results file
F<qtlcart.z3.css.1000>, use

    RCpermute.pl -R qtlcart.z3.rss -C qtlcart.z3.css.1000 > qtlcart.rcout

to compute the row/column means and standard deviations.   


=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>,  B<GetMaxLR.pl(1)>, B<Permute.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
