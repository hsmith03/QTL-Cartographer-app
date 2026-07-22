#
#               CSSupdate
#
local ($opt_C,$opt_f,$opt_h) = (4,"",0);

my $indata;
my @fields;
my @lratio;
my @sortedlratio;
my %sum;
my %ssum;
my $col;
my $i;
my $lr;
my $header;
my $window = 10.0;
my $nbp = 5;
my $Model = 3;
my $trait = 1;
my $cross = "BC1";
my $permutation = 0;
my %counts; my %lrsample;

$usage = $usage .  "$0: [-C column] [-f ZipermC.out ] [-h]   < input > output  \n";

getopts "C:f:h" ;
die $usage if ( $opt_h == 1 );

$ZMAPQTL = $opt_f;  #  This is the old results file
#
#   If a CSSupdate.out file is given, then process it.
#
if ( length($opt_f) > 0 ) {   
    $indata = 0;
	open ZMAPQTL  or die "Can't open CSSupdate.out results file $ZMAPQTL $!\n";
    while ( <ZMAPQTL> ) {
      chomp;
      ($first, $second, @rest) = split; 
       if ( $indata == 0 ) {
        if ( /-start/ ) {
          $indata = 1 ;
          $permutation = $rest[$#rest - 1];
          $permutation += 1;
        }
      } 
      else {
        $row = $first;
        $chromosome = $second;
        $marker = $rest[0];
        $position = $rest[1];
        $key =  $first;
        $sum{$key} = $second;
           
        $ssum{$key} = $rest[0];
        $counts{$key} = $rest[1];
        $indata = -1 if (/-end/ );
      }    
    }
	close(ZMAQTL);
#    $permutation++;	
}



$row = 0;
$indata = 0;
while (<>) {
      $indata = -1 if ( /^-e/ ) ;  # end of data
      chomp;
      ($first, $second, @rest) = split;
      if ( $indata == 0 ) {  #  we are in the header, so get relevant information
        $window = $second if ( $first eq "-window" );
        $nbp = $second if ( $first eq "-background" );
        $Model = $second if ( $first eq "-Model" );
        if ( $first eq "-trait" ) {
          $trait = $second;
          $traitname = $rest[$#rest];
        }
        $cross = $second if ( $first eq "-cross" );  
      } 
      elsif ($indata == 1)  {  # we are in the data, so update the line
        $row++;
        $chromosome = $first;
        $marker = $second;
        $position = $rest[0];
        $lr = $rest[$opt_C-3];
          $key =  $chromosome . ":" . $marker . ":" . $position ;   
          $counts{$key}++ if ( $lr > 0.0 ) ;
          $sum{$key} = $sum{$key} + $lr ;
          $ssum{$key} = $ssum{$key} + $lr * $lr ;
          print "\n  $key       $sum{$key}      $ssum{$key}       $counts{$key}";

      }
    
 
#
      if ( /^-s/ ) { # line prior to beginning of data
        $indata = 1; 
#
#   The next line starts the data.   At this point, print the header
#   
$header = "\#      123456789   -filetype CSSupdate.out
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
\#Position is from the left telomere on the chromosome
\#     Key            Sum            SumofSquares       Count      -perm $permutation   -start";

        print $header;
      }	 
}
#
#  Put an end token on the file
#

print "  -end\n";


#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

SSupdate.pl - Update the column sums and sums of squares

=head1 SYNOPSIS

   CSSupdate.pl [-f CWTfile] [-h] < input > output


=head1 DESCRIPTION

B<CSSupdate.pl> reads from the standard input and writes to the standard
output. It reads the results of a B<Zmapqtl> analysis and updates the column
sum and sum of squares file.   

=head1 OPTIONS

=over 4

=item  B<-f> 

This option requires an input filename that must exist.  It allows the user to 
specify the F<CSSupdate.out> file for processing.   If it is not given, then the 
script assumes that an initial file will be created.  

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=item B<-C>

requires an integer operand.  This should be the column from the F<Zmapqtl.out>
file that you want processed.

=back

=head1 EXAMPLE

See the example in the B<GetMaxLR.pl> manpage.  An alternate version
of the shell script is given below, but the example in B<GetMaxLR>
is cleaner and simpler to understand.   The one presented below
has the advantage of being able to input command line parameters
rather than having to edit the script.   This script has also been
rewritten as a B<Perl> script with its own man page (B<Permute.pl(1)>). 

B<CWTupdate.pl> was meant to be run in a shell script.  Here is an 
example of a B<c> shell program that allows calculate the experimentwise
threshold.


	#!/bin/csh
	#   Permute.csh
	#   Usage:  Permute.csh stem permutations email
	#  where stem is the filename stem.
	#        permutations is the number of permutations
	#  and   email is the user's email address
	#  Note:  This only works if you have set and used a filename stem.
	#
	if ( $1 == '-h' ) then
	echo "    Usage:  Permute.csh stem model permutations email"
	echo "Where"
	echo "          stem  = filename stem"
	echo "         model  = Zmapqtl Model"
	echo "  permutations  = number of permutations"
	echo "         email  = user's email address"
	echo " "
	echo "Now exiting"
	exit
	endif
	set templog=temp.log
	/usr/bin/rm -f $templog
	echo "Permutation test started " > $templog
	/usr/bin/date >>  $templog
	echo "Stem: " $1 >> $templog
	echo "Model: " $2 >> $templog
	echo "Reps: " $3 >> $templog
	echo "Email: " $4 >> $templog
	set bindir=/usr/local/bin
	set i=1
	/usr/bin/mv $1.log $1.logsave
	/usr/bin/mv $1.z $1.zsave
	/usr/bin/rm -f $1.z$2e
	$bindir/CWTupdate -C 4 < $1.z > $1.z$2c
	while ( $i < $3 )
	$bindir/Prune -A -V -i $1.cro -b 2  >>&  $templog
	nice $bindir/Zmapqtl -A -V -M $2 -i $1.crb  >>&  $templog
	$bindir/CWTupdate -f $1.z$2c -C 4 < $1.z  >> $1.z$2cc
	/usr/bin/mv -f $1.z$2cc $1.z$2c
	/usr/bin/rm -f $1.z
	@ i++
	end
	/usr/bin/mv $1.zsave $1.z
	/usr/bin/mv $1.logsave $1.log
	/usr/bin/date >>  $templog
	/usr/ucb/mail $4 <  $templog

Suppose you had a data set F<corn.cro> and a map file F<corn.map>.  To use the above
shell script, create a directory called F<cornperm> and copy the two files into it.
Run B<Qstats> on the files to initialize the F<qtlcart.rc> file, and B<SRmapqtl> to
rank a set of markers for use with composite interval mapping.  Make sure that the
B<QTL Cartographer> programs are installed in the F</usr/local/bin> subdirectory
(or change the 29th line above).   Then, to do a permutation test using interval
mapping with 1,000 repetitions, run 

	% Permute.csh corn 3 1000 your.email.address  &

(substituting your real email address above).  The script will email you a message when
it is complete. 

Note that the above example uses B<-C 4> for the B<CWTupdate> line in the loop.
If you want to use a different column of likelihood ratios, you can change 
that option.  You could create multiple files of the format F<ZipermC.out>
and collect the maximal likelihood ratio from different hypothesis tests by
having multiple instances of B<CWTupdate> in the loop.    



=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>,  B<GetMaxLR.pl(1)>, B<Permute.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
