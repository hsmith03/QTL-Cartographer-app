# 
#               GetMaxLR
#
local ($opt_C,$opt_h,$opt_i,$opt_j,$opt_r,$opt_t) = (4,0,0,0,1,0);

my $indata;
my @fields;
my @lratio;
my @sortedlratio;
my $col;
my $i;
my $header;
 
$usage = $usage . "$0: [-C column] [-h] [-i] [-j] [-r rep] ...  < input > output  \n";

getopts "C:hijr:t:" ;
$trait = $opt_t;

die $usage if ( $opt_h == 1 );

$indata = 0;
$i = 0;
while (<>) {
  $indata = -1 if ( /^-e/ ) ;  # end of data
#
  if ( $indata == 1 ) {       # if we are in the data, do something.
    chomp;
    @fields = split;
    $lratio[$i] = $fields[$opt_C -1];
    $i++;
  }
  if ( $indata == 0 ) {
    chomp;
    ($first, $second, @rest) = split;
    $window = $second if ( $first eq "-window" );
    $nbp = $second if ( $first eq "-background" );
    $Model = $second if ( $first eq "-Model" );
    if ( $first eq "-trait" ) {
      $trait = $second;
      $traitname = $rest[$#rest];
    }
    $cross = $second if ( $first eq "-cross" );    
  }
#
  if ( /^-s/ ) {# line prior to beginning of data
    $indata = 1;
    $traitname = "Joint Mapping" if ( $opt_j == 1);
$header = "\#      123456789   -filetype ZipermE.out
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
-trait          $trait     Analyzed trait $traitname
-cross          $cross     Cross
\#  Repetition    GlobalMax  -start";

	if ( $opt_i == 1 ) {
	  print $header;
	  exit;
	}
  }

}

@sortedlratio = sort by_number_descending @lratio;
print "\n $opt_r     $sortedlratio[0]";

#
sub by_number_descending {
        $b <=> $a;
}

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

GetMaxLR.pl - The the maximum likelihood ratio from a Zmapqtl output file

=head1 SYNOPSIS

   GetMaxLR.pl [-C column] [-i] [-h] [-j] [-r rep] [-t trait] < input > output


=head1 DESCRIPTION

B<GetMaxLR.pl> reads from the standard input and writes to the standard output.
It can do one of two things.   The first is to initialize a file of likelihood ratio
maxima so that the file would have the same output as a F<ZipermE.out> file.
The second is to read the results of a run of B<Zmapqtl> and pull out the
largest value in the user-specified column.  It is a B<Perl> script meant
to be run in a loop with B<Prune> and B<Zmapqtl> or B<JZmapqtl>.  

=head1 OPTIONS

=over 4

=item  B<-C> 

This option requires an integer value.  It allows the user to specify the F<Zmapqtl.out>
file column for processing.   

=item B<-i>

requires no operand.  If used, the script outputs a F<ZipermE.out> header to the
standard output and exits.  It does require an input file of the F<Zmapqtl.out> format.  
B<GetMaxLR> will get the model, cross, window size and number of background parameters
from the F<Zmapqtl.out> file.

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=item B<-j>

requires no operand.  It is a flag to indicate that the input file is
a joint-mapping file from B<JZmapqtl>.

=item B<-r>

requires the repetition number for the bootstrap or permutation test.

=item B<-t>

used with an integer is simply a way to pass on the trait value when using
B<JZmapqtl>.  It should only be used in conjuntion with the B<-j> and B<-i> 
options.   

=back

=head1 EXAMPLE

This B<Perl> program was meant to be run in a shell script.  Here is an 
example of a B<c> shell program that allows the user to calculate the experimentwise
threshold as well as comparisonwise values.


	#!/bin/csh
	#           Permute 
	#   Copyright (C) 2000 Christopher J. Basten 
	# Usage governed by the terms of the GNU  General Public License,  version 2 or higher
	#  See the file COPYING in this directory
	#
	#   This file was meant as an example.  You  will need to edit it 
	#   to work on your particular system with your data files.
	#
	#  Start by setting the variables needed.  
	# 
	set stem=corn                           #  filename stem
	set column=4                            #  LR column to process
	set model=3                             #  analysis model
	set reps=1000                           #  number of bootstraps
	set email=basten\@statgen.ncsu.edu      #  email address for notice
	set templog=temp.log                    #  temporary log file
	set qbin=/user/local/bin                #  where are the QTL Cart binaries
	set bin=/usr/bin                        #  where are the system programs
	#
	#   Should only need to change what is above.
	#
	$bin/rm -f $templog
	echo "Permutation test started " > $templog
	$bin/date >>  $templog
	$bin/echo "Stem: " $stem >> $templog
	$bin/echo "Model: " $model >> $templog
	$bin/echo "Reps: " $reps >> $templog
	$bin/echo "Email: " $email >> $templog
	$bin/mv $stem.log $stem.logsave
	$bin/mv $stem.z $stem.zsave
	$bin/rm -f $stem.z{$model}e
	set i=1
	$qbin/GetMaxLR -i < $stem.zsave > $stem.z{$model}.ewt
	$qbin/CWTupdate -C $column < $stem.zsave > $stem.z{$model}.cwt
	while ( $i <= $reps )
	$qbin/Prune -A -V -i $stem.cro -b 2  >>&  $templog
	$bin/nice $qbin/Zmapqtl -A -V -M $model -i $stem.crb  >>&  $templog
	$qbin/GetMaxLR -r $i -C $column < $stem.z  >> $stem.z{$model}.ewt
	$qbin/CWTupdate -f $stem.z{$model}.cwt -C $column < $stem.z  > $stem.z{$model}.newcwt
	$bin/mv $stem.z{$model}.newcwt $stem.z{$model}.cwt 
	$bin/rm -f $stem.z
	@ i++
	end
	$bin/echo "Now your can run EWThreshold on $stem.z$model.ewt" >> $templog
	$bin/mv $stem.zsave $stem.z
	$bin/mv $stem.logsave $stem.log
	$bin/date >>  $templog
	/usr/ucb/mail $email <  $templog

Suppose you had a data set F<corn.cro> and a map file F<corn.map>.  To use the above
shell script, create a directory called F<cornperm> and copy the two files into it.
Run B<Qstats> on the files to initialize the F<qtlcart.rc> file, and B<SRmapqtl> to
rank a set of markers for use with composite interval mapping.  Make sure that the
B<QTL Cartographer> programs and scripts are installed in the F</usr/local/bin> subdirectory
(or change the line setting F<qbin> above).   The script above is set to do 
1,000 permutations using interval mapping and restricting itself to column
four of the F<Zmapqtl.out> file.   Run it as follows:

	% Permute   &

The script will email you a message when it is complete. 

Note that the above example uses B<-C 4> for the B<GetMaxLR> line in the loop.
If you want to use a different column of likelihood ratios, you can change 
that option.  You could create multiple files of the format F<ZipermE.out>
and collect the maximal likelihood ratio from different hypothesis tests by
having multiple instances of B<GetMaxLR> in the loop.    

Once finished, you will have two files of interest.  The first will be F<corn.z3.e>
and will contain the maxima of the likelihoods for each of the permuted data sets.
Run B<EWThreshold> on this file to get your experimentwise thresholds.  The other
file, F<corn.z3.cwt>, will contain the comparisonwise thresholds.

The B<Permute> script above has been rewritten as a B<Perl> script and is
described in a man page B<Permute.pl(1)>.  

=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>, B<EWThreshold.pl(1)>, B<CWTupdate.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
