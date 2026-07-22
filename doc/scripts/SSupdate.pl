# 
#               SSupdate
#
local ($opt_l,$opt_a,$opt_d,$opt_I,$opt_f,$opt_h,$opt_c) = (4,0,0,0,"",0,0);

my $lrcol;
my $acol;
my $dcol;
my $indata;
my @fields;
my @lratio;
my $col;
my $i;
my $header;
my $window = 10.0;
my $nbp = 5;
my $Model = 3;
my $trait = 1;
my $cross = "BC1";
my $bootstrap = 0;
my %counts; my %lrsample;


$usage = $usage . "$0: [-l lrcol] [-a additive] [-d dominance] 
       [-I hypothesis] [-f Ziboot.out ] [-h] [-c]   < input > output  \n";

getopts "l:a:d:I:f:hc" ;
die $usage if ( $opt_h == 1 );

$lrcol = $opt_l;
$acol = $opt_a;
$dcol = $opt_d;
if ( $opt_I == 1 ) {
  $lrcol = 4;
  $acol = 7;
  $dcol = 0;

}
if ( $opt_I == 30 ) {
  $lrcol = 4;
  $acol = 8;
  $dcol = 10;

}
if ( $opt_I == 31 ) {
  $lrcol = 5;
  $acol = 8;
  $dcol = 10;
}
if ( $opt_I == 32 ) {
  $lrcol = 6;
  $acol = 8;
  $dcol = 10;
}
if ( $opt_I == 10 ) {
  $lrcol = 11;
  $acol = 7;
  $dcol = 0;
}
if ( $opt_I == 20 ) {
  $lrcol = 12;
  $acol = 0;
  $dcol = 9;
}

$ZMAPQTL = $opt_f;
#
#   If a Ziboot.out file is given, then process it.
#
if ( length($opt_f) > 0 ) {   
    $indata = 0;
	open ZMAPQTL  or die "Can't open Ziboot.out results file $ZMAPQTL $!\n";
    while ( <ZMAPQTL> ) {
      chomp;
      ($first, $second, @rest) = split; 
      $indata = -1 if (/-end/ );
       if ( $indata == 0 ) {
        $window = $second if ( $first eq "-window" );
        $nbp = $second if ( $first eq "-background" );
        $Model = $second if ( $first eq "-Model" );
        if ( $first eq "-trait" ) {
          $trait = $second;
          $traitname = $rest[$#rest];
        }
        $cross = $second if ( $first eq "-cross" );  
        if ( /-start/ ) {
          $indata = 1 ;
          $bootstrap = $rest[$#rest - 1];
          &print_header if ( $opt_c == 1 );
        }
      } 
      else {
        $chromosome = $first;
        $marker = $second;
        $position = $rest[0];
        $col = 1;
        $key =  $chromosome . ":" . $marker . ":" . $position ;    
        $lr{$key} = $rest[$col];  $col++;
        $lr2{$key} = $rest[$col];  $col++;
        $reps{$key} = $rest[$#rest];   
        if ( $opt_c == 1 ) {
          $lrmean = $lrvar = $lrstderr = -666;
          $lrmean = $lr{$key} / $reps{$key} if  $reps{$key} > 0;
          $lrvar = ($lr2{$key} - $lr{$key} * $lr{$key} / $reps{$key}) / ($reps{$key} - 1) if  $reps{$key} > 1; 
          $lrstderr = sqrt($lrvar)  if  $reps{$key} > 1 ;
          printf("\n%3d %4d %8.4f %9.3f %9.3f",$chromosome,   $marker,    $position,  $lrmean,   $lrstderr );       
        }
        if ( $acol > 0 ) {
          $a{$key} = $rest[$col];  $col++;
          $a2{$key} = $rest[$col];  $col++;
          if ( $opt_c == 1 ) {
            $amean = $avar = $astderr = -666;
            $amean = $a{$key} / $reps{$key} if  $reps{$key} > 0;
            $avar = ($a2{$key} - $a{$key} * $a{$key} / $reps{$key}) / ($reps{$key} - 1) if  $reps{$key} > 1;  
            $astderr = sqrt($avar) if  $reps{$key} > 1;      
            printf(" %9.3f %9.3f", $amean,   $astderr );       
          }
        }
        if ( $dcol > 0 ) {
          $d{$key} = $rest[$col];  $col++;
          $d2{$key} = $rest[$col];  $col++;
          if ( $opt_c == 1 ) {
            $dmean = $dvar = $dstderr = -666;
            $dmean = $d{$key} / $reps{$key} if  $reps{$key} > 0;
            $dvar = ($d2{$key} - $d{$key} * $d{$key} / $reps{$key}) / ($reps{$key} - 1) if  $reps{$key} > 1;    
            $dstderr = sqrt($dvar) if  $reps{$key} > 1;    
            printf(" %9.3f %9.3f", $dmean,   $dstderr );       
          }
          printf("  %d ",$reps{$key}) if $opt_c == 1 ;
        }
        $indata = -1 if (/-end/ );
      }    
    }
	close(ZMAQTL);
    $bootstrap++;	
    print "\n-end \n" if ( $opt_c == 1);
    exit if ( $opt_c == 1) ;
}
#
#  Process the standard input
#
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
        $chromosome = $first;
        $marker = $second;
        $position = $rest[0];
        $col = 1;
        if ( $opt_f eq "" ) {
           $lr = $lr2 = $a = $a2 = $d = $d2 = 0.0;
           $reps = 0;
        }
        else {
          $key =  $chromosome . ":" . $marker . ":" . $position ;   
          if ( $rest[$lrcol-3] > 0.0 ) {
              $reps = $reps{$key} +1;
			  $lr = $rest[$lrcol-3] + $lr{$key}; 
			  $lr2 = $rest[$lrcol-3] * $rest[$lrcol-3] + $lr2{$key} ;  
			  if ( $acol > 0 ) {
				$a = $rest[$acol-3] + $a{$key};  
				$a2 = $rest[$acol-3] * $rest[$acol-3] + $a2{$key};  
			  }
			  if ( $dcol > 0 ) {
				$d = $rest[$dcol-3]+ $d{$key};  
				$d2 = $rest[$dcol-3]* $rest[$dcol-3] + $d2{$key};  
			  }
		   }
		   else {
		      $reps = $reps{$key} ;
			  $lr =   $lr{$key}; 
			  $lr2 =  $lr2{$key} ;  
			  if ( $acol > 0 ) {
				$a =   $a{$key};  
				$a2 =   $a2{$key};  
			  }
			  if ( $dcol > 0 ) {
				$d = $d{$key};  
				$d2 =   $d2{$key};  
			  }
		   
		   
		   }
        }
        print "\n$chromosome    $marker    $position   $lr   $lr2";
        print " $a  $a2 " if ( $acol > 0 ) ;
        print " $d  $d2 " if ( $dcol > 0 ) ;     
        print " $reps";   
      } 
#
      if ( /^-s/ ) { # line prior to beginning of data
        $indata = 1; 
        &print_header;
      }	 
}
#
#  Put an end token on the file
#

print "\n-end\n";


sub print_header {

$filetype = "Ziboot.out" ;
$filetype = "Ziboots.out" if ( $opt_c == 1 ) ;
$colheads = "LR        LR2          a           a2             d              d2";
$colheads = "Mean(LR)  SD(LR)    Mean(a)    SD(a)     Mean(d)   SD(d)  " if ( $opt_c == 1);

$header = "\#      123456789   -filetype $filetype
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
\#Chrom Mark Position  $colheads      -boots $bootstrap   -start";

        print $header;
}

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

SSupdate.pl - Update the the sum and sum of squares for the bootstrap

=head1 SYNOPSIS

   SSupdate.pl [-l lrcol] [-a additive] [-d dominance] [-I hypothesis] 
         [-f SSfile] [-h] [-c] < input > output


=head1 DESCRIPTION

B<SSupdate.pl> reads from the standard input and writes to the standard
output. It can do one of two things.   The first is to initialize a file
of sums and sums of squares  in the F<Ziboot.out> format from a F<Zmapqtl.out>
file. The second is to read the results of a run of B<Zmapqtl> and the
current F<Ziboot.out> file and update the results.   It is a B<Perl>
script meant to be run in a loop with B<Prune> and B<Zmapqtl>.

=head1 OPTIONS

=over 4

=item  B<-f> 

This option requires an input filename that must exist.  It allows the user to specify 
the F<Ziboot.out> file for processing.   If it is not given, then the script assumes that
an initial file will be created.  

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=item B<-l>

requires an integer operand.  This should be the column from the F<Zmapqtl.out>
file with the likelihood ratio that you want processed.  By default, it is
4.

=item B<-a>

requires an integer operand.  This should be the column from the F<Zmapqtl.out>
file with the additive effect that you want processed.  By default, it is 0, which means
that the sum and sum of squares for the additive effect will not be updated.

=item B<-d>

requires an integer operand.  This should be the column from the F<Zmapqtl.out>
file with the dominance effect that you want processed.   By default, it is 0, which means
that the sum and sum of squares for the dominance effect will not be updated.


=item B<-I>

requires an integer operand.  This should be an hypothesis test code.  Possible values
are 1, 30, 31, 32, 10 and 20.   Using this option sets the proper values for
the B<-l>, B<-a> and B<-d> options, hence they are ignored in the presence of
the B<-I> option.

=item B<-h>

requires no operand.  If used, it processes the file specified with the B<-f> option.
This processing calculates the means and variances of the likelihood ratio and effects
for each position.   The program exits before reading from the standard input.

=back

=head1 EXAMPLE

This B<Perl> program was meant to be run in a shell script.  Here is an 
example of a B<c> shell program that allows calculate the sum and sum of 
squares for the likelihood ratio and additive effect in a bootstrap
experiment.

	#!/bin/csh
	#   Bootstrap
	#   Copyright (C) 2000 Christopher J. Basten 
	# Usage governed by the terms of the 
	# GNU  General Public License,  version 2 or higher
	#  See the file COPYING in this directory
	#
	#   This file was meant as an example.  You  will need to edit it 
	#   to work on your particular system with your data files.
	#
	#  Start by setting the variables needed.  
	# 
	set stem=corn                           #  filename stem
	set hypo=1                              #  hypothesis for SSupdate
	set model=3                             #  analysis model
	set reps=1000                           #  number of bootstraps
	set email=basten\@statgen.ncsu.edu      #  email address for notice
	set templog=temp.log                    #  temporary log file
	set qbin=/usr/local/bin                 #  where are the QTL Cart binaries
	set bin=/usr/bin                        #  where are the system programs
	#
	$bin/rm -f $templog
	echo "Bootstrap experiment started " > $templog
	$bin/date >>  $templog
	$bin/echo "Stem: " $stem >> $templog
	$bin/echo "Model: " $model >> $templog
	$bin/echo "Reps: " $reps >> $templog
	$bin/echo "Email: " $email >> $templog
	$bin/mv $stem.log $stem.logsave
	$bin/rm -f $stem.z${model}a
	$qbin/SSupdate -I $hypo < $stem.z > $stem.z$model.boot
	$bin/mv $stem.z $stem.zsave
	set i=1
	while ( $i <= $reps )
	$qbin/Prune -A -V -i $stem.cro -b 1  >>&  $templog
	$bin/nice $qbin/Zmapqtl -A -V -M $model -i $stem.crb >>&  $templog
	$qbin/SSupdate -I $hypo -f $stem.z$model.boot < $stem.z > $stem.z$model.new
	$bin/mv $stem.z$model.new $stem.z$model.boot
	$bin/rm $stem.z
	@ i++
	end
	$qbin/SSupdate -I $hypo -c -f $stem.z$model.boot > $stem.z$model.booted
	$bin/mv $stem.logsave $stem.log
	$bin/mv $stem.zsave $stem.z
	$bin/echo "Bootstrap experiment ended " >> $templog
	$bin/date >>  $templog
	/usr/ucb/mail $email <  $templog

Suppose you had a data set F<corn.cro> and a map file F<corn.map>.  To use the above
shell script, create a directory called F<cornboot> and copy the two files into it.
Run B<Qstats> on the files to initialize the F<qtlcart.rc> file, and B<SRmapqtl> to
rank a set of markers for use with composite interval mapping.  Make sure that the
B<QTL Cartographer> programs are installed in the F</usr/local/bin> subdirectory
(or change the F<qbin> line above).  Run the bootstrap with the following command: 

	% Bootstrap  &

The script will email you a message when it is complete.   The example above
uses interval mapping and does 1,000 bootstraps.   The script above has been
rewritten in B<Perl>:  Please look at the B<Bootstrap.pl> man page for more information.

Note that the above example uses B<-I 1> for the B<SSupdate> line in the loop
(this is set with the B<set hypo=1> line in the script).
This indicates that the dataset are the result of a backcross or
recombinant inbred line.  
If you want to use a different column of likelihood ratios, you can change 
that option.  You could create multiple files of the format F<Ziboot.out>
and collect the appropriate sums and sums of squares from different hypothesis tests by
having multiple instances of B<SSupdate> in the loop.    

=head1 HYPOTHESIS TESTS

Using the B<-I> option is an easier way to set the columns from the F<Zmapqtl.out> 
file that you want to process.  The following values are valid:

=over 4

=item   B<1> 

should be used with backcrosses or recombinant inbreds, that is 
only those crosses with two distinguishable marker types.    It will read the
likelihood ratio from column 4 and the additive effect from column 7.  The dominance
effect will be ignored. The likelihood ratio is for H1:H0.

=item B<10> 

can be used when more than three marker genotypes are distinguished.  Likelihood
ratios come from column 11, additive effects from column 7 and dominance 
effects are ignored.  The likelihood ratio is for H1:H0.

=item B<20> 

can be used when more than three marker genotypes are distinguished.  Likelihood
ratios come from column 12, dominance effects from column 9 and additive 
effects are ignored.    The likelihood ratio is for H2:H0.

=item B<30> 

can be used when more than three marker genotypes are distinguished.  Likelihood
ratios come from column 4, aditive effects from column 8 and dominance effects from column 10.
The likelihood ratio is for H3:H0.

=item B<31> 

can be used when more than three marker genotypes are distinguished.  Likelihood
ratios come from column 5, aditive effects from column 8 and dominance effects from column 10.
The likelihood ratio is for H3:H1.

=item B<32> 

can be used when more than three marker genotypes are distinguished.  Likelihood
ratios come from column 6, aditive effects from column 8 and dominance effects from column 10.
The likelihood ratio is for H3:H2.

=back

Recall that when you have two marker classes, there are two hypotheses:

=over 4

=item B<H0>

No QTL, that is the additive effect is zero

=item B<H1>

The additive effect is nonzero.

=back


In contrast, when you have three marker classes, then you have four hypotheses:

=over 4

=item B<H0>

No QTL, that is the additive and dominance effects are zero

=item B<H1>

The additive effect is nonzero, but the dominance effect is zero.

=item B<H2>

The dominance effect is nonzero, but the additive effect is zero.

=item B<H3>

Both the additive and dominance effects are nonzero.

=back


=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>, B<Qstats(1)>, B<SRmapqtl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
