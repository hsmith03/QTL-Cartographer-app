# 
#
#  Ttransform.pl
#
local ($opt_a, $opt_b, $opt_c, $opt_f,$opt_h) = (0.0,1.0,0.0,1,0);


my $i;
my $j;
my $n;
my $t;
my $o;
my $p;
my $header;
my @data;
my @chromosomes;
my $c;
my $indata;
my @fields;
my $dpts;
my $atfield;
my @firstmoment;
my @secondmoment;
my @mean;
my @variance;
my @samplesize;
my @carriages;
my $i1lines;



$usage = $usage . "$0: [-a a] [-b b] [-c c] [-f function]  [-h] <input >output \n";

getopts "a:b:c:f:h" ;

die $usage if ( $opt_h == 1 );

#  First, read in the Rcross.out formatted file.
#  pick up parameter values when found.
$atfield = 1;
$i1lines = $indata = 0;
while (<>) {
  if ( /^-q/ ) {
    $indata = -1;
  }
  if ( $indata == 1 ) {
    chomp;
    @fields = split;
    $id = $fields[0] if ( $atfield == 1) ;
    $i1lines +=1 if ( $id == 1 );
    for ( $i=0; $i<=$#fields; $i++ ) {
      $individuals[$id] = $individuals[$id]  . $fields[$i]   . "|";
    }
    $atfield = $atfield + $#fields + 1;
    $carriages[$i1lines] = $#fields + 1 if $id == 1;   #  keep track of carriage returns in  marker data
    if ( $atfield > $dpts ) {                          #  this will allow same formatting in output
      $atfield = 1 ; 
      chop( $individuals[$id] );  # get rid of trailing |      
    }
  }

  if ( /^-s/ ) {
    chomp;
    $dpts = $p + $t + $o + 1;
    $header = $header . "#\n# The traits in this file were transformed by Ttransform.pl\n# on $date_time#\n";
    $header = $header . $_ ;
    $indata = 1;
  }
  if ( $indata == 0 ) {
    $header = $header . $_ ;
    if ( /^-/ ) {
      @fields = split;
      $n = $fields[1] if ( $fields[0] =~ /-n/ ) ;      
      $o = $fields[1] if ( $fields[0] =~ /-otraits/ ) ;      
      $t = $fields[1] if ( $fields[0] =~ /-traits/ ) ;      
      $p = $fields[1] if ( $fields[0] =~ /-p/ ) ;      
    
    }
  
  }

}


# 
#  This bit of code calculates the sample mean and sample variance for each trait
#
for ( $j=1; $j<=$n ; $j++ ) {
  @fields = split /\|/, $individuals[$j] ;
  for ( $i= $p+1; $i<=$p+$t; $i++ ) {
    if ( !( $fields[$i] eq "." ) ) { 
      $firstmoment[$i - $p] = $firstmoment[$i - $p] + $fields[$i];
      $secondmoment[$i - $p] = $secondmoment[$i - $p] + $fields[$i] * $fields[$i] ;
      $samplesize[$i - $p] +=1;
    }
  }
}
for ( $i=1; $i<=$t; $i++ ) {
  $mean[$i] = $firstmoment[$i] / $samplesize[$i];
  $variance[$i] =   ( $secondmoment[$i] - $samplesize[$i] * $mean[$i]**2.0 )/ ( $samplesize[$i] - 1.0) ;
}
#
#  Now it is time to print out the new cro file, but do the transforms on the fly
#
print "$header\n";
for ( $j=1; $j<=$n ; $j++ ) {
  @fields = split /\|/, $individuals[$j] ;
  $i1lines = 1;
  $counter = 0;
  for ( $i=0; $i<= $p ; $i++ ) {
    $counter +=1;
    print "$fields[$i] " if ( $i <= 1 ) ;
    printf("%3s ", $fields[$i] ) if ( $i > 1 );
    if ( $counter == $carriages[$i1lines] ) {
      $counter = 0;
      $i1lines +=1;
      print "\n    ";    
    }  
  }
  for ( $i= $p+1; $i<=$p+$t; $i++ ) {
    $trait = $fields[$i] ;         # is no transformation
    if ( !( $trait eq "." ) ) {
		if ( $opt_f == 1 ) {
		  $trait = ($trait - $mean[$i-$p]) / sqrt( $variance[$i-$p] ) ; #  (y-bar(y))/ stdev(y)
		}
		elsif ( $opt_f == 2 ) {
		  $trait = ($trait - $mean[$i-$p])**2.0 / $variance[$i-$p]   ; #  (y-bar(y))^2/ var(y)
		}
		elsif ( $opt_f == 3 ) {
		  $trait = $opt_a + $opt_b * $trait  + $opt_c * $trait**2.0;   # linear/quadratic transform    
		}
		elsif ( $opt_f == 4 ) {
		  $trait = $opt_a + $opt_b * log($trait) if $trait > 0.0 ;     # is a natural-log transform
		  $trait = "." if $trait <= 0.0 ;   #  can't do non-positives with log.
		}
		elsif ( $opt_f == 5 ) {
		  $trait = $opt_a * exp($opt_b * $trait);     # ae^{by} to the power of the trait
		}
		elsif ( $opt_f == 6 ) {
	#    LOOK HERE  and put in any thing you want. 
		  $trait = $opt_a / ( $trait ** $opt_b);   #  a / y^b  (just an example)
		}
    }
    print "    " if ( $i > $p+1 ) ;
    print "     $trait\n";
  
  }
  for ( $i= $p+1+$t; $i<=$p+$t+$o; $i++ ) {
    print "        $fields[$i]\n";
  
  }
}

print "-q\n";

exit 0 ;

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Ttransform.pl - Shell of a script to allow programming of trait transformations

=head1 SYNOPSIS

  Ttransform.pl [-a a] [-b b] [-c c] [-f func] [-h]   < input > output


=head1 DESCRIPTION

B<Ttransform.pl> reads an F<Rcross.out> formatted file from the standard input
and writes same to the standard output.   You can customize the script to 
transform trait values.  

=head1 OPTIONS

You will need to customize the script to suit your needs.   Look for 
the phrase LOOK HERE in the script.     

=over 4

=item B<-a>

requires a real-valued operand.  It is 0.0 by default and used as explained 
below.

=item B<-b>

requires a real-valued operand.  It is 1.0 by default and used as explained 
below.

=item B<-c>

requires a real-valued operand.  It is 0.0 by default and used as explained 
below.

=item B<-f>

requires an integer-valued operand.  It is 1 by default tells B<Ttransform.pl>
which function to use.

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.



=back

=head1 TRANSFORMS

There are some predefined transforms available.   Let Y be the trait value,
a, b and c  real-valued  coefficients and Z the transformed trait.   Use the B<-f> option and
an integer from one to six to have the following transforms applied to all traits:

=over 4

=item 1

Subtract the mean and divide by the standard error.  This is the default.

Z = (Y - Mean(Y)) / STDERR(Y)

=item 2

Subtract the mean, square it and divide by the sample variance.

Z = (Y - Mean(Y))^2 / VAR(Y)

=item 3

This is a linear or quadratic tranform.  

Z = a + b Y + c Y^2 

=item 4

This is a natural log tranform.  

Z = a + b log(Y)  

=item 5

This is an expontial transform

Z = a exp(b Y)

=item 6

This is for you to code any transform that you like.  

=back

=head1 EXAMPLE

Suppose the file   F<mletest.cro> is in the current working
directory.   

    % Ttransform.pl < mletest.cro > mletest.cro.out

would transform all the traits in F<mletest.cro>.  The default transform is
to subtract the mean and divide by the standard error.  

=head1 CAVEATS

You should know some Perl to use this script.  


=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut

