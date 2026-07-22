# 
#               CWTupdate
#
local ($opt_H,$opt_r,$opt_t,$opt_s, $opt_l, $opt_w, $opt_h) = (1,0,0,0,0,80,0);

my $indata;
my @fields;
my $i;

my $header;

my $printflag;
my $atline;


$usage = $usage . "$0: [-H Hypothesis] [-l lines] [-w columns] [-r ] [-t ] 
                [-s ] [-h]   < input > output ";

getopts "H:l:w:rtsh" ;

die $usage if ( $opt_h == 1 );


$theline = 1;
$atline = $indata = 0;
while ( <> ) {   
  if ($indata == 0) {
    print ;
    $theline +=1;
    $theline = 1 if ( $theline == $opt_l );
       
  }
  $indata = -1 if ( /-cross/ )  ;  
  
  if ( /^-e/ ) {# end of data
    $theline = 1;
     $atline = $indata =  0  ;  
    printf("\n-e") if $opt_l == 0 ;
    print "\n";
  }
  if ( $indata == 1 ) {
    $atline += 1; 
    chomp;
    @fields = split ;
    printf("\n%3s  %3s  %8s  ",$fields[0], $fields[1], $fields[2]);
    printf("%12s  %10s  ",$fields[3], $fields[4] ) if ( $opt_H == 1  ) ;    
    printf("%12s  %10s  ",$fields[10], $fields[6] ) if ( $opt_H == 10 ) ;    
    printf("%12s  %10s  ",$fields[11], $fields[8] ) if ( $opt_H == 20 ) ;    
    printf("%12s  %10s  %10s  ",$fields[3], $fields[7], $fields[9] ) if ( $opt_H == 30 ) ;    
    printf("%12s  %10s  %10s  %10s  ",$fields[4], $fields[7], $fields[9], $fields[6]  ) if ( $opt_H == 31 ) ;    
    printf("%12s  %10s  %10s  %10s  ",$fields[5], $fields[7], $fields[9], $fields[8] ) if ( $opt_H == 32 ) ;    
    
    if ( $opt_r == 1 ) {
		printf("%s  ",$fields[5] ) if ( $opt_H == 1  ) ;    
		printf("%s  ",$fields[12]  ) if ( $opt_H == 10 ) ;    
		printf("%s  ",$fields[13] ) if ( $opt_H == 20 ) ;    
		printf("%s  ",$fields[14] ) if ( $opt_H == 30 ) ;    
		printf("%s  %s  ",$fields[12], $fields[14]   ) if ( $opt_H == 31 ) ;    
		printf("%s  %s  ",$fields[13], $fields[14]  ) if ( $opt_H == 32 ) ;    
    }
    if ( $opt_t == 1 ) {
		printf("%s  ",$fields[6] ) if ( $opt_H == 1  ) ;    
		printf("%s  ",$fields[15]  ) if ( $opt_H == 10 ) ;    
		printf("%s  ",$fields[16] ) if ( $opt_H == 20 ) ;    
		printf("%s  ",$fields[17] ) if ( $opt_H == 30 ) ;    
		printf("%s  %s  ",$fields[15], $fields[17]   ) if ( $opt_H == 31 ) ;    
		printf("%s  %s  ",$fields[16], $fields[17]  ) if ( $opt_H == 32 ) ;    
    }    
    if ( $opt_s == 1 ) {
		printf("%s  ",$fields[7] ) if ( $opt_H == 1  ) ;    
		printf("%s  ",$fields[18]  ) if ( $opt_H == 10 ) ;    
		printf("%s  ",$fields[19] ) if ( $opt_H == 20 ) ;    
		printf("%s  ",$fields[20] ) if ( $opt_H == 30 ) ;    
		printf("%s  %s  ",$fields[18], $fields[20]   ) if ( $opt_H == 31 ) ;    
		printf("%s  %s  ",$fields[19], $fields[20]  ) if ( $opt_H == 32 ) ;    
    }  
  }
  
  if ( /^-s/ ) {
    for ( $i=$theline; $opt_l > 0 && $i<= $opt_l+3    ; $i++ ) {
      print "\n";
    }
    $indata =  1 ;
    &print_header;
    printf("\n-s") if $opt_l == 0 ;
  }

  if ( $opt_l > 0 && $atline == $opt_l )  {
    print "\n";
    &print_header;
    $atline = 0;
  }

}

sub print_header {
  my $j;
  
    printf("#");
    for ( $j=1; $j < $opt_w ; $j ++ ) {
      printf("-");
    }
    printf("\n#Chrom Mark Position ");
    printf("  LR(H1:H0)      Est(a) ") if ( $opt_H == 1  ) ;    
    printf("  LR(H1:H0)     Est(a1)  " ) if ( $opt_H == 10 ) ;    
    printf("  LR(H2:H0)     Est(d2)  "  ) if ( $opt_H == 20 ) ;    
    printf("  LR(H3:H0)     Est(a3)     Est(d3)  " ) if ( $opt_H == 30 ) ;    
    printf("  LR(H3:H1)     Est(a3)     Est(d3)     Est(a1)  "   ) if ( $opt_H == 31 ) ;    
    printf("  LR(H3:H2)     Est(a3)     Est(d3)     Est(d2)  " ) if ( $opt_H == 32 ) ;    

    if ( $opt_r == 1 ) {
		printf("r2(H0:H1)  "  ) if ( $opt_H == 1  ) ;    
		printf("r2(H0:H1)  "  ) if ( $opt_H == 10 ) ;    
		printf("r2(H0:H2)  "  ) if ( $opt_H == 20 ) ;    
		printf("r2(H0:H3)  "  ) if ( $opt_H == 30 ) ;    
		printf("r2(H0:H1)    r2(H0:H3)  "    ) if ( $opt_H == 31 ) ;    
		printf("r2(H0:H2)    r2(H0:H3)  "  ) if ( $opt_H == 32 ) ;    
    }
    if ( $opt_t == 1 ) {
		printf("tr2(H0:H1) "  ) if ( $opt_H == 1  ) ;    
		printf("tr2(H0:H1) "  ) if ( $opt_H == 10 ) ;    
		printf("tr2(H0:H2) "  ) if ( $opt_H == 20 ) ;    
		printf("tr2(H0:H3) "  ) if ( $opt_H == 30 ) ;    
		printf("tr2(H0:H1)   tr2(H0:H3) "    ) if ( $opt_H == 31 ) ;    
		printf("tr2(H0:H2)   tr2(H0:H3) "  ) if ( $opt_H == 32 ) ;    
    }    
    if ( $opt_s == 1 ) {
		printf("   S1"  ) if ( $opt_H == 1  ) ;    
		printf("   S1"  ) if ( $opt_H == 10 ) ;    
		printf("   S2"  ) if ( $opt_H == 20 ) ;    
		printf("   S3"  ) if ( $opt_H == 30 ) ;    
		printf("   S1          S3"    ) if ( $opt_H == 31 ) ;    
		printf("   S2          S3"  ) if ( $opt_H == 32 ) ;    
    }  
    printf("\n#");
    for ( $j=1; $j < $opt_w ; $j ++ ) {
      printf("-");
    }
}
#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Ztrim.pl - Trim a Zmapqtl output file for terminal reading

=head1 SYNOPSIS

   Ztrim.pl  [-H Hypothesis]  [-l lines]  [-w width] [-t]  [-r]  
           [-s]  [-h] < input > output


=head1 DESCRIPTION

B<Ztrim.pl> reads from the standard input and writes to the standard
output. It reads the output of B<Zmapqtl> and eliminates columns based on
the value of I<Hypothesis>.  It is a quick way to make the output file fit
in a terminal window without having the text wrap.  

=head1 OPTIONS

=over 4

=item  B<-H> 

This option requires an integer value of 1, 10, 20, 30,31, 32 or 33.   The behavior
of the program for each value is explained below.

=item  B<-l> 

This option requires an integer value greater than 1.  It specifies how often the column headers should
be printed.   The default value of  0  means that headers will  be printed at the top of
each block of output.   If your terminal shows 25 lines, then it is useful to reprint the header every
21 lines so as to have a header on each page of the output.   (The header is 3 lines long...if you use 
B<more>, then you will need one more line for it.)

=item  B<-w> 

This option requires an integer value greater than 1.      It specifies character width of
 the output terminal.   The default is 80, which is a common value.   It only effects the
horizontal bars defining the headers.   

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=item B<-r>

Prints   R squared values 

=item B<-t>

Prints total R squared values

=item B<-s>

Prints residual test statistic values

=back

=head1 OUTPUT

Different values used with the I<-H> option allow the user to output
specific columns from the I<qtlcart.z> file.   All options print out 
columns for the chromosome, marker and postions.   All print at least one
likelihood ratio test statistic and at least one parameter estimate.   

For crosses with two marker classes (backcrosses and recombinant inbreds), 
we can estimate one parameter value, namely the additive effect.   The test is
H0 (no QTL) versus H1 (QTL with an additive effect I<a>).   

For crosses with three marker classes, we can estimate additive and dominance effects.
Let I<a> be the additive effect and I<d> be the dominance effect.  
We can   set up four hypotheses:  

=over 4

=item H0

No QTL:    I<a> = 0 and I<d> = 0.

=item H1

QTL with an additive effect:  I<a> not 0 and I<d> = 0.   We refer to the estimate of
the additive effect as I<a1> under this hypothesis.

=item H2

QTL with a dominance effect: I<a> = 0 and I<d> not 0.  We refer to the estimate of
the dominance effect as I<d2> under this hypothesis.

=item H3

QTL with both an additive and a dominance effect: I<a> not 0 and I<d> not 0.  
We refer to the estimate of
the additive effect as I<a3> and that of the dominance effect as I<d3> under this hypothesis.

=back


Using the following values with the I<-H> option  yields 

=over 4

=item B<1>

This is the default and is useful for crosses with two genotypic classes.
The output simply deletes the columns with zeros.    You are left with the
likelihood ratio of H1 to H0 and the estimate of the additive effect.

=item B<10>

This essentially prints the same information as 1 above except it does it for
crosses with three genotypic classes.  


=item B<20>

This prints the likelihood ratio for a test of dominance where the additive effect is zero.
You get the H2 to H0 likelihood ratio, the estimate of dominance, the R squared estimates and
the test statistic for normality of the residuals. 


=item B<30>

This prints the H3 to H0 likelihood ratio.
You get the hypothesis test likelihood ratio, the estimates of additivity and 
dominance, the R squared estimates and
the test statistic for normality of the residuals. 

=item B<31>  

This option prints everything that option 30 prints, but adds the test of H3 to H1 and
the I<a1> estimate.

=item B<32>

This option prints everything that option 30 prints, but adds the test of H3 to H2 and
the I<d2> estimate.

=back


=head1 EXAMPLE

  ./Ztrim.pl -l 20 -w 78 -H 30 -s < qtlcart.z | more

will print the test statistic H3:H0 and the estimates of the additive and
dominance effects for the results in  I<qtlcart.z>.  It will also reprint the 
column headers every 20 lines.   This works well for a terminal set to 24 lines
by 80 columns.   The headers take 3 lines each and the B<more> command chews up the
last line of each screen.   

=head1 SEE ALSO 

B<Zmapqtl(1)> 

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
