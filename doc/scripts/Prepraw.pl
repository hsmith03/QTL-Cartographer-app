#
#               Prepraw
#
use strict;
use vars qw(   $opt_h) ;
local (  $opt_h) = ( 0);

my $whichmarker; # which marker?    0,1,2,...m,m+1          0 means not a marker/trait
my $whichtrait;  # which trait?     0,1,2,...t,t+1          x+1 means finished with markers/traits
my %translator;  # translation table
my $nn; # sample size
my $mm; # number of markers
my $tt; # number of traits
my $counter; #  a counter
my $restoftheline;
my $previous;
my $case;
my @lfields; my @fields;
my $i;
my $temp;
my $tokencounter;


 
$usage = $usage . "$0:   [-h]  < input   \n";

getopts "s:h" ;
die $usage if ( $opt_h == 1 );

%translator = ("A", "A", 
               "B", "B", 
               "C", "C", 
               "D", "D",
               "H", "H",
               "-", "-"  );

$case = 0;               
$whichmarker = $whichtrait = 0;
$nn = $mm = $tt = 0;
$counter = 0;
$previous = "";

print "\#  -filetype mapmaker.raw\n\#  File processed by Prepraw on $date_time#\n";
while (<>) {
  $restoftheline = "";
  chomp;
  if ( /^\#/  || (length == 0 ) ) {
    print;
    print "\n";
  }
  else {
    @fields = split;
    if ( $fields[0] =~ /^\*/ ) {
      $tokencounter = 0;
      print "\n\#$previous had $counter tokens" if ( $counter > 0 ) ;
      print "  WARNING  "   if ( $counter > 0&&  $counter != $nn ) ;
      $counter = 0;
      $previous = $fields[0];
      printf("\n%-10s",$fields[0]);
      if ( $whichmarker <= $mm ) {
        $whichmarker++;
        $whichtrait++ if ( $whichmarker > $mm );
      }
      elsif ( $whichtrait <= $tt ) {
        $whichtrait++;
      }
      if ( $whichmarker <= $mm ) {
        for ( $i=1; $i<=$#fields; $i++ ) {
          $restoftheline = $restoftheline . $fields[$i] ;
        }        
        &finish_line;        
      }
      elsif ( $whichtrait <= $tt ) {
        for ( $i=1; $i<=$#fields; $i++ ) {
	          printf("%10s ", $fields[$i] );;
          $tokencounter++;
          if ( $tokencounter == 6 ) {
            print "\n          ";
            $tokencounter = 0;
          }
          $counter++;
        }            
      }    
    }
    else {
      $temp = $fields[0];
      $temp =~ tr/A-Z/a-z/;
      if ( $temp =~ /data/ ) {
        $previous = $temp;
        $_ =~ tr/A-Z/a-z/;
        print "$_\n";
      }
      elsif ( $previous =~ /data/ ) {
        $previous = "";
        $nn = shift(@fields);
        $mm = shift(@fields);
        $tt = shift(@fields);
        print "$nn   $mm   $tt\n";
        if ( defined($temp = shift(@fields)) ) {
			$temp =~ tr/A-Z/a-z/;
			if ( $temp =~ /case/ ) {
			  $case = 1;
			  $temp = shift(@fields);
			  $temp =~ tr/A-Z/a-z/;
			}
			if ( $temp =~ /symbols/ ) {  # take care of the symbols
			  
			  while ( defined($temp = shift(@fields)) ) {
				@lfields = split /=/, $temp ;  
				$translator{$lfields[0]} = $lfields[1];
			  
	#            $temp = shift(@fields);           
			  }               
			}
        }                        
      }
      else {
	     if ( $whichmarker <= $mm ) {
	        for ( $i=0; $i<=$#fields; $i++ ) {
	          $restoftheline = $restoftheline . $fields[$i] ;
	        }
	        
	        &finish_line;
	        
	        
	      }
	      elsif ( $whichtrait <= $tt ) {
	        for ( $i=0; $i<=$#fields; $i++ ) {
	          printf("%10s ", $fields[$i] );;
	          $tokencounter++;
              if ( $tokencounter == 6 ) {
                print "\n          ";
                $tokencounter = 0;
              }
	          $counter++;
	        }
	      
	      
	      }
      
      }
      
        
    }
    
    
    
  
  
  
  
  }

}

      print "\n\#$previous had $counter tokens" if ( $counter > 0 ) ;
      print "  WARNING  "   if ( $counter > 0&&  $counter != $nn ) ;
      print "\n";

sub finish_line {
  my $j; my $len; my $char;
  
  $len = length($restoftheline); 
  $counter = $counter + $len;
  for ( $j=0; $j<$len; $j++ ) {
    $char = substr($restoftheline, $j, 1)  ;
    print "$translator{$char} ";
    $tokencounter++;
    if ( $tokencounter == 30 ) {
      print "\n          ";
      $tokencounter = 0;
    }
  }
}

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Prepraw.pl - Input a mapmaker.raw file and prepare it for QTL Cartographer

=head1 SYNOPSIS

   Prepraw [-h]  < input  > output

=head1 DESCRIPTION

B<Prepraw.pl> reads from the standard input and writes to the standard output.
It reads a mapmaker.raw file and reformats it slightly.  It will use standard
tokens for each marker type.  It will count the number of individuals typed
for each marker, and print a warning if there are missing tokens.  


=head1 OPTIONS

=over 4


=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 EXAMPLE

Suppose you have a B<MAPMAKER> raw file called F<corn.raw>.
If you have trouble converting it with B<Rcross>, then 
try this.

	% Prepraw.pl < corn.raw > corn2.raw

and check if there are the proper number of data points for each
marker and trait.

=head1 SEE ALSO 

B<Rcross(1)> 

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
