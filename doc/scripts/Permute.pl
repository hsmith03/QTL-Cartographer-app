# 
#  Permute.pl
#
local ($opt_I,$opt_t,$opt_e,$opt_b,$opt_X,$opt_m,$opt_r,$opt_c,$opt_a,$opt_s,$opt_u,$opt_d,$opt_h) = (10,-1,"",$TheBinDir,"",0,1000,4,0,0,0,0,0);

my $stem;
my $bin;
my $model;
my $ihypo;
my $reps;
my $column;
my $i;
my $j;
my $traits;
my $whichtrait;
my $crofile;
my $templog;
my $output;
my $ewtfile;
my $sssfile;
my $doustfile;
my $dcolumn;
my $doutput;
my $oldfile;
my $newfile;
my $resource;
my $cssfile;
my $newcss;
my $currentz;
my $currentcro;
my %parameters;


$usage = $usage . "$0: [-t trait] [-e email] [-b bin] [-X stem] [-m model]  
 [-r reps] [-c column] [-I ihypo] [-a atrep] [-s] [-u] [-d] [-h]  
	The default values are:
		email      = $EmailAddress
		bin dir    = $TheBinDir
		stem       = $TheStem
		model      = $opt_m
		trait      = $opt_t
		reps       = $opt_r
		column     = $opt_c
		hypothesis = $opt_I
		atrep      = $opt_a
	  and the -s, -u and -d flags are off.  \n";

getopts "I:t:e:b:X:m:r:c:a:sudh" ;
die $usage if ( $opt_h == 1 );

if ( $opt_d == 1 ) {
   $opt_d = 0 unless ( $opt_I == 14 || $opt_I == 34 );
}

$traits = 0;
$whichtrait = -1;

$resource = "qtlcart.rc";
if ( -f $resource ) {
  open RESOURCE, $resource or die "Can't find resource file $resource : $!\n";
  while ( <RESOURCE> ) {
  
    if ( /^-/ ) {
      chomp;
      @fields = split;
      $parameters{$fields[0]} = $fields[1] ;     
    }
  
  }
  close RESOURCE;
  $model = $parameters{"-Model"} if length( $parameters{"-Model"}) > 0 ;
  $stem = $parameters{"-stem"} if length( $parameters{"-stem"}) > 0 ;
  $whichtrait = $parameters{"-whichtrait"} if length( $parameters{"-whichtrait"}) > 0 ;
  $traits = $parameters{"-traits"} if length( $parameters{"-traits"}) > 0 ;
  $crofile = "<" . $parameters{"-ifile"} if length( $parameters{"-ifile"}) > 0 ;

}
#  
#   If some command line options were used, then they override 
#   the values in qtlcart.rc.   
#
$bin    = $opt_b;         #  where are the  programs and scripts
#
$stem   = $opt_X if ( length($opt_X) > 0 );            #  filename stem
$stem   = $TheStem if ( length($stem) == 0 );         #  default is qtlcart
#
$model  = $opt_m if ( $opt_m > 0 );         #  analysis model
$model  = 3 if ( $model == 0 );             # set default to 3
#
$reps   = $opt_r;         #  number of bootstraps
$EmailAddress  = $opt_e if length($opt_e) > 0 ;         #  email address for notice
$column  = $opt_c;        #  LR column to process
$whichtrait = $opt_t unless ( $opt_t < 0 );     #  trait to analyze
$whichtrait = 1 if ( $whichtrait < 0 );         #  default is 1
$ihypo = $opt_I;          #  hypothesis for JZmapqtl
#
#
#     Determine the number of traits from the crofile
#
if ( $traits == 0 ) {
	$crofile = "<" . $stem . ".cro" unless length($crofile) > 0 ;
	open CROFILE, $crofile  or die "Can't open file $crofile : $!\n";
	while ( <CROFILE>   ) {
	  if ( /-traits/ ) {
		@fields = split;
		$traits = $fields[1];
	  }
	}
	close CROFILE;
}
#
#     First, delete any old temporary log file and then create a new
#     one with information about this run.
#
$templog =   $stem . ".templog" ;  #  Temporary log file
if ( $opt_a == 0 ) {
  unlink $templog if ( $opt_a == 0 );
  $templog = ">" . $stem . ".templog" ;  #  Want a new temporary log file
}
else {
  $templog = ">>" . $stem . ".templog"  ;  #  Temporary log file exists
}
open TEMPLOG, $templog or die "Can't open file $templog : $!\n";
print TEMPLOG  "Permutations  started $date_time\n" ;
print TEMPLOG  "restarting after $opt_a permutations\n" if ( $opt_a > 0 );
print TEMPLOG  "bin:     $bin\n" ;
print TEMPLOG  "Stem:    $stem\n" ;
print TEMPLOG  "Model:   $model\n" ;
print TEMPLOG  "Column:  $column\n" ;
print TEMPLOG  "Reps:    $reps\n" ;
print TEMPLOG  "Traits:  $traits";
print TEMPLOG  "  and you selected trait:   $whichtrait\n";
if ( $whichtrait > 0 && $whichtrait <= $traits ) {
  print TEMPLOG  "This is a single trait analysis with Zmapqtl\n"  ;
}
else {
  print TEMPLOG  "This is a multitrait analysis with JZmapqtl\n"  ;
}
print TEMPLOG  "Email:   $email\n" if length($email) > 0 ;
close TEMPLOG;
$templog = ">>" . $stem . ".templog" ;  #  make sure that we are appending to the log file
rename "$stem.log", "$stem.logsave";           #  rename the log file to save it

if ( $whichtrait > 0 && $whichtrait <= $traits ) {
  if ( $opt_u == 1 ) {
    $cssfile = $stem . ".z" . $model . ".css.$opt_a" ;    # column sums file
	$sssfile = $stem . ".z" . $model . ".rss" ;           # row sums file
    unlink     $sssfile if ( $opt_a == 0);        #  get rid of old  files if starting from scratch
	`$bin/SumLR.pl -i < $stem.z > $sssfile`  if ! -f $sssfile ;          # init the exp. wise file if none
  }
	$oldfile = $stem . ".z" . $model . ".cwt.$opt_a" ;    # comparison wise test results for previous iteration
	$ewtfile = $stem . ".z" . $model . ".ewt" ;           # exp. wise results file
    unlink   $oldfile, $ewtfile if ( $opt_a == 0);        #  get rid of old  files if starting from scratch
	`$bin/GetMaxLR.pl -i < $stem.z > $ewtfile`  if ! -f $ewtfile ;          # init the exp. wise file if none
	`$bin/CWTupdate.pl -C $column < $stem.z > $oldfile` if  ! -f $oldfile ; # init the comp. wise file
	for ( $i= $opt_a+1 ; $i <= $reps ; $i++ ) {
	    $newfile = $stem . ".z" . $model . ".cwt." . $i ;
	    $newcss  = $stem . ".z" . $model . ".css." . $i ;  #  new column sums file
		`$bin/Prune -A -V -i $stem.cro -b 2  `;            # permute the data
		$currentz = $stem . ".z." . $i ;
		$currentcro = $stem . ".cro." . $i; 
		rename "$stem.crb", $currentcro;
		`$bin/Zmapqtl -A -V  -t $whichtrait -M $model -i $currentcro -o $currentz `;         # analyze the data
		`$bin/GetMaxLR.pl -r $i -C $column < $currentz  >> $ewtfile`;      #  update the exp. wise file
		`$bin/SumLR.pl -r $i -C $column < $currentz  >> $sssfile` if ( $opt_u == 1 );      #  update the exp. wise file
		`$bin/CWTupdate.pl -f $oldfile -C $column < $currentz  > $newfile`;#  update the comp. wise file
		if ( $opt_u == 1 ) {
		  if ( -f $cssfile ) {
		    `$bin/CSSupdate.pl -f $cssfile -C $column < $currentz  > $newcss` ; #  update the column CSSfile
		  } else {
		    `$bin/CSSupdate.pl -C $column < $currentz  > $newcss`  ; #  Create the column CSSfile
		  }
		}
		if ( $opt_s == 0 ) {
		  unlink $currentz ;        # get rid of this analysis file
		  unlink $oldfile ;          # get rid of old cwt results
		  unlink $cssfile if ( $opt_u == 1 );
          unlink $currentcro ;      # get rid of permuted data set
        }
	    $oldfile = $newfile;         # set up for next iteration
	    $cssfile = $newcss if ( $opt_u == 1 );
	}
	$output = `$bin/EWThreshold.pl < $ewtfile`;      #  calculate the exp. wise test result for alpha = 0.05
    open TEMPLOG, $templog or die "Can't open file $templog : $!\n";
	print TEMPLOG $output ;
    close TEMPLOG;
}
else {
    $dcolumn = $column + 1;
	$ewtfile = $stem . ".z0"  . ".ewt" ;              # exp. wise results file
	$doustfile = $stem. ".z0" . ".gxe" ;
    for ( $j=0 ; $j <= $traits ; $j++ ) {             # Back up original analysis files
      rename "$stem.z$j",  "$stem.z$j.save" if -f "$stem.z$j" ;   
    }

	`$bin/GetMaxLR.pl -i -j -t $whichtrait < $stem.z0.save > $ewtfile`  if ! -f $ewtfile ;   # init the exp. wise file if none
	if ( $opt_u == 1 ) {
	  `$bin/SumLR.pl -i -j -t $whichtrait < $stem.z0.save > $sssfile`  if ! -f $sssfile ;   # init the sss file if none
	}
	if ( $opt_d == 1 ) {
	  `$bin/GetMaxLR.pl -i -j -t $whichtrait < $stem.z0.save > $doustfile`  if ! -f $doustfile ;   # GxE lr results
	}  
	for ( $i=$opt_a+1; $i <= $reps ; $i++ ) {
		`$bin/Prune -A -V -i $stem.cro -b 6  `;                            # permute the data
		rename "$stem.crb", "$stem.cro.$i";
		`$bin/JZmapqtl -A -V -I $ihypo -t $whichtrait -M $model -i $stem.cro.$i  `;   # analyze the data
		`$bin/GetMaxLR.pl -r $i -C $column < $stem.z0  >> $ewtfile`;       #  update the sss  file
		`$bin/GetMaxLR.pl -r $i -C $dcolumn < $stem.z0  >> $doustfile` if ( $opt_d == 1 );       #  update the doust file
		`$bin/SumLR.pl -r $i -C $column < $stem.z0  >> $sssfile` if ( $opt_u == 1 );       #  update the exp. wise file
		if ( $opt_s == 0 ) {
		  unlink "$stem.cro.$i";
          for ( $j=0 ; $j <= $traits ; $j++ ) {
            unlink "$stem.z$j"  if -f "$stem.z$j" ;    # get rid of analysis files of permuted data          
          }
         }
         else {
          for ( $j=0 ; $j <= $traits ; $j++ ) {
            rename "$stem.z$j", "$stem.z$j.$i" if -f "$stem.z$j" ;    # save analysis files        
          }
         
         
         }
	}
	$output = `$bin/EWThreshold.pl < $ewtfile`;       #  calculate the exp. wise test result for alpha = 0.05
	$doutput = `$bin/EWThreshold.pl < $doustfile` if ( $opt_d == 1 );       #  calculate the GxE  result for alpha = 0.05

    for ( $j=0 ; $j <= $traits ; $j++ ) {             # Copy back the original analysis files
      rename "$stem.z$j.save",  "$stem.z$j" if -f "$stem.z$j.save" ;   
    }
    open TEMPLOG, $templog or die "Can't open file $templog : $!\n";
	print TEMPLOG $output ;
	if ( $opt_d == 1 ) {
	  print TEMPLOG "\n  Here is the GxE output...\n";
	  print TEMPLOG $doutput ;
	}
    close TEMPLOG;
}

rename "$stem.logsave", "$stem.log";  # copy back the log file

$date_time = `date`;
open TEMPLOG, $templog or die "Can't open file $templog : $!\n";
print TEMPLOG "Permutations  ended $date_time\n" ;
close TEMPLOG;
$templog =   $stem . ".templog" ;  #  Temporary log file without the redirection

# `/usr/ucb/mail $email <  $templog` if length($email) > 0 ;   
# if this system has mail, you can mail a message of completion
`$MailProgram $EmailAddress <  $templog` if -x $MailProgram  ;

exit 0 ;


#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Permute.pl - Do a permutation analysis with Zmapqtl or JZmapqtl

=head1 SYNOPSIS

  Permute.pl  [-b bin] [-X stem] [-m model] [-r reps] [-e email] [-c column] 
      [-t trait] [-I ihypo] [-a atrep] [-s] [-u] [-h]   


=head1 DESCRIPTION

B<Permute.pl> iterates using B<Prune> and B<Zmapqtl> or B<JZmapqqtl>   to determine
significance thresholds based on permutation testing.  

=head1 OPTIONS

If the F<qtlcart.rc> file if it exists, B<Permute.pl> will first set its parameter values from
that file.   Any command line options will override the F<qtlcart.rc> values.   If a parameter
has not been set by either the F<qtlcart.rc> file or command line parameters, default values are
set.   

=over 4

=item  B<-b> 

This option requires the path to the B<QTL Cartographer> binaries and perl scripts.  
The default is   I<~/bin> directory.

=item B<-r>

requires an integer to control how many bootstrap iterations you want to do.  The default
is 1000.

=item  B<-X> 

This option allows you to specify the filename stem.  The default is I<qtlcart>.

=item  B<-e> 

This option requires an email address.   The temporary log file will be sent to
this address to indicate that the bootstrap is complete.   If blank, then no email
message will be sent.  The default is not to use this option.

=item  B<-m> 

This option allows you to specify the B<Zmapqtl> model to use.   The default is interval mapping (3).

=item  B<-c> 

Use this option to specify which column in the B<Zmapqtl> output file should be processed.  The default is
4.

=item  B<-t> 

Use this option to specify which trait to analyze.   If there are I<t> traits, then a value greater than
zero and less than or equal to I<t> will cause B<Zmapqtl> to be used.   Otherwise, B<JZmapqtl> will be used.
The default is 1.   Note that using this option with a negative value will set this option to 1.

=item  B<-I>

Use this option to specify the hypothesis test for  B<JZmapqtl>.  The default is 10.  It is ignored if
B<Zmapqtl> is used.

=item B<-s>

requires no operand.  This tells B<Permute.pl> to save the permuted datasets and their analytical results.
If used with a large number of permutations, a great deal of harddisk space will be used up.  It is mainly 
for debugging purposes.   

=item B<-a>

requires an integer operand indicating the last completed iteration.   Useful if your machine crashed during 
the permutation test.     

=item B<-u>

requires no operand.  This is the Unger flag.  If used, then the script will also calculate
the sum and sum of squares for each permutation over the entire genome, and for each site over
all permutations.   These values will be saved to files that can be processed with
B<RCpermute.pl> to calculate site and permutation means and variances.  

=item B<-d>

requires no operand.  This is the Doust flag.  If used, then the script will also calculate
the sum and sum of squares for the GxE scores.  It was put in at the request of Andrew Doust.
It only works with B<JZmapqtl> and models 14 or 34.  

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 EXAMPLE

Suppose the files F<mletest.map> and F<mletest.cro> are in the current working
directory.   

    % Zmapqtl -X mletest -M 3 -A -V 
    % Permute.pl  -b /usr/local/bin -c 4 -r 500

Will assume that the B<QTL Cartographer> programs are in F</usr/local/bin>. 
It will use interval mapping and do 500 permutations.  Note that you need to
do an initial B<Zmapqtl> run before beginning the permutation test.  This initial 
run will have created a F<qtlcart.rc> file that contains the model and stem information.

If your computer went down during the permutation run, you can pick up where you
left off.   Suppose that B<Permute.pl> had finished 323 permutations in the above 
example.   You would observe a file F<mletest.z3.cwt.323> in the current working
directory.   You could then  run

    % Permute.pl  -b /usr/local/bin -c 4 -r 500 -a 323

to continue with the permutation test starting at the 324th iteration and
finishing with iteration 500.   

Suppose we have another data set with multiple traits:  The map is in 
F<multitest.map> and the data in F<multitest.cro>.   Further suppose that this
data set has four traits and it is a backcross.   As above, assume that the
binaries are in F</usr/local/bin>.

    % JZmapqtl -X multitest -M 3 -t 5 -I 10 
    % Permute.pl -b /usr/local/bin  -c 5 -r 500 -I 10 -t 5

Will use all the traits in a multitrait analysis and permute the data 500 times.
The likelihood ratio in the F<multitest.z0> file will be the focus of the test.
We need to specify column 5 for use with B<JZmapqtl>.   


=head1 CAVEATS

The B<-s> option allows you to save the permuted datasets and analytical results.
For iteration I<i>, model I<m>, filename stem I<qtlcart> and single trait analysis,
there will be files F<qtlcart.cro.i>, F<qtlcart.zm.i> and F<qtlcart.zm.cwt.i>.
For multitrait analysis, the F<qtlcart.cro.i> will be saved along with the
set of trait analysis files F<qtlcart.zt.i>, where I<t> is
the trait.

This option is mainly for debugging purposes.
Be aware that a large number of permutations will use a great deal of disk space.  
You could modify the B<Permute.pl> script to compress these files to save disk
space.  

=head1 SEE ALSO 

B<Zmapqtl(1)>, B<JZmapqtl(1)>, B<Prune(1)>, B<CWTupdate.pl(1)>, B<GetMaxLR.pl(1)>,
B<EWThreshold.pl(1)>, B<RCpermute.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut

