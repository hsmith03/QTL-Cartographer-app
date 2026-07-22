#
#  Bootstrap.pl
#
#   Copyright (C) 2002 Christopher J. Basten 
#  Usage governed by the terms of the GNU  General Public License,  version 2 or higher
#  See the file COPYING in this directory
#
local ($opt_e,$opt_b,$opt_X,$opt_m,$opt_r, $opt_H,$opt_a,$opt_s,$opt_h) = ("",$TheBinDir,"",0,1000,1,0,0,0);

my $stem;
my $bin;
my $model;
my $reps;
my $ihypo;
my $i;
my $j;
my $templog;
my $output;
my $oldfile;
my $newfile;
my $resource;
my %parameters;


$usage = $usage . "$0: [-e email ] [-b bin ] [-X stem ] [-m model ] [-r iterations ] 
  [-H hypothesis ] [-a atrep] [-s] [-h] 
	The default values are:
		email      = $EmailAddress
		bin dir    = $TheBinDir
		stem       = $TheStem
		model      = $opt_m
		iterations = $opt_r
		hypothesis = $opt_H
		atrep      = $opt_a  \n";

getopts "e:b:X:m:r:H:a:sh" ;

die $usage if ( $opt_h == 1 );

$model = 0;
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
}
#  
#   If some command line options were used, then they override 
#   the values in qtlcart.rc.   
#
$bin    = $opt_b;         #  where are the  programs and scripts
#
$stem   = $opt_X if ( length($opt_X) > 0 );            #  filename stem
$stem   = $TheStem if ( length($stem) == 0 );          #  default is set in scripts.cfg (qtlcart)
#
$model  = $opt_m if ( $opt_m > 0 );         #  analysis model
$model  = 3 if ( $model == 0 );             # set default to 3
#
$reps   = $opt_r;         #  number of bootstraps
$EmailAddress  = $opt_e if length($opt_e) > 0 ;         #  email address for notice
$column  = $opt_c;        #  LR column to process
$ihypo = $opt_H;          #  hypothesis for JZmapqtl
#
if ( $opt_a == 0 ) {
  $templog = ">" . $stem . ".templog" if ( $opt_a == 0 );  #  Temporary log file
  unlink $templog if ($opt_a == 0 );
}
else {
  $templog = ">>" . $stem . ".templog" ;  #  Temporary log file
}
#
open TEMPLOG, $templog or die "Can't open file $templog : $!\n";
#
print TEMPLOG  "Bootstrap experiment started $date_time\n" ;
print TEMPLOG  "restarting after $opt_a iterations" if ( $opt_a > 0 );
print TEMPLOG  "bin:     $bin\n" ;
print TEMPLOG  "Stem:    $stem\n" ;
print TEMPLOG  "Model:   $model\n" ;
print TEMPLOG  "Hypo:    $ihypo\n" ;
print TEMPLOG  "Reps:    $reps\n" ;
print TEMPLOG  "Email:   $EmailAddress\n" if length($EmailAddress) > 0 ;
rename  "$stem.log", "$stem.logsave" if ( $opt_a == 0 );
$oldfile = $stem . ".z" . $model . ".b$opt_a" ;
`$bin/SSupdate.pl -I $ihypo < $stem.z > $oldfile` if ($opt_a == 0);
for ( $i= $opt_a+1; $i <= $reps ; $i++ ) {
    $newfile = $stem . ".z" . $model . ".b$i" ;
	`$bin/Prune -A -V -i $stem.cro -b 1`;   
	rename "$stem.crb", "$stem.cro.$i";
	`$bin/Zmapqtl -A -V -M $model -i $stem.cro.$i -o $stem.z.$i`;
	`$bin/SSupdate.pl -I $ihypo -f $oldfile < $stem.z.$i > $newfile`;
    if ( $opt_s == 0 ) {
	  unlink "$stem.z.$i" ;
	  unlink "$oldfile";
	  unlink "$stem.cro.$i";
	}
	$oldfile = $newfile;
}
$newfile = $stem . ".z" . $model . ".boot" ;
unlink "$newfile";
`$bin/SSupdate.pl -I $ihypo -c -f $oldfile > $newfile`;

rename "$stem.logsave", "$stem.log" if -f "$stem.logsave" ;
$date_time = `date`;
print TEMPLOG "Bootstrap experiment ended $date_time\n" ;

`$MailProgram $EmailAddress <  $templog` if -x $MailProgram  ;

close TEMPLOG;
exit 0 ;

#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Bootstrap.pl - Do a bootstrap analysis with Zmapqtl

=head1 SYNOPSIS

  Bootstrap.pl [-b bin] [ -X stem] [-m model] [-r iterations] [-e email] 
       [-H hypothesis] [-a rep] [-s] [-h]   


=head1 DESCRIPTION

B<Bootstrap.pl> iterates using B<Prune> and B<Zmapqtl>   to determine sampling variances 
based on a bootstrap resampling.   

=head1 OPTIONS

If the F<qtlcart.rc> file if it exists, B<Bootstrap.pl> will first set its parameter values from
that file.   Any command line options will override the F<qtlcart.rc> values.   If a parameter
has not been set by either the F<qtlcart.rc> file or command line parameters, default values are
set.   

=over 4

=item  B<-b> 

This option requires the path to the B<QTL Cartographer> binaries and perl scripts.  

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

=item  B<-H> 

Use this option to specify which hypothesis test you want to use.  The usual values are 1  or 30.  
The default is 1:  Unless there are three genotypic classes, do not change this option..   

=item B<-s>

requires no operand.  This tells B<Bootstrap.pl> to save the bootstrapped datasets and their analytical results.
If used with a large number of bootstraps, a great deal of harddisk space will be used up.  It is mainly 
for debugging purposes.   

=item B<-a>

requires an integer operand indicating the last completed iteration.   Useful if your machine crashed during 
the bootstrap.     

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.



=back

=head1 EXAMPLE

Suppose the files F<mletest.map> and F<mletest.cro> are in the current working
directory.   

    % Zmapqtl -X mletest -M 3 -A -V 
    % Bootstrap.pl  -b /home/basten/bin  -r 500

These commands assume that the B<QTL Cartographer> programs are in F</home/basten/bin>. 
It will use interval mapping and do 500 bootstrap iterations.  Note that you need to
do an initial B<Zmapqtl> run before beginning the bootstrap, and that initial run
will set the filename stem and the model for analysis.   

If your machine crashed during the bootstrap analysis, then you can restart where you
left off (provided that the crash  occurred during the B<Zmapqtl> run).   Suppose the
above completed 356 iterations and crashed during iteration 357.  Then

    % Bootstrap.pl  -b /home/basten/bin  -r 500 -a 356

would pick up the bootstrap at iteration 357 and complete at 500.  


=head1 CAVEATS

The B<-s> option allows you to save the bootstrapped datasets and analytical results.
For iteration I<i>, model I<m>, filename stem I<qtlcart> and single trait analysis,
there will be files F<qtlcart.cro.i>, F<qtlcart.zm.i> and F<qtlcart.zm.bi>.

This option is mainly for debugging purposes.
Be aware that a large number of iterations will use a great deal of disk space.  
You could modify the B<Bootstrap.pl> script to compress these files to save disk
space.  

=head1 SEE ALSO 

B<Zmapqtl(1)>, B<Prune(1)>, B<SSupdate.pl(1)>

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut

