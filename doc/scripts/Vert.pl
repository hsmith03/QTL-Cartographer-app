# 
#               Vert
#
local ($opt_i,$opt_o,$opt_h) = ("u","d",0);

$usage = $usage . "$0: [-i u||d||m ] [-o u||d||m]   [-h]   < input > output  \n";
getopts "i:o:h" ;
die $usage if ( $opt_h == 1 );

#  First, define the Input_Recod_Separator
$/ = "\r"   if ( $opt_i eq "m"  );    # Mac line-end is \r
$/ = "\n"   if ( $opt_i eq "u"  );    # Unix line-end is \n
$/ = "\r\n" if ( $opt_i eq "d"  );    # Dos line-end is \r\n

# Next, define the Output_Record_Separator and the Output directory
$\ = "\r"   if ( $opt_o eq "m"  );    # Mac line-end is \r
$\ = "\n"   if ( $opt_o eq "u"  );    # Unix line-end is \n
$\ = "\r\n" if ( $opt_o eq "d"  );    # Dos line-end is \r\n

while (<>) {
  chomp;
  print;
}



#
#  Pod documentation will follow.  Perl should ignore everything beyond this.
#

=head1 NAME

Vert.pl - Convert line endings

=head1 SYNOPSIS

   Vert.pl [-i u||d||m ] [-o u||d||m  ]  [-h] < input > output


=head1 DESCRIPTION

B<Vert.pl> reads from the standard input and writes to the standard
output. It converts the line endings of text files from one format
to another based on the optins B<-i> and B<-o>.   The default is to
convert UNIX files into DOS files.   

=head1 OPTIONS

=over 4

=item  B<-i> 

This option requires a one letter code indicating the type of 
input file.   Options are B<u> for UNIX, B<d> for DOS and 
B<m> for Macintosh.   

=item B<-o>

This option requires a one letter code indicating the type of 
output file.   Options are B<u> for UNIX, B<d> for DOS and 
B<m> for Macintosh.   

=item B<-h>

requires no operand.  If used, it prints a usage message and exits.

=back

=head1 EXAMPLE

Suppose F<corn.cro> and F<corn.map> are UNIX formatted files that you want
to convert to DOS.   The following sequence converts them to DOS format

	% Vert.pl -i u -o d < corn.cro > dcorn.cro
	% Vert.pl -i u -o d < corn.map > dcorn.map
	
No imagine that you have an entire directory of UNIX files and you want to 
convert them to Macintosh formatted files.   You could use 

	% foreach i ( * )
	? Vert.pl -i u -o m < $i > $i.mac
	? /usr/bin/mv $i.mac $i
	? end
	%

The above would replace all of the files with Macintosh formatted versions.  Alternatively,
you could copy them into a new directory:

	% mkdir macdir
	% foreach i ( * )
	? if ( -f $i ) Vert.pl -i u -o m < $i > macdir/$i 
	? end
	%

This script is useful for shared filesystems, that is when a fileserver can be
accessed by Windows, UNIX or Macintoshes.

=head1 AUTHORS

In general, it is best to contact us via email (basten@statgen.ncsu.edu).

	Christopher J. Basten, B. S. Weir and Z.-B. Zeng
	Department of Statistics, North Carolina State University
	Raleigh, NC 27695-7566, USA
	Phone: (919)515-1934

=cut
