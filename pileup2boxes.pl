#!/usr/bin/perl

=head1 NAME

 tab2boxes.pl
 Tool to count sequence ocurrences by box size 

=cut

=head1 SYPNOSIS

 tab2boxes.pl [-h] -i <inputtabfile> -s <id_sizes_file> [-b <box_size>]

=head2 I<Flags:>

=over


=item -i

B<inputtabfile>           tabular file (mandatory)

=item -s

B<id_sizes_file>          id_sizes file (two columns: ID, SIZE) (mandatory)

=item -b

B<box_size>               box size (1000000 by default)

=item -o 

B<out_file>             Output file for writing. Defaults to infile.bins

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script count ocurrences from a tab delimited file (generally sequence coordinates). These could be read from a .vcf file (the first 2 columns will be used)

 For example: pileup2boxes.pl -i myfile.tab -b 1000

   chr1       2           chr1   1000   2
   chr1     678           chr1   2000   1
   chr1    1230      =>   chr1   3000   0
   chr1    4536           chr1   4000   0
   chr1    4690           chr1   5000   2
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez <ab782@cornell.edu>
  Naama Menda <nm249@cornell.edu>

=cut

=head1 METHODS

 pileup2boxes.pl


=cut

use strict;
use warnings;
use Math::Round qw( nhimult );
use Getopt::Std;
use File::Slurp;

our ($opt_i, $opt_s, $opt_b, $opt_h, $opt_o, $opt_c);
getopts("i:s:b:o:c:h");
if (!$opt_i && !$opt_s && !$opt_b && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check them

my $in = $opt_i || 
    die("ARGUMENTS ERROR: -i <inputtabfile> option was not supplied.\n");

my $insz = $opt_s ||
    die("ARGUMENTS ERROR: -s <id_size_file> option was not supplied.\n");

my $box = $opt_b || 1000000;

my $out = $opt_o || $in . ".bins" ;

print STDERR "\n\n1) Parsing the ID sizes file.\n\n";

open my $sfh, '<', $insz or die("ERROR to open the file $insz: $!\n");

my %sizes = ();

while(<$sfh>) {
    chomp($_);
    my @data = split(/\t/, $_);
    $sizes{$data[0]} = $data[1];
    print STDERR "Size for $data[0] = $data[1] \n";
}


print STDERR "\n\n2) Parsing the $in file.\n\n";

open my $ifh, '<', $in or die("ERROR to open the file $in: $!\n");

my %id = ();
my $crr_box = $box;

my $l = 0;

if ($opt_c) { 
    print STDERR "opt_c : Counting sum of coverage values.\n\n";
}

while(<$ifh>) {
    chomp($_);
    $l++;

    ###print STDERR "\tParsing the line $l       \r";
    
    my @data = split(/\t/, $_);
    my $pos = $data[1];
    my $chr = $data[0];
    my $value = $data[2];
    #print STDERR "\n$chr \t $pos \t  box = $crr_box\n" ;
    if (exists $id{ $chr }) {
	if ( nhimult($box, $pos) == $crr_box) {
	    #print STDERR "INCREMENTING box $crr_box ($chr \t $pos ) \n";
	    if ($opt_c) { $id{ $chr }->{$crr_box} += ( $value -1 ) ; }
	    $id{ $chr }->{$crr_box} += 1;
	}
	else {
	    $crr_box = nhimult($box , $pos);
	    #print STDERR "NEW BOX: $crr_box ( $chr \t $pos) \n";
            if($opt_c) { $id{ $chr }->{$crr_box} =  $value  ; }
	    else { $id{ $chr }->{$crr_box} = 1; }
	}
    }
    else {
	$crr_box = nhimult($box , $pos);
	print STDERR  "new chromosome $chr. Box = $box , crr_box = $crr_box\n";
	print STDERR "NEW BOX: $crr_box ( $chr \t $pos) value = $value \n";

	if($opt_c) { $id{ $chr }->{$crr_box} = $value -1  ; }
	else { $id{ $chr }->{ $crr_box } =  1; }
    }
}
print STDERR "\n\n";

print STDERR "3) Producing output .\n\n";

## It will scan the sizes

foreach my $chr (sort keys %sizes) {
    print STDERR "Binning chromosome $chr \n";
    if (defined $id{$chr} ) {
	my $max = $sizes{$chr};
	my %feats = %{$id{$chr}};
	my $c = $box;

	while($c <= $max) {
	
	    if  ($c < $max ) {
		if (exists $feats{$c}) {
		    #print STDERR "Found box: $chr : $c\n";
		    write_file($out, { append => 1} , ( "$chr\t$c\t$feats{$c}\n")  );
		}
		else {
		    write_file($out, { append => 1} , ( "$chr\t$c\t0\n")  );
		}
	    }
	    if ($c == $max) {
		if (exists $feats{$c}) {
		    write_file($out, { append => 1} , ( "$chr\t$c\t$feats{$c}\n")  );
		}
		else {
		    write_file($out, { append => 1} , ( "$chr\t$c\t0\n")  );
		}
	    }
	    $c += $box;
	}
    }
}

print STDERR "\n\nDONE\n\n";


=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

      This script count ocurrences for a tab file (generally sequence coord.).

      For example: pileup2boxes.pl -i myfile.tab -b 1000

     chr1       2           chr1      0   2
     chr1     678           chr1   1000   1
     chr1    1230      =>   chr1   2000   0
     chr1    4536           chr1   3000   0
     chr1    4690           chr1   4000   2

    Usage:
       
      tab2boxes.pl [-h] -i <inputtabfile> -s <id_sizes_file> [-b <box_size>]

    Flags:

      -i <inputtabfile>           tabular file (mandatory)
      -s <id_sizes_file>          id_sizes file (two columns: ID, SIZE) (mand.)
      -b <box_size>               box size (1000000 by default)
      -o <out_file>               output file. Defaults to input_file.bins
      -h <help>                   print the help
     
EOF
exit (1);
}

