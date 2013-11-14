#!/usr/bin/perl

=head1 NAME

 introgres_stats.pl
  print introgression start, end, and number of gene models within that region

=cut

=head1 SYPNOSIS

 introgres_stats.pl [-h] -i <inputtabfile> -a <gff3 annotation file>

=head2 I<Flags:>

=over


=item -i

<input_tab_file>           tabular file (mandatory)

=item -a

B<gff3_annotation_file>      gff3 file with gene models for counting number of genes in each region


=item -o 

B<out_file>             Output file for writing. Defaults to infile.sum

=item -h

B<help>                   print the help

=back

=cut

=head1 AUTHORS

  Naama Menda <nm249@cornell.edu>

=cut



use strict;
use warnings;

use Getopt::Long;
use File::Slurp qw (write_file);
use Pod::Usage;
use List::MoreUtils qw( minmax );

my ( $in_file, $gff, $out, $help); 

GetOptions (
    "infile|i=s"  => \$in_file,
    "gff|g|a=s"   => \$gff,
    "out|o=s"     => \$out,
    "help"        => \$help)   # flag
    or  pod2usage(-verbose  => 2);

if ($help || !$in_file || !$gff  || !$out)  { pod2usage(-verbose  => 2); } 

# clear existing out file 
write_file( $out , '');

open (IN, "<$in_file") || die "Cannot open input file $in_file.\n";

my $err = $out . ".err" ;
my %starts;
my %ends;

#parse the gff3 file 
open ( GFF3 , "<$gff" ) || die "Cannot open gff3 file $gff \n";
my %gff_genes;

while ( <GFF3> ) { 
    my $line = $_;
    chomp($line) ;
    my @fields = split ("\t" , $line ) ;
    my $chr = $fields[0];
    my $type = $fields[2];
    no  warnings "uninitialized";
    if ( $type eq "gene" ) {
	push @{ $starts{ $chr } } , $fields[3];
	push @{ $ends{ $chr } } , $fields[4];
    }
}

my %stats;

while ( <IN> ) { 
    my $line = $_;
    chomp($line) ; 
    my @fields = split ("\t" , $line) ; 
    my $num = $fields[0];
    my $chr = $fields[1];
    my $pos = $fields[2];
    my $snps = $fields[3];
    $stats{$num}->{chr} = $chr;
    push ( @ { $stats{$num}->{pos} } , $pos ); # {pos} is a list of positions, later the start and end will be extracted
    $stats{$num}->{snps} += $snps;
}


foreach my $introg_num (sort { $a <=> $b } keys %stats ) {
    my @positions = @{$stats{$introg_num}->{pos} };
    my ($start, $end )  = minmax @positions;
    my $size = scalar( @positions ) ;
    my $real_start = $start - ( ($end - $start)/($size-1)) + 1 ; #subtract bin size from first position
    my $bin_size =  ( ( $end - $start ) / ( $size -1 ) ) * ( $size );

    my $chr = $stats{$introg_num}->{chr};
    my $snps = $stats{$introg_num}->{snps};
    my $gene_count;
    my @filtered_e =  ( grep { $_ < $end } @{ $starts{ $chr } } );
    my @filtered_s =  grep { $_ > $start } @filtered_e;
    $gene_count = scalar( @filtered_s) ;
    write_file( $out, { append => 1 } , join("\t" , ($introg_num, $chr, $bin_size, $real_start, $end, $snps, $gene_count, "\n") ) ) ;
}
