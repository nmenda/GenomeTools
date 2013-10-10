#!/usr/bin/perl

=head1 NAME

 extract_snp_freq.pl
 This script extracts possible snps between a target genome (mpileup file with position and nucleotide) and  a reference file with SNPs called using VarScan L<http://varscan.sourceforge.net/>. It also prints the allele frequency from the varscan file of the target genome.

=cut

=head1 SYPNOSIS

 extract_snp_freq.pl [-h] -target <input_mpileup_file> -var <input_varscan target file > -ref <varscan file of the reference > -o <output_file> [-c 10 -f 90 -m 10 ]

=head1 DESCRIPTION

 This script reads 3 files: mpileup and varscan files of a target genome ( a genome containing some SNPs of interest from a known or unknown origin ) 
 and a varscan file of a reference genome ( for example a possible donor for the target SNPs).
 It is a good idea to use the options for filtering  the varscan file of the reference genome, and extract SNPs with a certain coverage (-c )  and frequency ( -f )  cutoff, and filtering by coverage from the target mpileup file ( -m ).

 The output is a  subset mpileup file based on matching chromosome positions,  along with the alternative allele from the reference genomes and the allele frequency from the target varscan file in the last 2 columns of the new file.

 Next step would be for example:
 * Plotting the SNPs 
 * Filtering the output further by subtracting the SNPs from another genome (e.g. if you are trying to eliminate shared SNPs across several genomes) 


=head2 I<Flags:>

=over

=item -target

B<input mpileup file>             input mpileup file of the target genome (mandatory)

=item -var

B<input varscan file>           input varscan file of the target genome (mandatory)

=item -ref 

B<input varscan file2>           input varscan file of the reference genome (mandatory)
    
=item -o

B<output_file>            output file (mandatory)

=item -c 

B<coverage>   coverage cutoff for selecting SNPs from the reference varscan file (optional) 

=item -f

B<frequency>  Frequency cutoff for selecting SNPs from the reference varscan file (optional)

=item -m

B<mpileup coverage> coverage cutoff for filtering target positions with low genome coverage (optional)

=item -h

B<help>                   print this help doc

=back

=head1 AUTHOR

Naama Menda<nm249@cornell.edu>

=cut



use strict;
use warnings;
use File::Slurp;

use Getopt::Long;
use Pod::Usage;

my ( $target_mpileup, $target_var, $ref_var, $out, $help); 
my $coverage = 1;
my $frequency = 1;
my $mcov = 1;

GetOptions (
    "target=s" => \$target_mpileup,    # string
    "var=s"    => \$target_var,         # string
    "ref=s"    => \$ref_var,
    "out=s"    => \$out,
    "c|cov=i"  => \$coverage,
    "f|freq=i" => \$frequency,
    "mcov|m=i" => \$mcov,
    "help"  => \$help)   # flag
    or  pod2usage(-verbose  => 2);

if ($help || !$target_mpileup || !$target_var || !$ref_var || !$out)  { pod2usage(-verbose  => 2); } 

print STDERR " c = $coverage , f = $frequency , m = $mcov \n";
exit;
# list file is varscan with chr. in column 1 and position in clumn 2

open (MPILEUP, "<$target_mpileup") || die "Cannot open target mpileup file $target_mpileup.\n";
open (TVAR, "<$target_var") || die "Cannot open target varscan file $target_var.\n";
open (REFVAR, "<$ref_var") || die "Cannot open reference varscan file $ref_var.\n";

my %ref_varscan;

while ( <REFVAR> ) {
    my $var_line = $_;
    chomp($var_line);
    my @var = split ("\t", $var_line);
    my $chr = $var[0];
    my $pos = $var[1];
    my $var_allele = $var[3];
    my $string = $var[4];
    my @string_values = split(":", $string);
    my $cov = $string_values[1];
    my $freq = $string_values[4];
    chop($freq);
    if ( ( $coverage == 1 ) || ( $cov >= $coverage ) ) {
	if (  ( $frequency == 1 ) || ( $freq >= $frequency ) ) {
	    $ref_varscan{$chr . ":" . $pos }->{allele} = $var_allele;
	    print STDERR " Ref varscan: $chr \t $pos \t coverage = $cov \t frequency = $freq ($frequency) \n";
	}
    }
}

# this is the varcan file of the target genome. We need to extract from here the frequency of each allele 

while ( <TVAR> ) {
    my $var_line = $_;
    chomp($var_line);
    my @var = split ("\t", $var_line);
    my $chr = $var[0];
    my $pos = $var[1];
    my $string = $var[4];
    my @string_values = split(":", $string);
    my $freq = $string_values[4];
    if ( defined $ref_varscan{ $chr . ":" . $pos } ) {
	$ref_varscan{$chr . ":" . $pos }->{freq} =  $freq;
    }
    print STDERR "Target  varscan: $chr \t $pos \t $freq \n";
}


## Select only if SNP exists in the ref_varscan file. 
##Not all of them will have frequency from the target varscan file, since only shared SNPs will be in both files

while ( <MPILEUP> ) {  # my mpileup file. Column 1 is chromosome name, 2 is the position of the SNP, 4 is the coverage
    my $pileup_line = $_;
    chomp $pileup_line;
    my @cells= split (/\t/,$pileup_line);
    my $chr = $cells[0];
    my $pos = $cells[1];
    my $nuc = $cells[2];
    my $cov = $cells[3]; #optional filtering by mpileup coverage
    if ( ( exists $ref_varscan{$chr . ":" . $pos} ) && ( $cov >= $mcov )  ) {
	my $var_freq = $ref_varscan{$chr . ":" . $pos}->{freq};
	my $var_allele = $ref_varscan{$chr . ":" . $pos}->{allele};
	print STDERR "mpileup: $chr \t $pos \t $var_allele \t $cov \t" ;
	if ($var_freq) { print STDERR "$var_freq"; }
	print STDERR "\n" ;
	no warnings 'uninitialized';
	write_file($out, { append => 1} , ( $pileup_line, "\t", $var_allele , "\t", $var_freq ,"\n")  );
    }
    else {
	next;
    }
}
