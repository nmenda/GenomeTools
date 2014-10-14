#!/usr/bin/perl

=head1 NAME

 generate_r_snp_dist.pl
generate R commads for reading SNP and genome coverage data, and introgression regions. 
Output will be R graphs by chromosome , X axis is position , Y axis is SNP frequency, with introgression regions highlighted 
You can also mask introgressions using SNPs from another genome, with the -m option.
Second set of graphs can be generated for the genome coverage by chromosome 

=cut

=head1 SYPNOSIS

 generate_r_snp_dist.pl [-h] -s <base name for R table containing the SNP and coverage data>  -i <base name for R table containing the introgression regions> -o <output_file>
=head1 DESCRIPTION

 This script reads 2 R table names : SNP/Coverage data (read tab delimited file $V1=chromosome $V2=coordinate $V3=SNP count $V4=Coverage (sum of coverage for $V2 range}

SL2.50ch01  10000  8  181582
SL2.50ch01  20000  2  179344
.
.
.
SL2.50ch12 10000  35  177901


Variable 2 is the introgression data: $V1=introgression number  $V2=chromosome $V3=position $V4=SNP count $V5=coverage

1  SL2.50ch01  300000  48  170975
1  SL2.50ch01  310000  40  166692
1  SL2.50ch01  320000  22  17097
'   
'
'

Pass a 3rd table with introgressions for masking overlapping regions

1  SL2.50ch01  300000  33  110565                                                                                                                         
1  SL2.50ch01  310000  34  104331                                                                                                                          
2  SL2.50ch01  450000  20  259822  
2  SL2.50ch01  460000  74  321990
.
.

=head2 I<Flags:>

=over

=item -s | --snps

B<input R Table name>       input R table name for SNPs and genome coverage  (mandatory)

=item -i | --introg

B<input R table name>       input R table name for introgression regions in the genome (optional)

=item -m | --mask

B<input R table name>       input R table for introgression regions in a second genome for masking introgressions in genome 1.
      
=item -o

B<output_file>            output file (mandatory)


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
use Carp qw /croak/;

my ( $snps, $introg, $out, $mask, $help); 

GetOptions (
    "snps|s=s"   => \$snps,    # string
    "introg|i=s" => \$introg,
    "out|o=s"    => \$out,
    "mask|m=s"   => \$mask,
    "help"  => \$help)   # flag
    or  pod2usage(-verbose  => 2);

if ($help || !$snps || !$out)  { pod2usage(-verbose  => 2); } 




################################################################################

#first load the data into R. This has to be done before running the commands this script outputs 
# $snps <- read.table("bti87.snps.introg.10.tab", as.is=TRUE,header=FALSE, sep="\t")
# $introg <- read.table("8624h.introg.10.tab", as.is=TRUE,header=FALSE, sep="\t")

####################
#assign a variable name for each chromosome
my @chromosomes = qw | 01 02 03 04 05 06 07 08 09 10 11 12| ; 
my $chr_var;
my $int_var;
my $mask_var;

foreach my $chr ( @chromosomes) {
    $chr_var =  $snps . "_ch" . $chr ;
    write_file($out, { append =>1 } , $chr_var , " <- " , $snps , "[" , $snps , '$V1==\'SL2.50ch' , $chr , "',]" ,"\n\n" );
    
    $int_var = $introg . "_ch" . $chr ;
    $mask_var = $mask . "_ch" . $chr ;

    if ( $introg )  { 
	write_file($out, { append =>1 } , $int_var , " <- " , $introg , "[" , $introg , '$V2==\'SL2.50ch' , $chr , "',]" ,"\n\n" );
    }
    if ( $mask ) {
	write_file($out, {append =>1 } , $mask_var , " <- " , $mask, "[" , $mask , '$V2==\'SL2.50ch' , $chr , "',]" , "\n\n" );
    }
    #generate the coverage plot
    write_file($out , { append => 1 } ,  "pdf(file =", '"' , $chr_var , 'cov.pdf", width = 10, height = 4, )' , "\n\n" , 
	       "plot((" , $chr_var , '$V2)/1000000,' , $chr_var , '$V4/10000, type="h", ylim=rev(c(1,50)), col="black" , ylab="Coverage", xlab="coordinate (Mb)")' , "\n" , 
	       "dev.off()\n\n" );
    
    # generate the SNP plot with the introgression highlighted (if -i is used ) 
    write_file($out, {append =>1 } , 'pdf(file = "' , $chr_var, 'SNP_introg.pdf", width=10,height=6)' , "\n\n" ,
	       "plot((" , $chr_var , '$V2)/1000000,' , $chr_var , '$V3, type="h", ylim=c(1,100), col="black" , ylab="SNP Frequency", xlab="coordinate (Mb)",main = "Chromosome ' , $chr , '")' , "\n");
    if ( $introg) { 
	write_file($out, {append => 1 } , "lines((" , $int_var ,  '$V3)/1000000,' , $int_var , '$V4, type="h", col="red")' , "\n" );
    }
    if ( $mask ) {
	write_file( $out , {append => 1 } , "lines((" , $mask_var , '$V3)/1000000,' , $mask_var , '$V4, type="h", col="yellow")' , "\n" );
    }
    write_file($out, {append => 1 } , "dev.off()" , "\n\n" );

}
 




