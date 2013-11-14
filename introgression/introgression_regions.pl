#/usr/lib/perl 

use strict;
use File::Slurp;
use Getopt::Long;

my $tab_file = "my_data_file";
my $density =10;
my $bins = 5;
my $spaces = 4;
my $out_file = "introgressions.tab";
my $skip;

GetOptions (
    "infile=s"   => \$tab_file,   # string
    "density=i"  => \$density,    # numeric
    "bins=i"     => \$bins,       # flag
    "space=i"    => \$spaces,
    "out=s"      => \$out_file,
    "c_skip=s"     => \$skip # chromosome ID to skip
    );

print "Density = $density, bins = $bins, spaces = $spaces, out= $out_file \n\n";
my @lines = read_file($tab_file);
my $bin_count = 0;
my $introg_count = 1;
my $space_count;
my @first = split("\t" , $lines[0] ) ;
my $prev_chr  = $first[0]; 
my %data;

foreach my $line (@lines) {
    chomp $line;
    my @intersect  = split ("\t" , $line ); #  1 = pimpinellifolium, 2 = YP , 3 = BTI87, 4= pimpin.
    
    my $chr = $intersect[0];
    if ( $chr eq $skip ) { next; }
    my $loc = $intersect[1];
    my $snps = $intersect[2];
    my $coverage = $intersect[3];

    if (  $chr eq $prev_chr  ) {
	###print "snps = $snps ";
	if ( $snps < $density ) { $space_count ++ ; }
	##print " space count = $space_count ";
	if ( $snps >= $density  || $space_count <= $spaces  || ( $space_count > $spaces && $coverage < 50000 ) ) {
	    push( @{ $data{$introg_count} },  join ("\t" , @intersect ) );
	   ### print " BIN count = $bin_count \n";
	    if ($snps >= $density ) { 
		$space_count = 0 ; 
		$bin_count++;
	    }
	}
    }
    #no introgression when snps are less than $density, and number of bins is lower than $bins
    if ( ( $snps < $density ) && ( $bin_count < $bins )  && ( $space_count > $spaces ) ) {
	undef $data{ $introg_count } ;
	####print " NO INTROGRESSION. snps = $snps, bins = $bin_count , spaces = $spaces_count \n";
	$introg_count++;
	$bin_count = 0;
    } elsif (  ( $snps < $density ) && ( $bin_count >= $bins ) && ( ( $space_count > $spaces ) && ( $coverage > 50000 ) )  ) { #end of introgression when # of SNPs bellow threshold and number spaces has maxed, and coverage is sufficient 
	#####print "Moving to next introgression :  snps = $snps, bins = $bin_count , spaces = $spaces_count \n";
	print "End of introgression $introg_count. bins = $bin_count, spaces = $space_count\n";
	## chop the max. 4 bin tail of low SNP number 
	for (my $tail = 0; $tail < 4; $tail++) {
	    pop(@{ $data{$introg_count} });
	}
        ##
	$introg_count++;
	$bin_count = 0;
	$space_count = 0;
    }
    # change of chromosome
    if ( $chr ne $prev_chr ) { 
	print "Next chromosome : $chr\n\n";
	#check the last introgression of the previous chromosome
	if ( $bin_count < $bins ) { undef $data{ $introg_count } ;  }
	$prev_chr = $chr;
	$bin_count = 0;
	$space_count = 0;
	$introg_count++;
	push( @{ $data{$introg_count} }, join ("\t", @intersect )  );
	if ( $snps < $density ) { $space_count++; } 
	$bin_count++;
    }
}
my $i_num = 0;
foreach my $key (sort { $a <=> $b}  keys %data) {
    if (defined ($data{$key} ) ) { 
	$i_num++;
	print "Writing introgression $key ($i_num)...\n"; 
    }
    foreach my $bin ( @{ $data{$key} } ) {
	write_file( $out_file, {append => 1 }, ($i_num , "\t",  $bin, "\n") ) ;
    }
}

    



