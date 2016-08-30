#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV_XS;

my $dirname = '';
my $output = '';

GetOptions ('directory=s'=>\$dirname, 'output_file=s'=>\$output);

my $usage=<<INFO;
    Gather readcount files from panel of normals and merge them by variant.

    Example Usage: 
        perl pon_readcount_merge.pl --directory=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/BRCAnorms/ --output_file=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/pon_readcount_combined_BRCA.tsv

INFO

## Set up to parse TSV files
my $tsv = Text::CSV_XS->new ({
    binary => 1,        ## Allows non-ASCII characters including new lines
    eol => $/,          ## End of line string used
    auto_diag => 1,     ## Reports errors
    sep_char => "\t"    ## Defines character separator as tab
});

## Open output file
open my $out, '>', $output;

## Define hash
my %hash;

## Open the first file
#my $filename = join("/",$dirname,"H_LS-A1-A0SB-10B-01D-A142-09_readcounts.tsv");
#open my $file1, '<', $filename or die "Couldn't open File 1 ($filename)\n";

open my $file1, '<', join("/",$dirname,"H_LS-A1-A0SB-10B-01D-A142-09_readcounts.tsv") or die "Couldn't open File 1.\n";
open my $file2, '<', join("/",$dirname,"H_LS-A1-A0SD-10A-01D-A110-09_readcounts.tsv") or die "Couldn't open File 2.\n";

while (my $line = <$file1>) {
    chomp($line);
    my @row = split("\t", $line);
    my $ref = join("-",@row[0..4]);
    $hash{$ref} = $line;
}

## Create a hash to skip duplicates and counter to count them
my %dup;
my $counter = 0;

while (my $line2 = <$file2>) {
    chomp($line2);
    my @row2 = split("\t", $line2);
    my $check = join("-",@row2[0..4]);
    if($dup{$check}) {
        $counter++;
    } else {
        $hash{$check} = join("\t",$hash{$check},@row2[5..7]);
        $dup{$check} = "exists";
    }
}

print "$_ $hash{$_}\n" for (keys %hash);
print "$counter duplicates\n";

close $file1;
close $file2;
close $out;
