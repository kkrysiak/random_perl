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

#open my $file1, '<', join("/",$dirname,"H_LS-A1-A0SB-10B-01D-A142-09_readcounts.tsv") or die "Couldn't open File 1.\n";
#open my $file2, '<', join("/",$dirname,"H_LS-A1-A0SD-10A-01D-A110-09_readcounts.tsv") or die "Couldn't open File 2.\n";

my $firstfile = "";
#my $file1 = "";

opendir(my $dir, $dirname) or die "Couldn't open directory ($dirname).\n";

while (my $file = readdir($dir)) {
    
    ## check that it ends in .tsv
    next unless ($file =~ m/\.tsv$/);

    ## Get the variant list from the first file
    if($firstfile eq "") {
        ## Assign a value to first file to make this if statement false
        $firstfile = $file;
        ## Open the first file 
        open my $file1, '<', join("/",$dirname,$file) or die "Couldn't open file ($file).\n";
        while (my $line = <$file1>) {
            chomp($line);
            ## Split the first row into an array
            my @row = split("\t", $line);
            ## Create a hash reference out of the chromosome name, start, stop, ref, var joined by -
            my $ref = join("-",@row[0..4]);
            ## Assign the line to the value of the hash
            $hash{$ref} = $line;
        }

    ## Read in the rest of the readcount files
    } else {
        ## Create a hash to skip duplicates and counter to count them
        my %dup;
        my $counter = 0;
        ## Open the next file in the directory
        open my $file2, '<', join("/",$dirname,$file) or die "Couldn't open file ($file).\n";
        while (my $line2 = <$file2>) {
            chomp($line2);
            my @row2 = split("\t", $line2);
            ## Make the variant (chr, start, stop, ref, var) a key to check if it is in the hash
            my $check = join("-",@row2[0..4]);
            ## Check if that variant has been seen before (therefore a duplicate in this file
            if($dup{$check}) {
                ## Increment the duplicate counter
                $counter++;
            ## If it is not a duplicate add it to the hash
            } else {
                ## Appends the 3 readcount columns to the end of the growing hash
                $hash{$check} = join("\t",$hash{$check},@row2[5..7]);
                ## Add the variant to the duplicate check hash
                $dup{$check} = "exists";
            }
        }
        #print "$_ $hash{$_}\n" for (keys %hash);
        print "$counter duplicates $file\n";
    }
    close $file;
}

print $out "$hash{$_}\n" for (keys %hash);

#close $file;
#close $file1;
#close $file2;
close $out;
