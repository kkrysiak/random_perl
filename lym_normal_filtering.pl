#!/usr/bin/perl
use strict;
use warnings;

## Open readcount file
open(NORM, "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/lym_normal_readcounts.tsv") or die "Input file not found";
#open(NORM, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/test.tsv") or die "Test file not found";
## Open variant file
open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.tsv") or die "List of called variants not found";

## Create output files
open(OUT, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon_excluded.tsv");
open(OUT2, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.tsv");

## Declare cutoffs
my $VAF = 3;        ## minimum VAF
my $readcount = 2;  ## minimum readcount
my $recurrence = 5; ## minimum recurrence count

## Assign the variants being tested against to a hash
my %list = ();
my $row = 0;
while(my $line = <VARIANTS>) {
    chomp($line);
    $row+=1;
    if($row==1){
        print OUT "$line\tCount\n";
        print OUT2 "$line\n";
    } else {
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT)$/ && $vars[1] =~ /^\d+$/ && $vars[2] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
            ## Create a new string with the approved chromosome, start, reference base, variant base
            my $var_string = join("\t",@vars[0..4]);
            ## Assign the new string as a key in a hash
            $list{$var_string} = $line;
        }
    }
}
    

## Identify variants to be removed
while(my $line2 = <NORM>) {
    my @var = split("\t",$line2);
    ## Skip the header
    if($var[0] =~ /^chr/) {
    } else { 
        my @var = split("\t",$line2);
        my $variant = join("\t",@var[0..4]);
        my $count = 0;
        ## Iterate through the normal samples for each variant
        if($var[6] eq "NA") {
            print "WARNING: Data contains NA values\n";
        } else {
            for(my $v=7; $v<@var; $v+=3) {
                if($var[$v]>=$VAF && $var[$v-1]>=$readcount) {
                    $count+=1;
                }
            }
            ## Print the variants to the appropriate file (exclude or not)
            if($count>=$recurrence){
                print OUT "$list{$variant}\t$count\n";
            } else {
                if($var[6] eq "NA") {
                } else {
                    print OUT2 "$list{$variant}\n";
                }
            }
        }
    }
}

## Close files
close NORM;
close VARIANTS;
close OUT;
close OUT2;
