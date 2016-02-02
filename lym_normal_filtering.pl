#!/usr/bin/perl
use strict;
use warnings;

## Open readcount file
open(NORM, "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/lym_normal_indel_readcounts.tsv") or die "Input file not found";

## Create output file
open(OUT, "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/lym_normal_output.tsv");

## Declare cutoffs
my $VAF = 5; ## minimum VAF
my $readcount = 3; ## minimum readcount
my $recurrence = 3; ## minimum recurrence count



## Identify variants to be removed
while(my $line = <NORM>) {
    my @var = split("\t",$line);
    ## Skip the header
    if($var[0] =~ /^chr/) {
    } else { 
        my @var = split("\t",$line);
        my $variant = join("\t",@var[0..4]);
        my $count = 0;
        ## Iterate through the normal samples for each variant
        for($v=6; $v<@var; $v+=3) {
            ## Check the VAF and variant readcount for each normal sample, summarizing the number of samples
            if($var[$v]>=$VAF && $var[$v-1]>=$readcount) {
                $count += 1;
            }
        }
        if($count>=$recurrence){
            print $variant;
        }
    }
}

## Close files
close NORM;
close OUT;
