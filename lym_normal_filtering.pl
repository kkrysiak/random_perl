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

## Declare cutoffs for SNVs
my $VAF = 2.5;      ## minimum VAF
my $readcount = 3;  ## minimum readcount
my $recurrence = 5; ## minimum recurrence count

## Declare cutoffs for Indels
my $indel_VAF = 2.5;      ## minimum VAF
my $indel_readcount = 2;  ## minimum readcount
my $indel_recurrence = 7; ## minimum recurrence count

## Create a hashes of passed and failed variants
my %fail = ();
my %pass = ();

## Identify variants to be removed
while(my $line = <NORM>) {
    my @var = split("\t",$line);
    ## Skip the header
    if($var[0] =~ /^chr/) {
    } else { 
        my $variant = join("\t",@var[0..4]);
        my $count = 0;
        my $indel_count = 0;
        ## Iterate through the normal samples for each variant
        if($var[6] eq "NA") {
            print "WARNING: Data contains NA values\n";
        } else {
            if($var[3] =~ /-|0/ || $var[4] =~ /-|0/) {
                for(my $vi=7; $vi<@var; $vi+=3) {
                    if($var[$vi]>=$indel_VAF && $var[$vi-1]>=$indel_readcount) {
                        $indel_count+=1;
                    }
                }
                if($indel_count>=$indel_recurrence){
                    $fail{$variant} = $indel_count;
                } else {
                    $pass{$variant} = "PASSED";
                }
            } elsif($var[3] =~ /^(A|C|T|G)$/ && $var[4] =~ /^(A|C|T|G)$/) {
                for(my $v=7; $v<@var; $v+=3) {
                    if($var[$v]>=$VAF && $var[$v-1]>=$readcount) {
                        $count+=1;
                    }
                }
                if($count>=$recurrence){
                    $fail{$variant} = $count;
                } else {
                    $pass{$variant} = "PASSED";
                }
            } else {
                print "Variant format didn't match. Check the following is as expected:\n$variant\n";
            }
        }
    }
}

my $row = 0;

while(my $line2 = <VARIANTS>) {
    chomp($line2);
    $row+=1;
    if($row==1){
        print OUT "$line2\tCount\n";
        print OUT2 "$line2\n";
    } else {
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line2);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT)$/ && $vars[1] =~ /^\d+$/ && $vars[2] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
            ## Create a new string with the approved chromosome, start, reference base, variant base
            my $var_string = join("\t",@vars[0..4]);
            if (defined $fail{$var_string}) {
                print OUT "$line2\t$fail{$var_string}\n";
            } elsif(defined $pass{$var_string}) {
                print OUT2 "$line2\n";
            } else {
                print "Variant not found.\n$var_string\n";
            }
        } else {
            print "$line2\nCheck variant file formatting, indel Chr/Stop/Start/Ref/Var not properly formatted\n";
        }
    }
}


## Close files
close NORM;
close VARIANTS;
close OUT;
close OUT2;
