#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


## Declare cutoffs for SNVs
my $VAF = 2.5;
my $readcount = 3;
my $recurrence = 5;
my $indel_VAF = 2.5;
my $indel_readcount = 2;
my $indel_recurrence = 7;
my $normals = '';
my $variants = '';
my $prefix = 'PON';

GetOptions ('snv_min_vaf=f'=>\$VAF, 'snv_min_readcount=i'=>\$readcount, 'snv_min_recurrence=i'=>\$recurrence, 'indel_min_vaf=f'=>\$indel_VAF, 
            'indel_min_readcount=i'=>\$indel_readcount, 'indel_min_recurrence=i'=>\$indel_recurrence, 'readcount_file=s'=>\$normals, 'variant_file=s'=>\$variants, 
            'outfile_prefix=s'=>\$prefix);

my $usage=<<INFO;
    Example usage:
        PON_filtering.pl --readcount_file=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/lym_normal_readcounts.tsv --outfile_prefix='test' --variant_file=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.tsv --snv_min_vaf=5 --snv_min_readcount=3 --snv_min_recurrence=5 --indel_min_vaf=2

    Requires
        --readcount_file        add-reacount output file with only normal readcount output
        --variant_file          variant file with first 5 col (chr,start,stop,ref,var), additional columns will be retained

    Optional parameters
        --outfile_prefix        default="PON"   File name prefix for output files
        --snv_min_vaf           default=2.5     Minimum VAF allowed to count a SNV as detected in a normal sample
        --snv_min_readcount     default=3       Minimum number of variant reads to count a SNV as detected in a normal sample
        --snv_min_recurrence    default=5       Minimum number of normal samples where a SNV is detected for it to fail this filter
        --indel_min_vaf         default=2.5     Minimum VAF allowed to count an indel as detected in a normal sample
        --indel_min_readcount   default=2       Minimum number of variant reads to count an indel as detected in a normal sample
        --indel_min_recurrence  default=7       Minimum number of normal samples where an indel is detected for it to fail this filter

INFO

## Open input files
open my $norm, '<', $normals or die "Normal readcount file ($normals) not found.\n\n$usage";
open my $var, '<', $variants or die "Variant file ($variants) not found.\n\n$usage";

unless($VAF){
    print "Using default SNV minimum VAF: 2.5\n";
    $VAF = 2.5;
}
unless($readcount) {
    print "Using default SNV minimum readcount cutoff: 3\n";
    $readcount = 3;
}
unless($recurrence) {
    print "Using default SNV minimum recurrence count: 5\n";
    $recurrence = 5;
}
unless($indel_VAF) {
    print "Using default indel minimum VAF: 2.5\n";
    $indel_VAF = 2.5;
}
unless($indel_readcount) {
    print "Using default indel minimum readcount: 2\n";
    $indel_readcount = 2;
}
unless($indel_recurrence) {
    print "Using default indel minimum recurrence count: 7\n";
    $indel_recurrence = 7;
}

## Create output files
open my $out_pass, '>', join("", $prefix, "_PASS.tsv"); 
open my $out_fail, '>', join("", $prefix, "_FAIL.tsv");

#open(OUT, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon_excluded.tsv");
#open(OUT2, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.tsv");

## Create a hashes of passed and failed variants
my %fail = ();
my %pass = ();

## Identify variants to be removed
while(my $line = <$norm>) {
    my @var = split("\t",$line);
    ## Skip the header
    if($var[0] =~ /^chr/) {
    } else {
        ## Declare my variant as the first 5 columns, to be used as hash key 
        my $variant = join("\t",@var[0..4]);
        my $count = 0;
        my $indel_count = 0;
        ## Iterate through the normal samples for each variant
        if($var[6] eq "NA") {
            print "WARNING: Data contains NA values\n";
        } else {
            ## Check if the variant is an indel
            if($var[3] =~ /-|0/ || $var[4] =~ /-|0/) {
                ## Iterate through readcount columns, checking VAF and variant readcount
                for(my $vi=7; $vi<@var; $vi+=3) {
                    if($var[$vi]>=$indel_VAF && $var[$vi-1]>=$indel_readcount) {
                        $indel_count+=1;
                    }
                }
                ## Assign variant to pass or fail hash
                if($indel_count>=$indel_recurrence){
                    $fail{$variant} = $indel_count;
                } else {
                    $pass{$variant} = "PASSED";
                }
            ## Check if the variant is a SNV
            } elsif($var[3] =~ /^(A|C|T|G)$/ && $var[4] =~ /^(A|C|T|G)$/) {
                ## Iterate through readcount columns, checking VAF and variant readcount
                for(my $v=7; $v<@var; $v+=3) {
                    if($var[$v]>=$VAF && $var[$v-1]>=$readcount) {
                        $count+=1;
                    }
                }
                ## Assign variant to pass or fail hash
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

while(my $line2 = <$var>) {
    chomp($line2);
    $row+=1;
    if($row==1){
        print $out_fail "$line2\tCount\n";
        print $out_pass "$line2\n";
    } else {
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line2);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT)$/ && $vars[1] =~ /^\d+$/ && $vars[2] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
            ## Create a new string with the approved chromosome, start, reference base, variant base
            my $var_string = join("\t",@vars[0..4]);
            if (defined $fail{$var_string}) {
                print $out_fail "$line2\t$fail{$var_string}\n";
            } elsif(defined $pass{$var_string}) {
                print $out_pass "$line2\n";
            } else {
                print "Variant not found.\n$var_string\n";
            }
        } else {
            print "$line2\nCheck variant file formatting, indel Chr/Stop/Start/Ref/Var not properly formatted\n";
        }
    }
}


## Close files
close $norm;
close $var;
close $out_pass;
close $out_fail;
