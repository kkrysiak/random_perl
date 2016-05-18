#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

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
        PON_filtering.pl --readcount_file=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/lym_normal_readcounts.tsv 
        --outfile_prefix='test' --variant_file=/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.tsv --snv_min_vaf=5 
        --snv_min_readcount=3 --snv_min_recurrence=5 --indel_min_vaf=2

    Requires
        --readcount_file        add-reacount output file with ONLY normal readcount output. First 5 cols (chr,start,stop,ref,var) followed by ref_count, 
                                var_count, VAF for each normal sample
        --variant_file          Variant file with first 5 col (chr,start,stop,ref,var). Additional columns will be retained. All variants will be passed to
                                pass or fail file. Last column will indicate # of recurrences in normal samples or NA if not found in readcount file.

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

## Print parameters to user
print "Using SNV minimum VAF: $VAF\n";
print "Using SNV minimum readcount cutoff: $readcount\n";
print "Using SNV minimum recurrence count: $recurrence\n";
print "Using indel minimum VAF: $indel_VAF\n";
print "Using indel minimum readcount: $indel_readcount\n";
print "Using indel minimum recurrence count: $indel_recurrence\n";

## Create output files
open my $out_pass, '>', join("", $prefix, ".pass.tsv"); 
open my $out_fail, '>', join("", $prefix, ".fail.tsv");

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
                    $pass{$variant} = $indel_count;
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
                    $pass{$variant} = $count;
                }
            } else {
                print "Variant format didn't match. Check the following is as expected:\n$variant\n";
            }
        }
    }
}

## Initialize counters
my $row = 0;
my $passed = 0;
my $failed = 0;
my $unexpected = 0;
my $skipped = 0;

while(my $line2 = <$var>) {
    chomp($line2);
    $row+=1;
    if($row==1){
        print $out_fail "$line2\tnorm_recurrence_count\n";
        print $out_pass "$line2\tnorm_recurrence_count\n";
    } else {
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line2);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT|GL\S+)$/ && $vars[1] =~ /^\d+$/ && $vars[2] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
            ## Create a new string with the approved chromosome, start, reference base, variant base
            my $var_string = join("\t",@vars[0..4]);
            if (defined $fail{$var_string}) {
                print $out_fail "$line2\t$fail{$var_string}\n";
                $failed++;
            } elsif(defined $pass{$var_string}) {
                print $out_pass "$line2\t$pass{$var_string}\n";
                $passed++;
            } else {
                #print "\nWARNING: Variant not found in readcount file. Variant passed by default.\n$var_string\n";
                $unexpected++;
                print $out_pass "$line2\tNA\n";
            }
        } else {
            print "\nWARNING: Check variant file formatting, chr/stop/start/ref/var not properly formatted. Variant skipped.\n$line2\n";
            $skipped++;
        }
    }
}

my $total = $passed + $unexpected;

print "\n$total passed variants\n$failed failed variants\n";
if($unexpected>0){
    print "$unexpected variants not found in normal readcount file - Passed by default\n";
}
if($skipped>0){
    print "$skipped variants skipped (not in output files) due to unsupported variant formatting\n";
}

## Close files
close $norm;
close $var;
close $out_pass;
close $out_fail;
