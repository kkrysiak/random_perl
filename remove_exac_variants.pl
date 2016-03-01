#!/usr/bin/perl
use strict;
use warnings;

#### Open relavent input files
## Our variant file
open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.tsv") or die "Variant file not found";
#open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/All_Variants.lym_pon.tsv") or die "Variant file not found";
## Open Exac-annotated variants file
open(ANNOT, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/exac_matched_variants_output.tsv") or die "ExAc annotated file not found";
#open(ANNOT, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/matched_variants_output.tsv") or die "ExAc annotated file not found";

#### Create output files
open(KEEP, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.exac2.tsv");
open(FAILED, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.exac_excluded2.tsv");

## Declaire the exac allele frequency cutoff to separate the file
my $af_cutoff = 0.001;

## Create a hashes of passed and failed variants
my %fail = ();
my %pass = ();

## Initialize a variable for the exac allele frequency column
my $af_col = 0;

while(my $eline = <ANNOT>) {
    chomp($eline);
    ## If header, determine the allele frequency column
    if($eline =~ /^chr/) {
        my @colnames = split("\t",$eline);
        for my $col (@colnames) {
            last if($col =~ /exac_AF/);
            $af_col++;
        }
        if($colnames[$af_col] eq "exac_AF") {
            print "Using $colnames[$af_col] column for ExAC allele frequency filtering\n";
        } else {
            print "Can't find 'exac_AF' column for ExAC allele frequency filtering, check input file.\n";
        }
    } else {
        my @evars = split("\t", $eline);
        ## Create a key out of the chr,start,stop,ref,var
        my $e_string = join("\t",@evars[0..4]);
        ## Convert scientific notation of the allele frequency to something perl can handle
        my $af_value = sprintf("%.8f", $evars[$af_col]);
        ## Assign the line as the value to the fail or pass hash based on the allele frequency cutoff used
        if($af_value>$af_cutoff) {
            $fail{$e_string} = $af_value;
        } elsif($af_value <= $af_cutoff) {
            $pass{$e_string} = $af_value;
        } else {
            print "ExAc allele frequency not valid, check fromatting.\n";
        }
    }
}

while(my $fline = <VARIANTS>) {
    chomp($fline);
    ## skip the header line
    if($fline =~ /^chr/) {
        print KEEP "$fline\texac_AF\n";
        print FAILED "$fline\texac_AF\n";
    } else {
        ## Split the lines on tab
        my @fvars = split("\t", $fline);
        ## Create a key out of the chr,start,stop,ref,var to call the line from the correct hash
        my $f_string = join("\t",@fvars[0..4]);
        ## Print the line to the correct file
        if($fail{$f_string}) {
            print FAILED "$fline\t$fail{$f_string}\n";
        } elsif($pass{$f_string}) {
            print KEEP "$fline\t$pass{$f_string}\n";
        } else {
            print KEEP "$fline\tNA\n";
        }
    }
}

close VARIANTS;
close ANNOT;
close KEEP;
close FAILED;
