#!/usr/bin/perl
use strict;
use warnings;

#### Open relavent input files
## Our variant file
open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/lym_normal_filter/All_Variants.lym_pon.tsv") or die "Variant file not found";
## Open Exac-annotated variants file
open(ANNOT, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/matched_variants_output.tsv") or die "ExAc annotated file not found";

#### Create output files
open(KEEP, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.exac.tsv");
open(FAILED, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.lym_pon.exac_excluded.tsv");

## Declaire the exac allele frequency cutoff to separate the file
my $af_cutoff = 0.001;

## Create a hash of the variant file
my %filter = ();
my %fail = ();
my %pass = ();
## Pull header from the Exac annotated file
my $exac_header = "";


#while(my $fline = <VARIANTS>) {
#    chomp($fline);
#    ## Split the lines on tab
#    my @fvars = split("\t", $fline);
#    ## Create a key out of the chr,start,stop,ref,var and assign the line as the value in the hash
#    my $f_string = join("\t",@fvars[0..4];
#    $filter{$f_string} = $fline;
#}

while(my $eline = <ANNOT>) {
    chomp($eline);
    if($exac_header eq "") {
        print FAILED "$exac_header\n";
    }
    my @evars = split("\t", $eline);
    ## Create a key out of the chr,start,stop,ref,var
    my $e_string = join("\t",@evars[0..4]);
    ## Assign the line as the value to the fail or pass hash based on the allele frequency cutoff used
    if($evars[38]>$af_cutoff) {
        $fail{$e_string} = $eline;
    } elsif($evars[38]<=$af_cutoff) {
        $pass{$e_string} = join("\t",$evars[0..34],$evars[38]);
    } else {
        print "ExAc allele frequency not valid, check fromatting.\n";
    }
}

while(my $fline = <VARIANTS>) {
    chomp($fline);
    ## skip the header line
    if($fline =~ /^chr/) {
        print KEEP "$fline\texac_AF\n"
    } else {
        ## Split the lines on tab
        my @fvars = split("\t", $fline);
        ## Create a key out of the chr,start,stop,ref,var to call the line from the correct hash
        my $f_string = join("\t",@fvars[0..4]);
        ## Print the line to the correct file
        if($fail{$f_string}) {
            print FAILED "$fail{$f_string}\n";
        } elsif($pass{$f_string}) {
            print KEEP "$pass{$f_string}\n";
        } else {
            print KEEP "$fline\tNA\n";
        }
    }
}

close VARIANTS;
close ANNOT;
close KEEP;
close FAILED;
