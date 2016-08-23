#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV_XS;

## Define the top level directory where the summarized somatic validation results have been placed
my $dirname = "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/";
#opendir(DIR, $dirname) or die "Not able to open $dirname $!";

## Set up to parse TSV files
my $tsv = Text::CSV_XS->new ({
    binary => 1,        ## Allows non-ASCII characters including new lines
    eol => $/,          ## End of line string used
    auto_diag => 1,     ## Reports errors
    sep_char => "\t"    ## Defines character separator as tab
});

## Get the list of samples to include
open my $sample_list, '<', '/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/sample_list_f49450ad26db40e1b065b7a0a79d85e2.txt';
    chomp(my @samples = <$sample_list>);
close $sample_list;

my $file = join("",$dirname,$samples[0],'/snvs.indels.annotated');
#my $file = '/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/H_ML-1017096/snvs.indels.annotated';
open my $io, "<", $file or die "$file: $!";

my $header = $tsv->getline ($io);
print join("-", @$header), "\n\n";

while (my $row = $tsv->getline ($io)) {
    print join("-", @$row), "\n";
}
