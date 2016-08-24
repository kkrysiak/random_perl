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

## Define trv types to include
my @trv_keep = ("3_prime_untranslated_region","5_prime_untranslated_region","frame_shift_del","frame_shift_ins","in_frame_del","missense","nonsense","nonstop","splice_site","splice_site_del","splice_site_ins","in_frame_ins"); 
            ###### Keep RNA??


## Get the variant file header into an array
open my $header_file, '<', "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/header.txt";
my $header_line = <$header_file>;
chomp($header_line);
my @header = split(/\s+/, $header_line);
print "$header[0]\t$header[1]\t$header[12]\n";
close $header_file;

## Iterate through each file in the directory
foreach my $s (@samples) {
    my $file = join("",$dirname,$samples['$s'],'/snvs.indels.annotated');
    #my $file = '/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma_group2/collect_variants/H_ML-1017096/snvs.indels.annotated';
    open my $io, "<", $file or die "$file: $!";

    
    ## Pulls file in using column names as the key
    $tsv->column_names($tsv->getline ($io));
    while (my $row = $tsv->getline_hr ($io)) {

        ## Check if the trv type indicates the variant should be kept
        my $trv = $row->{trv_type};
        if( grep(/$trv/, @trv_keep) ) {    
            foreach my $h (@header) { 
                print "$row->{$h}\t";
            }
            print "$samples['$s']\n";
        }
    }
}
