#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV_XS;
use Getopt::Long;

my $dirname = '';
my $sample_list = '';
my $header_file = '';
my $outdir = '';
my $outpre = 'all_variants';

GetOptions ('directory=s'=>\$dirname, 'sample_list_file=s'=>\$sample_list, 'header_file=s'=>\$header_file, 'output_dir=s'=>\$outdir, 'output_prefix=s'=>\$outpre);

my $usage=<<INFO;
    Gather all variants in a Lymphoma subdirectory for filtering.

    Example usage:
        perl somatic_validation_summary_filter.pl --directory=/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/ --sample_list_file=/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/sample_list_f49450ad26db40e1b065b7a0a79d85e2.txt --header_file=/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/header.txt --output_dir=/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/ --output_prefix=all_variants

INFO

## Define the top level directory where the summarized somatic validation results have been placed
#my $dirname = "/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/";
#opendir(DIR, $dirname) or die "Not able to open $dirname $!";

## Set up to parse TSV files
my $tsv = Text::CSV_XS->new ({
    binary => 1,        ## Allows non-ASCII characters including new lines
    eol => $/,          ## End of line string used
    auto_diag => 1,     ## Reports errors
    sep_char => "\t"    ## Defines character separator as tab
});

## Get the list of samples to include
#open my $sample_list, '<', '/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/sample_list_f49450ad26db40e1b065b7a0a79d85e2.txt';
open my $sample_fh, '<', $sample_list or die "Sample list ($sample_list) not found.\n";
chomp(my @samples = <$sample_fh>);
close $sample_fh;

######### Define variables to keep variants
## Define chromosomes to keep
my @chr_keep = (1..22,"X","Y");
## Define transcript errors to allow
my @trans_keep = ("no_errors","-");
## Define minimum coverage, variant count and VAF
my $min_cov = 25;
my $min_var = 4;
my $min_vaf = 5;
## Define tiers to keep
my @tiers = ("tier1");
## Define trv types to include
my @trv_keep = ("3_prime_untranslated_region","5_prime_untranslated_region","frame_shift_del","frame_shift_ins","in_frame_del","missense","nonsense","nonstop","splice_site","splice_site_del","splice_site_ins","in_frame_ins"); 
            ###### Keep RNA??

######### Filter
## Get the variant file header into an array
#open my $header_file, '<', "/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/header.txt";
open my $header_fh, '<', $header_file or die "Header file ($header_file) not found.\n";
my $header_line = <$header_fh>;
chomp($header_line);
my @header = split(/\s+/, $header_line);
close $header_fh;

## Create output files
open my $all, '>', join("/",$outdir,join("",$outpre,".tsv"));
open my $chr, '>', join("/",$outdir,join("",$outpre,".chr.tsv"));
open my $trans, '>', join("/",$outdir,join("",$outpre,".chr.no_error.tsv"));
open my $cov, '>', join("/",$outdir,join("",$outpre,".chr.no_error.coverage.tsv"));
open my $tier_file, '>', join("/",$outdir,join("",$outpre,".chr.no_error.coverage.tier1.tsv"));
open my $coding, '>', join("/",$outdir,join("",$outpre,".chr.no_error.coverage.tier1.coding.tsv"));

## Print headers
print $all join("\t", @header), "\tsample\n";
print $chr join("\t", @header), "\tsample\n";
print $trans join("\t", @header), "\tsample\n";
print $cov join("\t", @header), "\tsample\n";
print $tier_file join("\t", @header), "\tsample\n";
print $coding join("\t", @header), "\tsample\n";

## Iterate through each file in the directory
foreach my $s (@samples) {
    my $file = join("",$dirname,$s,'/snvs.indels.annotated');
    #my $file = '/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/collect_variants/H_ML-1017096/snvs.indels.annotated';
    open my $io, "<", $file or die "$file: $!";
 
    ## Pulls file in using column names as the key
    $tsv->column_names($tsv->getline ($io));
    while (my $row = $tsv->getline_hr ($io)) {
        ## Print every line to the all variants file
        foreach my $g (@header) { 
            print $all "$row->{$g}\t";
        }
        print $all "$s\n";
        
        ## Remove unwanted chromosomes
        my $chrs = $row->{chromosome_name};
        if( grep(/$chrs/, @chr_keep) ) {
            foreach my $h (@header) {
                print $chr "$row->{$h}\t";
            }
            print $chr "$s\n";

            ## Remove variants with transcript errors
            my $tr_error = $row->{transcript_error};
            if( grep(/$tr_error/, @trans_keep) ) {
                foreach my $h (@header) {
                    print $trans "$row->{$h}\t";
                }
                print $trans "$s\n";

                ## Remove variants with insufficient coverage or VAF
                my $cov_tot = $row->{val_Tumor_ref_count} + $row->{val_Tumor_var_count};
                my $var_cov = $row->{val_Tumor_var_count};
                my $vaf = $row->{val_Tumor_VAF};
                if($cov_tot >= $min_cov && $var_cov >= $min_var && $vaf >= $min_vaf) {
                    foreach my $h (@header) {
                        print $cov "$row->{$h}\t";
                    }
                    print $cov "$s\n";

                    ## Check if the variant is in tiers to keep list
                    my $tier = $row->{tier};
                    if( grep(/$tier/, @tiers) ) {
                        foreach my $h (@header) { 
                            print $tier_file "$row->{$h}\t";
                        }
                        print $tier_file "$s\n";

                        ## Check if the trv type indicates the variant should be kept
                        my $trv = $row->{trv_type};
                        if( grep(/$trv/, @trv_keep) ) {    
                            foreach my $h (@header) { 
                                print $coding "$row->{$h}\t";
                            }
                            print $coding "$s\n";
                        }
                    }
                }
            }
        }
    }
}

close $all;
close $chr;
close $trans;
close $cov;
close $tier_file;
close $coding;
