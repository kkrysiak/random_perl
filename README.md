# random_perl
Collection of random perl scripts


#### Annotate my input list of lymphoma variants with the ExAC allele frequency
exac_filtering_lymphoma.pl
## Input:
    ## Modified MGI annotation format file with first 5 columns = chr/start/stop/ref/var
## Output:
    ## Original input file + 1 column of ExAC allele frequencies
    exac_matched_variants_output.tsv
    ## Original input file + 12 columns of ExAC data for additional review (chr,start,reference allele, alt allele, allele frequency, variant allele count (AC), adjusted allele count (AC_Adj), total allele count (AN), adjusted total allele count (AN_Adj), quality score, filter, info field)
    exac_matched_variants_output_fullannot.tsv
## Uses chr/start/ref/var to match the ExAC variants in the VCF file
## Several manipulations must occur to allow for match
    ## Separates multiple alternate alleles represented on a single line
    ## Trims prefixes/suffixes from ref/alt representations for indels and adjusts start site as necessary (Insertion - ExAC format: A/AAC, MGI format: -/AC)


#### Apply a panel of lymphoma normals filter to follicular lymphoma variants
lym_normal_filtering.pl
## Input:
    ## Modified MGI annotation format variant file with first 5 columns = chr/start/stop/ref/var
    ## Panel of normals add-readcount output file with first 5 columns as above followed by 3 columns (ref/alt/VAF) for each normal to be used in the filter
## Output:
    ## Variants that pass the VAF/Readcount/Recurrence cutoffs
        All_Variants.lym_pon.tsv
    ## Variants that fail the VAF/Readcount/Recurrence cutoffs with appended number of normals above the readcount/VAF cutoff
        All_Variants.lym_pon_excluded.tsv
## Uses chr/start/stop/ref/var to match variant to the panel of normals file
## Counts the number of normals with readcount/VAF >= set cutoffs and prints to the appropriate file

#### Apply an ExAC allele frequency filter to a variant file
remove_exac_variants.pl
## Input:
    ## Modified MGI annotation format file with first 5 columns = chr/start/stop/ref/var
    ## Output file from exac_filtering_lymphoma.pl (exac_matched_variants_output.tsv)
## Output:
    ## Variants that pass the ExAC allele frequency cutoff 
    All_Variants.lym_pon.exac.tsv
    ## Variants that fail the ExAC allele frequency cutoff
    All_Variants.lym_pon.exac_excluded.tsv
## Uses chr/start/stop/ref/var to match variants
## Exports variant file + either allele frequency or NA for passed variants, allele frequency for failed variants
