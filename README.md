random_perl
===========
Collection of random perl scripts. Use at your own discretion/risk and please let me know if there are any bugs.

exac_annotation_and_filter.pl
-----------------------------
Annotate variant file with ExAC adjusted allele frequencies and separate into passed or failed variants
- Uses chr/start/ref/var to match the ExAC variants in the VCF file
- Several manipulations to the ExAC variants must occur to allow for match
    - Separates multiple alternate alleles represented on a single line
    - Trims prefixes/suffixes from ref/alt representations for indels and adjusts start site as necessary (Insertion - ExAC format: A/AAC, MGI format: -/AC)
    - Uses the adjusted allele count divided by the total adjusted allele count (AC_Adj/AN_Adj) to calculate the adjusted allele count
    - AF is basic allele frequency based on AC/AN, adjusted figures only count alleles with DP >= 10 & GQ >= 20

*Input:*

    - Modified MGI annotation format file with first 5 columns = chr/start/stop/ref/var (User input), additional columns are ignored but retained
    - Accepts user-defined exac cutoff or uses default 0.001
    - ExAC release 2 VCF file

*Output:*
    
    2 files containing all original variants and annotation/extra columns + ExAC adjusted allele frequency column split into passed and failed variants


PON_filtering.pl
----------------
Filter variants using panel of normal BAM files 
- Uses chr/start/stop/ref/var to match variant to the panel of normals file
- Counts the number of normals with readcount/VAF >= set cutoffs and prints to the appropriate file
- SNV/indel cutoffs are distinct, defaults set for each

*Input:*

    Variant file with first 5 columns = chr/start/stop/ref/var, additional columns retained but ignored 
    Panel of normals add-readcount output file with first 5 columns as above followed by 3 columns (ref/alt/VAF) for each normal to be used in the filter
*Output:*
    
    Original input variants split into passed and failed + final column with normal recurrence count or PASS indicated


 

**User beware - Hard coded scripts**
====================================
sangerMGP_annotation_and_filter.pl
----------------------------------
Split mouse variants by their presence in the sanger mouse genomes project

*Input:* Hard coded input

    Modified MGI annotation format file with the first 5 columns = chr/start/stop/ref/var, specifically variants from the Stat1 project

    Sanger MGP Release version 2 (mm9) SNPs and indel files (REL-1211)
*Output* Hard coded, placed in the current working directory

    mgp_pass        All variants not in the MGP files or MGP variants with a quality score below a minimum (hard-coded at 1) cutoff + 1 col of MGP quality score

    mgp_fail        All variants that fail + MGP quality and MGP filter columns added on the end


lym_normal_filtering.pl
-----------------------
Apply a panel of lymphoma normals filter to follicular lymphoma variants
- Uses chr/start/stop/ref/var to match variant to the panel of normals file
- Counts the number of normals with readcount/VAF >= set cutoffs and prints to the appropriate file
- SNV/indel cutoffs are distinct, also hard coded
    - minimimum VAF 2.5%/2.5%
    - minimum readcount 3/2
    - minimum recurrence 5/7 

*Input:* Hard coded input

    Modified MGI annotation format variant file with first 5 columns = chr/start/stop/ref/var, specifically variants from the follicular lymphoma project

    Panel of normals add-readcount output file with first 5 columns as above followed by 3 columns (ref/alt/VAF) for each normal to be used in the filter,
        specifically variants from the follicular lymphoma project
*Output:* Hard coded, placed in the current working directory

    All_Variants.lym_pon.tsv                Variants that pass the VAF/Readcount/Recurrence cutoffs with new last column

    All_Variants.lym_pon_excluded.tsv       Variants that fail the VAF/Readcount/Recurrence cutoffs with new last column providing the number of normals 
                                            above the readcount/VAF cutoff


exac_filtering_lymphoma.pl
--------------------------
Annotate my input list of lymphoma variants with the ExAC allele frequency
- Uses chr/start/ref/var to match the ExAC variants in the VCF file
- Several manipulations must occur to allow for match
    - Separates multiple alternate alleles represented on a single line
    - Trims prefixes/suffixes from ref/alt representations for indels and adjusts start site as necessary (Insertion - ExAC format: A/AAC, MGI format: -/AC)
    - Uses the adjusted allele count divided by the total adjusted allele count (AC_Adj/AN_Adj) to calculate the adjusted allele count

*Input:* Hard coded input
    
    Modified MGI annotation format file with first 5 columns = chr/start/stop/ref/var, specifically variants from the follicular lymphoma project

    ExAC release 2 VCF file, GRCh37
*Output:* Hard coded, placed in the current working directory

    exac_matched_variants_output.tsv                Original input file + 1 column of ExAC adjusted allele frequencies (AC_adj/AN_adj)
    
    exac_matched_variants_output_fullannot.tsv      Original input file + 12 columns of ExAC data for additional review (chr,start,reference allele,
                                                    alt allele, allele frequency, variant allele count (AC), adjusted allele count (AC_Adj), total
                                                    allele count (AN), adjusted total allele count (AN_Adj), quality score, filter, info field

remove_exac_variants.pl
-----------------------
Apply an ExAC allele frequency filter to a variant file
- Uses chr/start/stop/ref/var to match variants
- Exports variant file + either allele frequency or NA for passed variants, allele frequency for failed variants
- Hard coded cutoff is <=0.001 pass

*Input:* Hard coded input
    
    Modified MGI annotation format file with first 5 columns = chr/start/stop/ref/var, currently exac annotated follicular lymphoma file
    
    Must have ExAC_adj_AF column to filter by
    
    Output file from exac_filtering_lymphoma.pl (exac_matched_variants_output.tsv)

*Output:* Hard coded, placed in the current working directory
    
    All_Variants.lym_pon.exac.tsv           Variants that pass the ExAC allele frequency cutoff 
    All_Variants.lym_pon.exac_excluded.tsv  Variants that fail the ExAC allele frequency cutoff

