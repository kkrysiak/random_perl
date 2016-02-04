#!/usr/bin/perl
use strict;
use warnings;

#### Open relavent input files
## Our variant file
open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/All_Variants.tsv") or die "Variant file not found";
## gz VCF file containing all ExAc variants
open(EXAC, "gunzip -c /gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/ExAC.r0.2.sites.vep.vcf.gz | ") or die "Can't open ExAc file";
#open(EXAC, "/gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/test.vcf") or die "Can't open test file";

#### Create output files
open(OUT, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/lymphoma/variant_files/exac_matched_variants_output.tsv");

## Create hash to hold final cleaned up test variants
my %variants = ();
## Pull header from test variant file
my $header = "";

#### Read in test variants line by line, clean them up and assign them to a hash
while(my $line = <VARIANTS>) {
    chomp($line);
    ## Print the first line of the file + 5 ExAc column headings to output file
    if($header eq "") {
        $header = join("\t",$line,"exac_chr","exac_start","exac_ref","exac_alt","exac_AF");
        print OUT "$header\n";
    }
    ## Check that the first character is a valid chromosome, otherwise print to the removed output file
    if($line =~ /^([0-9]|X|Y|M)/){
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT)$/ && $vars[1] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
            ## Change - to 0 to designate indels
            if($vars[3] =~ /-/){
                $vars[3] = 0;
            } elsif ($vars[4] =~ /-/){
                $vars[4] = 0;
            }
            ## Create a new string with the approved chromosome, start, reference base, variant base
            my $var_string = join("\t",$vars[0],$vars[1],@vars[3..4]);
            ## Assign the new string as a key in a hash
            $variants{$var_string} = $line;           
        } else {
            unless($line =~ /^chr/) {
                ## Print skipped lines
                print "Variant formatting problem:\n$line\n";
            }
        }
    } else {
        unless($line =~ /^chr/) {
            ## Print skipped lines
            print "Variant formatting problem:\n$line\n";
        }
    }
}

## Iterate through the ExAc VCF file line by line
while(my $line = <EXAC>){  ## read single line from the file
    chomp($line);
    ## Skip lines starting with #
    if($line =~ /^\#/){
#        print "$line\n";
    } else {
        my $ins = 0;
        ## create an array with the different columns of the line
        my @vars = split("\t", $line);
        
        ## Check that the chromosome and start site are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y|MT)$/ && $vars[1] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G)+/ && $vars[4] =~ /(A|C|T|G)+/) {
            ## Split the alt allele into a new array (some have more than 1 alternate allele)
            my @alt = split(",", $vars[4]);
            ## Initialize the final ExAc variable and a counter
            my $exac = "";
            my $i = 0;
            my $exac2 = "";
            
            ## Loop through all alternate alleles
            while(exists $alt[$i]){
                ## Initialize new ref/alt and start variables
                my $new_ref = $vars[3];
                my $new_var = $alt[$i];
                my $pos = $vars[1];

                ## Adjust ref/var to match our formatting as well as start position
                if(length($new_ref)>1 || length($new_var)>1) {
                    ## Create a set of strings for comparison
                    my $new_ref2 = $new_ref;
                    my $new_var2 = $new_var;
                    my $pos2 = $pos;

                    ## Remove matching suffixes
                    while(substr($new_ref,(length($new_ref)-1)) eq substr($new_var,(length($new_var)-1))) {
                        chop($new_ref);
                        chop($new_var);
                    }
                    ## Remove matching prefixes and iterate coordinate start position
                    while(substr($new_ref,0,1) eq substr($new_var,0,1)) {
                        ## Replace reference with 0 for insertions
                        if(length($new_ref)==1 && length($new_var)>1) {
                            $new_ref = 0;
                            $new_var = substr($new_var,1,length($new_var));
                        ## Replace variant with 0 for deletions 
                        } elsif(length($new_ref)>1 && length($new_var)==1) {
                            $new_ref = substr($new_ref,1,length($new_ref));
                            $new_var = 0;
                            $pos = $pos + 1;
                        ## Remove shared starting bases 
                        } elsif(length($new_ref)>1 && length($new_var)>1) {
                            $new_ref = substr($new_ref,1,length($new_ref));
                            $new_var = substr($new_var,1,length($new_var));
                            $pos = $pos + 1;
                        }
                    }
                    $exac = join("\t",$vars[0],$pos,$new_ref,$new_var);

                    ## Determine if there is another way to represent the variant
                    ## Remove matching prefixes and iterate coordinate start position as necessary (first this time)
                    while(substr($new_ref2,0,1) eq substr($new_var2,0,1)) {
                        ## Replace reference with 0 for insertions
                        if(length($new_ref2)==1 && length($new_var2)>1) {
                            $new_ref2 = 0;
                            $new_var2 = substr($new_var2,1,length($new_var2));
                        ## Replace variant with 0 for deletions 
                        } elsif(length($new_ref2)>1 && length($new_var2)==1) {
                            $new_ref2 = substr($new_ref2,1,length($new_ref2));
                            $new_var2 = 0;
                            $pos2 = $pos2 + 1;
                        ## Remove shared starting bases 
                        } elsif(length($new_ref2)>1 && length($new_var2)>1) {
                            $new_ref2 = substr($new_ref2,1,length($new_ref2));
                            $new_var2 = substr($new_var2,1,length($new_var2));
                            $pos2 = $pos2 + 1;
                        }
                        ## Increase the starting base position
                    }
                    ## Remove matching suffixes
                    while(substr($new_ref2,(length($new_ref2)-1)) eq substr($new_var2,(length($new_var2)-1))) {
                        if(length($new_ref2)==1 && length($new_var2)>1) {    
                            $new_ref2 = 0;
                            chop($new_var2);
                        } elsif(length($new_ref2)>1 && length($new_var2)==1) {
                            chop($new_ref2);
                            $new_var2 = 0;
                        } elsif(length($new_ref2)>1 && length($new_var2)>1) {
                            chop($new_ref2);
                            chop($new_var2);
                        }
                    }
                    ## Create a second exac string using the reverse order
                    $exac2 = join("\t",$vars[0],$pos2,$new_ref2,$new_var2);

                    ## Check if the 2 exac strings are different
                    if($exac eq $exac2) {
                        $exac2 = "";
                    }

                } else {
                    ## For variants with only 1 ref and 1 var base, create the exac string
                    $exac = join("\t",$vars[0],$vars[1],$vars[3],$alt[$i]);
                }

                ## Use the new ExAc string to see if it matches a key in the variant hash
                if($variants{$exac}) {
                    my @info = split(";",$vars[7]);
                    foreach my $l (@info) {
                        if ($l =~ /^AF/) {
                            ## Pull the allele frequency data (AF section)
#                            print "$l\n";
                            ## Like alternate alleles, allele frequences are separated by a , so pull the correct AF for the alt allele
                            $l =~ s/AF=//;
                            my @af = split(",",$l);

                            ## print the original test input line and the ExAc matched chr, ref, alt and allele frequency
                            print OUT "$variants{$exac}\t$exac\t$af[$i]\n";
                        }
                    }
                ## Use alternate exac string to query for the test variant                  
                } elsif($exac2 ne "") {
                    if($variants{$exac2}) {
                        my @info2 = split(";",$vars[7]);
                        foreach my $l2 (@info2) {
                            if ($l2 =~ /^AF/) {
                                ## Pull the allele frequency data (AF section)
                                ## Like alternate alleles, allele frequences are separated by a , so pull the correct AF for the alt allele
                                $l2 =~ s/AF=//;
                                my @af2 = split(",",$l2);

                                ## print the original test input line and the ExAc matched chr, ref, alt and allele frequency
                                print OUT "$variants{$exac2}\t$exac2\t$af2[$i]\n";
                            }
                        }
                    }    
                }

                $i++;        
            } 
        } else {
            print "$line\n";
        }
    }
}

## Close files
close VARIANTS;
close EXAC;
close OUT;
