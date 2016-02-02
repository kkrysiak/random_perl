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
open(OUT, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/matched_variants_output.tsv");
open(OUT2, ">/gscmnt/gc2547/mardiswilsonlab/kkrysiak/exac_release2/input_variants_removed.tsv");

## Create hash to hold final cleaned up test variants
my %variants = ();

#### Read in test variants line by line, clean them up and assign them to a hash
while(my $line = <VARIANTS>) {
    chomp($line);
    ## Check that the first character is a valid chromosome, otherwise print to the removed output file
    if($line =~ /^([0-9]|X|Y)/){
        ## Create an array with the different columns of the line
        my @vars = split("\t", $line);
        ## Check that the chromosome, start, reference base, variant base are formatted as expected
        if($vars[0] =~ /^(\d*|X|Y)$/ && $vars[1] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G|-|0)+/ && $vars[4] =~ /(A|C|T|G|-|0)+/) {
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
            ## Output skipped lines to a separate file
            print OUT2 "$line\n";
        }
    } else {
        ## Output skipped lines to a separate file
        print OUT2 "$line\n";
    }
}

#print scalar keys %variants;
#print "$_ $variants{$_}\n" for keys %variants;

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
        if($vars[0] =~ /^(\d*|X|Y)$/ && $vars[1] =~ /^\d+$/ && $vars[3] =~ /(A|C|T|G)+/ && $vars[4] =~ /(A|C|T|G)+/) {
            ## Split the alt allele into a new array (some have more than 1 alternate allele)
            my @alt = split(",", $vars[4]);
            ## Initialize the final ExAc variable and a counter
            my $exac = "";
            my $i = 0;
            
            ## Loop through all alternate alleles
            while(exists $alt[$i]){
                ## Initialize new ref/alt and start variables
                my $new_ref = $vars[3];
                my $new_var = $alt[$i];
                my $pos = $vars[1];

                ## Adjust ref/var to match our formatting as well as start position
                if(length($new_ref)>1 || length($new_var)>1) {
                    ## Remove matching suffixes
                    while(substr($new_ref,(length($new_ref)-1)) eq substr($new_var,(length($new_var)-1))) {
                        chop($new_ref);
                        chop($new_var);
                        my $a = substr($new_ref,(length($new_ref)-1));
                        my $b = substr($new_var,(length($new_var)-1));
#                        print "$vars[3]\t$alt[$i]\t$a\t$b\t$new_ref\t$new_var\n";
#                        print "$vars[3]\t$alt[$i]\t$vars[1]\tfirst:$new_ref\t$new_var\t$pos\n";
                    }
                    ## Remove matching prefixes and iterate coordinate start position
                    while(substr($new_ref,0,1) eq substr($new_var,0,1)) {
                        ## Replace reference with 0 for insertions and iterate start
                        if(length($new_ref)==1 && length($new_var)>1) {
                            $new_ref = 0;
                            $pos = $pos + 1;
                            $new_var = substr($new_var,1,length($new_var));
                        ## Replace variant with 0 for deletions 
                        } elsif(length($new_ref)>1 && length($new_var)==1) {
                            $new_ref = substr($new_ref,1,length($new_ref));
                            $pos = $pos + 1;   
                            $new_var = 0;
                        ## Remove shared starting bases 
                        } elsif(length($new_ref)>1 && length($new_var)>1) {
                            $new_ref = substr($new_ref,1,length($new_ref));
                            $pos = $pos + 1;
                            $new_var = substr($new_var,1,length($new_var));
                        }
#                        print "$vars[3]\t$alt[$i]\t$vars[1]\tsecond:$new_ref\t$new_var\t$pos\n";
                    }
#                    print "$vars[3]\t$alt[$i]\t$a\t$b\t$new_ref\t$new_var\n";
                    $exac = join("\t",$vars[0],$pos,$new_ref,$new_var);
#                    print "$vars[3]\t$alt[$i]\t$vars[1]\t\t$exac\n";



#                    if(length($new_ref)==1 && length($new_var)>1) {
#                        $new_ref = 0;
#                        my $c = substr($new_var,1,length($new_var));
#                        print "$vars[3]\t$alt[$i]\t$new_ref\t$c\n";
#                    }
                } else {
                    $exac = join("\t",$vars[0],$vars[1],$vars[3],$alt[$i]);
                }

#                if(length($vars[3])>1 || length($alt[$i])>1) {
#                    if(substr($vars[3],(length($vars[3])-1)) eq substr($alt[$i],(length($alt[$i])-1))) {
#                        chop($
#                        my $a = substr($vars[3],(length($vars[3])-1));
#                        my $b = substr($alt[$i],(length($alt[$i])-1));
#                        print "$vars[3]\t$alt[$i]\t$a\t$b\n";
#                    }
#                }
#                ## Above we made the ref/var alleles 0 for indels when creating the hash keys
#                ## Loop through all alternative alleles for a given variant
#                if(length($alt[$i])>1) {
#    ## TO DO: Handle 2bp reference alleles - for example chr 1:2234874 CG  C,CGG,TG 
#                    ## For alternate alleles that are more than 1 base, replace the reference allele with 0 (insertions)
#                    $exac = join("\t",$vars[0],$vars[1],$ins,$alt[$i]);
#                } elsif(length($vars[3])>1) {
#                    ## For reference alleles that are more than 1 base, replace the alternate allele with 0 (deletion)
#                    $exac = join("\t",$vars[0],$vars[1],$vars[3],$ins);
#                } else {
#                    $exac = join("\t",$vars[0],$vars[1],$vars[3],$alt[$i]);
#                }
#                
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
close OUT2;
