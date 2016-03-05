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
open(PASS, ">exac_pass.tsv");
open(FAIL, ">exac_fail.tsv");

## Declaire the exac allele frequency cutoff to separate the file
my $af_cutoff = 0.001;

## Create a hashes of passed and failed variants
my %fail = ();
my %pass = ();

## Initialize a variable to track progress
my $prog = "";

## Iterate through the ExAc VCF file line by line
while(my $line = <EXAC>){  ## read single line from the file
    chomp($line);
    ## Skip lines starting with #
    if($line =~ /^\#/){
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

                ## Print a progress message to the user
                if($prog eq $vars[0]) {
                } else {
                    print "Processing ExAC chromosome $vars[0]\n";
                    $prog = $vars[0];
                }

                ## Initialize allele frequency variable
                my $AF_adj = 0;
                ## Collect desired columns (adjusted allele count, adjusted total count) from INFO field
                my @info = split(";",$vars[7]);
                ## Set up to output allele counts and numbers
                my @AC = split(",",$info[0]);
                my @AC_adj = split(",",$info[3]);
                my $AN = $info[12];
                my $AN_adj = $info[15];
                $AC[$i] =~ s/AC=//;
                $AC_adj[$i] =~ s/AC_Adj=//;
                $AN =~ s/AN=//;
                $AN_adj =~ s/AN_Adj=//;
                ## Change the allele count to a number
                $AC_adj[$i] = sprintf("%.2f",$AC_adj[$i]);
                ## Create an adjusted allele frequency based on the adjusted allele counts as a decimal (not scientific notation)
                if($AC_adj[$i] == 0) {
                    $AF_adj = $AC_adj[$i];
                } else {
                    $AF_adj = sprintf("%.8f",$AC_adj[$i]/$AN_adj);
                }
#            print "$exac\t$AF_adj\n";
                ## Assign 
                if($AF_adj>$af_cutoff) {
                    $fail{$exac} = $AF_adj;
                    if($exac2 ne "") {
                        $fail{$exac2} = $AF_adj;
                    }
                } elsif($AF_adj <= $af_cutoff) {
                    $pass{$exac} = $AF_adj;
                    if($exac2 ne "") {
                        $pass{$exac2} = $AF_adj;
                    }
                } else {
                    print "ExAc allele frequency not valid, check formatting.\n";
                }                        
                $i++;        
            } 
        } else {
            print "$line\n";
        }
    }
}

## Initialize header
my $header = "";

## Read in variant file
while(my $line = <VARIANTS>) {
    chomp($line);
    ## Print the first line of the file + 5 ExAc column headings to output file
    if($header eq "") {
        $header = join("\t",$line,"ExAC_adj_AF");
        print PASS "$header\n";
        print FAIL "$header\n";
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
            ## Test to see if the variant is in the pass or fail hash and print to the appropriate file
            if($fail{$var_string}) {
                print FAIL "$line\t$fail{$var_string}\n";
            }elsif($pass{$var_string}) {
                print PASS "$line\t$pass{$var_string}\n";
            }else{
                print PASS "$line\tNA\n";
            }
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



## Close files
close VARIANTS;
close EXAC;
close PASS;
close FAIL;