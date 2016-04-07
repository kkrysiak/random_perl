#!/usr/bin/perl
use strict;
use warnings;

#### Open relavent input files
## Our variant file
open(VARIANTS, "</gscmnt/gc2547/mardiswilsonlab/kkrysiak/Stat1/maf_filtering/Table_S1.MAF_noquote.tsv") or die "Variant file not found";
## gz VCF file containing all ExAc variants
open(MGP_SNP, "gunzip -c /gscmnt/gc2547/mardiswilsonlab/kkrysiak/sanger_MGP/mgp.v2.snps.annot.reformat.vcf.gz | ") or die "Can't open ExAc file";
open(MGP_INDEL, "gunzip -c /gscmnt/gc2547/mardiswilsonlab/kkrysiak/sanger_MGP/mgp.v2.indels.annot.reformat.vcf.gz | ") or die "Can't open ExAc file";

#### Create output files
open(PASS, ">mgp_pass.tsv");
open(FAIL, ">mgp_fail.tsv");

## Declare QUAL cutoff to separate the file
my $qual_cutoff = 999;

## Create a hashes of passed and failed variants
my %fail = ();
my %pass = ();

## Initialize a variable to track progress
my $prog = "";

## Subroutine to extract the proper coordinates from the MGP files assign them to a pass or fail hash based on the qual cutoff
sub mgpTest {
    my $line = shift;
    $prog = shift;
    my $fail = shift;
    my $pass = shift;

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
            ## Initialize the final MGP variable and a counter
            my $mgp = "";
            my $i = 0;
            my $mgp2 = "";
            
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
                    $mgp = join("\t",$vars[0],$pos,$new_ref,$new_var);

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
                    ## Create a second mgp string using the reverse order
                    $mgp2 = join("\t",$vars[0],$pos2,$new_ref2,$new_var2);

                    ## Check if the 2 mgp strings are different
                    if($mgp eq $mgp2) {
                        $mgp2 = "";
                    }

                } else {
                    ## For variants with only 1 ref and 1 var base, create the mgp string
                    $mgp = join("\t",$vars[0],$vars[1],$vars[3],$alt[$i]);
                }

                ## Print a progress message to the user
                if($prog eq $vars[0]) {
                } else {
                    print "Processing MGP chromosome $vars[0]\n";
                    $prog = $vars[0];
                }

                ## Initialize and extract the qual variable
                my $qual = sprintf("%d", $vars[5]);

                ## Assign variants as passed or failed
                if($qual<$qual_cutoff) {
                    $$fail{$mgp} = $qual;
                    if($mgp2 ne "") {
                        $$fail{$mgp2} = $qual;
                    }
                } elsif($qual >= $qual_cutoff) {
                    $$pass{$mgp} = $qual;
                    if($mgp2 ne "") {
                        $$pass{$mgp2} = $qual;
                    }
                } else {
                    print "MGP quality score not valid, check formatting.\n";
                }                        
                $i++;        
            } 
        } else {
            print "$line\n";
        }
    }
}

## Iterate through the MGP SNP VCF file line by line
while(my $line = <MGP_SNP>){  ## read single line from the file
    chomp($line);
    mgpTest($line, $prog, \%fail, \%pass);
}
## Iterate through the MGP indel VCF file line by line
while(my $line = <MGP_INDEL>){  ## read single line from the file
    chomp($line);
    mgpTest($line, $prog, \%fail, \%pass);
}

## Initialize header
my $header = "";

## Read in variant file
while(my $line = <VARIANTS>) {
    chomp($line);
    ## Print the first line of the file + QUAL heading to output file
    if($header eq "") {
        $header = join("\t",$line,"MGP_QUAL_Score");
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
close MGP_SNP;
close MGP_INDEL;
close PASS;
close FAIL;
