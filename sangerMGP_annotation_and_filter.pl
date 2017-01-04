#!/usr/bin/perl
use strict;
use warnings;

#### Open relavent input files
## Our variant file
open(VARIANTS, "</gscmnt/gc2547/griffithlab/kkrysiak/Stat1/maf_filtering/Table_S1.MAF_noquote.tsv") or die "Variant file not found";
#open(VARIANTS, "</gscmnt/gc2547/griffithlab/kkrysiak/Stat1/maf_filtering/test.tsv") or die "Test file not found";

## gz VCF file containing all MGP variants
open(MGP_SNP, "gunzip -c /gscmnt/gc2547/griffithlab/kkrysiak/sanger_MGP/mgp.v2.snps.annot.reformat.vcf.gz | ") or die "Can't open ExAc file";
open(MGP_INDEL, "gunzip -c /gscmnt/gc2547/griffithlab/kkrysiak/sanger_MGP/mgp.v2.indels.annot.reformat.vcf.gz | ") or die "Can't open ExAc file";

#### Create output files
open(PASS, ">mgp_pass.tsv");
open(FAIL, ">mgp_fail.tsv");

## Declare QUAL cutoff to separate the file
my $qual_cutoff = 1;

## Initialize a hash of variants to be tested
my %variants = ();
## Initialize a hash of repeated variants to be tested
my %rep_vars = ();

## Initialize a variable to track progress
my $prog = "";
## Initialize header variable
my $header = "";

## Read in variant file and assign it to a hash
while(my $line = <VARIANTS>) {
    chomp($line);
    ## Print the first line of the file + QUAL heading to output file
    if($header eq "") {
        $header = join("\t",$line,"MGP_QUAL_Score");
        print PASS "$header\n";
        print FAIL "$header\tMGP_FILTER\n";
        print "Reading in variant file\n";
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
            ## Check if the variant is already a key in the hash
            if($variants{$var_string}) {
                push(@{$rep_vars{$var_string}}, $line);
            } else {
                ## assign new variant string as a key to the variant hash and the line to the value
                $variants{$var_string} = $line;        
            }
        } else {
            unless($line =~ /^chr/) {
                # Print skipped lines
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

## Subroutine to extract the proper coordinates from the MGP files assign them to a pass or fail hash based on the qual cutoff
sub mgpTest {
    my $line = shift;
    $prog = shift;
    my $variants = shift;
    ## %rep_vars is hash of arrays
    my $rep_vars_ref = shift;  
    my %rep_vars = %$rep_vars_ref;

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
                
                ## Check if the MGP variant is in the tested list of variants
                if($$variants{$mgp} || $$variants{$mgp2}) {
                    ## Assign variants as passed or failed and remove them from the hashes after printing to the appropriate file
                    if($vars[6] eq "PASS") {
                        print FAIL "$$variants{$mgp}\t$qual\t$vars[6]\n";
                        delete $$variants{$mgp};
                        ## Print variants with repeat entries (2+ lines in the input variant file)
                        if($rep_vars{$mgp}) {
                            foreach my $r (@{$rep_vars{$mgp}}) {
                                print FAIL "$r\t$qual\t$vars[6]\n";
                            }
                            delete $$rep_vars_ref{$mgp};
                        }
                        ## Repeat above for alterative representation of the variant if it exists
                        if($$variants{$mgp2}) {
                            print FAIL "$$variants{$mgp2}\t$qual\t$vars[6]\n";
                            delete $$variants{$mgp2};
                            if($rep_vars{$mgp2}) {
                                foreach my $r (@{$rep_vars{$mgp2}}) {
                                    print FAIL "$r\t$qual\t$vars[6]\n";
                                }
                                delete $$rep_vars_ref{$mgp2};
                            }
                        }
                    } else { 
                        if($qual >= $qual_cutoff) {
                            print FAIL "$$variants{$mgp}\t$qual\t$vars[6]\n";
                            delete $$variants{$mgp};
                            if($rep_vars{$mgp}) {
                                foreach my $r (@{$rep_vars{$mgp}}) {
                                    print FAIL "$r\t$qual\t$vars[6]\n";
                                }
                                delete $$rep_vars_ref{$mgp};
                            }
                            if($$variants{$mgp2}) {
                                print FAIL "$$variants{$mgp2}\t$qual\t$vars[6]\n";
                                delete $$variants{$mgp2};
                                if($rep_vars{$mgp2}) {
                                    foreach my $r (@{$rep_vars{$mgp2}}) {
                                        print FAIL "$r\t$qual\t$vars[6]\n";
                                    }                                    
                                    delete $$rep_vars_ref{$mgp2};
                                }
                            }
                        } elsif($qual < $qual_cutoff) {
                            print PASS "$$variants{$mgp}\t$qual\n";
                            delete $$variants{$mgp};
                            if($rep_vars{$mgp}) {
                                foreach my $r (@{$rep_vars{$mgp}}) {
                                    print PASS "$r\t$qual\n";
                                }
                                delete $$rep_vars_ref{$mgp};
                            }
                            if($$variants{$mgp2}) {
                                print PASS "$$variants{$mgp2}\t$qual\n";
                                delete $$variants{$mgp2};
                                if($rep_vars{$mgp2}) {
                                    foreach my $r (@{$rep_vars{$mgp2}}) {
                                        print PASS "$r\t$qual\n";
                                    }
                                    delete $$rep_vars_ref{$mgp2};
                                } 
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

## Iterate through the MGP SNP VCF file line by line
print "Testing SNPs\n";
while(my $line = <MGP_SNP>){  ## read single line from the file
    chomp($line);
    mgpTest($line, $prog, \%variants, \%rep_vars);
}
## Iterate through the MGP indel VCF file line by line
print "Testing indels\n";
while(my $line2 = <MGP_INDEL>){  ## read single line from the file
    chomp($line2);
    mgpTest($line2, $prog, \%variants, \%rep_vars);
}

print "Printing remaining variants\n";
print PASS "$variants{$_}\tNA\n" for (keys %variants);
## Print remaining repeated variants
foreach my $i (keys %rep_vars) {
    foreach my $n (@{$rep_vars{$i}}) {
        print PASS "$n\tNA\n";
    }
}

## Close files
close VARIANTS;
close MGP_SNP;
close MGP_INDEL;
close PASS;
close FAIL;
