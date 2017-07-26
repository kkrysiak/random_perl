## Modification of /gscmnt/gc2547/griffithlab/zskidmor/analysis_projects/HCC30/bam2xml.pl
## Originally written by Lee Trani
#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;

my ( $file,        $bed_file, $bed_dir, $help );
my ( %sample_hash, %out_hash, %id_sample );
my ( @samples,     @bed_dirs );

#file format:sample, sample, tissue, bam. bed flag is location of bed directories
GetOptions(
    "file:s" => \$file,
    "bed:s"  => \$bed_file,
    "help"   => \$help
) or die( print "error in command line arguments\n" );

help() if ($help);

open( FILE, $file ) || die $!;

#create hash from file
while (<FILE>) {
    chomp;
    my @data = split("\t");
    $sample_hash{ $data[0] }{ $data[2] } = $data[3];
}

my $count = 0;

#create sample track hash for xml file, add xml header
foreach my $k ( keys %sample_hash ) {
    chomp($k);
    my $out_file = "$k" . '.xml';

    if ($bed_file) {
        $bed_dir = $bed_file . "$k" . '.bed';
        print "$bed_dir\n";
        my $bed_result = bed_track( $bed_dir, $k );
        push( @bed_dirs, "$k\t$bed_result" );
    }

    push( @samples, $k );
    open( OUT, '>', $out_file );
    print OUT q(<!--<?xml version="1.0" encoding="UTF-8" standalone="no"?>-->
<Session genome="b37" hasGeneTrack="true" hasSequenceTrack="true" locus="All" version="8">
    <Resources>
);
    ## Check if bam files provided start with https for URL access
    my $url_check=0;
    my $url_pre = "https://gscweb.gsc.wustl.edu";

#add resource section to xml file
    foreach my $y ( keys %{ $sample_hash{$k} } ) {
        my $res = resource( $sample_hash{$k}{$y} );
        print OUT "\t\t$res\n";
        $count++;
        if($res =~ /http/) {
            $url_check=1;
        }
        my $id_samp_out = id_sample( $y, $k, $sample_hash{$k}{$y}, $count );
        $id_sample{$k}{$y} = $id_samp_out;
    }
    if ($bed_dir) {
        my $bed_res = resource($bed_dir);
        if($url_check=1) {
            $bed_res =~ s/(Resource path=\")/$1$url_pre/;
        }
        print OUT "\t\t$bed_res\n";
    }

    print OUT "\t</Resources>";
    close OUT;
}

#print track xml code to corresponding sample file
foreach my $panel_out ( keys %id_sample ) {

    my $out_file = "$panel_out" . '.xml';
    open( OUT, '>>', $out_file );
    foreach my $k ( keys %{ $id_sample{$panel_out} } ) {
        foreach my $sample (@samples) {
            if ( $panel_out eq $sample ) {
                print OUT "$id_sample{$panel_out}{$k}\n\t</Panel>";

            }
        }
    }
}

my $bed_pass = scalar(@bed_dirs);

#if there is a bed file add and close xml file, else just close xml file. 
if ( $bed_pass > 0 ) {
    foreach my $bed_print (@bed_dirs) {
        chomp($bed_print);
        my @data = split( "\t", $bed_print );

        my $out_file = "$data[0]" . '.xml';
        open( OUT, '>>', $out_file );
        print OUT "\n\t$data[1]";
        print OUT
qq(\n\t<PanelLayout dividerFractions="0.006550218340611353,0.22489082969432314,0.4497816593886463,0.6746724890829694,0.8995633187772926"/>
    <HiddenAttributes>
        <Attribute name="NAME"/>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
    </HiddenAttributes>
</Session>);

        close OUT;
    }
}
elsif ( $bed_pass == 0 ) {
    foreach my $sa (@samples) {
        my $out = "$sa" . '.xml';
        open( OUT, '>>', $out );
        print OUT
qq(\n\t<PanelLayout dividerFractions="0.006550218340611353,0.22489082969432314,0.4497816593886463,0.6746724890829694,0.8995633187772926"/>
    <HiddenAttributes>
        <Attribute name="NAME"/>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
    </HiddenAttributes>
</Session>);

        close OUT;
    }
}

close FILE;
exit;

###########################

#create resource pathways at the beginning of the xml file
sub resource {
    my $bam      = shift;
    my $resource = qq(<Resource path="$bam\"/>);
    return $resource;
}

#create bam tracks for igv
sub id_sample {
    my ( $id, $sample, $bam, $count ) = @_;
    my ($bam_ident) = $bam =~ /\w*.bam/g;
    my $print1 = qq(
    <Panel height="195" name="Panel$count" width="1244">
            <Track altColor="0,0,178" autoScale="true" color="175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" );
    my $print2 = "id=\"$bam\_coverage\" name=\"$sample\_coverage\"";
    my $print3 =
qq( showReference="false" snpThreshold="0.2" sortable="true" visible="true">
                 <DataRange baseline="0.0" drawBaseline="false" flipAxis="false" maximum="160.0" minimum="0.0" type="LINEAR"/>
             </Track>
             <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" );
    my $print4 = "id=\"$bam\" name=\"$sample     $id     $bam_ident\"";
    my $print5 = qq( showSpliceJunctions="false" sortable="true" visible="true">
                <RenderOptions colorByTag="" colorOption="READ_STRAND" flagUnmappedPairs="false" groupByTag="" maxInsertSize="1000" minInsertSize="50" shadeBasesOption="QUALITY" shadeCenters="true" showAllBases="false" sortByTag=""/>
            </Track>);
    my $print_final = "$print1 $print2" . "$print3 $print4 $print5";
    return $print_final;

}

#create bed track for igv
sub bed_track {
    my ( $bed_track, $sample ) = @_;
    my $track1 = qq(<Panel height="87" name="FeaturePanel" width="1244">
            <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;4.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" );
    my $track2 = "id=\"$bed_track\" name=\"$sample.bed\"";
    my $track3 =
qq( renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="4.0" minimum="0.0" type="LINEAR"/>
            </Track>
        </Panel>);

    my $final_track = "$track1 $track2 $track3";
    return $final_track;
}

sub help {
    print"
    usage: perl bam2xml.pl --file=HCC_bams.txt --bed=/gscmnt/gc2547/mardiswilsonlab/zskidmor/analysis_projects/HCC30/manual_review/
    --file=  column1     column2          column3     column4
             sample name full sample name wgs/capture bam_path

      H_MU_748892	H_MU_748892_1209081	wgs30_normal	/gscmnt/gc2000/info/build_merged_alignments/merged-alignment-blade14-2-15.gsc.wustl.edu-apipe-builder-24809-129159340/129159340.bam
      H_MU_753980	H_MU_753980_1209084	wgs30_normal	/gscmnt/gc7001/info/build_merged_alignments/merged-alignment-blade12-3-16.gsc.wustl.edu-mcordes-369-129345661/129345661.bam
      H_MU_5222	H_MU_5222_1209012	wgs30_normal	/gscmnt/gc7001/info/build_merged_alignments/merged-alignment-blade11-3-9.gsc.wustl.edu-apipe-builder-22433-135838105/135838105.bam
  
    --bed= optional, directory path to directories named after each sample that contain the bed file (sample.bed)\n             
";
    exit;

}
