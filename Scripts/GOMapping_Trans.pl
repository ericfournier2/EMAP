#! /usr/bin/perl

# Creates a mapping between GO term IDs and probe names, and prints it to
# standard output.


sub readGO {
    my $results = {};
    
    open(INPUT, "<", $_[0]);
    while(my $line=<INPUT>) {
        if(substr($line, 0, 1) ne "!") {
            my @entries = split("\t", $line);
            
            if(!exists($results->{$entries[$_[1]]})) {
                $results->{$entries[$_[1]]} = {};
            }
            $results->{$entries[$_[1]]}->{$entries[4]} = "";
        }
    }
    
    return $results;
}

# Load libraries.
use List::MoreUtils qw(uniq);
use ReadTable;

my $annotationFile = $ARGV[0];
# Cow epigenetic v1: "Annotations/v1/EDMA.Annotation.txt"
# Cow epigenetic v2: "Annotations/v2/EDMA.Annotation.txt"
# Pig epigenetic v1: "Annotations/pigv1/EDMA.Annotation.txt"
# Pig epigenetic v2: "Annotations/pigv2/EDMA.Annotation.txt"
# Cow transcriptomic: "Annotations/cow/Trans.annotation.txt"
# Pig transcriptomic: "Annotations/pig/Trans.annotation.txt"
my $annotationType = $ARGV[1];
if(($annotationType ne "Transcriptomic") && ($annotationType ne "Epigenetic")) {
    die "Invalid annotation type.";
}

my $goFile = $ARGV[2];
# Cow: "Annotations/cow/gene_association.goa"
# Pig: "Annotations/pig/gene_association.goa"

my $symbolMapFile = "";
if(scalar(@ARGV)>3) {
    $symbolMapFile = $ARGV[3];
}
# Cow: none
# Pig: "Annotations/RefSeq2Gene_Pig.txt"

# Load GO annotation.
# Source: http://geneontology.org/gene-associations/gene_association.goa_cow.gz
# http://geneontology.org/page/download-annotations for a list of all available annotations.
my @GOColNames = ("DB", "DBObjectID", "Symbol", "Qualifier", "GOID", "GOReference", "EvidenceCode", "WithOrFrom", "Aspect", "GeneName", "Synonym", "ObjectType", "Taxon", "Date", "AssignedBy", "Extension", "GeneProductFormID");
my $cowGOAssociations = readGO($goFile, 2);

# Read GO association files, but key them by database (Uniprot, SGD, MGI) IDs.
my $cowGOAssociationsDBID = readGO($goFile, 1);

# Load EDMA annotations.
my $annotations = readTable($annotationFile, "\t", "Probe");

my $symbolMap = {};
if($symbolMapFile ne "") {
    $symbolMap = readTable($symbolMapFile, "\t", "ID");
}

# Loop over all EMBV3 probes.
foreach my $probe (keys %{$annotations}) {
    print $probe . "\t";
    my @uniqGOIDs = ();
    my @egGOIDs = ();    
    my @symbols = ();
    
    if($annotationType eq "Epigenetic") {
        if(substr($probe, 0, 8) eq "EDMA_MET" ||
           substr($probe, 0, 5) eq "GT_HQ" || substr($probe, 0, 5) eq "GT_LQ" ||
           substr($probe, 0, 6) eq "GT_pig" || substr($probe, 0, 9) eq "GT_MET_ch" || substr($probe, 0, 9) eq "GT_MET_GL"
           ) {
            # Get a list of non-redundant gene symbols for this probe.
            my @symbolsAll = ( split(" ", $annotations->{$probe}->{"Proximal_Promoter"}),
                               split(" ", $annotations->{$probe}->{"Promoter"}), 
                               split(" ", $annotations->{$probe}->{"Exon"}),
                               split(" ", $annotations->{$probe}->{"Intron"}));
                               
            for(my $i=0; $i<scalar(@symbolsAll); ++$i) {
                $symbolsAll[$i] =~ s/-\d+//g;
            }
            
            @symbols = uniq(@symbolsAll);
            if($symbolMapFile ne "") {
                for(my $i=0; $i<scalar(@symbols); ++$i) {
                    if(exists($symbolMap->{$symbols[$i]})) {
                        $symbols[$i] = $symbolMap->{$symbols[$i]}->{"Symbol"};
                    } else {
                        print STDERR "No mapping for " . $symbols[$i] . "\n";
                    }
                }
            }
        }
    } else {
        @symbols = ($annotations->{$probe}->{"Gene_Symbol"})
    }
        
    # Loop over all symbols, and accumulate the GO IDs associated with them.
    my @GOIDs = ();
    foreach my $symbol (@symbols) {
        if(exists($cowGOAssociations->{$symbol})) {
            foreach my $GOTerm (keys %{$cowGOAssociations->{$symbol}}) {
                push(@GOIDs, $GOTerm);
            }
        }
    }
    @uniqGOIDs = uniq(@GOIDs);
    
    # Print out all unique symbols.
    print join(", ", @uniqGOIDs) . "\n";
}
