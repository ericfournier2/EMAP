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

my $version = $ARGV[0];

# Load GO annotation.
# Source: http://geneontology.org/gene-associations/gene_association.goa_cow.gz
# http://geneontology.org/page/download-annotations for a list of all available annotations.
my @GOColNames = ("DB", "DBObjectID", "Symbol", "Qualifier", "GOID", "GOReference", "EvidenceCode", "WithOrFrom", "Aspect", "GeneName", "Synonym", "ObjectType", "Taxon", "Date", "AssignedBy", "Extension", "GeneProductFormID");
my $cowGOAssociations = readGO("Annotations/gene_association.goa_cow", 2);

# Read GO association files, but key them by database (Uniprot, SGD, MGI) IDs.
my $cowGOAssociationsDBID = readGO("Annotations/gene_association.goa_cow", 1);

# Load EDMA annotations.
my $edmaAnnotations = readTable("Annotations/" . $version . "/EDMA.Annotation.txt", "\t", "Probe");

# Loop over all EMBV3 probes.
foreach my $probe (keys %{$edmaAnnotations}) {
    print $probe . "\t";
    
    my @uniqGOIDs = ();
    my @egGOIDs = ();    
    if(substr($probe, 0, 8) eq "EDMA_MET") {
        # Get a list of non-redundant gene symbols for this probe.
        my @symbolsAll = ( split(" ", $edmaAnnotations->{$probe}->{"Proximal_Promoter"}),
                           split(" ", $edmaAnnotations->{$probe}->{"Promoter"}), 
                           split(" ", $edmaAnnotations->{$probe}->{"Exon"}),
                           split(" ", $edmaAnnotations->{$probe}->{"Intron"}));
                           
        for(my $i=0; $i<scalar(@symbolsAll); ++$i) {
            $symbolsAll[$i] =~ s/-\d+//g;
        }
        
        my @symbols = uniq(@symbolsAll);
        
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
    }
    
    # Print out all unique symbols.
    print join(", ", @uniqGOIDs) . "\n";
}
