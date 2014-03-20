#!/bin/perl

# Takes the results of the limma analysis of an EDMA experiment, and
# generates a file ready for import within IPA. The gene -> p-value
# associations are made using a best-case scenario, IE, if a gene is
# associated with 10 probes within his exons, introns and promoter, then
# we associate that gene with the probe with the lowest p-value in the set.

use ReadTable;

# Read in annotation data and limma results.
my $annotationPath = $ARGV[0];
my $limmaPath = $ARGV[1];
my $annotation = readTable($annotationPath, "\t", "Probe");
my $limma = readTable($limmaPath, "\t", "Probe");

# Loop over all annotated probes.
my $results = {};
foreach my $probe (keys %{$annotation}) {
    # Make a list of all associated genes. Discount the "Distal promoter" category.
    my $allGeneString = $annotation->{$probe}->{"Exon"} . " " .
                        $annotation->{$probe}->{"Intron"} . " " .
                        $annotation->{$probe}->{"Proximal_Promoter"} . " " .
                        $annotation->{$probe}->{"Promoter"};
    
    # For all identified genes, check if the current probe's p-value is better than the p-value
    # currently associated with this gene. If it is, assign this probe's fold-change/p-value pair
    # to this gene.
    $allGeneString =~ s/-\d+\s/ /g;
    my @geneList = split(" ", $allGeneString);
    foreach my $gene (@geneList) {
        if(!exists($results->{$gene}) || ($results->{$gene}->{"P-value"} > $limma->{$probe}->{"P-value"})) {
            $results->{$gene} = $limma->{$probe};
            $results->{$gene}->{"Probe"} = $probe;
        }
    }
}

# Output the list of all surveyed genes and their associated "best" p-value.
foreach my $gene (keys %{$results}) {
    print $gene . "\t" . $results->{$gene}->{"Probe"} . "\t" . $results->{$gene}->{"Fold-change"} . "\t" . $results->{$gene}->{"P-value"} . "\n";
}
