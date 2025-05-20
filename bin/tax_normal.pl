#!/usr/bin/perl -w
use warnings;
use strict;

die "usage: perl $0 otu_rep.fasta otu_rep.temp.tax otu_rep.tax\n" unless(@ARGV == 3);
my ($rep_seq, $taxFile, $output) = @ARGV;

my %hash;

open INT, "$taxFile" or die "$!\n";
while (<INT>) {
    chomp;
    next if (/^name/ || /^\#/); 
    my @temp = split /\t/;
    $hash{$temp[0]} = $temp[1];  
}
close INT;

open INF, "$rep_seq" or die "$!\n";
open OUT, "> $output" or die "$!\n";
while (<INF>) {
    chomp;
    if (/^>/) {
        my @temp = split;
        $temp[0] =~ s/^>//;
        if (exists $hash{$temp[0]}) {
            print OUT "$temp[0]\t$hash{$temp[0]}\n";
        } else {
            print OUT "$temp[0]\tUnclassified\n";
        }
    }
}
close INF;
close OUT;

