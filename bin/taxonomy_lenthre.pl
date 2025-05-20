use strict;

die "Usage: perl $0 <blast_output_file> <tax_file> <output_file>\n" unless @ARGV == 3;

open my $blast_fh, '<', $ARGV[0] or die "Cannot open BLAST output file: $!\n";

open my $tax_fh, '<', $ARGV[1] or die "Cannot open tax file: $!\n";

open my $out_fh, '>', $ARGV[2] or die "Cannot open output file: $!\n";

# Read taxnomy file
my %tax_info;
while (<$tax_fh>) {
    chomp;
    my @fields = split;  
    my $acc = $fields[0]; 
    my $taxonomy = $fields[1];  
    $tax_info{$acc} = $taxonomy;  
}
close $tax_fh;

# Process blast out file
while (<$blast_fh>) {
    chomp;
    my @fields = split;
    my $query = $fields[0];
    my $acc = $fields[1];
    my $score = $fields[3];

    # Coverage length >100 as restrict 
    if ($score >= 100) {
        if (exists $tax_info{$acc}) {
            print $out_fh "$query\t$tax_info{$acc}\n";
        }
    } else {
        print $out_fh "$query\tunclassified\n";
    }
}

close $blast_fh;
close $out_fh;

