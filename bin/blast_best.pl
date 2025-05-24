use strict;

die "Usage: perl $0 <blast_output_file>\n" unless @ARGV == 1;

# 打开文件
open my $fh, '<', $ARGV[0] or die "Cannot open file: $!\n";

my %seen_queries;  # Record processed query sequences

while (<$fh>) {
    chomp;
    my @fields = split;  

    my $query = $fields[0];  
    
    if (!$seen_queries{$query}) {
        print "$_\n";  # Keep the top hit
        $seen_queries{$query} = 1;  
    }
}

close $fh;

