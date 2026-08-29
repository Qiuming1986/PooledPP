#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError);  # Import Gzip compression

# Subroutine to show usage/help information
sub usage {
    print <<END;
Usage: perl $0 --sv_file <SV signal file> --ref_genome <reference genome file> --sv_distance <distance between SV signals> --output_prefix <output file prefix>

Description:
This script generates pseudo-genomes based on the clustering and deduplication of SV signals from HiFi data. The input includes an SV signal file, a reference genome, and a specified distance between SV signals for clustering.

Options:
    --sv_file <file>       Input file with structural variant signals (required). Supports compressed (.gz) or uncompressed files.
    --ref_genome <file>    Reference genome file in FASTA format (required).
    --sv_distance <int>    Distance between SV signals for clustering (default: 250).
    --output_prefix <str>  Prefix for the output pseudo-genome files (required). The output files will be named <prefix>.Pseudogenome1.fa, <prefix>.Pseudogenome2.fa, etc.
    --help                 Show this help message.

Examples:
    perl $0 --sv_file sv_signals.txt --ref_genome reference.fa --sv_distance 250 --output_prefix pseudo

END
    exit;
}

# Parse command line options
my ($sv_file, $ref_genome, $sv_distance, $output_prefix, $help);

# Default value for sv_distance
$sv_distance = 250;

GetOptions(
    "sv_file=s"       => \$sv_file,
    "ref_genome=s"    => \$ref_genome,
    "sv_distance=i"   => \$sv_distance,
    "output_prefix=s" => \$output_prefix,
    "help"            => \$help,
) or usage();

# Show usage if help is requested or required parameters are missing
if ($help || !$sv_file || !$ref_genome || !$output_prefix) {
    usage();
}

# Function to extract sequence using samtools faidx
sub get_sequence {
    my ($chrom, $start, $end, $ref_genome) = @_;
    my $seq = `samtools faidx $ref_genome $chrom:$start-$end`;
    $seq =~ s/^>.*\n//;  # 移除FASTA的头信息
    $seq =~ s/\n//g;     # 移除换行符
    return $seq;
}

# Step to capture chromosome order from the reference genome .fai file
sub get_chromosome_order {
    my ($ref_genome) = @_;
    my @chromosome_order;
    my $fai_file = "$ref_genome.fai";
    open my $fai_fh, '<', $fai_file or die "Could not open $fai_file: $!";
    while (<$fai_fh>) {
        my ($chrom, $length) = split(/\t/);
        push @chromosome_order, $chrom;
    }
    close $fai_fh;
    return @chromosome_order;
}

# Open SV signal file (supports both compressed and uncompressed files)
my $in_fh;
if ($sv_file =~ /\.gz$/) {
    open($in_fh, "gzip -dc $sv_file |") or die "Could not open compressed input file: $sv_file\n";
} else {
    open($in_fh, '<', $sv_file) or die "Could not open input file: $sv_file\n";
}

# Skip header line
<$in_fh>;

# Get chromosome order from the reference genome
my @chromosome_order = get_chromosome_order($ref_genome);

# Initialize data structures to store SV info for each pseudo-genome and chromosome
my %sv_info;  # Stores SV signals for each pseudo-genome
my %genome_index_per_chr;  # Track how many pseudo-genomes per chromosome

# Process SV signals and store the coordinates and sequences
while (<$in_fh>) {
    chomp;
    next if (/^Chromosome/);  # Skip header line if encountered again

    my ($chrom, $seed_start, $seed_end, $ref_start, $ref_end, $sv_type, $support, $read_name, $seed_seq) = split(/\t/);

    # Initialize chromosome-specific genome index if not already set
    if (!exists $genome_index_per_chr{$chrom}) {
        $genome_index_per_chr{$chrom} = 0;
    }

    # Try to find a pseudo-genome where the current SV can fit
    my $found_genome = 0;
    for (my $i = 0; $i <= $genome_index_per_chr{$chrom}; $i++) {
        if (!$sv_info{$i}{$chrom}) {
            $sv_info{$i}{$chrom} = [];
        }

        # Check the last end position for the current chromosome in the current pseudo-genome
        my $last_sv = $sv_info{$i}{$chrom}->[-1];  # Get the last inserted SV
        my $last_end = $last_sv ? $last_sv->{end} : 0;

        if ($seed_start - $last_end >= $sv_distance) {
            # Add the SV to the pseudo-genome for this chromosome
            push @{$sv_info{$i}{$chrom}}, {start => $seed_start, end => $seed_end, seq => $seed_seq};
            $found_genome = 1;
            last;
        }
    }

    # If no suitable genome was found, create a new pseudo-genome entry for this chromosome
    if (!$found_genome) {
        my $new_index = ++$genome_index_per_chr{$chrom};
        $sv_info{$new_index}{$chrom} = [{start => $seed_start, end => $seed_end, seq => $seed_seq}];
    }
}

# Now generate the pseudo-genomes based on the stored SV info and output the corresponding SVs
my $counter = 1;  # Initialize the pseudo-genome counter
for my $i (sort { $a <=> $b } keys %sv_info) {  # Ensure keys are sorted numerically (by pseudo-genome ID)
    open my $fh, ">", "$output_prefix.Pseudogenome" . ($counter) . ".fa" or die "Could not open output file: $!";
    open my $sv_fh, ">", "$output_prefix.SV_info" . ($counter) . ".txt" or die "Could not open SV output file: $!";

    # Write chromosomes in the order provided by the reference genome
    for my $chrom (@chromosome_order) {
        next unless exists $sv_info{$i}{$chrom};  # Skip if the chromosome has no SVs in this pseudo-genome

        print $fh ">$chrom\n";

        # Use samtools to extract from the reference genome
        my $chr_length = `samtools faidx $ref_genome $chrom | grep -v ">" | wc -c`;  # Get chromosome length
        my $result_seq = "";

        # Keep track of the last position processed
        my $last_pos = 1;  # samtools faidx uses 1-based indexing

        # Process each SV in order
        for my $sv (@{$sv_info{$i}{$chrom}}) {
            # Add the reference sequence up to the current SV start
            $result_seq .= get_sequence($chrom, $last_pos, $sv->{start} - 1, $ref_genome);
            # Add the SV sequence
            $result_seq .= $sv->{seq};
            # Update the last processed position
            $last_pos = $sv->{end} + 1;  # Adjust for 1-based indexing

            # Output corresponding SV information
            print $sv_fh "$chrom\t$sv->{start}\t$sv->{end}\t$sv->{seq}\n";
        }

        # Add the remaining reference sequence after the last SV
        $result_seq .= get_sequence($chrom, $last_pos, $chr_length, $ref_genome);

        # Format the result sequence in 60-character lines
        my $formatted_seq = $result_seq =~ s/(.{1,60})/$1\n/gr;
        print $fh $formatted_seq;
    }

    close $fh;
    close $sv_fh;

    # Compress the SV signal file
    my $gzip_out = "$output_prefix.SV_info" . ($counter) . ".txt.gz";
    gzip "$output_prefix.SV_info" . ($counter) . ".txt" => $gzip_out
        or die "gzip failed: $GzipError\n";

    # Remove the original uncompressed SV file
    unlink "$output_prefix.SV_info" . ($counter) . ".txt";

    $counter++;  # Increment the pseudo-genome counter
}

# Generate the minigraph script
open my $script_fh, '>', "$output_prefix.minigraph.sh" or die "Could not open script file: $!";

# Print the base command for Minigraph
print $script_fh "minigraph -cxggs -t20 $ref_genome ";

# Loop over the pseudo-genomes in ascending order (based on the counter)
for (my $i = $counter - 1; $i >= 1; $i--) {
    print $script_fh "$output_prefix.Pseudogenome$i.fa ";
}

# Output the GFA file
print $script_fh "> $output_prefix.gfa\n";

# Close the script file
close $script_fh;

# Close input file
close $in_fh;

print "Pseudogenomes, SV signals, and Minigraph script have been successfully generated.\n";
