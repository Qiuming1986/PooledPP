#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Variables to store command line options
my ($infile, $outfile, $vg, $help);

# Get command line options
GetOptions(
    'input|i=s'  => \$infile,    # Input GFA file
    'output|o=s' => \$outfile,   # Output file prefix
    'vg|v=s'     => \$vg,        # Path to vg software
    'help|h'     => \$help       # Help option
) or die "Error in command line arguments\n";

# Display help message if requested or if required options are missing
if ($help) {
    print_help();
    exit(0);
}
if (!$infile || !$outfile || !$vg) {
    print "Error: Input, output prefix, and vg software path are required\n";
    print_help();
    exit(1);
}

# Remove '.vcf' suffix from the output prefix if it's mistakenly included
my $prefix = $outfile;
$prefix =~ s/\.vcf$//;

# Run vg commands: 
# 'vg view' converts the GFA file to VG format
system("$vg view -F -g $infile > $prefix.vg") == 0 or die "vg view command failed: $!";
# 'vg deconstruct' extracts variants from the VG file and writes to a VCF file
system("$vg deconstruct $prefix.vg -a > $prefix.vcf") == 0 or die "vg deconstruct command failed: $!";

# Open the generated VCF file and prepare to modify it
open(my $in, '<', "$prefix.vcf") or die "Cannot open intermediate VCF file: $!";
open(my $out, '>', "${prefix}.addGT.vcf") or die "Cannot open output VCF file: $!";

# Process VCF file to modify header and add GT field
while(<$in>){
    chomp;
    my @F = split(/\t/, $_);
    
    # Add GT field to the VCF header
    if (/^##INFO=\<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\"\>/) {
        print $out "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">\n";
        print $out "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    }
    # Copy other header lines as they are
    elsif (/^##/) {
        print $out "$_\n";
    }
    # Modify the #CHROM line to add sample name
    elsif (/^#C/) {
        print $out "$_\tsampleName\n";
    }
    # Add GT 1/1 to each data line
    else {
        print $out "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\tGT\t1/1\n";
    }
}

# Close file handles after processing the VCF file
close($in);
close($out);

# Compress the VCF file using bgzip
system("bgzip ${prefix}.addGT.vcf") == 0 or die "bgzip failed: $!";
# Index the compressed VCF file using tabix
system("tabix -p vcf ${prefix}.addGT.vcf.gz") == 0 or die "tabix failed: $!";
# Normalize multi-allelic records into bi-allelic records using bcftools
system("bcftools norm -m -both -o $outfile.split.vcf.gz -O z ${prefix}.addGT.vcf.gz") == 0 or die "bcftools norm failed: $!";

### Additional processing after bcftools normalization to modify genotype information ###

# Uncompress the split VCF file created by bcftools
system("gunzip $outfile.split.vcf.gz");

# Open the split VCF file and prepare to modify genotype (GT) field
open(my $in1, '<', "$outfile.split.vcf") or die "Cannot open intermediate VCF file: $!";
open(my $out1, '>', "$outfile.split.tmp.vcf") or die "Cannot open output VCF file: $!";

# Loop through the lines in the split VCF file
while(<$in1>){
    chomp;
    my @F = split(/\t/, $_);    
    # Keep all header lines unchanged
    if (/^#/) {
        print $out1 "$_\n";
    }
    # Modify the genotype field in each data line to '1/1'
    else {
        print $out1 "$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\tGT\t1/1\n";
    }
}

# Close file handles for split VCF
close($in1);
close($out1);

# Replace the original split VCF file with the modified one
system("mv $outfile.split.tmp.vcf $outfile.split.vcf");
# Compress the modified split VCF file
system("bgzip $outfile.split.vcf");
# Index the newly compressed split VCF file
system("tabix -p vcf $outfile.split.vcf.gz") == 0 or die "tabix failed: $!";

# Help message function
sub print_help {
    print <<EOF;
Usage: vcf_processing.pl [options]

Options:
  --input | -i     Input GFA file (required)
  --output | -o    Output file prefix (required)
  --vg | -v        Path to vg software (required)
  --help | -h      Display this help message

Description:
  This script processes a GFA file by first converting it to VG format using 'vg view',
  extracting variants using 'vg deconstruct', adding a GT (genotype) field to the header
  and each variant, and then compressing the VCF file using bgzip, creating an index using
  tabix, and normalizing multi-allelic records using bcftools norm.

  After normalization, the script adjusts the genotype (GT) field of all records to '1/1',
  compresses and indexes the resulting VCF file again.

  The output VCF files will have '.addGT.vcf.gz' and '.split.vcf.gz' as suffixes.
  For example, if the '--output' option is set to 'output', the final files will be:
  'output.addGT.vcf.gz' and 'output.split.vcf.gz'.
EOF
}

