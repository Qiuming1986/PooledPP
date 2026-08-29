#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError);  # Module to handle gzip compression
use Set::IntervalTree;
use File::Temp qw(tempfile);

$| = 1;  # Disable output buffering

# Declare variables for command-line options and set default values
my ($bam_file, $bamSortRead_file, $output_file, $min_mapq, $min_align_len, $min_sv_len, $max_sv_len, $mitochondria_chr,
    $anchor_len_insertion, $anchor_len_deletion, $anchor_len_CGR, $anchor_len_duplication, $anchor_len_complexInsertion, $anchor_len_complexDeletion, $anchor_len_SCGR,
    $anchor_minDepth, $anchor_maxDepth);

# Set default values for thresholds and parameters
$min_mapq = 20;                               # Minimum mapping quality (MAPQ)
$min_align_len = 1000;                        # Minimum alignment length (bp)
$min_sv_len = 50;                             # Minimum length for structural variants (SVs) like insertions/deletions
$max_sv_len = 10000;                          # Maximum length for SVs (bp)
$anchor_len_insertion = 230;                   # anchor length for insertion
$anchor_len_deletion = 290;                    # anchor length for deletion
$anchor_len_CGR = 260;                        # anchor length for CGR
$anchor_len_duplication  = 210;              # anchor length for duplication
$anchor_len_complexInsertion = 230;           # anchor length for complex insertion
$anchor_len_complexDeletion = 290;            # anchor length for complex deletion
$anchor_len_SCGR = 340;                       # anchor length for SCGR
$anchor_minDepth = 0.33;                      # minimum anchor depth ratio to overall mean
$anchor_maxDepth = 3.00;                      # Maximum anchor depth ratio to overall mean

# Parse command-line options
Getopt::Long::GetOptions(
    "bam=s" => \$bam_file,
    "bamSortRead=s" => \$bamSortRead_file,
    "output=s" => \$output_file,
    "min-mapq=i" => \$min_mapq,
    "min-align-len=i" => \$min_align_len,
    "min-sv-len=i" => \$min_sv_len,
    "max-sv-len=i" => \$max_sv_len,
    "mitochondria-chr=s" => \$mitochondria_chr,  # Option to skip mitochondrial reads
    "anchor-len-insertion=i" => \$anchor_len_insertion,        # Option to set anchor length
    "anchor-len-deletion=i" => \$anchor_len_deletion,
    "anchor-len-CGR=i" => \$anchor_len_CGR,
    "anchor-len-duplication=i" => \$anchor_len_duplication,
    "anchor-len-complexInsertion=i" => \$anchor_len_complexInsertion,
    "anchor-len-complexDeletion=i" => \$anchor_len_complexDeletion,
    "anchor-len-SCGR=i" => \$anchor_len_SCGR,
    "anchor-minDepth=f" => \$anchor_minDepth,        # Option to set anchor depth ratio to overall mean
    "anchor-maxDepth=f" => \$anchor_maxDepth,
    "help" => \&usage
) or die "Error in command-line arguments\n";

# Check for required arguments
if (!$bam_file || !$bamSortRead_file || !$output_file) {
    usage();
    exit;
}

# Open the input BAM file using samtools
open(IN, "-|", "samtools view $bamSortRead_file") or die "Error: Unable to open BAM file: $bam_file\n";

# Initialize a hash to store alignments grouped by read name
my %read_alignments;

my $step_size = 100;  # Step size is half of the window size

# Declare global variables to store depth information
my %depth_data;
my $total_depth = 0;
my $position_count = 0;

# Calculate average depth from BAM file
open(DEPTH1, "-|", "samtools depth -Q $min_mapq $bam_file") or die "Cannot open BAM file: $bam_file\n";
while (<DEPTH1>) {
    chomp;
    my ($chromosome, $position, $depth) = split /\t/;
    $total_depth += $depth;
    $position_count++;
}
my $mean_depth = $position_count ? ($total_depth / $position_count) : 0;
close(DEPTH1);

# Read BAM file and apply sliding window with multiple offsets
open(DEPTH2, "-|", "samtools depth -Q $min_mapq $bam_file") or die "Cannot open BAM file: $bam_file\n";

# Initialize sliding window data structure for multiple offsets
my %window_data;

# Traverse depth data to ensure each position is included in two adjacent windows
while (<DEPTH2>) {
    chomp;
    my ($chromosome, $position, $depth) = split /\t/;

    # Calculate start of the current window
    my $window_start = int(($position - 1) / $step_size) * $step_size + 1;
    $window_data{$chromosome}{$window_start}{sum} += $depth;
    $window_data{$chromosome}{$window_start}{count}++;

    # Calculate start of the previous window
    my $previous_window_start = $window_start - $step_size;
    if ($previous_window_start > 0) {
        $window_data{$chromosome}{$previous_window_start}{sum} += $depth;
        $window_data{$chromosome}{$previous_window_start}{count}++;
    }
}
close(DEPTH2);

print "Sequencing Depth processed complete!\n";


# Initialize variables for  graded intervals
my %graded_intervals;

# Analyze windows and mark abnormal intervals based on depth
foreach my $chromosome (keys %window_data) {
    foreach my $window_start (sort { $a <=> $b } keys %{$window_data{$chromosome}}) {
        my $window_end = $window_start + 2 * $step_size - 1;
        my $sum = $window_data{$chromosome}{$window_start}{sum};
        my $count = $window_data{$chromosome}{$window_start}{count};

        # Skip if the count is below 80% of the window size (missing data threshold)
        if ($count < 120) {
            next;
        }

        my $window_mean_depth = $count ? ($sum / $count) : 0;

        # Define grading based on fold-change ranges
        my $grade;
        if ($window_mean_depth >= $mean_depth * 0.8 && $window_mean_depth <= $mean_depth * 1.25) {
            $grade = "1st";
        } elsif (($window_mean_depth >= $mean_depth * 0.5 && $window_mean_depth < $mean_depth * 0.8) ||
                 ($window_mean_depth > $mean_depth * 1.25 && $window_mean_depth <= $mean_depth * 2)) {
            $grade = "2nd";
        } elsif ($window_mean_depth >= $mean_depth * $anchor_minDepth && $window_mean_depth < $mean_depth * 0.5 ||
                 $window_mean_depth > $mean_depth * 2 && $window_mean_depth <= $mean_depth * $anchor_maxDepth) {
            $grade = "3rd";
        }

        # Store graded windows if a grade was assigned
        if ($grade) {
            push @{$graded_intervals{$chromosome}{$grade}}, { start => $window_start, end => $window_end };
        }
    }
    # Remove processed chromosome data from %window_data to save memory
    delete $window_data{$chromosome};
}

print "Grading window identification complete!\n";

# Combine %graded_intervals for overlap removal
my (%first, %second, %third);
foreach my $chromosome (keys %graded_intervals) {

    # Process each grade's intervals
    foreach my $grade (keys %{$graded_intervals{$chromosome}}) {
        my ($fh_grade, $grade_filename) = tempfile(SUFFIX => ".bed");
        foreach my $interval (@{$graded_intervals{$chromosome}{$grade}}) {
            print $fh_grade join("\t", $chromosome, $interval->{start} - 1, $interval->{end}), "\n";
        }
        close($fh_grade);

        # Ensure sorted input and run bedtools merge
        my $merged_file = "$grade_filename.merged";
        system("sort -k1,1 -k2,2n $grade_filename | bedtools merge -i - > $merged_file") == 0
            or die "Error running bedtools merge: $!";

        # Read filtered results and store
        open my $in, '<', $merged_file or die "Cannot open merged file: $!";
        my @filtered_intervals;
        while (<$in>) {
            chomp;
            my ($chrom, $start, $end) = split("\t", $_);
            push @filtered_intervals, { start => $start + 1, end => $end };  # Convert to 1-based
        }
        close($in);

        # Store in respective hash based on grade
        if ($grade eq "1st") {
            $first{$chromosome} = \@filtered_intervals;
        } elsif ($grade eq "2nd") {
            $second{$chromosome} = \@filtered_intervals;
        } elsif ($grade eq "3rd") {
            $third{$chromosome} = \@filtered_intervals;
        }

        # Clean up temporary files
        unlink $grade_filename, $merged_file, $output_file;
    }
}

print "Graded intervals processing complete with overlap removal!\n";

#Log completion and build interval trees
my (%first_tree, %second_tree, %third_tree);
foreach my $chromosome (keys %first) {
    $first_tree{$chromosome} = Set::IntervalTree->new;
    $first_tree{$chromosome}->insert($_, $_->{start}, $_->{end}) for @{$first{$chromosome}};
}
foreach my $chromosome (keys %second) {
    $second_tree{$chromosome} = Set::IntervalTree->new;
    $second_tree{$chromosome}->insert($_, $_->{start}, $_->{end}) for @{$second{$chromosome}};
}
foreach my $chromosome (keys %third) {
    $third_tree{$chromosome} = Set::IntervalTree->new;
    $third_tree{$chromosome}->insert($_, $_->{start}, $_->{end}) for @{$third{$chromosome}};
}

print "Interval tree construction complete!\n";

# Clear entire hashes after usage to save memory
%window_data = ();
%first = ();
%second = ();
%third = ();

# Write the sorted output to a gzipped file
my $gz = IO::Compress::Gzip->new("$output_file.gz") or die "Error: Cannot open output file $output_file.gz: $GzipError\n";

# Print the header line for the output file
print $gz "Chromosome\tRefStart\tRefEnd\tRefStartBreakpoint\tRefEndBreakpoint\tReadName\tReadStart\tReadEnd\tReadStartBreakpoint\tReadEndBreakpoint\tType\tStrand\tIntegrity\tleftAnchorScore\tRightAnchorScore\treferencePath\tReadsPath\toriginSequence\tanchorSequence\n";

# 当前 read 的累积信息
my $prev_read = "";
my @alignments = ();

while (<IN>) {
    chomp;
    my @fields = split /\t/, $_;
    my $read_name = $fields[0];

    # Skip alignments with mapping quality below the minimum threshold
    next if ($fields[4] < $min_mapq);

    my $cigar_string = $fields[5];
    my $ref_chr = $fields[2];  # Reference chromosome/contig

    # If specified, skip alignments mapped to the mitochondrial chromosome
    if (defined $mitochondria_chr && $ref_chr eq $mitochondria_chr) {
        next;
    }

    # Calculate the alignment length and skip if it's below the defined threshold
    my $alignment_length = calculate_alignment_length($cigar_string);
    next if ($alignment_length < $min_align_len);
    
    # 如果当前 read 还在继续
    if ($prev_read eq "" || $read_name eq $prev_read) {
        push @alignments, parse_alignment(@fields);
        $prev_read = $read_name;
    } else {
        # 处理上一个 read
        process_read(\@alignments, $gz);
        @alignments = (parse_alignment(@fields));
        $prev_read = $read_name;
    }
}

# 处理最后一个 read
process_read(\@alignments, $gz);
$gz->close();
close(IN);
    

sub parse_alignment {
    my @f = @_;

    my $flag = $f[1];
    my $strand = ($flag & 0x10) ? "-" : "+";

    my ($read_start, $read_end, $ref_start, $ref_end) = calculate_read_positions($f[5], $f[9], $f[3]);

    return {
        read_name   => $f[0],
        ref_chr     => $f[2],
        ref_pos     => $f[3],
        mapq        => $f[4],
        cigar       => $f[5],
        seq         => $f[9],
        flag        => $flag,
        strand      => $strand,
        read_start  => $read_start,
        read_end    => $read_end,
        ref_start   => $ref_start,
        ref_end     => $ref_end
    };
}

sub process_read {
    my ($alignments_ref, $gz) = @_;
    my @sorted_alignments = sort { $a->{read_start} <=> $b->{read_start} } @$alignments_ref;
    return unless @sorted_alignments;

    my $read_name = $sorted_alignments[0]->{read_name};
    my %chromosomes_seen = map { $_->{ref_chr} => 1 } @sorted_alignments;
    my $num_segments = scalar(@sorted_alignments);
    my $read_length = length($sorted_alignments[0]->{seq});
    return if ($num_segments > 3 + 0.1 * ($read_length / 1000));

    if (@sorted_alignments >=3){
        for my $i (0 .. $#sorted_alignments - 2) {
            my $curr = $sorted_alignments[$i];
        
            # Search for split alignments matching SCGR conditions starting from index i+2
            for my $j ($i + 2 .. $#sorted_alignments) {
                my $next = $sorted_alignments[$j];
            
                # Ensure split alignment conditions are met
                if ($curr->{ref_chr} eq $next->{ref_chr} &&
                    $curr->{strand} eq $next->{strand} &&
                    $next->{read_start} > $curr->{read_end} &&
                    $next->{ref_start} > $curr->{ref_end}) {
                            
                    # Collect @new_alignments, including all split alignments from i to j
                    my @new_alignments = @sorted_alignments[$i .. $j];
        
                    my $svType = "SCGR";   # Define SV type as SCGR
                    my $leftIntegrity = "leftSpan"; # Left span integrity
                    my $rightIntegrity = "RightSpan"; # Right span integrity
        
                    # Process SCGR event, similar to how CGR (Complex Genomic Rearrangement) handles left and right breakpoints
                    my $left_breakpoint = $curr->{ref_end};   # Define left breakpoint
                    my $right_breakpoint = $next->{ref_start}; # Define right breakpoint
    
                    # Calculate the length of the SCGR event
                    my $SCGR_length_read = $next->{read_start} - $curr->{read_end} - 1;
                    my $SCGR_length = $next->{ref_start} - $curr->{ref_end};
    
                    # Ensure SCGR length is within allowed range
                    if ($SCGR_length <= $max_sv_len && $SCGR_length >= $min_sv_len && 
                        $SCGR_length_read <= $max_sv_len && $SCGR_length_read >= $min_sv_len) {
        
                    # Generate the output format, referencing ref_path and read_path, similar to how duplications are processed
                    my ($ref_path, $read_path) = find_SCGR_sv_paths(\@new_alignments);  # Call a new subroutine to generate paths with chromosome info 
        
                    # Identify anchors for both the current left and right alignments
                    my ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $left_anchor_score);
                    my ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $right_anchor_score);
            
                    # Find left anchor
                    ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $leftIntegrity, $left_anchor_score) = find_left_anchor(
                        $left_breakpoint, $curr->{cigar}, $curr->{ref_pos}, $svType, $curr->{ref_start}, $curr->{ref_end}, $curr->{read_start}, $curr->{read_end}, $leftIntegrity, $curr->{ref_chr}
                    );
        
                    # Find right anchor
                    ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $rightIntegrity, $right_anchor_score) = find_right_anchor(
                        $right_breakpoint, $next->{cigar}, $next->{ref_pos}, $svType, $next->{ref_start}, $next->{ref_end}, $next->{read_start}, $next->{read_end}, $rightIntegrity, $next->{ref_chr}
                    );
        
                    # Merge sequences between anchors
                    my $original_seq = substr($next->{seq}, $left_breakpoint_read, $right_breakpoint_read - ($left_breakpoint_read + 1) + 1);
                    my $merged_seq = substr($next->{seq}, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);
        
                    # Reverse complement the sequence if the strand is negative
                    if ($sorted_alignments[-1]->{strand} eq "-") {
                        # $merged_seq = revcomp($merged_seq);  # Uncomment if reverse complement is needed
                    }
        
                        # Define integrity and output the result
                        my $integrity = $leftIntegrity . $rightIntegrity;
        
                        # Format output for the merged SCGR
                        my $output = join("\t", $next->{ref_chr}, $left_anchor_ref - 1, $right_anchor_ref, $left_breakpoint_ref, 
                            $right_breakpoint_ref, $read_name, $left_anchor_read - 1, $right_anchor_read, $left_breakpoint_read + 1, $right_breakpoint_read, 
                            $svType, $next->{strand}, $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $original_seq, $merged_seq);
        
                        # Add output to the output lines
                        print $gz "$output\n";
                    }
                }
            }
        }
    }
    if (@sorted_alignments >= 2){
            
        #treat complex insertion and deletion
        for (my $i = 0; $i <= @sorted_alignments - 1; $i++) {
            my $current = $sorted_alignments[$i];
            my $type = "complex";
            my $output = process_indels($current, $read_name, $type);
            print $gz "$output\n" if $output;
        }
            
        #treat complex genomic arrangements
        my $control = 0;  # Variable to track valid CGRs
    
        for (my $i = 0; $i < @sorted_alignments - 1; $i++) {
            my $curr = $sorted_alignments[$i];
            my $next = $sorted_alignments[$i + 1];
                
            # Ensure split alignment conditions are met
            if ($curr->{ref_chr} eq $next->{ref_chr} &&
                $curr->{strand} eq $next->{strand} &&
                $next->{ref_start} + $min_sv_len > $curr->{ref_end}) {
    
                my $svType = "CGR";
                my $leftIntegrity = "leftSpan";
                my $rightIntegrity = "RightSpan";
    
                # Define breakpoints for the current and next alignments
                my ($curr_left_breakpoint,$curr_leftReads_breakpoint);
                my $next_right_breakpoint = $next->{ref_start};
                my $next_rightReads_breakpoint= $next->{read_start};
                if ($next->{read_start} < $curr->{read_end}){ #read position overlap
                    # Parse the CIGAR string into operations
                    my @cigar_operations;
                    while ($curr->{cigar} =~ /(\d+)([MIDNSHP=X])/g) {
                        my ($len, $op) = ($1, $2);
                        push @cigar_operations, { len => $len, op => $op };
                    }
                    my $read_pos = $curr->{read_end};
                    my $ref_pos = $curr->{ref_end};
                    # Traverse CIGAR operations to find the breakpoint
                    my $startIndex;
                    if ($cigar_operations[$#cigar_operations] eq "S"){
                        $startIndex = $#cigar_operations - 1;
                    }else{
                        $startIndex = $#cigar_operations;
                    }
                        
                    my $is_find = 0; ##check whether the left_breakpoint can be found
                    for (my $i = $startIndex; $i >= 0; $i--) {
                        my $op = $cigar_operations[$i]->{op};
                        my $len = $cigar_operations[$i]->{len};
                        
                        if ($read_pos < $next->{read_start}){
                            $curr_left_breakpoint = $ref_pos;
                            $curr_leftReads_breakpoint = $read_pos;
                            $is_find = 1;
                            last;
                        }
                        # Adjust current_ref_pos and current_read_pos during the loop
                        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
                            $ref_pos -= $len;  # Move left
                        }
                        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
                            $read_pos -= $len;  # Move left
                        }
                    }
                    if ($is_find == 0){
                        next;
                    }
                }else{#read position non-overlap
                    $curr_left_breakpoint = $curr->{ref_end};
                    $curr_leftReads_breakpoint = $curr->{read_end};
                }
    
                # Calculate CGR length based on reference and read breakpoints
                my $CGR_length_ref = $next_right_breakpoint - $curr_left_breakpoint;
                my $CGR_length_read = $next_rightReads_breakpoint - $curr_leftReads_breakpoint - 1;
    
                # Skip if the CGR length exceeds maximum SV length
                if ($CGR_length_ref > $max_sv_len || $CGR_length_read > $max_sv_len) {
                    #if (($i == @sorted_alignments - 2) && ($control == 0)) {
                    #}
                    next;
                } else {
                    $control++;
                }
    
                # Identify anchors for both current and next alignments
                my ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $left_anchor_score);
                my ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $right_anchor_score);
    
                ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $leftIntegrity, $left_anchor_score) = find_left_anchor(
                    $curr_left_breakpoint, $curr->{cigar}, $curr->{ref_pos}, $svType, $curr->{ref_start}, $curr->{ref_end}, $curr->{read_start}, $curr->{read_end}, $leftIntegrity, $curr->{ref_chr}
                );
    
                ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $rightIntegrity, $right_anchor_score) = find_right_anchor(
                    $next_right_breakpoint, $next->{cigar}, $next->{ref_pos}, $svType, $next->{ref_start}, $next->{ref_end}, $next->{read_start}, $next->{read_end}, $rightIntegrity, $next->{ref_chr}
                );
    
                if ($left_breakpoint_read == $right_breakpoint_read){
                    # Merge sequences between anchors
                    my $original_seq =  ".";
                    my $merged_seq = substr($curr->{seq}, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);
                                
                    # Define integrity and output the result
                    my $integrity = $leftIntegrity . $rightIntegrity;
    
                    # Format output for the merged CGR
                    my $output = join("\t", $curr->{ref_chr}, $left_anchor_ref, $right_anchor_ref, $left_breakpoint_ref, $right_breakpoint_ref, 
                        $read_name, $left_anchor_read - 1, $right_anchor_read, $left_breakpoint_read + 1, $right_breakpoint_read + 1, 
                        $svType, $curr->{strand}, $integrity, $left_anchor_score, $right_anchor_score, ".", ".", $original_seq, $merged_seq);
    
                    print $gz "$output\n";                            
                }else{
                    # Merge sequences between anchors
                    my $original_seq = substr($curr->{seq}, $left_breakpoint_read, $right_breakpoint_read - ($left_breakpoint_read + 1) + 1);
                    my $merged_seq = substr($curr->{seq}, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);
    
                    # Reverse complement the sequence if the strand is negative
                    if ($curr->{strand} eq "-") {
                        # $merged_seq = revcomp($merged_seq);  # Uncomment if reverse complement is needed
                    }
    
                    # Define integrity and output the result
                    my $integrity = $leftIntegrity . $rightIntegrity;
    
                    # Format output for the merged CGR
                    my $output = join("\t", $curr->{ref_chr}, $left_anchor_ref - 1, $right_anchor_ref, $left_breakpoint_ref, $right_breakpoint_ref, 
                        $read_name, $left_anchor_read - 1, $right_anchor_read, $left_breakpoint_read + 1, $right_breakpoint_read, 
                        $svType, $curr->{strand}, $integrity, $left_anchor_score, $right_anchor_score, ".", ".", $original_seq, $merged_seq);
    
                    print $gz "$output\n";
                }
            }
        }
    }
    if (@sorted_alignments == 1){
        my $type = "simple";
        my $output = process_indels($sorted_alignments[0], $read_name, $type);
        print $gz "$output\n" if $output;
    }
        
    if ((scalar(keys %chromosomes_seen) == 1) && (@sorted_alignments > 1)) {
            
        # All alignments are on the same strand
        my $all_on_same_strand = 1;
        my $strand_pattern = $sorted_alignments[0]->{strand};  # Start with the first strand
        for my $i (1 .. $#sorted_alignments) {
            if ($sorted_alignments[$i]->{strand} ne $strand_pattern) {
                $all_on_same_strand = 0;
                last;
            }
        }
                        
        if ($all_on_same_strand) {
    
            #check reference duplication
            my $is_refDup = 0;
            my $min_index;
            my $leftBreakpoint;
            my $max_index;
            my $rightBreakpoint;
            my $is_read_overlap = 0;
            # check reference position overlap
            for (my $i = 0; $i < @sorted_alignments - 1; $i++) {
                my $curr = $sorted_alignments[$i];
                
                for (my $j = $i + 1; $j < @sorted_alignments; $j++) {
                    my $next = $sorted_alignments[$j];
                            
                    if ($next->{read_start} + $min_sv_len < $curr->{read_end}){
                        $is_read_overlap = 1;
                    }
                            
                    # Calculate the overlap start (maximum of the two start positions)
                    my $overlap_start = ($curr->{ref_start} > $next->{ref_start}) ? $curr->{ref_start} : $next->{ref_start};
    
                    # Calculate the overlap end (minimum of the two end positions)
                    my $overlap_end = ($curr->{ref_end} < $next->{ref_end}) ? $curr->{ref_end} : $next->{ref_end};
    
                    # If overlap exists, calculate the overlap length
                    if ($overlap_start < $overlap_end) {
                        my $overlap_length = $overlap_end - $overlap_start;
                        if ($overlap_length >= $min_sv_len && $overlap_length <= $max_sv_len) {
                            $is_refDup = 1;
                            if (defined $min_index){
                                if ($i< $min_index){
                                    $min_index = $i;
                                    $leftBreakpoint = $next->{ref_start};
                                }
                                if ($j > $max_index){
                                    $max_index = $j;
                                    $rightBreakpoint = $curr->{ref_end};
                                }
                            }else{
                                $min_index = $i;
                                $leftBreakpoint = $next->{ref_start};
                                $max_index = $j;
                                $rightBreakpoint = $curr->{ref_end};
                            }
                        }
                    }
                }
            }
            if (($is_refDup == 1) && ($is_read_overlap == 0)){
            
                # Collect @new_alignments, including all split alignments from $min_index to $max_index
                my @new_alignments = @sorted_alignments[$min_index .. $max_index];
     
                my $svType = "Duplication";
                my $first_alignment = $sorted_alignments[$min_index];
                my $leftIntegrity = ($first_alignment->{ref_start} < $leftBreakpoint) ? "leftSpan" : (($first_alignment->{ref_start} == $leftBreakpoint) ? "leftCover" : "leftUncover");
    
                my $last_alignment = $sorted_alignments[$max_index];
                my $rightIntegrity = ($last_alignment->{ref_end} > $rightBreakpoint) ? "RightSpan" : (($last_alignment->{ref_end} == $rightBreakpoint) ? "RightCover" : "RightUncover");
    
                # Find left and right anchors
                my ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $left_anchor_score);
                my ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $right_anchor_score);    
    
                if ($rightIntegrity eq "RightCover" || $rightIntegrity eq "RightUncover") {
                    $right_breakpoint_ref = $rightBreakpoint;
                    $right_anchor_ref = $rightBreakpoint;
                    $right_breakpoint_read = $last_alignment->{read_end};
                    $right_anchor_read = $last_alignment->{read_end};
                    $right_anchor_score = ".";
                } elsif ($rightIntegrity eq "RightSpan") {
                    ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $rightIntegrity, $right_anchor_score) =
                        find_right_anchor($rightBreakpoint, $last_alignment->{cigar}, $last_alignment->{ref_pos}, $svType, $last_alignment->{ref_start}, $last_alignment->{ref_end}, $last_alignment->{read_start}, $last_alignment->{read_end}, $rightIntegrity, $last_alignment->{ref_chr});
                }
    
                if ($leftIntegrity eq "leftCover" || $leftIntegrity eq "leftUncover") {
                    $left_breakpoint_ref = $leftBreakpoint;
                    $left_anchor_ref = $leftBreakpoint;
                    $left_breakpoint_read = $first_alignment->{read_start};
                    $left_anchor_read = $first_alignment->{read_start};
                    $left_anchor_score = ".";
                } elsif ($leftIntegrity eq "leftSpan") {
                    ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $leftIntegrity, $left_anchor_score) =
                        find_left_anchor($leftBreakpoint, $first_alignment->{cigar}, $first_alignment->{ref_pos}, $svType, $first_alignment->{ref_start}, $first_alignment->{ref_end}, $first_alignment->{read_start}, $first_alignment->{read_end}, $leftIntegrity, $first_alignment->{ref_chr});
                    }
    
                my $integrity = $leftIntegrity . $rightIntegrity;
    
                # Find reference and read paths
                my ($ref_path, $read_path) = find_sv_paths(\@new_alignments);
    
                # Prepare the final output based on integrity
                if ($integrity eq "leftSpanRightSpan") {
                    my $origin_seq = substr($sorted_alignments[1]->{seq}, $left_breakpoint_read - 1, $right_breakpoint_read - $left_breakpoint_read + 1);
                    my $merged_seq = substr($sorted_alignments[1]->{seq}, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);
    
                    my $output = join("\t", $sorted_alignments[1]->{ref_chr}, $left_anchor_ref, $right_anchor_ref,
                        $left_breakpoint_ref, $right_breakpoint_ref, $read_name, $left_anchor_read - 1, $right_anchor_read,
                        $left_breakpoint_read, $right_breakpoint_read, $svType, $sorted_alignments[0]->{strand}, $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $origin_seq, $merged_seq);
                    print $gz "$output\n";
                } elsif ($integrity eq "leftUncoverRightSpan" || $integrity eq "leftCoverRightSpan") {
                    my $output = join("\t", $sorted_alignments[1]->{ref_chr}, $left_anchor_ref, $right_anchor_ref,
                        $left_breakpoint_ref, $right_breakpoint_ref, $read_name, 1, $right_anchor_read,
                        1, $right_breakpoint_read, $svType, $sorted_alignments[0]->{strand}, $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $sorted_alignments[1]->{seq}, ".");
                    print $gz "$output\n";
                } elsif ($integrity eq "leftSpanRightCover" || $integrity eq "leftSpanRightUncover") {
                    my $output = join("\t", $sorted_alignments[1]->{ref_chr}, $left_anchor_ref, $right_anchor_ref,
                        $left_breakpoint_ref, $right_breakpoint_ref, $read_name, $left_anchor_read - 1, $right_anchor_read,
                        $left_breakpoint_read, $right_breakpoint_read, $svType, $sorted_alignments[0]->{strand}, $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $sorted_alignments[1]->{seq}, ".");
                    print $gz "$output\n";
                } else {
                    my $output = join("\t", $sorted_alignments[1]->{ref_chr}, $left_anchor_ref, $right_anchor_ref,
                        $left_breakpoint_ref, $right_breakpoint_ref, $read_name, 1, $right_anchor_read,
                        1, $right_breakpoint_read, $svType, $sorted_alignments[0]->{strand}, $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $sorted_alignments[1]->{seq}, ".");
                    print $gz "$output\n";
                }
            }
            if ($is_read_overlap == 1){
                for (my $i = 0; $i < @sorted_alignments - 1; $i++) {
                    print "$sorted_alignments[$i]->{ref_chr}\t$sorted_alignments[$i]->{ref_start}\t$sorted_alignments[$i]->{ref_end}\t$sorted_alignments[$i]->{read_start}\t$sorted_alignments[$i]->{read_end}\n";
                }    
            }
        }
    }   
}

# Subroutine to process insertions and deletions using anchor length
sub process_indels {
    my ($alignment, $read_name, $type) = @_;
    my $cigar = $alignment->{cigar};
    my $ref_pos = $alignment->{ref_pos};
    my $seq = $alignment->{seq};
    my $read_pos = 1;
    my $strand = $alignment->{strand};
    my $output;

    # Parse CIGAR string into operations and store them in an array with index tracking
    my @cigar_operations;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);
        push @cigar_operations, { len => $len, op => $op };
    }
    
    # Iterate through the cigar_operations array to process each operation
    for (my $cigar_index = 0; $cigar_index < @cigar_operations; $cigar_index++) {
        my $op = $cigar_operations[$cigar_index]->{op};
        my $len = $cigar_operations[$cigar_index]->{len};

        if ($op eq 'D') {  # Deletion
            if ($len >= $min_sv_len && $len <= $max_sv_len) {
                # Reference genome left and right breakpoints for deletion
                my $left_breakpoint = $ref_pos;
                my $right_breakpoint = $ref_pos + $len;

                my $svType;
                if ($type eq "simple") {
                    $svType = "Deletion";
                }else{
                    $svType = "ComplexDeletion";
                }
                my $leftIntegrity = "leftSpan";
                my $rightIntegrity = "RightSpan";

                # Find the left and right anchors for the reference and read
                my ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $left_anchor_score);
                my ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $right_anchor_score);
                
                ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $leftIntegrity, $left_anchor_score) = find_left_anchor(
                    $left_breakpoint, $cigar, $alignment->{ref_pos}, $svType, $alignment->{ref_start}, $alignment->{ref_end}, $alignment->{read_start}, $alignment->{read_end}, $leftIntegrity, $alignment->{ref_chr}
                );
                
                ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $rightIntegrity, $right_anchor_score) = find_right_anchor(
                    $right_breakpoint, $cigar, $alignment->{ref_pos}, $svType, $alignment->{ref_start}, $alignment->{ref_end}, $alignment->{read_start}, $alignment->{read_end}, $rightIntegrity, $alignment->{ref_chr}
                );

                # No original sequence for deletion, use a placeholder '.'
                my $orgin_seq = ".";
                my $deletion_seq = substr($seq, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);

                # Reverse complement the sequence if necessary (for negative strand)
                #$deletion_seq = revcomp($deletion_seq) if $strand eq "-";

                my $integrity = $leftIntegrity . $rightIntegrity;

                # Generate the output for deletion adjusted by anchors
                $output = join("\t", $alignment->{ref_chr}, $left_anchor_ref, $right_anchor_ref, $left_breakpoint_ref, $right_breakpoint_ref, $read_name,
                               $left_anchor_read - 1, $right_anchor_read, $left_breakpoint_read, $right_breakpoint_read + 1, $svType, $strand, $integrity, $left_anchor_score, $right_anchor_score, ".", ".", $orgin_seq, $deletion_seq);
            }
            $ref_pos += $len;  # Move the reference position after deletion
        }
        elsif ($op eq 'I') {  # Insertion
            if ($len >= $min_sv_len && $len <= $max_sv_len) {
                # Insertion points are based on read position; reference position remains unchanged
                my $insertion_start_ref = $ref_pos;
                my $insertion_end_ref = $ref_pos;  # No change in reference genome for insertion

                my $svType;
                if ($type eq "simple") {
                    $svType = "Insertion";
                }else{
                    $svType = "ComplexInsertion";
                }

                my $leftIntegrity = "leftSpan";
                my $rightIntegrity = "RightSpan";

                # Find the anchor positions in the reference genome and read
                my ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $left_anchor_score);
                my ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $right_anchor_score);
                
                ($left_breakpoint_ref, $left_anchor_ref, $left_breakpoint_read, $left_anchor_read, $leftIntegrity, $left_anchor_score) = find_left_anchor(
                    $insertion_start_ref, $cigar, $alignment->{ref_pos}, $svType, $alignment->{ref_start}, $alignment->{ref_end}, $alignment->{read_start}, $alignment->{read_end}, $leftIntegrity, $alignment->{ref_chr}
                );
                
                ($right_breakpoint_ref, $right_anchor_ref, $right_breakpoint_read, $right_anchor_read, $rightIntegrity, $right_anchor_score) = find_right_anchor(
                    $insertion_end_ref, $cigar, $alignment->{ref_pos}, $svType, $alignment->{ref_start}, $alignment->{ref_end}, $alignment->{read_start}, $alignment->{read_end}, $rightIntegrity, $alignment->{ref_chr}
                );

                # Extract insertion sequence from read based on anchor boundaries
                my $orgin_seq = substr($seq, $left_breakpoint_read, $right_breakpoint_read - ($left_breakpoint_read + 1) + 1);
                my $insertion_seq = substr($seq, $left_anchor_read - 2, $right_anchor_read - ($left_anchor_read - 1) + 1);

                # Reverse complement if necessary (for negative strand)
                #$insertion_seq = revcomp($insertion_seq) if $strand eq "-";

                my $integrity = $leftIntegrity . $rightIntegrity;

                # Generate output for insertion
                $output = join("\t", $alignment->{ref_chr}, $left_anchor_ref, $right_anchor_ref, $left_breakpoint_ref, $right_breakpoint_ref, $read_name,
                               $left_anchor_read - 1, $right_anchor_read, $left_breakpoint_read + 1, $right_breakpoint_read, $svType, $strand, $integrity, $left_anchor_score, $right_anchor_score, ".", ".", $orgin_seq, $insertion_seq);
            }
            $read_pos += $len;  # Move the read position after insertion
        }
        elsif ($op eq 'M' || $op eq '=' || $op eq 'X') {  # Match or Mismatch
            # Update both reference and read positions
            $ref_pos += $len;
            $read_pos += $len;
        }
        elsif ($op eq 'S') {  # Soft clip
            # Update read position
            $read_pos += $len;
        }
    }

    return $output;
}

# Subroutine to reverse complement a DNA sequence
sub revcomp {
    my ($seq) = @_;
    # Reverse the sequence
    $seq = reverse($seq);
    # Translate ACGT to their complements TGCA
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}

# Subroutine to calculate alignment length from CIGAR string
sub calculate_alignment_length {
    my ($cigar) = @_;
    my $length = 0;

    # Parse the CIGAR string and calculate the alignment length
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);

        # Add length for operations that consume reference bases
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'N' || $op eq 'D') {
            $length += $len;
        }
    }
    return $length;
}

# Subroutine to calculate read_start, read_end, ref_start, and ref_end based on CIGAR string
sub calculate_read_positions {
    my ($cigar, $seq, $ref_pos) = @_;  # $ref_pos is the starting reference genome position

    my $read_start = -1;   # Initialize as -1, to be updated later
    my $read_end = -1;     # Initialize as -1, to be updated later
    my $ref_start = -1;    # Initialize as -1, to be updated later
    my $ref_end = -1;      # Initialize as -1, to be updated later

    my $current_read_pos = 1;          # Current read position, starts from 1
    my $current_ref_pos = $ref_pos;    # Current reference genome position, starts from the passed reference position

    # Parse the CIGAR string and track positions of M, =, and X operations
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($length, $op) = ($1, $2);

        # Only handle M, =, X, D, N operations as they affect alignment
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            if ($read_start == -1) {
                # Set the start position for both read and reference genome on the first occurrence of M, =, or X
                $read_start= $current_read_pos;
                $ref_start = $current_ref_pos;
            }
            # Update the end position for both read and reference genome on each occurrence of M, =, or X
            $read_end = $current_read_pos + $length - 1;
            $ref_end = $current_ref_pos + $length;
        }

        # Move the read position based on CIGAR operations that consume read bases
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
            $current_read_pos += $length;  # Move read position
        }
        # Move the reference genome position based on CIGAR operations that consume reference bases
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            $current_ref_pos += $length;   # Move reference position
        }
    }

    # Return the start and end duplication positions for both read and reference genome
    return ($read_start, $read_end, $ref_start, $ref_end);
}

# Subroutine to find the leftmost anchor
sub find_left_anchor {
    my ($breakpoint, $cigar, $ref_pos, $svType, $ref_start, $ref_end, $read_start, $read_end, $integrity, $ref_chr) = @_;

    my $anchor_len;
    if ($svType eq "Insertion"){
        $anchor_len = $anchor_len_insertion;
    }elsif ($svType eq "Deletion"){
        $anchor_len = $anchor_len_deletion;
    }elsif ($svType eq "CGR"){
        $anchor_len = $anchor_len_CGR;
    }elsif ($svType eq "Duplication"){
        $anchor_len = $anchor_len_duplication;
    }elsif ($svType eq "ComplexInsertion"){
        $anchor_len = $anchor_len_complexInsertion;
    }elsif ($svType eq "ComplexDeletion"){
        $anchor_len = $anchor_len_complexDeletion;
    }elsif ($svType eq "SCGR"){
        $anchor_len = $anchor_len_SCGR;
    }
    
    # Cache for previous overlap region query
    my ($cached_start, $cached_end, $cached_overlap_result) = (-1, -1, undef);

    # Parse the CIGAR string into operations
    my @cigar_operations;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);
        push @cigar_operations, { len => $len, op => $op };
    }

    # Initialize current reference and read positions
    my $current_ref_pos = $ref_pos;
    my $current_read_pos = 1;
    my $cigar_index;
    my $breakpoint_read;
    my $breakpoint_ref = $breakpoint;
    my $last_len;

    # Traverse CIGAR operations to find the breakpoint
    for (my $i = 0; $i <= $#cigar_operations; $i++) {
        my $op = $cigar_operations[$i]->{op};
        my $len = $cigar_operations[$i]->{len};

        # Check for matching positions and adjust accordingly
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            if ($current_ref_pos == $breakpoint) {
                if ($svType eq "Insertion") {
                    $cigar_index = $i;
                    $breakpoint_read = $current_read_pos - $last_len - 1;
                    last;
                } else {
                    $cigar_index = $i;
                    $breakpoint_read = $current_read_pos;
                    last;
                }
            } elsif ($current_ref_pos > $breakpoint) {
                $cigar_index = $i - 1;
                $breakpoint_read = $current_read_pos + $len - ($current_ref_pos + $len - $breakpoint) - 1;
                $current_read_pos -= $last_len;
                $current_ref_pos -= $last_len;
                last;
            }
            $current_ref_pos += $len;
        }

        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
            $current_read_pos += $len;
            $last_len = $len;
        }
    }

    # The cigar index for CGR was the last index;
    if (!defined $cigar_index && ($svType eq "CGR" || $svType eq "SCGR" || $svType eq "Duplication")) {
        $cigar_index = $#cigar_operations;  # Set to last index
        $breakpoint_ref = $ref_end;
        $breakpoint_read = $read_end;
        $current_read_pos -= $cigar_operations[$#cigar_operations]->{len};
    }

    # Find the maximum_score;
    my $score = 0;
    my ($anchor_ref_pos, $anchor_read_pos);
    my $curr_len = 0; #compare length

    # Traverse CIGAR operations backward until identity and length thresholds are met
    for (my $i = $cigar_index - 1; $i >= 0; $i--) {
        my $op = $cigar_operations[$i]->{op};
        my $len = $cigar_operations[$i]->{len};

        if (($op eq '=') && ($len >= $anchor_len)){
            
            my $char_ref_end = $current_ref_pos; #end position for candidate anchor
            my $char_ref_start = $current_ref_pos - $len; #start position for candidate anchor

            if ($score > 2){
                if ($len > $curr_len){
                    if ($first_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                            $overlaps = $cached_overlap_result;
                        } else {
                            $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                            if (@$overlaps) {
                                my $region = $overlaps->[0];
                                $cached_overlap_result = $overlaps;
                                $cached_start = $region->{start};
                                $cached_end = $region->{end};
                            }
                        }
   
                        if (@$overlaps) {
                            $score = $len/100000 + 2;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos - $len;
                            $anchor_read_pos = $current_read_pos - $len;
                        }
                    }
                }
            } elsif ($score < 2 && $score > 1){
                
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
   
                        #undate score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos - $len;
                        $anchor_read_pos = $current_read_pos - $len;                        
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    if ($len > $curr_len){
                        if ($second_tree{$ref_chr}) {
                            # Check cache or fetch new overlaps
                            my $overlaps;
                            if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                                $overlaps = $cached_overlap_result;
                            } else {
                                $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                                if (@$overlaps) {
                                    my $region = $overlaps->[0];
                                    $cached_overlap_result = $overlaps;
                                    $cached_start = $region->{start};
                                    $cached_end = $region->{end};
                                }
                            }
   
                            if (@$overlaps) {
                                $score = $len/100000 + 1;
                                $curr_len = $len;
                                $anchor_ref_pos = $current_ref_pos - $len;
                                $anchor_read_pos = $current_read_pos - $len;
                            }
                        }
                    }
                }    
            } elsif ($score < 1 && $score > 0){
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
   
                        #update score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos - $len;
                        $anchor_read_pos = $current_read_pos - $len;
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    my $is_second = 0; #check whether overlapped the second grade
                    
                    if ($second_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                        if (@$overlaps) {
                            #update cache
                            my $region = $overlaps->[0];
                            $cached_overlap_result = $overlaps;
                            $cached_start = $region->{start};
                            $cached_end = $region->{end};

                            #update score
                            $score = $len/100000 + 1;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos - $len;
                            $anchor_read_pos = $current_read_pos - $len;
                            $is_second = 1;
                        }
                    }
                    
                    if ($is_second == 0){
                        if ($len > $curr_len){
                            if ($third_tree{$ref_chr}) {
                                # Check cache or fetch new overlaps
                                my $overlaps;
                                if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                                    $overlaps = $cached_overlap_result;
                                } else {
                                    $overlaps = $third_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                                    if (@$overlaps) {
                                        my $region = $overlaps->[0];
                                        $cached_overlap_result = $overlaps;
                                        $cached_start = $region->{start};
                                        $cached_end = $region->{end};
                                    }
                                }
   
                                if (@$overlaps) {
                                    $score = $len/100000;
                                    $curr_len = $len;
                                    $anchor_ref_pos = $current_ref_pos - $len;
                                    $anchor_read_pos = $current_read_pos - $len;
                                }
                            }
                        }
                    }
                }
            } elsif ($score == 0){
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache 
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
    
                        #update score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos - $len;
                        $anchor_read_pos = $current_read_pos - $len;
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    my $is_second = 0; #check whether overlapped the second grade
                    
                    if ($second_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                        if (@$overlaps) {
                            #update cache
                            my $region = $overlaps->[0];
                            $cached_overlap_result = $overlaps;
                            $cached_start = $region->{start};
                            $cached_end = $region->{end};

                            #update score
                            $score = $len/100000 + 1;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos - $len;
                            $anchor_read_pos = $current_read_pos - $len;
                            $is_second = 1;
                        }
                    }
                    
                    if ($is_second == 0){
                        if ($third_tree{$ref_chr}) {
                            # Check cache or fetch new overlaps
                            my $overlaps;
                            $overlaps = $third_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                            if (@$overlaps) {
                                #update cache
                                my $region = $overlaps->[0];
                                $cached_overlap_result = $overlaps;
                                $cached_start = $region->{start};
                                $cached_end = $region->{end};

                                #update score
                                $score = $len/100000;
                                $curr_len = $len;
                                $anchor_ref_pos = $current_ref_pos - $len;
                                $anchor_read_pos = $current_read_pos - $len;
                            }
                        }
                    }
                }
            }
        }

        # Adjust current_ref_pos and current_read_pos during the loop
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            $current_ref_pos -= $len;  # Move left
        }
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
            $current_read_pos -= $len;  # Move left
        }
    }
    
    if ($score > 0){
        return ($breakpoint_ref, $anchor_ref_pos + 1, $breakpoint_read, $anchor_read_pos + 1, $integrity, $score);
    }else{
        return ($breakpoint_ref, $breakpoint_ref, $breakpoint_read, $breakpoint_read, "leftCover", ".");
    }
}

# Subroutine to find the rightmost anchor
sub find_right_anchor {
    my ($breakpoint, $cigar, $ref_pos, $svType, $ref_start, $ref_end, $read_start, $read_end, $integrity, $ref_chr) = @_;

    my $anchor_len;
    if ($svType eq "Insertion"){
        $anchor_len = $anchor_len_insertion;
    }elsif ($svType eq "Deletion"){
        $anchor_len = $anchor_len_deletion;
    }elsif ($svType eq "CGR"){
        $anchor_len = $anchor_len_CGR;
    }elsif ($svType eq "Duplication"){
        $anchor_len = $anchor_len_duplication;
    }elsif ($svType eq "ComplexInsertion"){
        $anchor_len = $anchor_len_complexInsertion;
    }elsif ($svType eq "ComplexDeletion"){
        $anchor_len = $anchor_len_complexDeletion;
    }elsif ($svType eq "SCGR"){
        $anchor_len = $anchor_len_SCGR;
    }

    # Cache for previous overlap region query
    my ($cached_start, $cached_end, $cached_overlap_result) = (-1, -1, undef);

    # Parse the CIGAR string into operations
    my @cigar_operations;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($len, $op) = ($1, $2);
        push @cigar_operations, { len => $len, op => $op };
    }

    # Initialize current reference and read positions
    my $current_ref_pos = $ref_pos;
    my $current_read_pos = 0;  # Start at 0
    my $cigar_index;
    my $breakpoint_read;
    my $breakpoint_ref = $breakpoint;

    my $last_operation;
    my $last_len;
    
    # Traverse CIGAR operations to find the breakpoint
    for (my $i = 0; $i <= $#cigar_operations; $i++) {
        my $op = $cigar_operations[$i]->{op};
        my $len = $cigar_operations[$i]->{len};

        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            if ($current_ref_pos == $breakpoint) {
                $cigar_index = $i;
                $breakpoint_read = $current_read_pos;
                last;
            } elsif ($current_ref_pos > $breakpoint) {
                $cigar_index = $i;
                if ($last_operation eq "I") {
                    $breakpoint_read = $current_read_pos + $len - ($current_ref_pos + $len - $breakpoint) - $last_len;
                } else {
                    $breakpoint_read = $current_read_pos + $len - ($current_ref_pos + $len - $breakpoint);
                }
                last;
            }
            $current_ref_pos += $len;
        }

        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
            $current_read_pos += $len;
        }
        $last_operation = $op;
        $last_len = $len;
    }

    # If no appropriate breakpoint is found, use default values
    if (!defined $cigar_index && ($svType eq "Duplication" || $svType eq "SCGR" || $svType eq "CGR")) {
        $cigar_index = $#cigar_operations;
        $breakpoint_ref = $ref_end;
        $breakpoint_read = $read_end;
    }

    # Find the maximum_score;
    my $score = 0;
    my ($anchor_ref_pos, $anchor_read_pos);
    my $curr_len = 0; #compare length

    # Traverse CIGAR operations forward until identity and length thresholds are met
    for (my $i = $cigar_index; $i <= $#cigar_operations; $i++) {
        my $op = $cigar_operations[$i]->{op};
        my $len = $cigar_operations[$i]->{len};

        if (($op eq '=') && ($len >= $anchor_len)) {
            my $char_ref_start = $current_ref_pos;
            my $char_ref_end = $char_ref_start + $len;
            
            if ($score > 2){
                if ($len > $curr_len){
                    if ($first_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                            $overlaps = $cached_overlap_result;
                        } else {
                            $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                            if (@$overlaps) {
                                my $region = $overlaps->[0];
                                $cached_overlap_result = $overlaps;
                                $cached_start = $region->{start};
                                $cached_end = $region->{end};
                            }
                        }
   
                        if (@$overlaps) {
                            $score = $len/100000 + 2;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos + $len;
                            $anchor_read_pos = $current_read_pos + $len;
                        }
                    }
                }
            } elsif ($score < 2 && $score > 1){
                
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
   
                        #undate score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos + $len;
                        $anchor_read_pos = $current_read_pos + $len;                        
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    if ($len > $curr_len){
                        if ($second_tree{$ref_chr}) {
                            # Check cache or fetch new overlaps
                            my $overlaps;
                            if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                                $overlaps = $cached_overlap_result;
                            } else {
                                $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                                if (@$overlaps) {
                                    my $region = $overlaps->[0];
                                    $cached_overlap_result = $overlaps;
                                    $cached_start = $region->{start};
                                    $cached_end = $region->{end};
                                }
                            }
   
                            if (@$overlaps) {
                                $score = $len/100000 + 1;
                                $curr_len = $len;
                                $anchor_ref_pos = $current_ref_pos + $len;
                                $anchor_read_pos = $current_read_pos + $len;
                            }
                        }
                    }
                }    
            }elsif ($score < 1 && $score > 0){
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
   
                        #update score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos + $len;
                        $anchor_read_pos = $current_read_pos + $len;
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    my $is_second = 0; #check whether overlapped the second grade
                    
                    if ($second_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                        if (@$overlaps) {
                            #update cache
                            my $region = $overlaps->[0];
                            $cached_overlap_result = $overlaps;
                            $cached_start = $region->{start};
                            $cached_end = $region->{end};

                            #update score
                            $score = $len/100000 + 1;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos + $len;
                            $anchor_read_pos = $current_read_pos + $len;
                            $is_second = 1;
                        }
                    }
                    
                    if ($is_second == 0){
                        if ($len > $curr_len){
                            if ($third_tree{$ref_chr}) {
                                # Check cache or fetch new overlaps
                                my $overlaps;
                                if ($char_ref_start >= $cached_start && $char_ref_end <= $cached_end) {
                                    $overlaps = $cached_overlap_result;
                                } else {
                                    $overlaps = $third_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                                    if (@$overlaps) {
                                        my $region = $overlaps->[0];
                                        $cached_overlap_result = $overlaps;
                                        $cached_start = $region->{start};
                                        $cached_end = $region->{end};
                                    }
                                }
   
                                if (@$overlaps) {
                                    $score = $len/100000;
                                    $curr_len = $len;
                                    $anchor_ref_pos = $current_ref_pos + $len;
                                    $anchor_read_pos = $current_read_pos + $len;
                                }
                            }
                        }
                    }
                }
            } elsif ($score == 0){
                my $is_first = 0; #check whether overlapped the first grade 
                
                if ($first_tree{$ref_chr}) {
                    # Check cache or fetch new overlaps
                    my $overlaps;
                    $overlaps = $first_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                    if (@$overlaps) {
                        #update cache 
                        my $region = $overlaps->[0];
                        $cached_overlap_result = $overlaps;
                        $cached_start = $region->{start};
                        $cached_end = $region->{end};
    
                        #update score
                        $score = $len/100000 + 2;
                        $curr_len = $len;
                        $anchor_ref_pos = $current_ref_pos + $len;
                        $anchor_read_pos = $current_read_pos + $len;
                        $is_first = 1;
                    }
                }
                
                if ($is_first == 0){
                    my $is_second = 0; #check whether overlapped the second grade
                    
                    if ($second_tree{$ref_chr}) {
                        # Check cache or fetch new overlaps
                        my $overlaps;
                        $overlaps = $second_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                        if (@$overlaps) {
                            #update cache
                            my $region = $overlaps->[0];
                            $cached_overlap_result = $overlaps;
                            $cached_start = $region->{start};
                            $cached_end = $region->{end};

                            #update score
                            $score = $len/100000 + 1;
                            $curr_len = $len;
                            $anchor_ref_pos = $current_ref_pos + $len;
                            $anchor_read_pos = $current_read_pos + $len;
                            $is_second = 1;
                        }
                    }
                    
                    if ($is_second == 0){
                        if ($third_tree{$ref_chr}) {
                            # Check cache or fetch new overlaps
                            my $overlaps;
                            $overlaps = $third_tree{$ref_chr}->fetch($char_ref_start, $char_ref_end);
                            if (@$overlaps) {
                                #update cache
                                my $region = $overlaps->[0];
                                $cached_overlap_result = $overlaps;
                                $cached_start = $region->{start};
                                $cached_end = $region->{end};

                                #update score
                                $score = $len/100000;
                                $curr_len = $len;
                                $anchor_ref_pos = $current_ref_pos + $len;
                                $anchor_read_pos = $current_read_pos + $len;
                            }
                        }
                    }
                }
            }
        }

        # Adjust current_ref_pos and current_read_pos as we move forward
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'D' || $op eq 'N') {
            $current_ref_pos += $len;
        }
        if ($op eq 'M' || $op eq '=' || $op eq 'X' || $op eq 'I' || $op eq 'S') {
            $current_read_pos += $len;
        }
    }
    
    if ($score > 0){
        return ($breakpoint_ref, $anchor_ref_pos + 1, $breakpoint_read, $anchor_read_pos + 1, $integrity, $score);
    }else{
        return ($breakpoint_ref, $breakpoint_ref, $breakpoint_read, $breakpoint_read, "RightCover", ".");
    }    

}

# calculate path for duplication
sub find_sv_paths {
    my ($alignments) = @_;

    my (@ref_path, @read_path);

    # Iterate through all alignments and collect ref and read paths
    for my $align (@$alignments) {
        my $ref_segment  = $align->{ref_start} . "-" . $align->{ref_end};
        my $read_segment = $align->{read_start} . "-" . $align->{read_end};

        push @ref_path, $ref_segment;
        push @read_path, $read_segment;
    }

    # Join the paths with '->'
    my $ref_path_str  = join "->", @ref_path;
    my $read_path_str = join "->", @read_path;

    return ($ref_path_str, $read_path_str);
}

#calculate path for SCGR
sub find_SCGR_sv_paths {
    my ($alignments) = @_;

    my (@ref_path, @read_path);

    # Iterate through all alignments and collect ref and read paths
    for my $align (@$alignments) {
        my $ref_segment  = $align->{ref_chr} . ":" . $align->{ref_start} . "-" . $align->{ref_end};  # Include chromosome info for reference
        my $read_segment = $align->{read_start} . "-" . $align->{read_end};  # Read segment does not need chromosome

        push @ref_path, $ref_segment;
        push @read_path, $read_segment;
    }

    # Join the paths with '->'
    my $ref_path_str  = join "->", @ref_path;
    my $read_path_str = join "->", @read_path;

    return ($ref_path_str, $read_path_str);
}

# Usage help function
sub usage {
    print <<END;
Usage: perl $0 --bam <sorted BAM file by genome coordinate> --bamSortRead <sorted BAM file by read name> --output <SV signature table> [options]

Options:
    --bam                                   Input BAM file by genome coordinate (required)
    --bamSortRead                           Input BAM file by read name (required)
    --output                                Output SV signature table file, and be compressed (required)
    --min-mapq                              Minimum mapping quality (default: 20)
    --min-align-len                         Minimum alignment length (default: 1000 bp)
    --min-sv-len                            Minimum length for structural variants (SVs) (default: 50 bp)
    --max-sv-len                            Maximum length for structural variants (default: 10,000 bp)
    --mitochondria-chr                      Sequence name for mitochondrial chromosome to exclude from analysis (optional)
    --anchor-len-insertion                  Minimum anchor length for insertion (default: 230 bp)
    --anchor-len-deletion                   Minimum anchor length for deletion (default: 290 bp)
    --anchor-len-CGR                        Minimum anchor length for CGRs (complex genomic rearrangements) (default: 260 bp)
    --anchor-len-duplication                Minimum anchor length for duplication (default: 210 bp)
    --anchor-len-complexInsertion           Minimum anchor length for complex insertion in which split alignment exists (default: 230 bp)
    --anchor-len-complexDeletion            Minimum anchor length for complex deletion in which split alignment exists (default: 290 bp)
    --anchor-len-SCGR                       Minimum anchor length for SCGRs (super CGRs) (default: 340 bp)
    --anchor-minDepth                       Minimum anchor depth rattio to overall mean (default: 0.33)
    --anchor-maxDepth                       Maximum anchor depth rate to overall mean for insertion (default: 3.00)       
    --help                                  Show this help message and exit

Examples:
    perl $0 --bam sample.bam --bamSortRead sample.sortRead.bam --output sv_signature_table --min-mapq 20 --min-sv-len 50 --max-sv-len 10000

END
}
