#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError);
use Data::Dumper;

# Define variables and default values
my ($input_file, $output_file, 
    $max_distance_insertion, $max_distance_deletion, $max_distance_CGR, $max_distance_duplication, $max_distance_SCGR,
    $min_reads_support_insertion, $min_reads_support_deletion, $min_reads_support_CGR, $min_reads_support_SCGR, $min_reads_support_duplication,
    $weight_complexINSDEL, $weight_unCover, $block_gap, $node_max_distance, $intermediate_output_file, $help);

# Default values for parameters (Insertion, Deletion, CGR, Duplication, SCGR)
$max_distance_insertion      = 130;
$max_distance_deletion       = 70;
$max_distance_CGR            = 490;  # Example default for CGR
$max_distance_duplication    = 220;  # Example default for Duplication
$max_distance_SCGR           = 500;  # Example default for SCGR

$min_reads_support_insertion      = 2;
$min_reads_support_deletion       = 2;
$min_reads_support_CGR            = 13; # Example default for CGR
$min_reads_support_duplication    = 18; # Example default for Duplication
$min_reads_support_SCGR           = 18; # Example default for SCGR

$weight_complexINSDEL = 1;
$weight_unCover = 0.5;
$block_gap = 500;
$node_max_distance = 10;

# Get command line options including help
GetOptions(
    "input=s"                         => \$input_file,
    "output=s"                        => \$output_file,
    "max-distance-insertion=i"        => \$max_distance_insertion,
    "max-distance-deletion=i"         => \$max_distance_deletion,
    "max-distance-CGR=i"              => \$max_distance_CGR,
    "max-distance-duplication=i"      => \$max_distance_duplication,
    "max-distance-SCGR=i"             => \$max_distance_SCGR,
    "min-reads-support-insertion=i"   => \$min_reads_support_insertion,
    "min-reads-support-deletion=i"    => \$min_reads_support_deletion,
    "min-reads-support-CGR=i"         => \$min_reads_support_CGR,
    "min-reads-support-duplication=i" => \$min_reads_support_duplication,
    "min-reads-support-SCGR=i"        => \$min_reads_support_SCGR,
    "weight-complexINSDEL=f"          => \$weight_complexINSDEL,
    "weight-unCover=f"                => \$weight_unCover,
    "block-gap=i"                     => \$block_gap,
    "node-max-distance=i"             => \$node_max_distance, 
    "intermediate-output=s"           => \$intermediate_output_file,
    "help"                            => \$help,
) or usage();

# Show usage if --help is provided or required options are missing
if ($help || !$input_file || !$output_file) {
    usage();
    exit;
}

# Open input file
my $in_fh;
if ($input_file =~ /\.gz$/) {
    open($in_fh, "gzip -dc $input_file |") or die "Could not open compressed input file: $input_file\n";
} else {
    open($in_fh, '<', $input_file) or die "Could not open input file: $input_file\n";
}

# Hash to track the min and max path points for each chromosome, ref_start, and ref_end
my %node_bounds = ();

# Define separate hash tables for each SV type
my %insertion_sv_signals;
my %deletion_sv_signals;
my %CGR_sv_signals;
my %duplication_sv_signals;
my %SCGR_sv_signals;

# Hash to store already processed SVs based on ref_start, ref_end, and sv_type
my %seen_sv = ();

my %countSig; #count SVs signatures for each type

# Read the input file and distribute SV signals based on their type
while (<$in_fh>) {
    chomp;  # Remove newline character
    next if (/Chromosome/);  # Skip the header line
    my @fields = split /\t/;  # Split the line into fields based on tab character

    # Unpack fields into variables
    my ($chromosome, $anchor_ref_start, $anchor_ref_end, $ref_start, $ref_end, $read_name, 
        $anchor_read_start, $anchor_read_end, $read_start, $read_end, $sv_type, $strand, 
        $integrity, $left_anchor_score, $right_anchor_score, $ref_path, $read_path, $origin_seq, $anchor_seq) = @fields;

    my $support = 1;  # Initialize support count to 1 for each line

    # Adjust support based on SV type and integrity
    if ($sv_type eq "ComplexDeletion") {
        $sv_type = "Deletion";  # Simplify SV type for processing
        $support *= $weight_complexINSDEL;  # Adjust support value
    } elsif ($sv_type eq "ComplexInsertion") {
        $sv_type = "Insertion";  # Simplify SV type for processing
        $support *= $weight_complexINSDEL;  # Adjust support value
    }

    # Modify support based on integrity values
    if ($integrity =~ /leftUncoverRightUncover/) {
        $support *= $weight_unCover * $weight_unCover;  # Increase support based on conditions
    } elsif ($integrity =~ /leftUncover/ || $integrity =~ /RightUncover/) {
        $support *= $weight_unCover;  # Increase support based on conditions
    }
    
    $countSig{$sv_type}++;
    
    my $score;
    if ($integrity =~ /leftSpanRightSpan/){
        $score = $left_anchor_score + $right_anchor_score;
    }elsif ($integrity =~ /leftSpan/) {
        $score = $left_anchor_score;
    }elsif ($integrity =~ /RightSpan/){
        $score = $right_anchor_score;
    }else{
        $score = ".";
    }

    # Duplication-specific merging logic
    if ($sv_type eq "Duplication") {
        my $merged = 0;  # Control variable to track whether a merge happened
        
        my $curr_sv = { # Initialize the current SV as the current line being processed (structure is correct)
            chromosome      => $chromosome,
            sv_type         => $sv_type,
            anchor_ref_start  => $anchor_ref_start,
            anchor_ref_end    => $anchor_ref_end,
            ref_start       => $ref_start,
            ref_end         => $ref_end,
            support         => $support,
            read_name       => $read_name,
            anchor_read_start => $anchor_read_start,
            anchor_read_end   => $anchor_read_end,
            read_start      => $read_start,
            read_end        => $read_end,
            strand          => $strand,
            integrity       => $integrity,
            ref_path        => $ref_path,
            read_path       => $read_path,
            origin_seq      => $origin_seq,
            anchor_seq        => $anchor_seq,
            score           => $score,
        };
    
        # Create initial sv_key
        my $sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name});
    
        my $continue_merging = 1;  # Control variable for the merge loop
    
        my $control = 0; # Initialize control counter for loop iterations
        
        # For the first merge, initialize bounds using initialize_bounds_for_nodes
        my @nodes = split /->/, $curr_sv->{ref_path};
        my @min_max_merge_bounds = initialize_bounds_for_nodes(\@nodes);
        
        #print "Starting merging loop...\n";
    
        while ($continue_merging) {
            $continue_merging = 0;  # Assume no further merges will happen
    
            # 每次循环使用最新的 curr_sv 更新 nodes 和 bounds
            $control++;
            #print "Control: $control\n";
    
            if ($control > 1) {
                # For subsequent merges, use the existing node bounds from curr_sv
                if (exists $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}) {
                    @min_max_merge_bounds = @{$node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}};
                    $sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name}); 
                } else {
                    print "Warning: No bounds found for current SV: $curr_sv->{read_name}\n";         
                }
            }
    
            # Check node bounds and perform merging if possible
            if (exists $node_bounds{$chromosome}{$ref_start}{$ref_end}) {
                #my $bound_count = scalar keys %{$node_bounds{$chromosome}{$ref_start}{$ref_end}};
                #print "Number of node bounds for $chromosome:$ref_start-$ref_end: $bound_count\n";
                foreach my $existing_read_name (keys %{$node_bounds{$chromosome}{$ref_start}{$ref_end}}) {
                    #print "Checking existing path: $existing_read_name\t$seen_sv{$sv_key}\n";
                    
                    my $existing_sv_key = join("_", $chromosome, $ref_start, $ref_end, $existing_read_name);
                    my $existing_sv = $seen_sv{$existing_sv_key};  
                    next unless (exists $seen_sv{$existing_sv_key});
                    next if ($existing_sv_key eq $sv_key);
                    my @min_max_target_bounds = @{$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}};
    
                    # Ensure node counts match
                    if (@nodes != @min_max_merge_bounds || @nodes != @min_max_target_bounds) {
                        next;  # Skip this iteration if counts don't match
                    }
    
                    #print "当前合并节点范围: " . Dumper(\@min_max_merge_bounds) . "\n";
                    #print "目标节点范围: " . Dumper(\@min_max_target_bounds) . "\n";
    
                    # Check if nodes can be merged
                    my $merge;
                    if (@nodes == 2) {
                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                    } elsif (@nodes > 2) {
                        $merge = can_merge_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance) if $merge;
                    }
    
                    if ($merge) {
                        #print "Merging successful.\n";
                        # 更新节点边界
                        update_node_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                        #print "$existing_sv->{integrity}\t$curr_sv->{integrity}\n";
    
                        # 删除旧信号前进行检查，存在即删除
                        #print "Before deletion, duplication signals: ", scalar(@{$duplication_sv_signals{$chromosome}}), "\n";
                        @{$duplication_sv_signals{$chromosome}} = grep { 
                            $_->{read_name} ne $existing_read_name || $_->{ref_start} != $ref_start || $_->{ref_end} != $ref_end 
                        } @{$duplication_sv_signals{$chromosome}};
                        #print "After deletion of existing, duplication signals: ", scalar(@{$duplication_sv_signals{$chromosome}}), "\n";

                        @{$duplication_sv_signals{$chromosome}} = grep { 
                            $_->{read_name} ne $curr_sv->{read_name} || $_->{ref_start} != $ref_start || $_->{ref_end} != $ref_end 
                        } @{$duplication_sv_signals{$chromosome}};
                        #print "After deletion of current, duplication signals: ", scalar(@{$duplication_sv_signals{$chromosome}}), "\n";
    
                        # 更新合并逻辑
                        if ($existing_sv->{integrity} eq "leftSpanRightSpan" && $curr_sv->{integrity} eq "leftSpanRightSpan") {
                            # Case 1: Both have "leftSpanRightSpan"
                            my $curr_score = $curr_sv->{score};
                            my $existing_score = $existing_sv->{score};
                            if ($curr_score > $existing_score) {
                                $curr_sv->{support} += $existing_sv->{support};
                                delete $seen_sv{$existing_sv_key};
                                # 更新 node_bounds，删除旧的 existing_read_name
                                #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};        
                                $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}} = \@min_max_merge_bounds;
                            } else {
                                $existing_sv->{support} += $curr_sv->{support};
                                delete $seen_sv{$sv_key};
                                $curr_sv = $existing_sv; # 将 curr_sv 更新为 existing_sv
                                #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};  
                                $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name} = \@min_max_merge_bounds;
                            }
                        } elsif ($existing_sv->{integrity} eq "leftSpanRightSpan" || $curr_sv->{integrity} eq "leftSpanRightSpan") {
                            # Case 2: One has "leftSpanRightSpan"
                            if ($curr_sv->{integrity} eq "leftSpanRightSpan") {
                                $curr_sv->{support} += $existing_sv->{support};
                                delete $seen_sv{$existing_sv_key};
                                # 删除旧的 existing_read_name
                                #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};        
                                $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}} = \@min_max_merge_bounds;
                            } else {
                                $existing_sv->{support} += $curr_sv->{support};
                                delete $seen_sv{$sv_key};
                                $curr_sv = $existing_sv; # 将 curr_sv 更新为 existing_sv
                                #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}; 
                                $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name} = \@min_max_merge_bounds;
                            }
                        } else {
                            # Case 3: Neither has "leftSpanRightSpan"
                            $existing_sv->{support} += $curr_sv->{support}; 
                            delete $seen_sv{$sv_key};
                            $curr_sv = $existing_sv; # 更新 curr_sv 为 existing_sv
                            #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};                            
                            $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name} = \@min_max_merge_bounds;
                        }
    
                        # 将合并后的 curr_sv 推回信号列表
                        push @{$duplication_sv_signals{$chromosome}}, $curr_sv;
                        #print "After push, duplication signals: ", scalar(@{$duplication_sv_signals{$chromosome}}), "\n";
                        $sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name});
                        $seen_sv{$sv_key} = $curr_sv;
    
                        # **移除已经合并的节点**
                        #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};                 
                        #delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};              
                        $continue_merging = 1;
                        $merged = 1;
                        last;  # 退出 foreach，回到 while 循环
                    }
                }
            }
    
            # Break the loop if no further merging is needed
            if ($continue_merging == 0) {
                #print "No further merging possible, exiting loop.\n";
                last;  # Exit the loop for further merging
            }
        }
    
        # If no merging occurred, store as new duplication SV
        if (!$merged) {
            $node_bounds{$chromosome}{$ref_start}{$ref_end}{$read_name} = \@min_max_merge_bounds;
            $seen_sv{$sv_key} = $curr_sv;
            #print "$seen_sv{$sv_key}\n";
            push @{$duplication_sv_signals{$chromosome}}, $seen_sv{$sv_key};
        }
    # 处理 SCGR 类型的 SV 合并
    } elsif ($sv_type eq "SCGR") {
        my $merged = 0;  # 记录是否发生合并
    
        my $curr_sv = {
            chromosome      => $chromosome,
            sv_type         => $sv_type,
            anchor_ref_start  => $anchor_ref_start,
            anchor_ref_end    => $anchor_ref_end,
            ref_start       => $ref_start,
            ref_end         => $ref_end,
            support         => $support,
            read_name       => $read_name,
            anchor_read_start => $anchor_read_start,
            anchor_read_end   => $anchor_read_end,
            read_start      => $read_start,
            read_end        => $read_end,
            strand          => $strand,
            integrity       => $integrity,
            ref_path        => $ref_path,
            read_path       => $read_path,
            origin_seq      => $origin_seq,
            anchor_seq        => $anchor_seq,
            score           => $score,
        };
        
        # Create initial sv_key
        my $sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name},$curr_sv->{ref_path});
    
        # 提取 ref_path 中的节点 (chr:pos-range)，并将 pos-range 部分用于节点的边界初始化
        my @nodes = split /->/, $curr_sv->{ref_path};
    
        foreach my $index (0..$#nodes) {
            my ($chr, $path) = split/:/,$nodes[$index];
            my ($start, $end) = split /-/, $path;
            if ($start == 1) { #a weired bug: if start is equal 1, the warnning in line 358 will ocuur and the corresponding loop will not stop
                $start = 2;
            }
            my $new_path = $chr.":".$start."-".$end;
        }
    
        my $continue_merging = 1;  # Control variable for the merge loop
        
        my $control = 0; # Initialize control counter for loop iterations
    
        # 基于 range 初始化节点边界
        my @min_max_merge_bounds = initialize_bounds_SCGR_for_nodes(\@nodes); 
    
        
        # SCGR 合并逻辑
        while ($continue_merging) {
            $continue_merging = 0;  # 假设没有进一步的合并
    
            $control++;
            if ($control > 1) {
                # For subsequent merges, use the existing node bounds from curr_sv
                if (exists $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}{$curr_sv->{ref_path}}) {
                    @min_max_merge_bounds = @{$node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}{$curr_sv->{ref_path}}};
                    $sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name}, $curr_sv->{ref_path}); 
                } else {
					print "$control\n";
                    print "Warning: No bounds found for current SV: $curr_sv->{read_name}\n";         
                }
            }
        
            # Check node bounds and perform merging if possible
            if (exists $node_bounds{$chromosome}{$ref_start}{$ref_end}){
                foreach my $existing_read_name (keys %{$node_bounds{$chromosome}{$ref_start}{$ref_end}}) {
					foreach my $existing_ref_path (keys %{$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}}){
         
						my $existing_sv_key = join("_", $chromosome, $ref_start, $ref_end, $existing_read_name, $existing_ref_path);
						my $existing_sv = $seen_sv{$existing_sv_key};  
						next unless (exists $seen_sv{$existing_sv_key});
						next if ($existing_sv_key eq $sv_key);
						my @min_max_target_bounds = @{$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}{$existing_ref_path}};
        
						# Ensure node counts match
						if (@nodes != @min_max_merge_bounds || @nodes != @min_max_target_bounds) {
							next;  # Skip this iteration if counts don't match
						}
    
						# Check if nodes can be merged
						my $merge = can_merge_SCGR_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
    
        
						if ($merge) {
							update_node_SCGR_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
        
							#print "Before deletion, SCGR signals: ", scalar(@{$SCGR_sv_signals{$chromosome}}), "\n";
							@{$SCGR_sv_signals{$chromosome}} = grep { 
								$_->{read_name} ne $existing_read_name && $_->{ref_path} ne $existing_ref_path && $_->{ref_start} != $ref_start && $_->{ref_end} != $ref_end 
							} @{$SCGR_sv_signals{$chromosome}};
							#print "After deletion of existing, SCGR signals: ", scalar(@{$SCGR_sv_signals{$chromosome}}), "\n";
    
							@{$SCGR_sv_signals{$chromosome}} = grep { 
								$_->{read_name} ne $curr_sv->{read_name} && $_->{ref_path} ne $curr_sv->{ref_path} && $_->{ref_start} != $ref_start && $_->{ref_end} != $ref_end 
							} @{$SCGR_sv_signals{$chromosome}};
							#print "After deletion of current, SCGR signals: ", scalar(@{$SCGR_sv_signals{$chromosome}}), "\n";
        
							# 更新合并逻辑
							if ($existing_sv->{integrity} eq "leftSpanRightSpan" && $curr_sv->{integrity} eq "leftSpanRightSpan") {
								# Case 1: Both have "leftSpanRightSpan"
								my $curr_score = $curr_sv->{score};
								my $existing_score = $existing_sv->{score};
								if ($curr_score > $existing_score) {
									$curr_sv->{support} += $existing_sv->{support};
									delete $seen_sv{$existing_sv_key};
									# 更新 node_bounds，删除旧的 existing_read_name
									#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};        
									$node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}{$curr_sv->{ref_path}} = \@min_max_merge_bounds;
								} else {
									$existing_sv->{support} += $curr_sv->{support};
									delete $seen_sv{$sv_key};
									$curr_sv = $existing_sv; # 将 curr_sv 更新为 existing_sv
									#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};  
									$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}{$existing_ref_path} = \@min_max_merge_bounds;
								}
							} elsif ($existing_sv->{integrity} eq "leftSpanRightSpan" || $curr_sv->{integrity} eq "leftSpanRightSpan") {
								# Case 2: One has "leftSpanRightSpan"
								if ($curr_sv->{integrity} eq "leftSpanRightSpan") {
									$curr_sv->{support} += $existing_sv->{support};
									delete $seen_sv{$existing_sv_key};
									# 删除旧的 existing_read_name
									#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};        
									$node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}{$curr_sv->{ref_path}} = \@min_max_merge_bounds;
								} else {
									$existing_sv->{support} += $curr_sv->{support};
									delete $seen_sv{$sv_key};
									$curr_sv = $existing_sv; # 将 curr_sv 更新为 existing_sv
									#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}}; 
									$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}{$existing_ref_path} = \@min_max_merge_bounds;
								}
							}else {
								# Case 3: Neither has "leftSpanRightSpan"
								$existing_sv->{support} += $curr_sv->{support}; 
								delete $seen_sv{$sv_key};
								$curr_sv = $existing_sv; # 更新 curr_sv 为 existing_sv
								#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};                            
								$node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name}{$existing_ref_path} = \@min_max_merge_bounds;
							}
        
							# 将合并后的 curr_sv 推回信号列表
							push @{$SCGR_sv_signals{$chromosome}}, $curr_sv;
							#print "After push, SCGR signals: ", scalar(@{$SCGR_sv_signals{$chromosome}}), "\n";
							$sv_key = join("_", $chromosome, $ref_start, $ref_end, $curr_sv->{read_name}, $curr_sv->{ref_path});
							$seen_sv{$sv_key} = $curr_sv;
        
							# **移除已经合并的节点**
							#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$existing_read_name};                 
							#delete $node_bounds{$chromosome}{$ref_start}{$ref_end}{$curr_sv->{read_name}};              
							$continue_merging = 1;
							$merged = 1;
							last;  # 退出 foreach，回到 while 循环
						}
                    }
                }
            }
        
            # Break the loop if no further merging is needed
            if ($continue_merging == 0) {
                #print "No further merging possible, exiting loop.\n";
                last;  # Exit the loop for further merging
            }
        }
        
        # If no merging occurred, store as new SCGR SV
        if (!$merged) {
            $node_bounds{$chromosome}{$ref_start}{$ref_end}{$read_name}{$ref_path} = \@min_max_merge_bounds;
            $seen_sv{$sv_key} = $curr_sv;
            #print "$seen_sv{$sv_key}\n";
            push @{$SCGR_sv_signals{$chromosome}}, $seen_sv{$sv_key};
        }
    }else {
        # General strategy for non-Duplication types
        my $sv_key = join("_", $chromosome, $ref_start, $ref_end, $sv_type);

        if (exists $seen_sv{$sv_key}) {
            my $existing_sv = $seen_sv{$sv_key};

            # Handle SVs with leftSpanRightSpan integrity
            if ($existing_sv->{integrity} eq "leftSpanRightSpan" && $integrity eq "leftSpanRightSpan") {
                my $curr_score = $score;
                my $existing_score = $existing_sv->{score};
                if ($curr_score > $existing_score) {
                    # Keep the longer anchor_seq and update other fields
                    $support += $existing_sv->{support};
                    $seen_sv{$sv_key} = {
                        chromosome      => $chromosome,
                        sv_type         => $sv_type,
                        anchor_ref_start  => $anchor_ref_start,
                        anchor_ref_end    => $anchor_ref_end,
                        ref_start       => $ref_start,
                        ref_end         => $ref_end,
                        support         => $support,
                        read_name       => $read_name,
                        anchor_read_start => $anchor_read_start,
                        anchor_read_end   => $anchor_read_end,
                        read_start      => $read_start,
                        read_end        => $read_end,
                        strand          => $strand,
                        integrity       => $integrity,
                        ref_path        => $ref_path,
                        read_path       => $read_path,
                        origin_seq      => $origin_seq,
                        anchor_seq        => $anchor_seq,  # Keep the longest anchor_seq
                        score           => $score,
                    };
                } else {
                    # Only update support if the existing anchor_seq is longer
                    $existing_sv->{support} += $support;
                    $seen_sv{$sv_key} = $existing_sv;
                }
            } elsif ($existing_sv->{integrity} eq "leftSpanRightSpan" || $integrity eq "leftSpanRightSpan") {
                if ($integrity eq "leftSpanRightSpan") {
                    # Prefer the current SV with leftSpanRightSpan
                    $support += $existing_sv->{support};
                    $seen_sv{$sv_key} = {
                        chromosome      => $chromosome,
                        sv_type         => $sv_type,
                        anchor_ref_start  => $anchor_ref_start,
                        anchor_ref_end    => $anchor_ref_end,
                        ref_start       => $ref_start,
                        ref_end         => $ref_end,
                        support         => $support,
                        read_name       => $read_name,
                        anchor_read_start => $anchor_read_start,
                        anchor_read_end   => $anchor_read_end,
                        read_start      => $read_start,
                        read_end        => $read_end,
                        strand          => $strand,
                        integrity       => $integrity,
                        ref_path        => $ref_path,
                        read_path       => $read_path,
                        origin_seq      => $origin_seq,
                        anchor_seq        => $anchor_seq,
                        score            => $score,
                    };
                } else {
                    # Update support only
                    $existing_sv->{support} += $support;
                    $seen_sv{$sv_key} = $existing_sv;
                }
            } else {
                # Combine support if neither has leftSpanRightSpan
                $existing_sv->{support} += $support;
                $seen_sv{$sv_key} = $existing_sv;
            }

            # Update final SV data based on SV type
            if ($sv_type eq "Insertion") {
                @{$insertion_sv_signals{$chromosome}} = grep { $_->{ref_start} ne $ref_start || $_->{ref_end} ne $ref_end } @{$insertion_sv_signals{$chromosome}};
                push @{$insertion_sv_signals{$chromosome}}, $seen_sv{$sv_key};
            } elsif ($sv_type eq "Deletion") {
                @{$deletion_sv_signals{$chromosome}} = grep { $_->{ref_start} ne $ref_start || $_->{ref_end} ne $ref_end } @{$deletion_sv_signals{$chromosome}};
                push @{$deletion_sv_signals{$chromosome}}, $seen_sv{$sv_key};
            } elsif ($sv_type eq "CGR") {
                @{$CGR_sv_signals{$chromosome}} = grep { $_->{ref_start} ne $ref_start || $_->{ref_end} ne $ref_end } @{$CGR_sv_signals{$chromosome}};
                push @{$CGR_sv_signals{$chromosome}}, $seen_sv{$sv_key};
            }
        } else {
            # If this SV is new, store it in the appropriate SV type hash
            my $sv_data = {
                chromosome      => $chromosome,
                sv_type         => $sv_type,
                anchor_ref_start  => $anchor_ref_start,
                anchor_ref_end    => $anchor_ref_end,
                ref_start       => $ref_start,
                ref_end         => $ref_end,
                support         => $support,
                read_name       => $read_name,
                anchor_read_start => $anchor_read_start,
                anchor_read_end   => $anchor_read_end,
                read_start      => $read_start,
                read_end        => $read_end,
                strand          => $strand,
                integrity       => $integrity,
                ref_path        => $ref_path,
                read_path       => $read_path,
                origin_seq      => $origin_seq,
                anchor_seq        => $anchor_seq,
                score           => $score,
            };

            $seen_sv{$sv_key} = $sv_data;  # Store SV

            # Store in the appropriate SV type array
            if ($sv_type eq "Insertion") {
                push @{$insertion_sv_signals{$chromosome}}, $sv_data;
            } elsif ($sv_type eq "Deletion") {
                push @{$deletion_sv_signals{$chromosome}}, $sv_data;
            } elsif ($sv_type eq "CGR") {
                push @{$CGR_sv_signals{$chromosome}}, $sv_data;
            }
        }
    }
}

# Process each SV type separately
my @sv_types = (
    { type => 'Insertion', signals => \%insertion_sv_signals },
    { type => 'Deletion', signals => \%deletion_sv_signals },
    { type => 'CGR', signals => \%CGR_sv_signals },
    { type => 'SCGR', signals => \%SCGR_sv_signals },
    { type => 'Duplication', signals => \%duplication_sv_signals }
);

# 初始化存储所有 SV 数据的数组
my @all_sv_data;

my %countIntermediate; #count the SVs for removing redundancy the same chromosome, start, end, type and path
# 处理每个 SV 类型，汇总数据
foreach my $sv_type_ref (@sv_types) {
    my $sv_signals = $sv_type_ref->{signals};
    my $sv_type = $sv_type_ref->{type};
    
    $countIntermediate{$sv_type} = 0;
    
    foreach my $chromosome (sort keys %$sv_signals) {
        foreach my $sv (@{$sv_signals->{$chromosome}}) {
            # 将 SV 数据加入到全局数组中，用于排序
            push @all_sv_data, {
                chromosome      => $sv->{chromosome},
                anchor_ref_start  => $sv->{anchor_ref_start},
                anchor_ref_end    => $sv->{anchor_ref_end},
                ref_start       => $sv->{ref_start},
                ref_end         => $sv->{ref_end},
                read_name       => $sv->{read_name},
                sv_type         => $sv->{sv_type},
                integrity       => $sv->{integrity},
                support         => $sv->{support},
                strand          => $sv->{strand},
                ref_path        => $sv->{ref_path},
                read_path       => $sv->{read_path},
                anchor_seq      => $sv->{anchor_seq},
                score           => $sv->{score},
            };
            $countIntermediate{$sv_type}++;
        }
    }
}

# 对所有 SV 数据按排序规则排序
@all_sv_data = sort {
    my ($a_chr, $a_start) = ($a->{chromosome}, $a->{ref_start});
    my ($b_chr, $b_start) = ($b->{chromosome}, $b->{ref_start});
    my ($a_num) = $a_chr =~ /(\d+)/ ? $1 : 999;
    my ($b_num) = $b_chr =~ /(\d+)/ ? $1 : 999;
    my $cmp = $a_num <=> $b_num;
    return $cmp || $a_start <=> $b_start;
} @all_sv_data;

if ($intermediate_output_file){
    # 打开中间结果文件为压缩文件
    my $intermediate_gz = IO::Compress::Gzip->new("$intermediate_output_file.gz") 
        or die "Error: Cannot open intermediate output file $intermediate_output_file.gz: $GzipError\n";

    # 中间文件输出头
    print $intermediate_gz "Chromosome\tanchor_Ref_Start\tanchor_Ref_End\tRef_Start\tRef_End\tRead_Name\tSV_Type\tSupport\tStrand\tRef_Path\tRead_Path\tanchor_Sequence\n";

    # 将排序后的数据输出到中间文件
    foreach my $sv (@all_sv_data) {
        print $intermediate_gz join("\t", 
            $sv->{chromosome},
            $sv->{anchor_ref_start},
            $sv->{anchor_ref_end},
            $sv->{ref_start},
            $sv->{ref_end},
            $sv->{read_name},
            $sv->{sv_type},
            $sv->{integrity},
            $sv->{support},
            $sv->{strand},
            $sv->{score},
            $sv->{ref_path},
            $sv->{read_path},
            $sv->{anchor_seq},
        ), "\n";
    }

    # 关闭中间文件句柄
    close($intermediate_gz);
}

# Open the gzip output file for writing
my $gz = IO::Compress::Gzip->new("$output_file.gz") or die "Error: Cannot open output file $output_file.gz: $GzipError\n";

# Store all output lines in an array for later sorting
my @output_lines;
print $gz "Chromosome\tanchor_Ref_Start\tanchor_Ref_End\tRef_Start\tRef_End\tSV_Type\tSupport\tRead_Name\tanchor_Sequence\n";

my %countSVcluster; #count the SVs that were finally clustered.

print "input file read complete!\n";
# Process each SV type
foreach my $sv_type_ref (@sv_types) {
    my $sv_type_name = $sv_type_ref->{type};  # Get the SV type name
    my $sv_count = 0;
   
    my ($max_distance, $min_reads_support);
    if ($sv_type_name eq "Insertion"){
        $max_distance = $max_distance_insertion;
        $min_reads_support = $min_reads_support_insertion;
    } elsif ($sv_type_name eq "Deletion"){
        $max_distance = $max_distance_deletion;
        $min_reads_support = $min_reads_support_deletion;
    } elsif ($sv_type_name eq "CGR"){
        $max_distance = $max_distance_CGR;
        $min_reads_support = $min_reads_support_CGR;
    } elsif ($sv_type_name eq "Duplication"){
        $max_distance = $max_distance_duplication;
        $min_reads_support = $min_reads_support_duplication;
    }elsif ($sv_type_name eq "SCGR"){
        $max_distance = $max_distance_SCGR;
        $min_reads_support = $min_reads_support_SCGR;
    }

    # Iterate over each chromosome in the SV signals
    for my $chromosome (sort custom_chromosome_sort keys %{$sv_type_ref->{signals}}) {
        my @sv_list = @{$sv_type_ref->{signals}{$chromosome}};  # Get the SV list for the chromosome
        @sv_list = sort { 
            $a->{ref_start} <=> $b->{ref_start} ||
            $a->{ref_end}   <=> $b->{ref_end} ||
            custom_chromosome_sort($a->{read_name}, $b->{read_name})
        } @sv_list;

        my @first_cluster;
        # Check if sv_list is not empty
        if (@sv_list) {
            # First-level clustering based on ref_start
            my @first_cluster = ($sv_list[0]) if @sv_list;                           
            my $first_cluster_max_ref_start = $sv_list[0]->{ref_start};  # Max ref_start of the current cluster

            for my $i (1 .. $#sv_list) {
                my $sv = $sv_list[$i];
				#print "$sv->{ref_start}\t$first_cluster_max_ref_start\n";

                if ($sv->{ref_start} - $first_cluster_max_ref_start > $block_gap) {
					#print "first cluster detected!\n";
                    if (any_left_span_right_span(@first_cluster)) {
                        if (sum_support(\@first_cluster) >= $min_reads_support) {

                            # Sort the large cluster by ref_end to define smaller clusters
                            @first_cluster = sort { 
                            $a->{ref_end} <=> $b->{ref_end}      # First, sort by ref_end
                                ||
                            $a->{ref_start}   <=> $b->{ref_start}        # If ref_start is the same, sort by ref_start
                                ||
                            custom_chromosome_sort($a->{read_name}, $b->{read_name})  # Sort by read_name using custom_chromosome_sort
                            } @first_cluster;    
							#print "first cluster order!\n";
                            my @second_cluster = ($first_cluster[0]);                    
                            my $second_cluster_max_ref_end = $first_cluster[0]->{ref_end};

                            for my $i (1 .. $#first_cluster) {
                                my $sv = $first_cluster[$i]; 
                                
                                if ($sv->{ref_end} - $second_cluster_max_ref_end > $block_gap) {
                                    if (any_left_span_right_span(@second_cluster)) {
                                        if (sum_support(\@second_cluster) >= $min_reads_support) {
                                            
                                            my @third_cluster = ();
											#print "Process third cluster\n";
                                            # Function to cluster SV signals within a block and collect results for sorting

                                            # First, filter SV list to only include those with leftSpanRightSpan integrity
                                            my @filtered_sv_list = grep { $_->{integrity} =~ /leftSpanRightSpan/ && $_->{score} > 1 }@second_cluster;

                                            # Skip if no SV with leftSpanRightSpan integrity is found
                                            if (!@filtered_sv_list) {
                                                print "No SVs with 'leftSpanRightSpan' integrity found. Skipping clustering.\n";
                                            }

                                            # Sort SVs by support, then by length of anchor_seq, ref_start, ref_end, and read_name
                                            @filtered_sv_list = sort {
                                                $b->{support} <=> $a->{support} ||
                                                $b->{score} <=> $a->{score} ||
                                                $a->{ref_start} <=> $b->{ref_start} ||
                                                $a->{ref_end} <=> $b->{ref_end} ||
                                                $a->{read_name} cmp $b->{read_name}
                                            } @filtered_sv_list;

                                            # Loop through the filtered SVs to cluster them
                                            while (@filtered_sv_list) {
                                                # Take the first SV as the initial candidate for clustering
                                                my $peak_sv = shift @filtered_sv_list;

                                                # If $peak_sv is not defined, skip processing
                                                next unless defined $peak_sv;

                                                # Initialize peak_start and peak_end correctly after selecting peak_sv
                                                my $peak_start = $peak_sv->{ref_start};
                                                my $peak_end = $peak_sv->{ref_end};

                                                # Further refine the cluster to include only those within max_distance of both start and end
                                                my @final_cluster = grep {
                                                    abs($_->{ref_start} - $peak_start) <= $max_distance &&
                                                    abs($_->{ref_end} - $peak_end) <= $max_distance
                                                } @second_cluster;

                                                # Skip empty clusters
                                                next unless @final_cluster;

                                                 # Duplication-specific logic
                                                if ($sv_type_name eq "Duplication") {
                                                    #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                                    # Initialize the final duplication cluster to include peak_sv
                                                    my @final_duplication_cluster = ($peak_sv);
                                                    my $continue_merging = 1;

                                                    while ($continue_merging) {
                                                        $continue_merging = 0;  # Assume no further merges will happen

                                                        foreach my $sv (@final_cluster) {
                                                            next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                            # Check if the SV can be merged with peak_sv based on node bounds
                                                            if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}) {                        
                                                                my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}};

                                                                # Ensure node counts are the same for merging
                                                                if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}) {
                                                                    my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}};

                                                                    # Skip if node counts are not equal
                                                                    #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                                    #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                                    # Ensure node counts match
                                                                    if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                        next;  # Skip this iteration if counts don't match
                                                                    }

                                                                    # 如果通过，继续输出
                                                                    #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                                    # Check if nodes can be merged
                                                                    my $merge;
                                                                    if (@min_max_merge_bounds == 2) {
                                                                        # Case for exactly two nodes: only check first and last nodes
                                                                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                    } elsif (@min_max_merge_bounds > 2) {
                                                                        # Case for more than two nodes: check intermediate and first-last nodes
                                                                        $merge = can_merge_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance) if $merge;
                                                                    }

                                                                    if ($merge) {
                                                                        # Update bounds after merging
                                                                        update_node_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                        @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}} = \@min_max_merge_bounds;
                                                                        #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                        # Add this SV to the final duplication cluster
                                                                        push @final_duplication_cluster, $sv;

                                                                        # Remove this SV from @final_cluster to avoid double-checking
                                                                        @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                         # Set flag to continue merging in the next iteration
                                                                        $continue_merging = 1;
                                                                        last;  # Exit loop and reprocess after merge
                                                                    }
                                                                }
                                                            }
                                                        }

                                                        # Ensure @final_duplication_cluster is not empty
                                                        if (!@final_duplication_cluster) {
                                                            print "Final duplication cluster is empty, skipping further processing.\n";
                                                            next;
                                                        }

                                                        # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                        @final_cluster = @final_duplication_cluster;
                                                    }

                                                    # Check if final_cluster is valid and non-empty
                                                    if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                        print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                        next;
                                                    }
                                                }

                                                # SCGR-specific logic
                                                if ($sv_type_name eq "SCGR") {
                                                    #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                                    # Initialize the final duplication cluster to include peak_sv
                                                    my @final_SCGR_cluster = ($peak_sv);
                                                    my $continue_merging = 1;

                                                    while ($continue_merging) {
                                                        $continue_merging = 0;  # Assume no further merges will happen

                                                        foreach my $sv (@final_cluster) {
                                                            next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                            # Check if the SV can be merged with peak_sv based on node bounds
                                                            if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}) {                        
                                                                my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}};

                                                                # Ensure node counts are the same for merging
                                                                if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}) {
                                                                    my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}};

                                                                    # Skip if node counts are not equal
                                                                    #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                                    #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                                    # Ensure node counts match
                                                                    if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                        next;  # Skip this iteration if counts don't match
                                                                    }

                                                                    # 如果通过，继续输出
                                                                    #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                                    # Check if nodes can be merged
                                                                    my $merge = can_merge_SCGR_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                       

                                                                    if ($merge) {
                                                                        # Update bounds after merging
                                                                        update_node_SCGR_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                        @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}} = \@min_max_merge_bounds;
                                                                        #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                        # Add this SV to the final duplication cluster
                                                                        push @final_SCGR_cluster, $sv;

                                                                        # Remove this SV from @final_cluster to avoid double-checking
                                                                        @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                         # Set flag to continue merging in the next iteration
                                                                        $continue_merging = 1;
                                                                        last;  # Exit loop and reprocess after merge
                                                                    }
                                                                }
                                                            }
                                                        }

                                                        # Ensure @final_duplication_cluster is not empty
                                                        if (!@final_SCGR_cluster) {
                                                            print "Final SCGR cluster is empty, skipping further processing.\n";
                                                            next;
                                                        }

                                                        # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                        @final_cluster = @final_SCGR_cluster;
                                                    }

                                                    # Check if final_cluster is valid and non-empty
                                                    if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                        print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                        next;
                                                    }
                                                }
                                                
                                                # Calculate the total support for the cluster
                                                 my $total_support = sum_support(\@final_cluster);

                                                # If the total support meets the minimum threshold, finalize the cluster
                                                if ($total_support >= $min_reads_support) {
                                                    # Update peak_sv with the total support
                                                     $peak_sv->{support} = $total_support;

                                                    # Use peak_sv as the representative SV
                                                    my $representative_sv = $peak_sv;

                                                    # Store the cluster's SVs and the representative SV
                                                    push @third_cluster, { representative_sv => $representative_sv, svs => \@final_cluster };

                                                    # Remove clustered SVs from the list (only remove those in @final_cluster)
                                                    @second_cluster = grep {
                                                        !in_cluster($_, \@final_cluster)
                                                    } @second_cluster;

                                                    # Also remove overlapping SVs from @filtered_sv_list
                                                    @filtered_sv_list = grep {
                                                        !in_cluster($_, \@final_cluster)
                                                    } @filtered_sv_list;

                                                    # After removing clustered SVs, restart the loop to process remaining SVs
                                                    redo if @filtered_sv_list;
                                                } else {
                                                    # If the total support does not meet the threshold, remove peak_sv from @filtered_sv_list
                                                    @filtered_sv_list = grep { $_ ne $peak_sv } @filtered_sv_list;

                                                    # Continue to process other SVs from @filtered_sv_list
                                                    redo if @filtered_sv_list;  # Restart loop if there are still SVs to process             
                                                    last;  # If no more SVs in filtered list, exit loop
                                                }
                                            }

                                            foreach my $cluster (@third_cluster) {
                                                my $representative_sv = $cluster->{representative_sv};

                                                # Ensure the representative SV contains the necessary fields
                                                my $chromosome = $representative_sv->{chromosome} // '';
                                                my $sv_type = $representative_sv->{sv_type} // '';
                                                my $anchor_ref_start = $representative_sv->{anchor_ref_start} // '';
                                                my $anchor_ref_end = $representative_sv->{anchor_ref_end} // '';
                                                my $ref_start = $representative_sv->{ref_start} // '';
                                                my $ref_end = $representative_sv->{ref_end} // '';
                                                my $read_name = $representative_sv->{read_name} // '';
                                                my $anchor_seq = $representative_sv->{anchor_seq} // '';
                                                my $support = $representative_sv->{support} // '';
                                                
                                                # Collect output line for later sorting
                                                push @output_lines, join("\t",
                                                    $chromosome,
                                                    $anchor_ref_start,
                                                    $anchor_ref_end,
                                                    $ref_start,
                                                    $ref_end,
                                                    $sv_type,
                                                    $support,
                                                    $read_name,
                                                    $anchor_seq
                                                );
                                                $countSVcluster{$sv_type}++;
                                            }
                                        }
                                    }
                                    @second_cluster = ();
                                    $second_cluster_max_ref_end = $sv->{ref_end};
                                }
                                push @second_cluster, $sv;
                            }
                            if (any_left_span_right_span(@second_cluster)) {
                                if (sum_support(\@second_cluster) >= $min_reads_support) {
                                    my @third_cluster = ();
                                    # Function to cluster SV signals within a block and collect results for sorting

                                    # First, filter SV list to only include those with leftSpanRightSpan integrity
                                    my @filtered_sv_list = grep { $_->{integrity} =~ /leftSpanRightSpan/ && $_->{score} > 1 }@second_cluster;

                                    # Skip if no SV with leftSpanRightSpan integrity is found
                                    if (!@filtered_sv_list) {
                                        print "No SVs with 'leftSpanRightSpan' integrity found. Skipping clustering.\n";
                                    }

                                    # Sort SVs by support, then by length of anchor_seq, ref_start, ref_end, and read_name
                                    @filtered_sv_list = sort {
                                        $b->{support} <=> $a->{support} ||
                                        $b->{score} <=> $a->{score} ||
                                        $a->{ref_start} <=> $b->{ref_start} ||
                                        $a->{ref_end} <=> $b->{ref_end} ||
                                        $a->{read_name} cmp $b->{read_name}
                                    } @filtered_sv_list;

                                    # Loop through the filtered SVs to cluster them
                                    while (@filtered_sv_list) {
                                        # Take the first SV as the initial candidate for clustering
                                        my $peak_sv = shift @filtered_sv_list;

                                        # If $peak_sv is not defined, skip processing
                                        next unless defined $peak_sv;

                                        # Initialize peak_start and peak_end correctly after selecting peak_sv
                                        my $peak_start = $peak_sv->{ref_start};
                                        my $peak_end = $peak_sv->{ref_end};

                                        # Further refine the cluster to include only those within max_distance of both start and end
                                        my @final_cluster = grep {
                                            abs($_->{ref_start} - $peak_start) <= $max_distance &&
                                            abs($_->{ref_end} - $peak_end) <= $max_distance
                                        } @second_cluster;

                                        # Skip empty clusters
                                        next unless @final_cluster;

                                         # Duplication-specific logic
                                        if ($sv_type_name eq "Duplication") {
                                            #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                            # Initialize the final duplication cluster to include peak_sv
                                            my @final_duplication_cluster = ($peak_sv);
                                            my $continue_merging = 1;

                                            while ($continue_merging) {
                                                $continue_merging = 0;  # Assume no further merges will happen

                                                foreach my $sv (@final_cluster) {
                                                    next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                    # Check if the SV can be merged with peak_sv based on node bounds
                                                    if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}) {                        
                                                        my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}};

                                                        # Ensure node counts are the same for merging
                                                        if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}) {
                                                            my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}};
                                                            # Skip if node counts are not equal

                                                            #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                            #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                            # Ensure node counts match
                                                            if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                next;  # Skip this iteration if counts don't match
                                                            }

                                                            # 如果通过，继续输出
                                                            #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                            # Check if nodes can be merged
                                                            my $merge;
                                                            if (@min_max_merge_bounds == 2) {
                                                                # Case for exactly two nodes: only check first and last nodes
                                                                $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                            } elsif (@min_max_merge_bounds > 2) {
                                                                # Case for more than two nodes: check intermediate and first-last nodes
                                                                $merge = can_merge_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance) if $merge;
                                                            }

                                                            if ($merge) {
                                                                # Update bounds after merging
                                                                update_node_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}} = \@min_max_merge_bounds;
                                                                #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                # Add this SV to the final duplication cluster
                                                                push @final_duplication_cluster, $sv;

                                                                # Remove this SV from @final_cluster to avoid double-checking
                                                                @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                 # Set flag to continue merging in the next iteration
                                                                $continue_merging = 1;
                                                                last;  # Exit loop and reprocess after merg3
                                                            }
                                                        }
                                                    }
                                                }

                                                # Ensure @final_duplication_cluster is not empty
                                                if (!@final_duplication_cluster) {
                                                    print "Final duplication cluster is empty, skipping further processing.\n";
                                                    next;
                                                }

                                                # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                @final_cluster = @final_duplication_cluster;
                                            }

                                            # Check if final_cluster is valid and non-empty
                                            if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                next;
                                            }
                                        }

                                        # SCGR-specific logic
                                        if ($sv_type_name eq "SCGR") {
                                            #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                            # Initialize the final duplication cluster to include peak_sv
                                            my @final_SCGR_cluster = ($peak_sv);
                                            my $continue_merging = 1;

                                            while ($continue_merging) {
                                                $continue_merging = 0;  # Assume no further merges will happen

                                                foreach my $sv (@final_cluster) {
                                                    next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                    # Check if the SV can be merged with peak_sv based on node bounds
                                                    if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}) {                        
                                                        my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}};

                                                        # Ensure node counts are the same for merging
                                                        if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}) {
                                                            my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}};

                                                            # Skip if node counts are not equal
                                                            #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                            #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                            # Ensure node counts match
                                                            if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                next;  # Skip this iteration if counts don't match
                                                            }

                                                            # 如果通过，继续输出
                                                            #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                            # Check if nodes can be merged
                                                            my $merge = can_merge_SCGR_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                       

                                                            if ($merge) {
                                                                # Update bounds after merging
                                                                update_node_SCGR_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}} = \@min_max_merge_bounds;
                                                                #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                # Add this SV to the final duplication cluster
                                                                push @final_SCGR_cluster, $sv;

                                                                # Remove this SV from @final_cluster to avoid double-checking
                                                                @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                # Set flag to continue merging in the next iteration
                                                                $continue_merging = 1;
                                                                last;  # Exit loop and reprocess after merge
                                                            }
                                                        }
                                                    }
                                                }

                                                # Ensure @final_duplication_cluster is not empty
                                                if (!@final_SCGR_cluster) {
                                                    print "Final SCGR cluster is empty, skipping further processing.\n";
                                                    next;
                                                }

                                                # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                @final_cluster = @final_SCGR_cluster;
                                            }

                                            # Check if final_cluster is valid and non-empty
                                            if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                next;
                                            }
                                        }

                                        # Calculate the total support for the cluster
                                         my $total_support = sum_support(\@final_cluster);

                                        # If the total support meets the minimum threshold, finalize the cluster
                                        if ($total_support >= $min_reads_support) {
                                            # Update peak_sv with the total support
                                             $peak_sv->{support} = $total_support;

                                            # Use peak_sv as the representative SV
                                            my $representative_sv = $peak_sv;

                                            # Store the cluster's SVs and the representative SV
                                            push @third_cluster, { representative_sv => $representative_sv, svs => \@final_cluster };

                                            # Remove clustered SVs from the list (only remove those in @final_cluster)
                                            @second_cluster = grep {
                                                !in_cluster($_, \@final_cluster)
                                            } @second_cluster;

                                            # Also remove overlapping SVs from @filtered_sv_list
                                            @filtered_sv_list = grep {
                                                !in_cluster($_, \@final_cluster)
                                            } @filtered_sv_list;

                                            # After removing clustered SVs, restart the loop to process remaining SVs
                                            redo if @filtered_sv_list;
                                        } else {
                                            # If the total support does not meet the threshold, remove peak_sv from @filtered_sv_list
                                            @filtered_sv_list = grep { $_ ne $peak_sv } @filtered_sv_list;

                                            # Continue to process other SVs from @filtered_sv_list
                                            redo if @filtered_sv_list;  # Restart loop if there are still SVs to process             
                                            last;  # If no more SVs in filtered list, exit loop
                                        }
                                    }

                                    foreach my $cluster (@third_cluster) {
                                        my $representative_sv = $cluster->{representative_sv};

                                        # Ensure the representative SV contains the necessary fields
                                        my $chromosome = $representative_sv->{chromosome} // '';
                                        my $sv_type = $representative_sv->{sv_type} // '';
                                        my $anchor_ref_start = $representative_sv->{anchor_ref_start} // '';
                                        my $anchor_ref_end = $representative_sv->{anchor_ref_end} // '';
                                        my $ref_start = $representative_sv->{ref_start} // '';
                                        my $ref_end = $representative_sv->{ref_end} // '';
                                        my $read_name = $representative_sv->{read_name} // '';
                                        my $anchor_seq = $representative_sv->{anchor_seq} // '';
                                        my $support = $representative_sv->{support} // '';
                                                

                                        # Calculate total support by summing up the support of all SVs in the cluster
                                        #my $total_support = sum_support($cluster->{svs});

                                        # Collect output line for later sorting
                                        push @output_lines, join("\t",
                                            $chromosome,
                                            $anchor_ref_start,
                                            $anchor_ref_end,
                                            $ref_start,
                                            $ref_end,
                                            $sv_type,
                                            $support,
                                            $read_name,
                                            $anchor_seq
                                        );
                                        $countSVcluster{$sv_type}++;
                                    }
                                }    
                            }
                        }    
                    }
                    @first_cluster = ();
                    $first_cluster_max_ref_start = $sv->{ref_start};  # Update max ref_start for the new cluster
                }
                push @first_cluster, $sv;
            }

            if (any_left_span_right_span(@first_cluster)) {
                if (sum_support(\@first_cluster) >= $min_reads_support) {
                    # Sort the large cluster by ref_end to define smaller clusters

                    @first_cluster = sort { 
                    $a->{ref_end} <=> $b->{ref_end}      # First, sort by ref_end
                                ||
                    $a->{ref_start}   <=> $b->{ref_start}        # If ref_start is the same, sort by ref_start
                                ||
                    custom_chromosome_sort($a->{read_name}, $b->{read_name})  # Sort by read_name using custom_chromosome_sort
                    } @first_cluster; 
					#print "last first cluster started\n";
    
                    my @second_cluster = ($first_cluster[0]);
                    my $second_cluster_max_ref_end = $first_cluster[0]->{ref_end};

                    for my $i (1 .. $#first_cluster) {
                        my $sv = $first_cluster[$i];

                        if ($sv->{ref_end} - $second_cluster_max_ref_end > $block_gap) {
                            if (any_left_span_right_span(@second_cluster)) {
                                if (sum_support(\@second_cluster) >= $min_reads_support) {
                                            
                                    my @third_cluster = ();
                                    # Function to cluster SV signals within a block and collect results for sorting

                                    # First, filter SV list to only include those with leftSpanRightSpan integrity
                                    my @filtered_sv_list = grep { $_->{integrity} =~ /leftSpanRightSpan/ && $_->{score} > 1 }@second_cluster;

                                    # Skip if no SV with leftSpanRightSpan integrity is found
                                    if (!@filtered_sv_list) {
                                        print "No SVs with 'leftSpanRightSpan' integrity found. Skipping clustering.\n";
                                    }

                                    # Sort SVs by support, then by length of anchor_seq, ref_start, ref_end, and read_name
                                    @filtered_sv_list = sort {
                                        $b->{support} <=> $a->{support} ||
                                        $b->{score} <=> $a->{score} ||
                                        $a->{ref_start} <=> $b->{ref_start} ||
                                        $a->{ref_end} <=> $b->{ref_end} ||
                                        $a->{read_name} cmp $b->{read_name}
                                    } @filtered_sv_list;

                                    # Loop through the filtered SVs to cluster them
                                    while (@filtered_sv_list) {
                                        # Take the first SV as the initial candidate for clustering
                                        my $peak_sv = shift @filtered_sv_list;

                                        # If $peak_sv is not defined, skip processing
                                         next unless defined $peak_sv;

                                        # Initialize peak_start and peak_end correctly after selecting peak_sv
                                        my $peak_start = $peak_sv->{ref_start};
                                        my $peak_end = $peak_sv->{ref_end};

                                        # Further refine the cluster to include only those within max_distance of both start and end
                                        my @final_cluster = grep {
                                            abs($_->{ref_start} - $peak_start) <= $max_distance &&
                                            abs($_->{ref_end} - $peak_end) <= $max_distance
                                        } @second_cluster;

                                        # Skip empty clusters
                                        next unless @final_cluster;

                                         # Duplication-specific logic
                                        if ($sv_type_name eq "Duplication") {
                                            #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                            # Initialize the final duplication cluster to include peak_sv
                                            my @final_duplication_cluster = ($peak_sv);
                                            my $continue_merging = 1;

                                            while ($continue_merging) {
                                                $continue_merging = 0;  # Assume no further merges will happen

                                                foreach my $sv (@final_cluster) {
                                                    next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                    # Check if the SV can be merged with peak_sv based on node bounds
                                                    if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}) {                        
                                                        my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}};

                                                        # Ensure node counts are the same for merging
                                                        if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}) {
                                                            my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}};
                                                            
                                                            #print "$peak_sv->{read_name}: Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);
                                                            # Skip if node counts are not equal
                                                            #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                            #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                            # Ensure node counts match
                                                            if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                next;  # Skip this iteration if counts don't match
                                                            }

                                                            # 如果通过，继续输出
                                                            #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                            # Check if nodes can be merged
                                                            my $merge;
                                                            if (@min_max_merge_bounds == 2) {
                                                                # Case for exactly two nodes: only check first and last nodes
                                                                $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                            } elsif (@min_max_merge_bounds > 2) {
                                                                # Case for more than two nodes: check intermediate and first-last nodes
                                                                $merge = can_merge_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance) if $merge;
                                                            }

                                                            if ($merge) {
                                                                # Update bounds after merging
                                                                update_node_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}} = \@min_max_merge_bounds;
                                                                #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                # Add this SV to the final duplication cluster
                                                                push @final_duplication_cluster, $sv;

                                                                # Remove this SV from @final_cluster to avoid double-checking
                                                                @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                 # Set flag to continue merging in the next iteration
                                                                $continue_merging = 1;
                                                                last;  # Exit loop and reprocess after merg3
                                                            }
                                                        }
                                                    }
                                                }

                                                # Ensure @final_duplication_cluster is not empty
                                                if (!@final_duplication_cluster) {
                                                    print "Final duplication cluster is empty, skipping further processing.\n";
                                                    next;
                                                }

                                                # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                @final_cluster = @final_duplication_cluster;
                                            }

                                            # Check if final_cluster is valid and non-empty
                                            if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                next;
                                            }
                                        }

                                        # SCGR-specific logic
										#print "last SCGR started!\n";
                                        if ($sv_type_name eq "SCGR") {
                                            #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                            # Initialize the final duplication cluster to include peak_sv
                                            my @final_SCGR_cluster = ($peak_sv);
                                            my $continue_merging = 1;

                                            while ($continue_merging) {
                                                $continue_merging = 0;  # Assume no further merges will happen

                                                foreach my $sv (@final_cluster) {
                                                    next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                                    # Check if the SV can be merged with peak_sv based on node bounds
                                                    if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}) {                        
                                                        my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}};

                                                        # Ensure node counts are the same for merging
                                                        if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}) {
                                                            my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}};

                                                            # Skip if node counts are not equal
                                                            #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                            #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                            # Ensure node counts match
                                                            if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                                next;  # Skip this iteration if counts don't match
                                                            }

                                                            # 如果通过，继续输出
                                                            #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                            # Check if nodes can be merged
                                                            my $merge = can_merge_SCGR_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                       

                                                            if ($merge) {
                                                                # Update bounds after merging
                                                                update_node_SCGR_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                                @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}} = \@min_max_merge_bounds;
                                                                #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                                # Add this SV to the final duplication cluster
                                                                push @final_SCGR_cluster, $sv;

                                                                # Remove this SV from @final_cluster to avoid double-checking
                                                                @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                                # Set flag to continue merging in the next iteration
                                                                $continue_merging = 1;
                                                                last;  # Exit loop and reprocess after merge
                                                            }
                                                        }
                                                    }
                                                }

                                                # Ensure @final_duplication_cluster is not empty
                                                if (!@final_SCGR_cluster) {
                                                    print "Final SCGR cluster is empty, skipping further processing.\n";
                                                    next;
                                                }

                                                # Replace the @final_cluster with @final_duplication_cluster for further processing
                                                @final_cluster = @final_SCGR_cluster;
                                            }

                                            # Check if final_cluster is valid and non-empty
                                            if (!@final_cluster || scalar(@final_cluster) == 0) {
                                                print "Warning: final_cluster is undefined or empty, skipping.\n";
                                                next;
                                            }
                                        }

                                        # Calculate the total support for the cluster
                                         my $total_support = sum_support(\@final_cluster);

                                        # If the total support meets the minimum threshold, finalize the cluster
                                        if ($total_support >= $min_reads_support) {
                                            # Update peak_sv with the total support
                                             $peak_sv->{support} = $total_support;

                                            # Use peak_sv as the representative SV
                                            my $representative_sv = $peak_sv;

                                            # Store the cluster's SVs and the representative SV
                                            push @third_cluster, { representative_sv => $representative_sv, svs => \@final_cluster };

                                            # Remove clustered SVs from the list (only remove those in @final_cluster)
                                            @second_cluster = grep {
                                                !in_cluster($_, \@final_cluster)
                                            } @second_cluster;

                                            # Also remove overlapping SVs from @filtered_sv_list
                                             @filtered_sv_list = grep {
                                                !in_cluster($_, \@final_cluster)
                                            } @filtered_sv_list;

                                            # After removing clustered SVs, restart the loop to process remaining SVs
                                            redo if @filtered_sv_list;
                                        } else {
                                            # If the total support does not meet the threshold, remove peak_sv from @filtered_sv_list
                                            @filtered_sv_list = grep { $_ ne $peak_sv } @filtered_sv_list;

                                            # Continue to process other SVs from @filtered_sv_list
                                            redo if @filtered_sv_list;  # Restart loop if there are still SVs to process             
                                            last;  # If no more SVs in filtered list, exit loop
                                        }
                                    }

                                    foreach my $cluster (@third_cluster) {
                                        my $representative_sv = $cluster->{representative_sv};

                                        # Ensure the representative SV contains the necessary fields
                                        my $chromosome = $representative_sv->{chromosome} // '';
                                        my $sv_type = $representative_sv->{sv_type} // '';
                                        my $anchor_ref_start = $representative_sv->{anchor_ref_start} // '';
                                        my $anchor_ref_end = $representative_sv->{anchor_ref_end} // '';
                                        my $ref_start = $representative_sv->{ref_start} // '';
                                        my $ref_end = $representative_sv->{ref_end} // '';
                                        my $read_name = $representative_sv->{read_name} // '';
                                        my $anchor_seq = $representative_sv->{anchor_seq} // '';
                                        my $support = $representative_sv->{support} // '';
                                                

                                        # Calculate total support by summing up the support of all SVs in the cluster
                                        #my $total_support = sum_support($cluster->{svs});

                                        # Collect output line for later sorting
                                        push @output_lines, join("\t",
                                            $chromosome,
                                            $anchor_ref_start,
                                            $anchor_ref_end,
                                            $ref_start,
                                            $ref_end,
                                            $sv_type,
                                            $support,
                                            $read_name,
                                            $anchor_seq
                                        );
                                        $countSVcluster{$sv_type}++;
                                    }
                                }
                            }
                            @second_cluster = ();
                            $second_cluster_max_ref_end = $sv->{ref_end};
                        }
                        push @second_cluster, $sv;
                    }
                    if (any_left_span_right_span(@second_cluster)) {
                        if (sum_support(\@second_cluster) >= $min_reads_support) {
                            my @third_cluster = ();
                            # Function to cluster SV signals within a block and collect results for sorting

                            # First, filter SV list to only include those with leftSpanRightSpan integrity
                            my @filtered_sv_list = grep { $_->{integrity} =~ /leftSpanRightSpan/ && $_->{score} > 1 }@second_cluster;

                            # Skip if no SV with leftSpanRightSpan integrity is found
                            if (!@filtered_sv_list) {
                                print "No SVs with 'leftSpanRightSpan' integrity found. Skipping clustering.\n";
                            }

                            # Sort SVs by support, then by length of anchor_seq, ref_start, ref_end, and read_name
                            @filtered_sv_list = sort {
                                $b->{support} <=> $a->{support} ||
                                $b->{score} <=> $a->{score} ||
                                $a->{ref_start} <=> $b->{ref_start} ||
                                $a->{ref_end} <=> $b->{ref_end} ||
                                $a->{read_name} cmp $b->{read_name}
                            } @filtered_sv_list;
							#my $c = @filtered_sv_list;
							#print "$filtered_sv_list[0]->{ref_path}\t$filtered_sv_list[1]->{ref_path}\n";

                            # Loop through the filtered SVs to cluster them
                            while (@filtered_sv_list) {
                                # Take the first SV as the initial candidate for clustering
                                my $peak_sv = shift @filtered_sv_list;

                                # If $peak_sv is not defined, skip processing
                                next unless defined $peak_sv;

                                # Initialize peak_start and peak_end correctly after selecting peak_sv
                                my $peak_start = $peak_sv->{ref_start};
                                my $peak_end = $peak_sv->{ref_end};
								#print "$peak_sv->{ref_path}\n";

                                # Further refine the cluster to include only those within max_distance of both start and end
                                my @final_cluster = grep {
                                    abs($_->{ref_start} - $peak_start) <= $max_distance &&
                                    abs($_->{ref_end} - $peak_end) <= $max_distance
                                } @second_cluster;

                                # Skip empty clusters
                                next unless @final_cluster;

                                 # Duplication-specific logic
                                if ($sv_type_name eq "Duplication") {
                                    #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                    # Initialize the final duplication cluster to include peak_sv
                                    my @final_duplication_cluster = ($peak_sv);
                                    my $continue_merging = 1;

                                    while ($continue_merging) {
                                        $continue_merging = 0;  # Assume no further merges will happen

                                        foreach my $sv (@final_cluster) {
                                            next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                            # Check if the SV can be merged with peak_sv based on node bounds
                                            if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}) {                        
                                                my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}};

                                                # Ensure node counts are the same for merging
                                                if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}) {
                                                    my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}};
                                                        
                                                    # Skip if node counts are not equal
                                                    #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                    #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";

                                                    # Ensure node counts match
                                                    if (@min_max_merge_bounds != @min_max_target_bounds) {
                                                        next;  # Skip this iteration if counts don't match
                                                    }

                                                    # 如果通过，继续输出
                                                    #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                    # Check if nodes can be merged
                                                    my $merge;
                                                    if (@min_max_merge_bounds == 2) {
                                                        # Case for exactly two nodes: only check first and last nodes
                                                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                    } elsif (@min_max_merge_bounds > 2) {
                                                        # Case for more than two nodes: check intermediate and first-last nodes
                                                        $merge = can_merge_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                        $merge = can_merge_first_last_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance) if $merge;
                                                    }

                                                    if ($merge) {
                                                        # Update bounds after merging
                                                        update_node_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                        @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}} = \@min_max_merge_bounds;

                                                        # Add this SV to the final duplication cluster
                                                        push @final_duplication_cluster, $sv;

                                                        # Remove this SV from @final_cluster to avoid double-checking
                                                        @final_cluster = grep { $_ ne $sv } @final_cluster;
                                
                                                         # Set flag to continue merging in the next iteration
                                                        $continue_merging = 1;
                                                        last;  # Exit loop and reprocess after merg3
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    # Ensure @final_duplication_cluster is not empty
                                    if (!@final_duplication_cluster) {
                                        print "Final duplication cluster is empty, skipping further processing.\n";
                                        next;
                                    }

                                    # Replace the @final_cluster with @final_duplication_cluster for further processing
                                    @final_cluster = @final_duplication_cluster;
                                }

                                # SCGR-specific logic
                                if ($sv_type_name eq "SCGR") {
									#print "last SCGR started!\n";
                                    #print "Processing final_cluster with " . scalar(@final_cluster) . " elements for Duplication.\n";
                                    # Initialize the final duplication cluster to include peak_sv
                                    my @final_SCGR_cluster = ($peak_sv);
                                    my $continue_merging = 1;

                                    while ($continue_merging) {
                                        $continue_merging = 0;  # Assume no further merges will happen
										#print "$peak_sv->{ref_path}\n";

                                        foreach my $sv (@final_cluster) {
                                            next if $sv eq $peak_sv;  # Skip the peak_sv itself

                                            # Check if the SV can be merged with peak_sv based on node bounds
                                            if (exists $node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}) {                        
                                                my @min_max_target_bounds = @{$node_bounds{$sv->{chromosome}}{$sv->{ref_start}}{$sv->{ref_end}}{$sv->{read_name}}{$sv->{ref_path}}};

                                                # Ensure node counts are the same for merging
                                                if (exists $node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}) {
                                                    my @min_max_merge_bounds = @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}};

                                                    # Skip if node counts are not equal
                                                    #print "Number of min_max_merge_bounds: " . scalar(@min_max_merge_bounds) . "\n";
                                                    #print "Number of min_max_target_bounds: " . scalar(@min_max_target_bounds) . "\n";
													
													#my $c1 = @min_max_merge_bounds;
													#my $c2 = @min_max_target_bounds;
													#print "$sv->{ref_path}\t$c2\t$peak_sv->{ref_path}\t$c1\n";
                                                    # Ensure node counts match
                                                    if (@min_max_merge_bounds != @min_max_target_bounds) {
														#print "$sv->{ref_path}\t$c2\t$peak_sv->{ref_path}\t$c1\n";
                                                        next;  # Skip this iteration if counts don't match
                                                    }

                                                    # 如果通过，继续输出
                                                    #print "$peak_sv->{read_name}\t$peak_sv->{support}\t$sv->{read_name}\t$sv->{support}\n";

                                                    # Check if nodes can be merged
                                                    my $merge = can_merge_SCGR_nodes(\@min_max_merge_bounds, \@min_max_target_bounds, $node_max_distance);
                                                                       

                                                    if ($merge) {
                                                        # Update bounds after merging
                                                        update_node_SCGR_bounds(\@min_max_merge_bounds, \@min_max_target_bounds);
                                                        @{$node_bounds{$peak_sv->{chromosome}}{$peak_sv->{ref_start}}{$peak_sv->{ref_end}}{$peak_sv->{read_name}}{$peak_sv->{ref_path}}} = \@min_max_merge_bounds;
                                                        #print "Updated min_max_merge_bounds: " . Dumper(\@min_max_merge_bounds);

                                                        # Add this SV to the final duplication cluster
                                                        push @final_SCGR_cluster, $sv;

                                                        # Remove this SV from @final_cluster to avoid double-checking
                                                        @final_cluster = grep { $_ ne $sv } @final_cluster;
														#print "$sv->{ref_path}\t$peak_sv->{ref_path}\tmerged\n";
                                
                                                        # Set flag to continue merging in the next iteration
                                                        $continue_merging = 1;
                                                        last;  # Exit loop and reprocess after merge
                                                    }
                                                }
                                            }
                                        }

                                        # Ensure @final_duplication_cluster is not empty
                                        if (!@final_SCGR_cluster) {
                                            print "Final SCGR cluster is empty, skipping further processing.\n";
                                            next;
                                        }

                                        # Replace the @final_cluster with @final_duplication_cluster for further processing
                                        @final_cluster = @final_SCGR_cluster;
                                    }

                                    # Check if final_cluster is valid and non-empty
                                    if (!@final_cluster || scalar(@final_cluster) == 0) {
                                        print "Warning: final_cluster is undefined or empty, skipping.\n";
                                        next;
                                    }
                                }
                                

                                # Calculate the total support for the cluster
                                 my $total_support = sum_support(\@final_cluster);

                                # If the total support meets the minimum threshold, finalize the cluster
                                if ($total_support >= $min_reads_support) {
                                    # Update peak_sv with the total support
                                     $peak_sv->{support} = $total_support;

                                    # Use peak_sv as the representative SV
                                    my $representative_sv = $peak_sv;

                                    # Store the cluster's SVs and the representative SV
                                    push @third_cluster, { representative_sv => $representative_sv, svs => \@final_cluster };

                                    # Remove clustered SVs from the list (only remove those in @final_cluster)
                                    @second_cluster = grep {
                                        !in_cluster($_, \@final_cluster)
                                    } @second_cluster;

                                    # Also remove overlapping SVs from @filtered_sv_list
                                    @filtered_sv_list = grep {
                                        !in_cluster($_, \@final_cluster)
                                    } @filtered_sv_list;

                                    # After removing clustered SVs, restart the loop to process remaining SVs
                                    redo if @filtered_sv_list;
                                } else {
                                    # If the total support does not meet the threshold, remove peak_sv from @filtered_sv_list
                                    @filtered_sv_list = grep { $_ ne $peak_sv } @filtered_sv_list;

                                    # Continue to process other SVs from @filtered_sv_list
                                    redo if @filtered_sv_list;  # Restart loop if there are still SVs to process             
                                    last;  # If no more SVs in filtered list, exit loop
                                }
                            }

                            foreach my $cluster (@third_cluster) {
                                my $representative_sv = $cluster->{representative_sv};

                                # Ensure the representative SV contains the necessary fields
                                my $chromosome = $representative_sv->{chromosome} // '';
                                my $sv_type = $representative_sv->{sv_type} // '';
                                my $anchor_ref_start = $representative_sv->{anchor_ref_start} // '';
                                my $anchor_ref_end = $representative_sv->{anchor_ref_end} // '';
                                my $ref_start = $representative_sv->{ref_start} // '';
                                my $ref_end = $representative_sv->{ref_end} // '';
                                my $read_name = $representative_sv->{read_name} // '';
                                my $anchor_seq = $representative_sv->{anchor_seq} // '';
                                my $support = $representative_sv->{support} // '';
                                                

                                # Calculate total support by summing up the support of all SVs in the cluster
                                #my $total_support = sum_support($cluster->{svs});

                                # Collect output line for later sorting
                                push @output_lines, join("\t",
                                    $chromosome,
                                    $anchor_ref_start,
                                    $anchor_ref_end,
                                    $ref_start,
                                    $ref_end,
                                    $sv_type,
                                    $support,
                                    $read_name,
                                    $anchor_seq
                                );
                                $countSVcluster{$sv_type}++;
                            }    
                        }
                    }    
                }
            }
        }
    }
}


# Now, after processing all blocks, sort the output lines
@output_lines = sort {
    my ($a_chr, $a_start) = (split /\t/, $a)[0, 3];  # Extract chromosome and ref_start for line A
    my ($b_chr, $b_start) = (split /\t/, $b)[0, 3];  # Extract chromosome and ref_start for line B

    # Extract numeric part of chromosome for comparison, default to 999 for non-numeric chromosomes
    my ($a_num) = $a_chr =~ /(\d+)/ ? $1 : 999;
    my ($b_num) = $b_chr =~ /(\d+)/ ? $1 : 999;

    # Compare numeric chromosome values first
    my $cmp = $a_num <=> $b_num;

    # If chromosomes are the same, compare by start position (ref_start)
    return $cmp || $a_start <=> $b_start;
} @output_lines;

# Finally, write the sorted output lines to the gzip file
foreach my $line (@output_lines) {
    print $gz "$line\n";  # 将行写入 gzip 文件
}

# Subroutine to initialize node bounds
sub initialize_bounds_for_nodes {
    my ($nodes) = @_;
    my @min_max_bounds;
    foreach my $node (@$nodes) {
        my ($start, $end) = split /-/, $node;
        push @min_max_bounds, [$start, $start, $end, $end];  # Initialize bounds for each node
    }
    return @min_max_bounds;
}

# Subroutine to check if nodes (from the second to the second-last) can be merged
sub can_merge_nodes {
    my ($merge_bounds, $target_bounds, $node_max_distance) = @_;

    # Initialize a variable to track if merging is possible
    my $can_merge = 1;  # Assume merging is possible

    # Check the intermediate nodes (from second to second-last)
    for my $i (1 .. $#$merge_bounds - 1) {
        # Get the boundary information of the merging node
        my ($start_merge_min, $start_merge_max, $end_merge_min, $end_merge_max) = @{$merge_bounds->[$i]};
        
        # Get the boundary information of the node being merged
        my ($start_target_min, $start_target_max, $end_target_min, $end_target_max) = @{$target_bounds->[$i]};

        # Check if the start and end boundaries of the merging node and target node overlap
        my $start_overlap = 0;
        my $end_overlap = 0;

        # Check start overlap (expand merging signal boundaries, check if overlaps with the target)
        if (($start_merge_min - $node_max_distance <= $start_target_max) && 
            ($start_merge_max + $node_max_distance >= $start_target_min)) {
            $start_overlap = 1;
        }

        # Check end overlap (expand merging signal boundaries, check if overlaps with the target)
        if (($end_merge_min - $node_max_distance <= $end_target_max) && 
            ($end_merge_max + $node_max_distance >= $end_target_min)) {
            $end_overlap = 1;
        }

        # If neither start nor end overlap, merging is not possible
        if (!$start_overlap || !$end_overlap) {
            $can_merge = 0;  # Set the control variable to 0, indicating that merging is not possible
            last;  # Exit the loop early since we found a case where merging fails
        }
    }

    return $can_merge;  # Return 1 if merging is possible, otherwise 0
}

# Subroutine to check if first and last nodes can be merged
sub can_merge_first_last_nodes {
    my ($merge_bounds, $target_bounds, $node_max_distance) = @_;

    # Initialize variables to track if merging is possible
    my $start_overlap = 0;
    my $end_overlap = 0;

    # Check first node bounds
    my ($start_merge_min_first, $start_merge_max_first, $end_merge_min_first, $end_merge_max_first) = @{$merge_bounds->[0]};
    my ($start_target_min_first, $start_target_max_first, $end_target_min_first, $end_target_max_first) = @{$target_bounds->[0]};

    # Check start overlap (expand merging signal boundaries, check if overlaps with the target)
    if (($end_merge_min_first - $node_max_distance <= $end_target_max_first) && 
        ($end_merge_max_first + $node_max_distance >= $end_target_min_first)) {
        $end_overlap = 1;
    }

    # Check last node bounds
    my ($start_merge_min_last, $start_merge_max_last, $end_merge_min_last, $end_merge_max_last) = @{$merge_bounds->[-1]};
    my ($start_target_min_last, $start_target_max_last, $end_target_min_last, $end_target_max_last) = @{$target_bounds->[-1]};

    # Check end overlap (expand merging signal boundaries, check if overlaps with the target)
    if (($start_merge_min_last - $node_max_distance <= $start_target_max_last) && 
        ($start_merge_max_last + $node_max_distance >= $start_target_min_last)) {
        $start_overlap = 1;
    }
    
    # If either start or end does not overlap, merging is not possible
    if (!$start_overlap || !$end_overlap) {
        return 0;  # Return 0 if no overlap is found
    }

    return 1;  # Return 1 if overlap is found
}

# Subroutine to update node bounds after merging
sub update_node_bounds {
    my ($merge_bounds, $target_bounds) = @_;

    for my $i (0 .. $#$merge_bounds) {
        # Get the boundary information of the merging node
        my ($start_merge_min, $start_merge_max, $end_merge_min, $end_merge_max) = @{$merge_bounds->[$i]};
        
        # Get the boundary information of the node being merged
        my ($start_target_min, $start_target_max, $end_target_min, $end_target_max) = @{$target_bounds->[$i]};

        # Update the min and max bounds for start and end
        my $start_min_bound = ($start_merge_min < $start_target_min) ? $start_merge_min : $start_target_min;
        my $start_max_bound = ($start_merge_max > $start_target_max) ? $start_merge_max : $start_target_max;
        my $end_min_bound   = ($end_merge_min < $end_target_min) ? $end_merge_min : $end_target_min;
        my $end_max_bound   = ($end_merge_max > $end_target_max) ? $end_merge_max : $end_target_max;

        # Store the updated bounds
        $merge_bounds->[$i] = [$start_min_bound, $start_max_bound, $end_min_bound, $end_max_bound];
    }
}

# Custom sorting subroutine for chromosomes
sub custom_chromosome_sort {
    my ($a_prefix, $a_suffix) = split_chromosome($a);
    my ($b_prefix, $b_suffix) = split_chromosome($b);

    # First, sort by the prefix (numerical first, then lexicographical)
    if ($a_prefix ne $b_prefix) {
        # If both are numbers, compare numerically
        if ($a_prefix =~ /^\d+$/ && $b_prefix =~ /^\d+$/) {
            return $a_prefix <=> $b_prefix;
        }
        # Otherwise, compare lexicographically
        return $a_prefix cmp $b_prefix;
    }

    # If prefixes are equal, compare by suffix (numerically if both are digits)
    return $a_suffix <=> $b_suffix if $a_suffix =~ /^\d+$/ && $b_suffix =~ /^\d+$/;

    # Otherwise, compare lexicographically
    return $a_suffix cmp $b_suffix;
}

# Helper subroutine to split a chromosome name into prefix and suffix
sub split_chromosome {
    my ($chr) = @_;

    # Match the prefix (non-digit part) and suffix (digit part if it exists)
    if ($chr =~ /^([^\d]+)(\d*)$/) {
        return ($1, $2 // 0);  # Return prefix and suffix, with suffix defaulting to 0
    }

    # If it's purely numeric, return the number as the prefix, and an empty suffix
    return ($chr, 0);
}

# Helper function to check if any SV in the cluster has integrity "leftSpanRightSpan"
sub any_left_span_right_span {
    my @cluster = @_;
    foreach my $sv (@cluster) {
        if ($sv->{integrity} eq 'leftSpanRightSpan' && $sv->{score} > 1) {
            return 1;  # Return true if found
        }
    }
    return 0;  # Return false if none found
}

# Helper function to calculate total support
sub sum_support {
    my ($svs) = @_;
    my $total_support = 0;
    $total_support += $_->{support} for @$svs;
    return $total_support;
}

# Helper function to check if an SV is in a cluster
sub in_cluster {
    my ($sv, $cluster) = @_;
    for my $c (@$cluster) {
        return 1 if $sv->{ref_start} == $c->{ref_start} && $sv->{ref_end} == $c->{ref_end} && $sv->{read_name} eq $c->{read_name} && $sv->{ref_path} eq $c->{ref_path};
    }
    return 0;
}

# Subroutine to initialize node bounds
sub initialize_bounds_SCGR_for_nodes {
    
    my ($nodes) = @_;
    my @min_max_bounds;
    foreach my $node (@$nodes) {
        my ($chr, $path) = split/:/,$node;
        my ($start, $end) = split /-/, $path;
        push @min_max_bounds, [$chr, $start, $start, $end, $end];  # Initialize bounds for each node
    }
    return @min_max_bounds;
}

# Subroutine to check if nodes (from the second to the second-last) can be merged
sub can_merge_SCGR_nodes {
    my ($merge_bounds, $target_bounds, $node_max_distance) = @_;

    # Initialize a variable to track if merging is possible
    my $can_merge = 1;  # Assume merging is possible

    # Check the intermediate nodes (from second to second-last)
    for my $i (1 .. $#$merge_bounds - 1) {
        # Get the boundary information of the merging node
        my ($chr_merge, $start_merge_min, $start_merge_max, $end_merge_min, $end_merge_max) = @{$merge_bounds->[$i]};
        
        # Get the boundary information of the node being merged
        my ($chr_target, $start_target_min, $start_target_max, $end_target_min, $end_target_max) = @{$target_bounds->[$i]};

        # Check if the start and end boundaries of the merging node and target node overlap
        my $start_overlap = 0;
        my $end_overlap = 0;
        my $chr_overlap = 0;
        
        # Check start overlap (expand merging signal boundaries, check if overlaps with the target)
        if (($start_merge_min - $node_max_distance <= $start_target_max) && 
            ($start_merge_max + $node_max_distance >= $start_target_min)) {
            $start_overlap = 1;
        }

        # Check end overlap (expand merging signal boundaries, check if overlaps with the target)
        if (($end_merge_min - $node_max_distance <= $end_target_max) && 
            ($end_merge_max + $node_max_distance >= $end_target_min)) {
            $end_overlap = 1;
        }
        
        #check the chromosomes are same
        if ($chr_merge eq $chr_target){
            $chr_overlap = 1;
        }

        # If neither start nor end overlap, merging is not possible
        if (!$start_overlap || !$end_overlap || !$chr_overlap) {
            $can_merge = 0;  # Set the control variable to 0, indicating that merging is not possible
            last;  # Exit the loop early since we found a case where merging fails
        }
    }

    return $can_merge;  # Return 1 if merging is possible, otherwise 0
}

# Subroutine to update node bounds after merging
sub update_node_SCGR_bounds {
    my ($merge_bounds, $target_bounds) = @_;

    for my $i (1 .. $#$merge_bounds - 1) {
        # Get the boundary information of the merging node
        my ($chr_merge, $start_merge_min, $start_merge_max, $end_merge_min, $end_merge_max) = @{$merge_bounds->[$i]};
        
        # Get the boundary information of the node being merged
        my ($chr_target, $start_target_min, $start_target_max, $end_target_min, $end_target_max) = @{$target_bounds->[$i]};

        # Update the min and max bounds for start and end
        my $start_min_bound = ($start_merge_min < $start_target_min) ? $start_merge_min : $start_target_min;
        my $start_max_bound = ($start_merge_max > $start_target_max) ? $start_merge_max : $start_target_max;
        my $end_min_bound   = ($end_merge_min < $end_target_min) ? $end_merge_min : $end_target_min;
        my $end_max_bound   = ($end_merge_max > $end_target_max) ? $end_merge_max : $end_target_max;

        # Store the updated bounds
        $merge_bounds->[$i] = [$chr_merge, $start_min_bound, $start_max_bound, $end_min_bound, $end_max_bound];
    }
}

my $total_written  = 0;
foreach my $key (keys %countIntermediate){
    unless (exists $countSVcluster{$key}){
        $countSVcluster{$key} = 0;
    }
    unless (exists $countSig{$key}){
        $countSig{$key} = 0;
    }
    unless (exists $countIntermediate{$key}){
        $countIntermediate{$key} = 0;
    }
    $total_written += $countSVcluster{$key};
    print "SV Type: $key, Number of SV signatures: $countSig{$key}, Number of SV signatures for redundancy: $countIntermediate{$key}, Number written: $countSVcluster{$key}, completed\n";
}
print "SV clustering complete. $total_written SVs were written to $output_file.gz\n";

# Subroutine to show usage
sub usage {
    print <<END;
Usage: perl $0 --input <SV input file> --output <Clustered SV output file> [options]

Options:
    --input                                  Input file with structural variant data (required). Can be plain text or gzipped (.gz).
    --output                                 Output file for clustered structural variant data (required).
    --max-distance-insertion <int>           Maximum distance for clustering Insertion SVs (default: 130).
    --max-distance-deletion <int>            Maximum distance for clustering Deletion SVs (default: 70).
    --max-distance-CGR <int>                 Maximum distance for clustering CGR (complex genomic rearrangement) SVs (default: 490).
    --max-distance-duplication <int>         Maximum distance for clustering Duplication SVs (default: 220).
    --max-distance-SCGR <int>                Maximum distance for clustering SCGR SVs (default: 500).
    --min-reads-support-insertion <int>      Minimum reads support for Insertion SVs (default: 2).
    --min-reads-support-deletion <int>       Minimum reads support for Deletion SVs (default: 2).
    --min-reads-support-CGR <int>            Minimum reads support for CGR SVs (default: 13).
    --min-reads-support-duplication <int>    Minimum reads support for Duplication SVs (default: 18).
    --min-reads-support-SCGR <int>           Minimum reads support for SCGR SVs (default: 18).
    --weight-complexINSDEL <float>           Weight for complex insertion/deletion reads (default: 1).
    --weight-unCover <float>                 Weight for uncovered reads at ends (default: 0.5).
    --block-gap <int>                        Gap size between SVs to split blocks (default: 500).
    --node-max-distance <int>                Maximum distance between split alignment points for Duplication SV merging (default: 10).
    --intermediate-output                    Output file for merging the SVs signature with same chromosome, start, end, and ref_path (default: not out).
    --help                                   Show this help message.

Examples:
    perl $0 --input sv_data.txt --output clustered_sv_output.txt --min-reads-support-insertion 5 --max-distance-deletion 10
    perl $0 --input sv_data.txt.gz --output clustered_sv_output.txt --min-reads-support-duplication 5 --max-distance-deletion 10

END
}
