#!/usr/bin/perl -w

use lib('/home/munzmatt/lib/perl5');

use strict;
use warnings;
use List::Util qw(sum);
use Statistics::Distributions;

#################################################
# Author: Matthias Munz
# Contact: matthias<DOT>munz<AT>gmx<DOT>de
#################################################


# This program runs the Friedman-Rufsky (FR) Test
# Input: Two datasets consisting of n-dimensional data points
# Workflow:
# 1. Compute distance matrix
# 2. Create minimal spanning tree MST using Prim
# 3. Compute number of runs
# 4. Compute mean, variance, permutation parameter and quantity W
# 5. Compute similarity measure


if (!$ARGV[0] || !$ARGV[1] || !$ARGV[2] || !$ARGV[3]){
	print STDOUT "Input: <dataset1> <dataset2> <header_flag> <print_dist_flag>
	<dataset1>\t\t\ttab-delimited
	<dataset2>\t\t\ttab-delimited
	<header_flag>\t\t{-1,1}, files contain header
	<print_dist_flag>\t{-1,1}, print distance matrix
	";
	exit(0);
}

my $data_file1 = $ARGV[0];
my $data_file2 = $ARGV[1];
my $header = $ARGV[2]; $header = 0 if $ARGV[2] == -1;
my $dist = $ARGV[3]; $dist = 0 if $ARGV[3] == -1;

# Parse data (tsv files)
my @dataset1 = &_parse_data($data_file1, $header);
my @dataset2 = &_parse_data($data_file2, $header);

#################################################
# Workflow

# Step 1
# Compute distance table
my $distance_table = &_create_distance_table([@dataset1, @dataset2]);
if ($dist){
	&_print_matrix($distance_table);
}

# Step 2
# Compute MST Prim Algo
my $tree = &_prim($distance_table);

# Step 3
# Compute number of runs
my $runs = &_count_FR_runs($tree, @dataset1 - 1);
print "Number of FR-Runs: ".$runs."\n";

# Step 4
# Compute permutation parameter C
my $C = &_compute_FR_permutation_par($tree);
print "FR-Permutation parameter C: $C\n";

# Compute variance
my $variance = &_compute_FR_variance($C, scalar @dataset1, scalar @dataset2);
print "FR-Variance: $variance\n";

# Compute mean
my $mean = &_compute_FR_mean(scalar @dataset1, scalar @dataset2);
print "FR-Mean: $mean\n";

# Compute quantity W (value of random variable)
my $W = &_compute_FR_quantity($variance, $mean);
print "FR-Quantity W: $W\n";

# Step 5
# Compute similarity measure (area under density curve left of W)
# 1 is the best similarity
my $uprob=Statistics::Distributions::uprob($W);
print "Datasets 'A' and 'B' have a similarity score of ".((1-$uprob)*100)."\n";


#################################################
# Routines

# Parse data (tsv files)
# Row: date, Column: attribute
sub _parse_data {
	my ($data_file, $header) = @_;
	
	my $first = 1;
	my @result;
	open(IN, "<$data_file") || die "Can't open file $data_file: $!";
	while(<IN>){
		
		if ($header && $first){
			$first = 0;
			next;
		}
		chomp($_);
		my @date = split(/\t/, $_);
		push(@result, \@date);
	}
	
	return @result;
}

# Compute quantity W
# Source of formular: Shao, Jie; Distribution-based Similarity Measures for Multi-dimensional Point Set Retrieval Applications; 2008
sub _compute_FR_quantity {
	my ($variance, $mean) = @_;
	return ($runs - $mean)/sqrt($variance);
}

# Compute mean
# Source of formular: Shao, Jie; Distribution-based Similarity Measures for Multi-dimensional Point Set Retrieval Applications; 2008
sub _compute_FR_mean {
	my ($m, $n) = @_;
	my $N = $m + $n;
	
	return ((2*$m*$n)/$N)+1;
}

# Compute variance
# Source of formular: Shao, Jie; Distribution-based Similarity Measures for Multi-dimensional Point Set Retrieval Applications; 2008
sub _compute_FR_variance {
	my ($C, $m, $n) = @_;
	my $N = $m + $n;
	
	my $variance = ((2*$m*$n)/($N*($N-1))) * 
					( (((2*$m*$n)-$N)/$N) + 
					((($C-$N+2)/(($N-2)*($N-3))) * (($N*($N-1))-(4*$m*$n)+2)) );
					
	return $variance;
}

# Compute permutation parameter C
# Source of formular: Shao, Jie; Distribution-based Similarity Measures for Multi-dimensional Point Set Retrieval Applications; 2008
sub _compute_FR_permutation_par {
	my $tree = shift;
	
	my %vertice2n_edges;
	
	foreach (@{$tree}){
		if (exists($vertice2n_edges{${$_}[0]})){
			$vertice2n_edges{${$_}[0]}++;
		}
		else {
			$vertice2n_edges{${$_}[0]} = 1;
		}
		
		if (exists($vertice2n_edges{${$_}[1]})){
			$vertice2n_edges{${$_}[1]}++;
		}
		else {
			$vertice2n_edges{${$_}[1]} = 1;
		}
	}
	
	my $C = 0;
	foreach (keys(%vertice2n_edges)){
		my $degree = $vertice2n_edges{$_};
		$C += $degree * ($degree - 1);
	}
	
	return $C/2;
}

# Count number of runs
# in the context of the FR-Test
sub _count_FR_runs {
	my ($tree, $s1_max_idx) = @_;
	
	my $runs = 1;
	
	foreach (@{$tree}){
		if ((${$_}[0] <= $s1_max_idx && ${$_}[1] > $s1_max_idx) ||
			(${$_}[0] > $s1_max_idx && ${$_}[1] <= $s1_max_idx)){
			$runs++;
		}
	}
	
	return $runs;
}

# Compute euclidian distance
sub _get_distance_between {
  if (@_ != 2)
  {
    die "You must supply 2 coordinates";
  }

  my $a_ref = shift;
  my $b_ref = shift;

  if (@{$a_ref} != @{$b_ref}){
    die "Coordinates do not have the same number of elements";
  }
  else {
    my @x;
    for (my $i = 0; $i < @{$a_ref}; $i++){
    	push(@x, (${$a_ref}[$i] - ${$b_ref}[$i])**2);
    }
    return sqrt(sum(@x));
  }
}

# Create distance table for data points
sub _create_distance_table {
	my $data_points = shift;
	
	
	my $cols = 0;
	my @result;
	for (my $i = 0; $i < scalar @{$data_points}; $i++){
		$cols++;
		
		for (my $j = 0; $j < $cols; $j++){
			if ($i == $j){
				$result[$i][$j] = 0;
			}
			else {
				$result[$i][$j] = _get_distance_between(${$data_points}[$i], ${$data_points}[$j]);
			}
		}
	}
	
	return \@result;
}

# Visualize any matrix e.g. distance matrix
sub _print_matrix {
    my $D = shift;
    my $size = @{$D};
    for(my $a=0;$a<$size;$a++){
    	my $cols = @{${$D}[$a]};
        for(my $b=0;$b<$cols;$b++){
        	my $n = sprintf "%.2f", ${$D}[$a][$b];
        	print $n."\t";
           	# printf "%4d", $n;
        }
        print "\n";
    }
}

# Create MST with Prim algorithm
sub _prim {

    my ($i,$total_cost);
    my (@visited, @cost, @pi, @tree);
    my $D = shift;
    my $size = @{$D};

    for($i=0; $i<$size ;$i++){
        $visited[$i] = 0; # No nodes visited yet
        $cost[$i] = "inf"; # All edges have infinite costs
    }

    $cost[0] = 0; # Initial node has cost 0
    $pi[0] = -1; # ...and no parents


	my $counter = 0;
    while(1){
    	
        my $mincost = "inf";
        my $u;


        # Select node that has minimum weight
        for($i=0; $i<$size ;$i++){
        	
        	if ($mincost eq "inf"){
        		if ($visited[$i] == 0 && !($cost[$i] eq "inf")){
        			$mincost = $cost[$i];
                	$u = $i; # Select node with min costs
        		}
        	}
        	elsif (!($cost[$i] eq "inf")){
        		if ($visited[$i] == 0 && $cost[$i] < $mincost){
        			$mincost = $cost[$i];
                	$u = $i; # Select node with min costs
        		}
        	}
        }


        # All visited
        if($mincost eq "inf"){
            print "Total weight of MST: ".$total_cost."\n";
            for($i=0; $i<$size ;$i++){
        	 	if ($visited[$i] == 0){
        	 		print STDERR "Node ".{$i+1}." unreachable\n";
        	 	}
        	 }
            last;
        }
        # Else add new costs to total costs
        else{
            $total_cost += $mincost;
        }
        $visited[$u] = 1;
        
        
        # Update costs and parents
        for($i=0; $i<$size ;$i++){
        	
        	my $row_idx = $u;
        	my $col_idx = $i;
        	
        	# Change row and col idx
        	if (@{${$D}[$u]} < $i + 1){
        		$row_idx = $i;
        		$col_idx = $u;
        	}
        	
        	# If node has already been visited or has infinite distance to selected node u
            if($visited[$i] || (${$D}[$row_idx][$col_idx]) eq "inf"){
                next;
            }
            # If existing costs are higher or equal than costs to selected node u
            elsif(!(${$D}[$row_idx][$col_idx] eq "inf") && ($cost[$i] eq "inf" || ${$D}[$row_idx][$col_idx] <= $cost[$i])){
                $cost[$i] = ${$D}[$row_idx][$col_idx];
                $pi[$i] = $u;
            }
        }
        
        
        if($pi[$u] == -1){
            # print "Add first node: ".($u+1)."\n";
        }
        else{
            # print "Edge added: e<".($pi[$u]+1).",".($u+1).">\n";
            push(@tree, [$pi[$u], $u]);
        }
    }
    
    return \@tree;
    
}