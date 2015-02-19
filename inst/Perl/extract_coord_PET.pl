#!/usr/bin/env perl;

# a script to extract coordinates of fragments from PET read file

use warnings;
use strict;
use Cwd;

my ($infile, $outfile, $summaryfile, $format) = @ARGV;

# process PET read file

open IN, "$infile" or die "Cannot open input file $infile\n";

my %seen_pos =();		# position (leftmost position, starting from 1)
#my %seen_pair =();		# pair number (1 or 2)
my %seen_str = ();		# strand ('F' or 'R')

my $start;
my $end;

my %nline = ();
my %start_pos = ();
my %end_pos = ();

while(<IN>){
	# parse
	
	chomp;	
	my ($id, $seq, $map, $t3, $t4, $t5, $chrt, $pos, $str, @rest) = split /\t/, $_;
	my ($id_only, $pair) = split /\//, $id;
	my $read_length = length $seq;
	
	# process coordinates
	
	if ( exists $seen_pos{$id_only} ) {		
		# if other pair was already observed (w.r.t. id), 
		# then write the coordinates of both reads in the pair
		
		#if ( $seen_pair{$id_only} ne $pair && $seen_str{$id_only} ne $str ) {
		if ( $seen_str{$id_only} ne $str ) {
			# pos1 & pos2 should have different strands and pair numbers
			
			my $pos1 = $seen_pos{$id_only};
			my $pos2 = $pos;
			
			# in the pair, write upstream pos first, 
			# then downstream pos + read length - 1.
			# [Note] eland_result: pos = leftmost position, starting from 1.
			
			if ( $pos1 < $pos2 && $str eq "R" ) {
				$pos2 += $read_length - 1;
				$start = $pos1;
				$end = $pos2;
			} elsif ( $pos1 > $pos2 && $seen_str{$id_only} eq "R" ) {
				$pos1 += $read_length - 1;
				$start = $pos2;
				$end = $pos1;
			} else {
				print "check id $id!\n";
			}
			
			# process coordinates
		    
		    if ( exists $start_pos{$chrt} ) {
		        # if there is already a matrix for chromosome, update it
		        
		        push @{$start_pos{$chrt}}, $start;
		        push @{$end_pos{$chrt}}, $end;
		        
		        # count proper output lines
		        
		        ${$nline{$chrt}}++;
		    } else {
		        # if there is no matrix for chromosome yet, construct one
		        
		        #@{$start_pos{$chrt}} = $start;
		        #@{$end_pos{$chrt}} = $end;
		        
		        @{$start_pos{$chrt}} = ();
		        @{$end_pos{$chrt}} = ();
		        
		        push @{$start_pos{$chrt}}, $start;
		        push @{$end_pos{$chrt}}, $end;
		        
		        # count proper output lines
		        
		        ${$nline{$chrt}} = 1;
		    }	
		} else {
			print "inappropriate pair for id $id.\n";
		}
	} else {
		# if other pair was not observed yet, record it
		
		$seen_pos{$id_only} = $pos;
		#$seen_pair{$id_only} = $pair;
		#$seen_pair{$id_only} = 1;
		$seen_str{$id_only} = $str;
	}
}

close( IN );

# write processed reads

foreach my $chr_id (keys %start_pos) {
    my $outfile_chr = $outfile."_".$chr_id;
    open OUT, ">$outfile_chr" or die "Cannot open $outfile_chr\n";
     
    my @start_pos_chr = @{$start_pos{$chr_id}};
    my @end_pos_chr = @{$end_pos{$chr_id}};
    
    for( my $i = 0; $i < scalar(@start_pos_chr); $i++ ) {
        print OUT "$chr_id\t$start_pos_chr[$i]\t$end_pos_chr[$i]\n";
    }
    
    close OUT;
}

# summary of processed file (chrID, # lines)

open OUT, ">$summaryfile" or die "Cannot open $summaryfile\n";
foreach my $chr_id (keys %nline) {
    print OUT "$chr_id\t${$nline{$chr_id}}\n";
}
close( OUT ); 
