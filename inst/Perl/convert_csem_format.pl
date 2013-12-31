# convert output from stand-alone csem to CSEM BED format

#!/usr/bin/env perl;
use warnings;
use strict;

# parse command arguments

if ( scalar(@ARGV) < 6 ) {
	print "\nUsage:\n";
	print "perl convert_csem_format.pl [csem_file] [outfile] [outfile_format] [ref_file]\n";
	print "\t[read_length] [pseudo_tags?]\n";
	print "\nExample:\n";
	print "perl convert_csem_format.pl example.csem example_final.csem bed hg19.ref 36 N\n\n";
	exit;
}

my ( $csem_file, $outfile_name, $outfile_format, $ref_file, $read_length, $pseudo_tags ) = @ARGV;

# chromosome info

open IN, "$ref_file" or die "Cannot open $ref_file\n";

<IN>;	# skip the first line

my $size_vec = <IN>;
chomp($size_vec);
my @size_list = split( /\s+/, $size_vec );

my $chr_vec = <IN>;
chomp($chr_vec);
my @chr_list = split( /\s+/, $chr_vec );

close IN;

# post-process chromosome & position 

open( IN, $csem_file ) or die "Cannot open csem file!\n";
open( OUT, ">", $outfile_name ) or die "Cannot open output file!\n";

<IN>;	# skip the first line (number of uni-reads and multi-reads)
	
while (my $line = <IN>) {
	chomp($line);
	my @element = split( /\s/, $line );
		# assume columns are separated by some white space
	
	# check for invalid line: may cause error in exporting step
	
	if ( scalar(@element) < 5 ) {
		next;
	}
	
	# post-process lines
	
	my ($id, $chr, $str, $loc, $prob) = @element;
	if ( $outfile_format ne "bed" ) {
		# first base is 0 in bowtie or BED
		# first base is 1 in table or GFF
		$loc++;
	}
	my $chrname = $chr_list[$chr];	# translate chromosome
	
	# write down processed lines
	# - generate pseudo-tags, if necessary (adapted from "round_tag_to_integer.pl")
	
	if ( $pseudo_tags eq "Y" ) {	
		# if we want to generate pseudo-tags,
		# then threshold prob at 0.5 & round prob to integer (i.e., set to one)
		# (exclude prob = 0.5 as well in order to avoid a read appears more than once)
		
		if ( $prob <= 0.5 ) {
			next;
		} else {
			$prob = 1;
		}
	}
			
	# write down results
	
	my $start;
	my $end;
	my $score;
	
	#my $id_final = $readID[$id];
	my $id_final = $id;
	
	if ( $outfile_format eq "bed" ) {
		# BED 
		# - name: read ID
		# - score = prob * 1000
		
		$start = $loc;
		$end = $start + $read_length - 1;
		my $name = $id_final;
		$score = $prob * 1000;
		
		print OUT "$chrname\t$start\t$end\t$name\t$score\t$str\n";
	} elsif ( $outfile_format eq "gff" ) {
		# GFF
		# - source: "CSEM"
		# - feature: read ID
		# - score = prob * 1000
		 
		$start = $loc;
		$end = $start + $read_length - 1;
		my $source = "CSEM";
		my $feature = $id_final;
		$score = $prob * 1000;
		
		print OUT "$chrname\t$source\t$feature\t$start\t$end\t$score\t$str\t.\t.\n";	
	} else {
		print "Inappropriate output file format!\n";
		exit 1;
	}
}

close( IN );
close( OUT );
