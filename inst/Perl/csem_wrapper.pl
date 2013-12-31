# Wrapper for CSEM with Bowtie
# Written by Dongjun Chung, Sep. 8, 2011

#!/usr/bin/env perl;
use warnings;
use strict;
use File::Temp qw/tempfile/;
use File::Temp qw/tmpnam/;

# parse command arguments

die "Usage: perl csem_wrapper.pl [infile_name] [infile_format] [outfile_name] [outfile_format] [ref_genome] [pseudo_tags] [n_mismatch] [maxpos] [window_size] [n_iter] [n_core]" unless @ARGV == 11;

my ( $infile_name, $infile_format, $outfile_name, $outfile_format, $ref_genome, $pseudo_tags, $n_mismatch, $maxpos, $window_size, $n_iter, $n_core ) = @ARGV;

# construct ref genome file (adapted from "genRef.pl")

open ( IN, "bowtie-inspect -s $ref_genome |" ) or die "Cannot run bowtie-inspect!\n";

my $line;

my $size = 0;
my (@names, @lens) = ();
	# $size: # of chromosomes
	# @lens: chromosome size
	# @names: chromosome name

for (my $i = 0; $i < 3; $i++) {
    # skip unnecessary lines
    $line = <IN>;
}

while ( $line = <IN> ) {
    ++$size;
    chomp($line);
    my ($seqn, $name, $len) = split(/[ \t]+/, $line);
    push(@names, $name);
    push(@lens, $len);
}
close(IN);

my ($fh, $temp_reffile) = tempfile();
print $fh "$size\n";
print $fh "@lens\n";
print $fh "@names\n";
close($fh);

# extract read length from FASTA/FASTQ files

open( IN, $infile_name ) or die "Cannot open tag file!\n";

$line = <IN>;	
if ( $infile_format eq "fasta" ) {
	while ( $line =~ /^>/ ) {
		$line = <IN>;	
	}
} elsif ( $infile_format eq "fastq" ) {
	while ( $line =~ /^@/ ) {
		$line = <IN>;	
	}
} else {
	print "Inappropriate aligned read file format!\n";
	exit 1;
}
chomp($line);
my $read_length = length $line;

close( IN );

# extract read ID

open( IN, $infile_name ) or die "Cannot open tag file!\n";

my @readID = ();	
if ( $infile_format eq "fasta" ) {
	foreach $line (<IN>) { 
		chomp($line);
		if ( $line =~ /^>(\S+)/ ) {
			push @readID, $1;
		}	
	}
} elsif ( $infile_format eq "fastq" ) {
	foreach $line (<IN>) {
		chomp($line);
		if ( $line =~ /^@(\S+)/ ) {
			push @readID, $1;
		}
	}
} else {
	print "Inappropriate aligned read file format!\n";
	exit 1;
}

close( IN );

# run bowtie & csem

my $outfile_temp = tmpnam();

if ( $infile_format eq "fasta" ) {
	system( "bowtie -f -v $n_mismatch -a -m $maxpos -p $n_core --quiet --concise $ref_genome $infile_name | csem $temp_reffile $window_size $n_iter $outfile_temp > /dev/null" ) == 0 or die "Error occurs while running either bowtie or csem!"
} elsif ( $infile_format eq "fastq" ) {
	system( "bowtie -q -v $n_mismatch -a -m $maxpos -p $n_core --quiet --concise $ref_genome $infile_name | csem $temp_reffile $window_size $n_iter $outfile_temp > /dev/null" ) == 0 or die "Error occurs while running either bowtie or csem!"
} else {
	print "Inappropriate aligned read file format!\n";
	exit 1;
}

# post-process chromosome & position 

open( IN, $outfile_temp ) or die "Cannot open csem file!\n";
open( OUT, ">", $outfile_name ) or die "Cannot open output file!\n";
	
foreach $line (<IN>) {
	chomp($line);
	my @element = split( /\s/, $line );
		# assume columns are separated by some white space
		
	# check for invalid line: may cause error in exporting step
	
	if ( scalar(@element)<5 ) {
		next;
	}
	
	# post-process lines
	
	my ($id, $chr, $str, $loc, $prob) = @element;
	if ( $outfile_format ne "bed" ) {
		# first base is 0 in bowtie or BED
		# first base is 1 in table or GFF
		$loc++;
	}
	my $chrname = $names[$chr];	# translate chromosome
	
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
	
	my $id_final = $readID[$id];
	#my $id_final = $id;
	
	if ( $outfile_format eq "table" ) {
		print OUT "$id_final\t$chrname $str $loc $prob\n";
	} elsif ( $outfile_format eq "bed" ) {
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
