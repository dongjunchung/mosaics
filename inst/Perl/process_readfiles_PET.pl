###################################################################
#
#   Process read files into bin-level files (PET)
#
#   Command arguments: 
#   - infile: Input file (directory + file name)
#   - outdir: Directory of output files
#   - format: File format
#   - binsize: Bin size
#   - collapse: Maximum # of reads allowed at each position (skip if collapse=0)
#   - bychr: Construct bin-level files by chromosome? (Y or N)
#	- chrinfo: Is the file for chromosome info provided?
#	- chrfile: File name for chromosome info (chr size)
#	- @excludeChr: Chromosomes to be excluded (vector)
#
#   Supported file format: 
#   - eland_result, sam
#
#   Note
#   - chromosomes are extracted from read files, after excluding invalid lines
#   - uni-reads are assumed
#
###################################################################

#!/usr/bin/env perl;
use warnings;
use strict;
use Cwd;
use File::Basename;

my ($infile, $outdir, $format, $binsize, $collapse, $bychr, 
	$chrinfo, $chrfile, @exclude_chr) = @ARGV;

# extract only filename from $infile

my @pr = fileparse( $infile );
my $filename = $pr[0];

# remember current working directory

my $pwd = cwd();

# construct bin-level data based on chromsome info, if provided

my %bin_count = ();
my $bin_start = 0;
my $bin_stop = 0;

if ( $chrinfo eq "Y" ) {
	open CHR, "$chrfile" or die "Cannot open $chrfile\n";
	
	while (<CHR>) {
		chomp;
		my ( $chrname, $chrsize ) = split /\s+/, $_;
		
		$bin_start = 0;
		$bin_stop = int($chrsize/$binsize);
	        
		for (my $i = $bin_start; $i <= $bin_stop; $i++) {
			${$bin_count{$chrname}}[$i] = 0;
		}
	}
	
	close CHR;
}

# chromosomes to be excluded

my $ecyn;
my %ec_hash;

if ( scalar(@exclude_chr) == 0 ) {
	$ecyn = "N";
} else {
	$ecyn = "Y";
	@ec_hash{ @exclude_chr } = 0;	
}

# process PET read file

open IN, "$infile" or die "Cannot open input file $infile\n";

my %seen =();

my %seen_pos =();		# position (leftmost position, starting from 1)
my %seen_str = ();		# strand ('F' or 'R')

my $start;
my $end;

while(<IN>){
	# parse
	
	chomp;	
    
	# procee read file, based on "format" option
    
	my @parsed;
    
	if ( $format eq "eland_result" ) {
		@parsed = eland_result( $_ );
	} elsif ( $format eq "sam" ) {
		@parsed = sam( $_ );
	} else {
        	# Unsupported file format -> exit and return 1 to environment
        
        	exit 1;
	}
    
	# skip if invalid line
    
	next if ( $parsed[0] == 0 );
    
	# otherwise, process it
    
	my ($status, $id, $chrt, $pos, $str, $read_length, $prob ) = @parsed;
	my ($id_only, $pair) = split /\//, $id;
    
	# skip if chrt \in @exclude_chr
    
	if ( $ecyn eq "Y" && exists $ec_hash{ $chrt } ) {
		next;
	}
	
	# process coordinates
	
	if ( exists $seen_pos{$id_only} ) {		
		# if other pair was already observed (w.r.t. id), 
		# then write the coordinates of both reads in the pair
		
		if ( $seen_str{$id_only} ne $str ) {
			# pos1 & pos2 should have different strands and pair numbers
			
			my $pos1 = $seen_pos{$id_only};
			my $pos2 = $pos;
			
			# in the pair, write upstream pos first, 
			# then downstream pos + read length - 1.
			# [Note] eland_result: pos = leftmost position, starting from 1.
			
			#if ( $pos1 < $pos2 && $str eq "R" ) {
			#	$pos2 += $read_length - 1;
			#	$start = $pos1;
			#	$end = $pos2;
			#} elsif ( $pos1 > $pos2 && $seen_str{$id_only} eq "R" ) {
			#	$pos1 += $read_length - 1;
			#	$start = $pos2;
			#	$end = $pos1;
			#} else {
			#	print "check id $id!\n";
			#}
			if ( $str eq "R" ) {
				$pos2 += $read_length - 1;
				$start = $pos1;
				$end = $pos2;
			} elsif ( $seen_str{$id_only} eq "R" ) {
				$pos1 += $read_length - 1;
				$start = $pos2;
				$end = $pos1;
			} else {
				print "check id $id!\n";
			}
    
		    # process to bin-level files if collapse condition is satisfied
		    
		    my $id_collapse = join("",$chrt,$start,$end);
		    $seen{$id_collapse}++;
		    
		    if ( $collapse > 0 && $seen{$id_collapse} > $collapse ) {
		        next;   
		    }
    
			# update bin count
    
			if ( exists $bin_count{$chrt} ) {
				# if there is already a matrix for chromosome, update it
        
				$bin_start = int($start/$binsize) ;
				$bin_stop = int($end/$binsize) ;
				for (my $i = $bin_start; $i <= $bin_stop; $i++) {
					${$bin_count{$chrt}}[$i] += $prob;
				}
			} else {
				# if there is no matrix for chromosome yet, 
        
				if ( $chrinfo eq "Y" ) {
					# if chr info is provided, do not construct new one
	        
					next;
				} else {
					# if no chr info is provided, construct one
	        
					@{$bin_count{$chrt}} = ();
					$bin_start = int($start/$binsize) ;
					$bin_stop = int($end/$binsize) ;
					for (my $i = $bin_start; $i <= $bin_stop; $i++) {
						${$bin_count{$chrt}}[$i] += $prob;
					}
				}
			}
			
			# delete the key to save memory
			
			delete $seen_pos{$id_only};
			delete $seen_str{$id_only};
		} else {
			print "inappropriate pair for id $id.\n";
		}
	} else {
		# if other pair was not observed yet, record it
		
		$seen_pos{$id_only} = $pos;
		$seen_str{$id_only} = $str;
	}
}

close( IN );

# move to output directory

chdir($outdir);

# write bin-level files

if ( $bychr eq "N" ) {
    # genome-wide version: all chromosome in one file
    
    my $outfile = $filename."_bin".$binsize.".txt";
    open OUT, ">$outfile" or die "Cannot open $outfile\n";
    
    foreach my $chr_id (keys %bin_count) {      
        my @bin_count_chr = @{$bin_count{$chr_id}};
        
        for( my $i = 0; $i< scalar(@bin_count_chr); $i++ ){
            my $coord = $i*$binsize;
            if ( $bin_count_chr[$i] ) {
                print OUT "$chr_id\t$coord\t$bin_count_chr[$i]\n";
            } else {
                print OUT "$chr_id\t$coord\t0\n";
            }
        }       
    }
    
    close OUT;
} else {
    # chromosome version: one chromosome in each file
    
    foreach my $chr_id (keys %bin_count) {
        #my $outfile = $chr_id."_".$filename."_bin".$binsize.".txt";
        my $outfile = $filename."_bin".$binsize."_".$chr_id.".txt";
        open OUT, ">$outfile" or die "Cannot open $outfile\n";
        
        my @bin_count_chr = @{$bin_count{$chr_id}};
        
        for( my $i = 0; $i< scalar(@bin_count_chr); $i++ ){
            my $coord = $i*$binsize;
            if ( $bin_count_chr[$i] ) {
                print OUT "$chr_id\t$coord\t$bin_count_chr[$i]\n";
            }
            else {
                print OUT "$chr_id\t$coord\t0\n";
            }
        }
        
        close OUT;
    }
}


################################################################################
#                                                                              #
#                              subroutines                                     #
#                                                                              #
################################################################################

# -------------------------- subroutine: eland result -------------------------- 

sub eland_result {
    my $line = shift;  
    
    # parse line    
    
    my ($id, $seq, $map, $t3, $t4, $t5, $chrt, $pos, $str, @rest) = split /\t/, $line;
    my $read_length = length $seq;
    
    # exclude invalid lines
    
    if ( $map eq "QC" || $map eq "NM" ) {
        my @status = 0;
        return @status;
    }
    
    if ( $chrt eq "QC" || $chrt eq "NM" ) {
        my @status = 0;
        return @status;
    }
    
    # exclude multi-reads
    
    if ( $map eq "R0" || $map eq "R1" || $map eq "R2" ) {
        my @status = 0;
        return @status;
    }
    
    if ( $chrt eq "R0" || $chrt eq "R1" || $chrt eq "R2" ) {
        my @status = 0;
        return @status;
    }
    
    my @pos_sets = split( /,/, $pos );
    if ( scalar(@pos_sets) > 1 ) {
        my @status = 0;
        return @status;
    }
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $id, $chrt, $pos, $str, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: SAM -----------------------------------

sub sam {
    my $line = shift;   
    
    # parse line
    
    if ( $line =~ /^[^@].+/ ) { 
        # exclude lines starting with "@" (comments)
        
        my ($id, $bwflag, $chrt, $pos, $t2, $t3, $t4, $t5, $t6, $seq, $t7 ) = split /\t/, $line;
        $pos = int($pos);
        my $read_length = length $seq;
        
        # strand
        
        my $str;
        
        if ( $bwflag & 4 or $bwflag & 512 or $bwflag & 1024 ) {
            # exclude invalid lines
                
            my @status = 0;
            return @status;
        } elsif ( $bwflag & 16 ) {
            # reverse strand
            
            $str = "R";
        } else {
            $str = "F";
        }   
    
        # exclude multi-reads?
        
        # return
        
        my $status = 1;
    	my $prob = 1;
        my @parsed = ( $status, $id, $chrt, $pos, $str, $read_length, $prob );
        return @parsed;         
    } else {
        my @status = 0;
        return @status; 
    }
}
