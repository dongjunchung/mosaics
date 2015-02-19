#!/usr/bin/env perl;

# a script to extract coordinates of fragments from SET read file

use warnings;
use strict;
use Cwd;

my ($infile, $outfile, $summaryfile, $format, $frag_length) = @ARGV;

# process SET read file

open IN, "$infile" or die "Cannot open input file $infile\n";

my $start;
my $end;

my %nline = ();
my %start_pos = ();
my %end_pos = ();
my %str_vec = ();

while(<IN>){
	# parse
	
	chomp;
	
	# procee read file, based on "format" option
	
	my @parsed;
	
	if ( $format eq "eland_result" ) {
		@parsed = eland_result( $_, $frag_length );
	} elsif ( $format eq "eland_extended" ) {
		@parsed = eland_extended( $_, $frag_length );
	} elsif ( $format eq "eland_export" ) {
		@parsed = eland_export( $_, $frag_length );
	} elsif ( $format eq "bowtie" ) {
		@parsed = bowtie( $_, $frag_length );
	} elsif ( $format eq "sam" ) {
		@parsed = sam( $_, $frag_length );
    } elsif ( $format eq "bed" ) {
        @parsed = bed( $_, $frag_length );
	} else {
		# Unsupported file format -> exit and return 1 to R environment
		
		exit 1;
	}
	
	# skip if invalid line
	
	next if ( $parsed[0] == 0 );
	
	# otherwise, process it
	
	my ($status, $chrt, $pos, $str, $L_tmp, $read_length ) = @parsed;
	
	if ( $str eq "F" ) {
		$start = $pos;
		$end = $pos + $frag_length - 1;
	} elsif ( $str eq "R" ) {
	    $start = $pos + $read_length - $frag_length;
		$end = $pos + $read_length - 1;
		if ( $start <= 0 ) { $start = 1; }
	} else {
	    print "Inappropriate strand!\n";
	}
	
	# process coordinates
    
    if ( exists $start_pos{$chrt} ) {
        # if there is already a matrix for chromosome, update it
        
        push @{$start_pos{$chrt}}, $start;
        push @{$end_pos{$chrt}}, $end;
        push @{$str_vec{$chrt}}, $str;
        
        # count proper output lines
        
        ${$nline{$chrt}}++;
    } else {
        # if there is no matrix for chromosome yet, construct one
        
        #@{$start_pos{$chrt}} = $start;
        #@{$end_pos{$chrt}} = $end;
        #@{$str_vec{$chrt}} = $str;
        
        @{$start_pos{$chrt}} = ();
        @{$end_pos{$chrt}} = ();
        @{$str_vec{$chrt}} = ();
        
        push @{$start_pos{$chrt}}, $start;
        push @{$end_pos{$chrt}}, $end;
        push @{$str_vec{$chrt}}, $str;
        
        # count proper output lines
        
        ${$nline{$chrt}} = 1;
    }	
}

close( IN );

# write processed reads

foreach my $chr_id (keys %start_pos) {
    my $outfile_chr = $outfile."_".$chr_id;
    open OUT, ">$outfile_chr" or die "Cannot open $outfile_chr\n";
     
    my @start_pos_chr = @{$start_pos{$chr_id}};
    my @end_pos_chr = @{$end_pos{$chr_id}};
    my @str_vec_chr = @{$str_vec{$chr_id}};
    
    for( my $i = 0; $i < scalar(@start_pos_chr); $i++ ) {
        print OUT "$chr_id\t$start_pos_chr[$i]\t$end_pos_chr[$i]\t$str_vec_chr[$i]\n";
    }
    
    close OUT;
}

# summary of processed file (chrID, # lines)

open OUT, ">$summaryfile" or die "Cannot open $summaryfile\n";
foreach my $chr_id (keys %nline) {
    print OUT "$chr_id\t${$nline{$chr_id}}\n";
}
close( OUT ); 

################################################################################
#                                                                              #
#                   subroutines (from mosaics package)                         #
#                                                                              #
################################################################################

# -------------------------- subroutine: eland result -------------------------- 

sub eland_result {
	my $line = shift;	
	my $L = shift;
	
	# parse line	
	
	my ($t1, $seq, $map, $t3, $t4, $t5, $chrt, $pos, $str, @rest) = split /\t/, $line;
	
	# exclude invalid lines
	
	if ( $chrt eq "QC" || $chrt eq "NM" ) {
		my @status = 0;
		return @status;
	}
	
	# exclude multi-reads
	
	if ( $chrt eq "R0" || $chrt eq "R1" || $chrt eq "R2" ) {
		my @status = 0;
		return @status;
	}
	
	my @pos_sets = split( /,/, $pos );
	if ( scalar(@pos_sets) > 1 ) {
		my @status = 0;
		return @status;
	}
	
	# fragment length adjustment
	
	my $read_length = length $seq;
	my $L_tmp = $L;  
	if ( $L < $read_length ) {
		$L_tmp = $read_length;
	}
	
	# return
	
	my $status = 1;
	my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length );
	return @parsed;
}

# -------------------------- subroutine: eland extended ------------------------ 

sub eland_extended {
	my $line = shift;	
	my $L = shift;
	
	# parse line	
	
	my ($t1, $seq, $map, $map_result) = split /\t/, $line;
	
	# exclude invalid lines
	
	if ( $map eq "QC" || $map eq "NM" || $map eq "RM" ) {
		my @status = 0;
		return @status;
	}
	
	# exclude multi-reads
	
	if ( $map ne "1:0:0" && $map ne "0:1:0" && $map ne "0:0:1" ) {
		my @status = 0;
		return @status;
	}
    
	# process chromosome, position, and strand
    
	$map_result =~ m/(chr\w+.fa):(\d+)([FR])/;
		# use "chr\w+.fa" in order to cover chrX & chrY as well
	my $chrt = $1;
	my $pos = $2;
	my $str = $3;
	
	# fragment length adjustment
	
	my $read_length = length $seq;
	my $L_tmp = $L;  
	if ( $L < $read_length ) {
		$L_tmp = $read_length;
	}
	
	# return
	
	my $status = 1;
	my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length );
	return @parsed;
}

# -------------------------- subroutine: eland export -------------------------- 

sub eland_export {
	my $line = shift;	
	my $L = shift;
	
	# parse line	
	
	my ($t1, $t2, $t3, $t4, $t5, $t6, $t7, $t8, $seq, $t9, $chrt, $pos, $str, @rest) = split /\t/, $line;
	
	# exclude invalid lines
	
	if ( $chrt eq "QC" || $chrt eq "NM" ) {
		my @status = 0;
		return @status;
	}
	
	# exclude multi-reads
	
	my @chrt_split = split( /\:/, $chrt );
	if ( scalar(@chrt_split) > 1 ) {
		my @status = 0;
		return @status;
	}
	
	# fragment length adjustment
	
	my $read_length = length $seq;
	my $L_tmp = $L;  
	if ( $L < $read_length ) {
		$L_tmp = $read_length;
	}
	
	# return
	
	my $status = 1;
	my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length );
	return @parsed;
}

# -------------------------- subroutine: bowtie default ------------------------ 

sub bowtie {
	my $line = shift;	
	my $L = shift;
	
	# parse line	
	
	my ( $read_idx, $str_org, $chrt, $pos, $seq, @rest ) = split /\t/, $line;
	
	# exclude invalid lines?
	# exclude multi-reads?
	
	# pos: 0-based offset -> adjust
	
	$pos++;
	
	# str: "+" and "-" -> convert to "F" and "R", respectively
	
	my $str;	
	if ( $str_org eq "+" ) {
		$str = "F";
	} elsif ( $str_org eq "-" ) {
		$str = "R";
	}
	
	# fragment length adjustment
	
	my $read_length = length $seq;
	my $L_tmp = $L;  
	if ( $L < $read_length ) {
		$L_tmp = $read_length;
	}
	
	# return
	
	my $status = 1;
	my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length );
	return @parsed;
}

# -------------------------- subroutine: SAM -----------------------------------

sub sam {
	my $line = shift;	
	my $L = shift;
	
	# parse line
	
	if ( $line =~ /^[^@].+/ ) {	
		# exclude lines starting with "@" (comments)
	    
		my ($t1, $bwflag, $chrt, $pos, $t2, $t3, $t4, $t5, $t6, $seq, $t7 ) = split /\t/, $line;
		$pos = int($pos);
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
		
		# fragment length adjustment
		
		my $read_length = length $seq;
		my $L_tmp = $L;  
		if ( $L < $read_length ) {
			$L_tmp = $read_length;
		}
		
		# return
		
		my $status = 1;
		my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length );
		return @parsed;       	
	} else {
		my @status = 0;
		return @status;	
	}
}

# -------------------------- subroutine: BED -------------------------- 

sub bed {
    my $line = shift;   
    my $L = shift;
    
    # skip track line
    
    my ($element1, @element2) = split /\t/, $line;
    
    if ( $element1 eq "track" ) {
        my @status = 0;
        return @status;
    }
    
    # parse line    
    
    my ($chrt, $start, $end, $t1, $t2, $str_org, @rest) = split /\t/, $line;
    
    # pos: 0-based offset -> adjust
    
    my $pos = $start;
    $pos++;
    
    # str: "+" and "-" -> convert to "F" and "R", respectively
    
    my $str;
    if ( $str_org eq "+" ) {
        $str = "F";
    } elsif ( $str_org eq "-" ) {
        $str = "R";
    }
    
    # fragment length adjustment (fixed from ver. 0.9.9)
    # [Note] BED format: The chromEnd base is not included in the display of the feature.
    # [Note] BED format covers tagAlign (BED6) format of ENCODE.
    
    #my $read_length = $end - $start + 1;
    my $read_length = $end - $start;
    my $L_tmp = $L;  
    if ( $L < $read_length ) {
        $L_tmp = $read_length;
    }
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

