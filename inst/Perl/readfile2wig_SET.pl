###################################################################
#
#   Process read files into WIG files
#
#   Command arguments: 
#   - infile: Input file (directory + file name)
#   - outdir: Directory of output files
#   - format: File format
#   - span: wig span
#	- norm_const: normalizing constant
#   - L: Expected fragment length
#   - collapse: Maximum # of reads allowed at each position (skip if collapse=0)
#   - bychr: Construct bin-level files by chromosome? (Y or N)
#	- chrinfo: Is the file for chromosome info provided?
#	- chrfile: File name for chromosome info (chr size)
#	- @excludeChr: Chromosomes to be excluded (vector)
#
#   Supported file format: 
#   - eland_result, eland_extended, eland_export, bowtie, sam, bed, csem
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

my ($infile, $outdir, $format, $span, $norm_const, $L, $collapse, $bychr, 
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
		$bin_stop = int($chrsize/$span);
	        
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

# load read file and process it
# (chromosome information is extracted from read file)

open IN, "$infile" or die "Cannot open $infile\n";

my %seen =();

while(<IN>){
    chomp;
    
    # procee read file, based on "format" option
    
    my @parsed;
    
    if ( $format eq "eland_result" ) {
        @parsed = eland_result( $_, $L );
    } elsif ( $format eq "eland_extended" ) {
        @parsed = eland_extended( $_, $L );
    } elsif ( $format eq "eland_export" ) {
        @parsed = eland_export( $_, $L );
    } elsif ( $format eq "bowtie" ) {
        @parsed = bowtie( $_, $L );
    } elsif ( $format eq "sam" ) {
        @parsed = sam( $_, $L );
    } elsif ( $format eq "bed" ) {
        @parsed = bed( $_, $L );
    } elsif ( $format eq "csem" ) {
        @parsed = csem( $_, $L );
    } else {
        # Unsupported file format -> exit and return 1 to R environment
        
        exit 1;
    }
    
    # skip if invalid line
    
    next if ( $parsed[0] == 0 );
    
    # otherwise, process it
    
    my ($status, $chrt, $pos, $str, $L_tmp, $read_length, $prob ) = @parsed;
    
    # skip if chrt \in @exclude_chr
    
    if ( $ecyn eq "Y" && exists $ec_hash{ $chrt } ) {
	    next;
    }
    
    # process to bin-level files if collapse condition is satisfied
    
    my $id_collapse = join("",$chrt,$pos,$str);
    $seen{$id_collapse}++;
    
    if ( $collapse > 0 && $seen{$id_collapse} > $collapse ) {
        next;   
    }    
    
    # adjust position if it is reverse strand
        
    $pos = $pos + $read_length - $L if ( $str eq "R" );
    if ( $pos <= 0 ) {
	    $pos = 1;
    }
    
    # update bin count
    
    if ( exists $bin_count{$chrt} ) {
        # if there is already a matrix for chromosome, update it
        
        $bin_start = int(($pos-1)/$span) ;
        $bin_stop = int(($pos + $L_tmp - 1 - 1 )/$span) ;
        for (my $i = $bin_start; $i <= $bin_stop; $i++) {
            ${$bin_count{$chrt}}[$i] += $prob;
        }
    } else {
        # if there is no matrix for chromosome yet, construct one
        
        @{$bin_count{$chrt}} = ();
        $bin_start = int(($pos-1)/$span) ;
        $bin_stop = int(($pos + $L_tmp - 1 - 1 )/$span) ;
        for (my $i = $bin_start; $i <= $bin_stop; $i++) {
            ${$bin_count{$chrt}}[$i] += $prob;
        }
    }
}
    
close IN;

# move to output directory

chdir($outdir);

# write bin-level files

if ( $bychr eq "N" ) {
    # genome-wide version: all chromosome in one file
    
	my $outfile = $filename."_fragL".$L."_span".$span.".wig";
	open OUT, ">$outfile" or die "Cannot open $outfile\n";
	
	print OUT "track type=wiggle_0 name=\"".$outfile."\" ";
	print OUT "description=\"".$outfile."\"\n";
	
	foreach my $chr_id (keys %bin_count) {      
	    my @bin_count_chr = @{$bin_count{$chr_id}};
	    
	    print OUT "variableStep chrom=$chr_id span=$span\n";
	
	    for( my $i = 0; $i < scalar(@bin_count_chr); $i++ ){
	        my $coord = $i*$span + 1;	        
	        if ( $bin_count_chr[$i] ) {
		    	my $count_final = $bin_count_chr[$i] * $norm_const;
		    	print OUT "$coord $count_final\n";
	        } else {
	            print OUT "$coord 0\n";
	        }
	    }
	}
	
	close OUT;
} else {
    # chromosome version: one chromosome in each file
	
	foreach my $chr_id (keys %bin_count) {      	    
		my $outfile = $filename."_fragL".$L."_span".$span."_".$chr_id.".wig";
		open OUT, ">$outfile" or die "Cannot open $outfile\n";
		
		print OUT "track type=wiggle_0 name=\"".$outfile."\" ";
		print OUT "description=\"".$outfile."\"\n";
		
	    my @bin_count_chr = @{$bin_count{$chr_id}};
	    
	    print OUT "variableStep chrom=$chr_id span=$span\n";
	
	    for( my $i = 0; $i < scalar(@bin_count_chr); $i++ ){
	        my $coord = $i*$span + 1;	        
	        if ( $bin_count_chr[$i] ) {
		    	my $count_final = $bin_count_chr[$i] * $norm_const;
		    	print OUT "$coord $count_final\n";
	        } else {
	            print OUT "$coord 0\n";
	        }
	    }
	}
	
	close OUT;
}


################################################################################
#                                                                              #
#                              subroutines                                     #
#                                                                              #
################################################################################

# -------------------------- subroutine: eland result -------------------------- 

sub eland_result {
    my $line = shift;   
    my $L = shift;
    
    # parse line    
    
    my ($t1, $seq, $map, $t3, $t4, $t5, $chrt, $pos, $str, @rest) = split /\t/, $line;
    
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
    
    # fragment length adjustment
    
    my $read_length = length $seq;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
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
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
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
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
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
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
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
        #if ( $L < $read_length ) {
        #    $L_tmp = $read_length;
        #}
        
        # return
        
        my $status = 1;
    	my $prob = 1;
        my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
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
    
    # fragment length adjustment
    
    my $read_length = $end - $start + 1;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = 1;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}

# -------------------------- subroutine: CSEM -------------------------- 

sub csem {
    my $line = shift;   
    my $L = shift;
    
    # skip track line
    
    my ($element1, @element2) = split /\s+/, $line;
    
    if ( $element1 eq "track" ) {
        my @status = 0;
        return @status;
    }
    
    # parse line    
    
    my ($chrt, $start, $end, $t1, $score, $str_org, @rest) = split /\t/, $line;
    
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
    
    # fragment length adjustment
    
    my $read_length = $end - $start + 1;
    my $L_tmp = $L;  
    #if ( $L < $read_length ) {
    #    $L_tmp = $read_length;
    #}
    
    # return
    
    my $status = 1;
    my $prob = $score / 1000;
    my @parsed = ( $status, $chrt, $pos, $str, $L_tmp, $read_length, $prob );
    return @parsed;
}
