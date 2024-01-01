#!/usr/bin/perl
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Perl script to calculate diploid reads on the infile.
# The infile should be in the merged_nodups form: no duplicates, >= 14 fields,
# laid out as:
#
# str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 cigar1 seq1 mapq2 cigar2 seq2
#
# The script also requires two versions of a SNP file:
# - The "chr_pos" site file lists on each line the sorted locations of the SNPs
# - The "paternal_maternal" file lists each SNP as chr:pos paternal_SNP maternal_SNP
#
# These files can be created using the script vcftotxt.awk from a phased VCF file
#
# Usage:  diploid.pl -s [chr_pos site file] -o [paternal_maternal SNP file] [infile or stream]
# Juicer version 1.5
use File::Basename;
use POSIX;
use List::Util qw[min max];
use Getopt::Std;
use vars qw/ $opt_s $opt_l $opt_d $opt_o $opt_h /;

# Check arguments
getopts('s:o:hl');

my $site_file;
my $phased_file;
my $star=0;

$" = "\t";

if ($opt_h) {
  print STDERR "Usage: diploid.pl <infile>\n";
  print STDERR " <infile>: file in intermediate format to calculate statistics on, can be stream\n";
  exit;
}
if ($opt_s) {  
  $site_file = $opt_s;
}

if ($opt_o) {
  $phased_file = $opt_o;
}

if ($opt_l) {
  $star=1;
}

if (scalar(@ARGV)==0) {
  # print STDERR "No input file specified, reading from input stream\n";
}

# Global variables for calculating statistics
my %chromosomes;
my %maternal_snps;
my %paternal_snps;

# read in SNP site file and store as multidimensional array
open FILE, $site_file or die $!;
while (<FILE>) {
  my @locs = split;
  my $key = shift(@locs);
  my $ref = \@locs;
  $chromosomes{$key} = $ref;
}
close FILE;

# read in SNP definition file and store in hashtable
open FILE, $phased_file or die $!;
while (<FILE>) {
  my @locs = split;
  my $key = $locs[0];
  # $paternal_snps{$key} = $locs[1];
  # $maternal_snps{$key} = $locs[2];
  $paternal_snps{$key} = $locs[1];
  $maternal_snps{$key} = $locs[2];  
}
close FILE;

# read in infile and find SNPs
my $checkedfile=0;
while (<>) {
  my @record = split;
  my $oldrecord = join(' ',@record);

  # holds the read assignments
  my @read1 = ();
  my @read2 = ();

  # holds the SNP assignments (nucleotides)
  my @snp1 = ();
  my @snp2 = ();

  # holds the SNP positions
  my @snppos1 = ();
  my @snppos2 = ();

  my $orgpos1 = -100;
  my $orgpos2 = -100;

  if ($checkedfile == 0) {
    # check that the file has the format expected
    # right now just check the sequence strings are in proper place
    if ($record[9] !~ /\A[acgtn]+\z/i) {
      # print STDERR "Expected DNA strings in fields 11 and 14, instead see " . $record[10] . " and " . $record[13] . "\n";
      print STDERR "Exiting.";
      exit(1);
    }
    else {
      $checkedfile = 1;
    }
  }

  ## Here may be how it deal with the changes I saw
  # set "original" position for both reads. then we can process cigar/sequence string forward.
  # ignore mitochrondria
  # First read:

  # count Ms,Ds,Ns,Xs,=s for sequence length 


  ## SHOULD CHECK CAREFULLY WHETHER I need to do $orgpos1 = $record[3] - $seqlength1 + 1;
  ## NOT DONE YET!! 26-04-2021

  my $seqlength1=0; 
  my $currstr=$record[5];

  my $where = $currstr =~ /[0-9]+[M|D|N|X|=|S|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
  while ($where > 0) {
    $seqlength1 += substr($currstr, ($where)-1, $RLENGTH - 1) + 0;
    $currstr = substr($currstr, ($where + $RLENGTH)-1);
    $where = $currstr =~ /[0-9]+[M|D|N|X|=|S|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
  }

  if (($record[1] == 0 && $record[2] ne "MT") || $star == 1) {
    $orgpos1 = $record[3];
  }
  elsif ($record[2] ne "MT") {
    # reverse strand, find original position
    $orgpos1 = $record[3];
    # $orgpos1 = $record[3] - $seqlength1 + 1;
  }  

  # count Ms,Ds,Ns,Xs,=s for sequence length 
  # my $seqlength2=0; 
  # my $currstr=$record[12];
  # my $where = $currstr =~ /[0-9]+[M|D|N|X|=|S|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
  # while ($where > 0) {
  #   $seqlength2 += substr($currstr, ($where)-1, $RLENGTH - 1) + 0;
  #   $currstr = substr($currstr, ($where + $RLENGTH)-1);
  #   $where = $currstr =~ /[0-9]+[M|D|N|X|=|S|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
  # }

  # Second read:
  # if (($record[4] == 0 && $record[5] ne "MT") || $star == 1) {
  #   $orgpos2 = $record[6];
  # }
  # elsif ($record[5] ne "MT") {
  #   # reverse strand, find original position
  #   $orgpos2 = $record[6] - $seqlength2 + 1;
  # }

  # find first read position in the SNP site array
  my $ind1 = &bsearch($orgpos1,$chromosomes{$record[2]});
  my $orgpos11 = $orgpos1;
  my $orgpos22 = $orgpos2;
  # first read might land on SNP: if on forward strand, position + length overlaps, if on
  # reverse strand, position - length overlaps.
  # ind1 is first hit; keep incrementing it until no longer on read.
  # shouldn't matter reverse or forward strand at this point, just use orgpos and cigar and
  # count forward
  while ($orgpos11 < $chromosomes{$record[2]}->[$ind1] && $orgpos11+$seqlength1 >= $chromosomes{$record[2]}->[$ind1]) {
    # need to count forward from orgpos1.
    my $currstr=$record[5];
    my $where = $currstr =~ /[0-9]+[M|D|S|I|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
    my $seqstr=$record[9];

    $orgpos1=$orgpos11;
    while ($where > 0) {
      my $char = substr($currstr, ($where)-1 + $RLENGTH - 1,1);
      my $len  = substr($currstr, ($where)-1, $RLENGTH - 1) + 0;
      my $ind;

      # match. check if match spans SNP, if so, take it. if not, advance cigar string
      # and sequence string and keep parsing cigar
      if ($char eq "M") {
        if ($orgpos1 + $len > $chromosomes{$record[2]}->[$ind1]) {
          # matching spans SNP
          $ind = $chromosomes{$record[2]}->[$ind1] - $orgpos1;
          push @snp1,substr $seqstr, $ind, 1; 
          push @snppos1,$chromosomes{$record[2]}->[$ind1];
          $where = 0;
        }
        else {
          # advance cigar and sequence and keep going.
          $seqstr = substr($seqstr, $len);
          $orgpos1 = $orgpos1 + $len;
          $currstr = substr($currstr, ($where + $RLENGTH)-1);
          $where = $currstr =~ /[0-9]+[M|D|S|I|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
        }
      }
      elsif ($char eq "D" || $char eq "N") {
        # delete. if deletion spans SNP, we don't have it in our read.
        if ($orgpos1 + $len > $chromosomes{$record[2]}->[$ind1]) {
          $where = 0;
        }
        else {
          # doesn't span SNP, update genomic position and continue.
          $orgpos1 = $orgpos1 + $len;
          $currstr = substr($currstr, ($where + $RLENGTH)-1);
          $where = $currstr =~ /[0-9]+[M|D|S|I|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
        }
      }
      elsif ($char eq "I") {
        # insertion, advance sequence string, doesn't affect genomic position
        $currstr = substr($currstr, ($where + $RLENGTH)-1);
        $where = $currstr =~ /[0-9]+[M|D|S|I|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
        $seqstr = substr($seqstr, $len);
      }
      elsif ($char eq "S" || $char eq "H") {
        if ($orgpos1 + $len > $chromosomes{$record[2]}->[$ind1]) {
          # skip spans SNP
          $where = 0;
        }
        else {
          # skip does not span SNP, advance cigar and sequence and keep going.
          # hard clipped bases do not appear in seqstr
          if ($char eq "S") {
            $seqstr = substr($seqstr, $len);
          }
          $orgpos1 = $orgpos1 + $len;
          $currstr = substr($currstr, ($where + $RLENGTH)-1);
          $where = $currstr =~ /[0-9]+[M|D|S|I|H|N]/ ? scalar($RLENGTH = length($&), $RSTART = length($`)+1) : 0;
        }
      }
    }
    $ind1++;
  }
  if (scalar @snp1 > 0) {
    for my $i (0 .. $#snp1) {
      my $key = $record[2] . ":" . $snppos1[$i];
      if ($snp1[$i] eq $paternal_snps{$key}) {
        $read1[$i] = "paternal";
      }
      elsif ($snp1[$i] eq $maternal_snps{$key}) {
        $read1[$i] = "maternal";
      }
      else {
        $read1[$i] = "mismatch";
      }
    }
  }

  if (scalar @read1 > 0) {
    print STDOUT "SNP";
    for my $i (0 .. $#snp1) {
      print STDOUT ":" . $snp1[$i] . ":" . $snppos1[$i] . ":" . $read1[$i];
    }
    print STDOUT "| ";       
    print "@record"; 
    print STDOUT "\n";
  }
}

# Binary search, array passed by reference
# search array of integers a for given integer x
# return index where found or upper index if not found
sub bsearch {
  my ($x, $a) = @_;    # search for x in array a
  my ($l, $u) = (0, @$a - 1);  # lower, upper end of search interval
  my $i;             # index of probe
  while ($l <= $u) {
    $i = int(($l + $u)/2);
    if ($a->[$i] < $x) {
      $l = $i+1;
    }
    elsif ($a->[$i] > $x) {
      $u = $i-1;
    }
    else {
      return $i; # found
    }
  }
  return $l;         # not found, return upper
}
