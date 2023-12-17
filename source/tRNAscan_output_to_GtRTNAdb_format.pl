#!/usr/bin/perl

use strict;
use warnings;

if ($#ARGV < 0) {
    die "\nConvert tRNAscan secondary structure file to fasta format\n",
	"\nUsage: tRNAscan_output_to_GtRTNAdb_format.pl <sec struct file(.ss)> <genome name (e.g.: Mus_musculus)> <generate mature tRNA: 0 or 1> <mito tRNA: 0 or 1> <mask anticodons with Ns: 0 or 1> \n\n";
}				

my $ss_file = $ARGV[0];

my $org_name = "";
my $mature = 0;
my $mito = 0;
my $maskAnticodon = 0;

if (scalar(@ARGV) > 1) {
    $org_name = $ARGV[1];
}
if (scalar(@ARGV) > 2) {
    $mature = $ARGV[2];
}
if (scalar(@ARGV) > 3) {
    $mito = $ARGV[3];
}
if (scalar(@ARGV) > 4) {
    $maskAnticodon = $ARGV[4];
}

open (SSFILE, $ss_file) || die "Couldn't find $ss_file\n";

my $line;
my $SeqName;
my $SeqStart = 0;
my $SeqEnd = 0;
my $SeqLen;
my $chr;
my $aaType;
my $anticodon;
my $anticodonStart = 0;
my $anticodonEnd = 0;
my $score;
my $pseudoTag = "";
my $SeqDescription;
my $intronStart = 0;
my $intronEnd = 0;
my $Seq;
my $orgSeq;
my $strand;
my $isotype;

my %isotypeCount;
my %seqCount;


while ($line = <SSFILE>) {
    if ($line =~ /^(\S+)\s+\((\d+)\-(\d+)\)\s+Length:\s(\d+)\sbp/) {
		$SeqName = $1; #chr3.trna350
		$SeqStart = $2; #1983884
		$SeqEnd = $3; #1983811
		$SeqLen = $4; #74
		$chr = (split /\./, $SeqName)[0]; #chr3
    }	
	elsif ($line =~ /^Type:\s(\S+)\s+Anticodon:\s(\S+) at (\d+)-(\d+).+Score:\s(\S+)/) {
		$aaType = $1; #Pro
		$anticodon = $2; #TGG
		$anticodonStart = $3; #33
		$anticodonEnd = $4; #35
		$score = $5; #59.8
    }
    elsif ($line =~ /pseudogene/) {
		$pseudoTag = "Pseudo";
    }
	elsif ($line =~ /^Possible intron: (\d+)-(\d+)/) {
		$intronStart = $1; #38
		$intronEnd = $2; #50
	}
	elsif ($line =~ /^Seq:\s(\S+)$/) {
		
		$Seq = $1;
		
		$orgSeq = $Seq;
		
		$line = <SSFILE>;
		$line = <SSFILE>;
		
		$isotype = $aaType . "-" . $anticodon;
		if (exists $isotypeCount{$isotype}) {
			$isotypeCount{$isotype} += 1;
		} else {
			$isotypeCount{$isotype} = 1;
		}
		
		if (exists $seqCount{$orgSeq}) {
			$seqCount{$orgSeq} += 1;
		} else {
			$seqCount{$orgSeq} = 1;
		}
		
		if ($SeqStart < $SeqEnd){
			$strand = "+";
		}else {
			$strand = "-";
			($SeqStart, $SeqEnd) = ($SeqEnd, $SeqStart);
		}
		
		if ($maskAnticodon) {
			$Seq = substr($Seq, 0, $anticodonStart - 1) . "NNN" . substr($Seq, $anticodonEnd);
		}
		else {
			if ($mature) {
				$Seq =~ s/T/U/g;
				if ($intronStart > 0) {
					$Seq = substr($Seq, 0, $intronStart - 1) . substr($Seq, $intronEnd);
					$SeqLen = length($Seq);
				}

			}
			
		}
		
		if ($mature) {
			$SeqDescription = "tRNA-$aaType-$anticodon-$isotypeCount{$isotype}-$seqCount{$orgSeq} (tRNAscan-SE ID: $SeqName) $aaType ($anticodon) $SeqLen bp mature sequence Sc: $score $chr:$SeqStart-$SeqEnd ($strand)";
		}else{
			$SeqDescription = "tRNA-$aaType-$anticodon-$isotypeCount{$isotype}-$seqCount{$orgSeq} (tRNAscan-SE ID: $SeqName) $aaType ($anticodon) $SeqLen bp Sc: $score $chr:$SeqStart-$SeqEnd ($strand)";
		}

		if ($mito) {
			$SeqDescription = "mt_".$SeqDescription;
		}
		if ($org_name ne "") {
			$SeqDescription = $org_name."_".$SeqDescription;
		}
		if ($pseudoTag ne ""){
			$SeqDescription = $SeqDescription." ".$pseudoTag;
		}
		
		print ">$SeqDescription\n"; 
		for (my $pos = 0; $pos < length($Seq); $pos += 60){
			my $tempSeq = substr(uc($Seq),$pos,60);
			print $tempSeq, "\n";
		}
		
		$SeqName = "";
		$SeqStart = 0;
		$SeqEnd = 0;
		$chr = "";
		$pseudoTag = "";
		$intronStart = 0;
		$intronEnd = 0;
		$SeqLen = 0;
		$isotype = "";
		$anticodon = "";
		$anticodonStart = 0;
		$anticodonEnd = 0;
		$score = 0.0;
		$SeqDescription = "";
		$Seq = "";
		$orgSeq = "";
		$strand = "";		
		
	}
}