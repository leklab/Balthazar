#!/bin/perl


@R1_files = glob("R1.23_*.fastq");


foreach $f (@R1_files){
	$r2 = $f;
	$r2 =~ s/R1/R2/;

	print "./find_R2_23.pl $f > $r2\n";
	`./find_R2_23.pl $f > $r2`;


}
