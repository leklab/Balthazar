#!/usr/bin/perl



while(<>){
	if($i%4 == 0){
		$line = $_;
		chomp $line;
		#print "key: $line\n";
		@fields = split(/\s+/,$line);
		$read_id{$fields[0]} = 1;
	}

	$i++;

	#if($i > 28){
	#	last;
	#}
}

#foreach $k (keys %read_id){
	#print "key: $k\n";

#}


open(F, "<2pool3pool_R2_001.fastq");
$i = 0;
$print_read = 0;

while(<F>){

	if($i%4 == 0){
		$line = $_;
		#chomp $line;
		@fields = split(/\s+/,$line);

		if(defined($read_id{$fields[0]})){
			print $line;
			$print_read = 1;
		}
		else{
			$print_read = 0;
		}
	}
	elsif($print_read == 1){
		print $_;
	}


	$i++;

	#if($i > 1000){
	#	last;
	#}

}

close(F);
