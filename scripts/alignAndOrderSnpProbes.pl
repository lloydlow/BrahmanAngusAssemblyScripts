#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --mem=10000
#SBATCH --ntasks-per-node=1
# This script is provided by Derek M Bickhart. It is free to use for academic/non-profit research in the hopes it will be useful!
# Output: base.tab <- snpnames and mapping coordinates
# Output: base.segs <- identified chromosome segments
# Output: base.conflicts <- conflicting segments
# Output: base.stats <- general statistics
# Output: (if circos) base(folder) -> basefilenames
# version 2: 11/16/2017	added conflict file generation
# version 3: 2/2/2018 added mashmap alignment format parsing
# version 4: 3/1/2018 modifying probe alignment sequence algorithm. Tracking variant site rather than forward 5' alignment position
# version 5: 10/9/2018 generate circos data files for plotting

use strict;
use Getopt::Std;

my $usage = "perl $0 (-a <assembly fasta> -p <probe fasta>) || (-n <nucmer aligns>) || (-g <general alignment format> -m <median algn length to filter>) -c <circos file flag> ((if c) -f <query fasta fai file> && -r <reference fasta fai file>) -o <output file basename>\n";
my %opts;
getopt('apongmfr', \%opts);

unless(((defined($opts{'a'}) && defined($opts{'p'})) || defined($opts{'n'}) || defined($opts{'g'})) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $circosOut;
my %qchrs; # {chr} => length
my %rchrs; # {chr} => length
if(defined($opts{'c'})){
	unless(defined($opts{'f'}) && defined($opts{'r'})){
		print "Need option f and r defined for circos plot!\n";
		print $usage;
		exit;
	}
	open(my $IN, "< $opts{f}") || die "Could not open query fai index file!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$qchrs{$segs[0]} = $segs[1];
	}
	close $IN;
	
	open(my $IN, "< $opts{r}") || die "Could not open reference fai index file!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$rchrs{$segs[0]} = $segs[1];
	}
	close $IN;
	
	$circosOut = $opts{'o'} . "_circos";
	mkdir($circosOut) || print "$!";
}


my $unmaps = 0; my $maps = 0; 
my %aligns; # {ochr}->{opos} = [] ->[probe, chr, pos, orient]
my %qcounts; # {chr}->{ochr} = count
my %qcvals; #(only for g option) {scaffold}->{ochr}->[cum perc len, algn count]
my %scaflens; #(only for g option) {scaffold} = length

## READ INPUT ##

if(defined($opts{'a'})){
	open(my $IN, "module load bwa; bwa mem $opts{a} $opts{p} |") || die "Could not begin BWA alignments!\n";
	while(my $line = <$IN>){
		if($line =~ /^@/){
			next;
		}

		chomp $line;
		my @segs = split(/\t/, $line);
		my @rnsegs = split(/\./, $segs[0]);
		if($segs[1] & 2048){
			next;
		}
		#if($rnsegs[1] == 0){next;} # Takes care of probes without prior chromosome alignments
		if($segs[2] eq "*"){
			$unmaps++; # Count unmapped probes
		}else{
			$maps++;
		}
		$segs[2] =~ s/[\|\;]/_/g; # Convert bad characters to underscores
		my $orient = ($segs[1] & 16)? "-" : "+";
		my $alen = calcAlignLen($segs[5]);
		my $pos = ($orient eq "+")? $alen + $segs[3] + 1 : $segs[3];
		$qcounts{$segs[2]}->{$rnsegs[1]} += 1;
		push(@{$aligns{$rnsegs[1]}->{$rnsegs[2]}} ,[$rnsegs[0], $segs[2], $pos, $orient]);
	}
	close $IN;
}elsif(defined($opts{'n'})){
	open(my $IN, "< $opts{n}") || die "Could not open nucmer aligns!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$maps++;
		my $orient = ($segs[2] > $segs[3])? "-" : "+";
		$qcounts{$segs[12]}->{$segs[11]} += 1;
		push(@{$aligns{$segs[11]}->{$segs[0]}} , ["none", $segs[12], $segs[2], $orient]);
	}
	close $IN;
}elsif(defined($opts{'g'})){
	open(my $IN, "< $opts{g}") || die "could not open general alignment file!\n";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\s+/, $line);
		$maps++;
		my $orient = $segs[4];
		if($segs[3] - $segs[2] < $opts{'m'}){next;} # Skip alignments less than the designated median
		push(@{$aligns{$segs[5]}->{$segs[7]}}, ["none", $segs[0], $segs[2], $orient]);
		push(@{$aligns{$segs[5]}->{$segs[8]}}, ["none", $segs[0], $segs[3], $orient]);
		$qcvals{$segs[0]}->{$segs[5]}->[0] += ($segs[3] - $segs[2]) / $segs[1];
		$qcvals{$segs[0]}->{$segs[5]}->[1] += 1;
		$qcounts{$segs[0]}->{$segs[5]} += 1;
		$scaflens{$segs[0]} = $segs[1];
	}
	close $IN;
}

## OUTPUT routines ##

if(defined($opts{'g'})){
	# print out alignment percentage stats
	open(my $OUT, "> $opts{o}.alnstats");
	foreach my $scaff (sort {scalar(keys(%{$qcvals{$b}})) <=> scalar(keys(%{$qcvals{$a}}))} keys(%qcvals)){
		print {$OUT} "$scaff\t$scaflens{$scaff}";
		foreach my $k (sort{$qcvals{$scaff}->{$b}->[0] <=> $qcvals{$scaff}->{$a}->[0]} keys(%{$qcvals{$scaff}})){
			my $acount = $qcvals{$scaff}->{$k}->[1];
			my $perc = sprintf("%.3f", $qcvals{$scaff}->{$k}->[0]);
			print {$OUT} "\t$k:$acount:$perc";
		}
		print {$OUT} "\n";
	}
	close $OUT;
}

if(defined($opts{'c'})){
	printKaryotype(\%qchrs, \%rchrs, $circosOut);
}


open(my $OUT, "> $opts{o}.tab");
open(my $STATS, "> $opts{o}.stats");
open(my $SEGS, "> $opts{o}.segs");
open(my $CON, "> $opts{o}.conflicts");
print {$STATS} "Mapping probes: $maps\tUnmapped probes: $unmaps\n";
my %qconsensus; # {refchr}->{qchr} = refchr
# Conflict mapping
foreach my $chr (sort {scalar(keys(%{$qcounts{$b}})) <=> scalar(keys(%{$qcounts{$a}}))} keys(%qcounts)){
	# Sorted by number of alternative chromosome alignments
	if(scalar(keys(%{$qcounts{$chr}})) > 1){
		my @multAligns;
		my $totCount = 0;
		foreach my $c (sort{$qcounts{$chr}->{$b} <=> $qcounts{$chr}->{$a}}keys(%{$qcounts{$chr}})){
			# Sorted again by number of mapped segments
			if($qcounts{$chr}->{$c} > 1){
				push(@multAligns, $c);
			}
			$totCount += $qcounts{$chr}->{$c};
		}
		if($totCount == 0){$totCount = 1;}
		$qconsensus{$multAligns[0]}->{$chr} = $multAligns[0];
		if(scalar(@multAligns) > 1){
			print {$CON} "$chr";
			foreach my $c (@multAligns){
				my $ratio = sprintf("%.3f", $qcounts{$chr}->{$c} / $totCount);
				print {$CON} "\t$c:" . $qcounts{$chr}->{$c} . ":$ratio";
			}
			print {$CON} "\n";
		}
	}else{
		my @k = keys(%{$qcounts{$chr}});
		$qconsensus{$k[0]}->{$chr} = $k[0];
	}
}
close $CON;


foreach my $chr (sort{$a <=> $b} keys(%aligns)){
	my ($consensus, $values) = determineConsensus($aligns{$chr});
	print {$STATS} "Ref $chr consensus:";
	for (my $x = 0; $x < scalar(@{$consensus}); $x++){
		print {$STATS} " $consensus->[$x]:$values->[$x]";
	}
	print {$STATS} "\n";
	if(defined($opts{'g'})){
		my ($refblocks, $qblocks) = condenseAlignments($aligns{$chr}, 1000000);
		for(my $x = 0; $x < scalar(@{$refblocks}); $x++){
			my $ref = $refblocks->[$x];
			my $query = $qblocks->[$x];
			my $rlen = abs($ref->[1] - $ref->[0]);
			my $qlen = abs($query->[0] - $query->[1]);
			my $orient = ($query->[0] < $query->[1])? "+" : "-";
			print {$SEGS} "$chr\t$ref->[0]\t$ref->[1]\t$query->[2]\t$query->[0]\t$query->[1]\t$orient\t$rlen\t$qlen\n";
		}
		
		if(defined($opts{'c'})){
			printLinks(\%qchrs, \%rchrs, $chr, $refblocks, $qblocks, 1000000, $circosOut);
		}
	}else{

		my ($refblocks, $qblocks) = identifyAndCondenseSegs($aligns{$chr}, $consensus->[0], $qconsensus{$chr});
		for(my $x = 0; $x < scalar(@{$refblocks}); $x++){
			my $ref = $refblocks->[$x];
			my $query = $qblocks->[$x];
			my $rlen = abs($ref->[1] - $ref->[0]);
			my $qlen = abs($query->[0] - $query->[1]);
			my $orient = ($query->[0] < $query->[1])? "+" : "-";
			print {$SEGS} "$chr\t$ref->[0]\t$ref->[1]\t$query->[2]\t$query->[0]\t$query->[1]\t$orient\t$rlen\t$qlen\n";
		}
		
		if(defined($opts{'c'})){
			printLinks(\%qchrs, \%rchrs, $chr, $refblocks, $qblocks, 1000000, $circosOut);
		}
	}
	
	foreach my $pos (sort{$a <=> $b} keys(%{$aligns{$chr}})){
		foreach my $arrayref (@{$aligns{$chr}->{$pos}}){
			print {$OUT} join("\t", @{$arrayref});
			print {$OUT} "\t$chr\t$pos\n";
		}
	}
}
close $OUT;
close $STATS;
close $SEGS;


exit;

sub condenseAlignments{
	my ($hashref, $thresh) = @_;
	# Logic: because these are placed alignments, I can sort them based on the original reference coordinates
	# No need for complex error tolerance, but I will discard condensed alignments below a certain bp length
	#my %aligns; # {ochr}->{opos} = [] ->[probe, chr, pos, orient]
	my @refblock; # [start, end]
	my @queryblock; # [start, end, chr, orient]
	
	my @buff; my $count = 0; my $segs = 0;
	foreach my $pos (sort{$a <=> $b} keys(%{$hashref})){
		#my $query = $hashref->{$pos};
		foreach my $query (@{$hashref->{$pos}}){
			if($segs == 0){
				#initializing segments
				push(@refblock, [$pos, $pos]);
				push(@queryblock, [$query->[2], $query->[2], $query->[1], $query->[3]]);
				$segs++;
			}else{
				my $last = $queryblock[-1];
				if($last->[2] eq $query->[1] && $last->[3] eq $query->[3] 
					&& min(abs($last->[0] - $query->[2]), abs($last->[1] - $query->[2])) < $thresh){
					# This is an overlapping alignment
					$queryblock[-1]->[1] = $query->[2];
					$refblock[-1]->[1] = $pos;
					$count++;
				}else{
					# starting a new segment
					push(@refblock, [$pos, $pos]);
					push(@queryblock, [$query->[2], $query->[2], $query->[1], $query->[3]]);
					$segs++;
				}
			}
		}
	}
	return \@refblock, \@queryblock;
}

sub min{
	my ($a, $b) = @_;
	return ($a < $b)? $a : $b;
}

sub identifyAndCondenseSegs{
	my ($hashref, $consensus, $qconsensus) = @_;
	# Logic: tolerate one deviation in consensus, otherwise condense region into a block
	my @refblock; # [start, end]
	my @queryblock; # [start, end, chr, testbit]

	my @buff; my $count = 0; my $segs = 0; my $skip = 0;
	foreach my $pos (sort{$a <=> $b} keys(%{$hashref})){
		#my $query = $hashref->{$pos};
		foreach my $query (@{$hashref->{$pos}}){
			if($query->[1] eq "*"){next;} # skip unmapped segs
			push(@buff, [$pos, $query]);
			if($count < 2){
				$count++;
				# Fill the initial container buffer
				next;
			}	
			
			if(scalar(@refblock) - 1 < $segs){
				# starting a new segment
				push(@refblock, [$buff[0]->[0], $buff[0]->[0]]);
				push(@queryblock, [$buff[0]->[1]->[2], $buff[0]->[1]->[2], $buff[0]->[1]->[1]]);
			}

			# Test two consecutive probes in the middle of the window to see if they match expectations
			my $comparator = $buff[0]->[1];
			my $test1 = $buff[1]->[1];
			my $test2 = $buff[2]->[1];
			# avg pairwise distance between reference probes in this view
			my $refDist = (($buff[1]->[0] - $buff[0]->[0]) + ($buff[2]->[0] - $buff[1]->[0])) / 2;
			my $t1dist = abs($test1->[2] - $comparator->[2]);
			my $t2dist = abs($test2->[2] - $comparator->[2]);
			
			if($skip){
				$skip = 0;
			}else{
				# passed the test bit for singleton deviations in consensus
				if(($t1dist > 5 * $refDist && $t2dist > 5 * $refDist) ||
					(($test1->[1] ne $consensus && $test2->[1] ne $consensus)
					&& (!exists($qconsensus->{$test1->[1]}) && !exists($qconsensus->{$test2->[1]})))){
					$segs++; # The conditional now knows to start a new segment
				}elsif($t1dist > 5 * $refDist || ($test1->[1] ne $consensus && !exists($qconsensus->{$test1->[1]}))){
					# We don't want singletons to screw up our segments
					$skip = 1;
				}
				# Update the current segments
				$refblock[-1]->[1] = $buff[0]->[0];
				$queryblock[-1]->[1] = $buff[0]->[1]->[2];
			}

			shift(@buff); # Remove the preceeding buffer item
		}
	}
	if(scalar(@buff) < 3 && scalar(@refblock) < 1){
		# For alignments with fewer lines than chromosomes
		push(@refblock, [$buff[0]->[0], $buff[0]->[0]]);
		push(@queryblock, [$buff[0]->[1]->[2], $buff[0]->[1]->[2], $buff[0]->[1]->[1]]);
		if(scalar(@buff) > 1){
			if($buff[1]->[1]->[2] eq $consensus || exists($qconsensus->{$buff[1]->[1]->[2]})){
				$refblock[0]->[1] = $buff[1]->[0];
				$queryblock[0]->[1] = $buff[1]->[1]->[2];
			}else{
				push(@refblock, [$buff[1]->[0], $buff[1]->[0]]);
				push(@queryblock, [$buff[1]->[1]->[2], $buff[1]->[1]->[2], $buff[1]->[1]->[1]]);
			}
		}
	}else{
	# Update the final segments for this chr
	$refblock[-1]->[1] = $buff[-1]->[0];
	$queryblock[-1]->[1] = $buff[-1]->[1]->[2];
	}

	return \@refblock, \@queryblock;
}

		

sub determineConsensus{
	my ($hashref) = @_;
	# The input is all of the mapped probes from a reference chr
	# All we need to do is to determine the highest mapping percentile chr from the mapping chr
	my %chrs;
	foreach my $pos (keys(%{$hashref})){
		foreach my $query (@{$hashref->{$pos}}){
		$chrs{$query->[1]} += 1;
		}
	}
	my @consensus = sort{$chrs{$b} <=> $chrs{$a}} keys(%chrs);
	my @values = map{$chrs{$_}} @consensus;
	return \@consensus, \@values;
}

sub calcAlignLen{
	my ($cigar) = @_;
	
	my $len = 0;
	while($cigar =~ /(\d+)(\D{1})/g){
		my $count = $1;
		my $val = $2;
		if($val =~ /[MIS=X]/){
			$len += $count;
		}
	}
	return $len;
}
sub getSortedChrs{
	my ($cref) = @_;
	return sort {
		my ($c1, $x) = $a =~ /(chr)*(.+)/;
		my ($c2, $y) = $b =~ /(chr)*(.+)/;
		if($x eq "X"){
			$x = 500;
		}elsif($x eq "Y"){
			$x = 501;
		}elsif($x eq "M" || $x eq "MT"){
			$x = 502;
		}

		if($y eq "X"){
			$y = 500;
		}elsif($y eq "Y"){
			$y = 501;
		}elsif($y eq "M" || $x eq "MT"){
			$y = 502;
		}
		$x <=> $y}
	keys(%{$cref});
}

sub printKaryotype{
	my ($qref, $rref, $outbase) = @_;
	open(my $OUT, "> $outbase/$outbase.karyotype.txt");
	my @qchrs = getSortedChrs($qref);
	my @rchrs = getSortedChrs($rref);
	for(my $x = 0; $x < scalar(@qchrs); $x++){
		print {$OUT} "chr - q$x $qchrs[$x] 0 " . $qref->{$qchrs[$x]} . " grey\n";
	}
	for(my $x = 0; $x < scalar(@rchrs); $x++){
		print {$OUT} "chr - r$x $rchrs[$x] 0 " . $rref->{$rchrs[$x]} . " black\n";
	}
	close $OUT;
	return \@qchrs, \@rchrs;
}

sub printLinks{
	my ($qref, $rref, $rchr, $refblocks, $qblocks, $thresh, $outbase) = @_;
	#refblocks [start, end]
	#queryblocks [start, end, chr, testbit]
	my @qchrs = getSortedChrs($qref);
	my @rchrs = getSortedChrs($rref);
	my %qlookup = map { $_ => $qchrs[$_] } 0..$#qchrs;
	my %rlookup = map { $_ => $rchrs[$_] } 0..$#rchrs;
	
	
	my $linkNum = 0;
	my $ridxchr = "r" . $rlookup{$rchr};
	open(my $OUT, ">> $outbase/$outbase.links.txt");
	for(my $x = 0; $x < scalar(@{$refblocks}); $x++){
		my $rlen = $refblocks->[$x]->[1] - $refblocks->[$x]->[0];
		my $qlen = $qblocks->[$x]->[1] - $qblocks->[$x]->[0];
		if($rlen < $thresh || $qlen < $thresh){next;}
		
		my $qidxchr = "q" . $qlookup{$qblocks->[$x]->[2]};
		
		print {$OUT} "link$linkNum $ridxchr " . $refblocks->[$x]->[0] . " " . $refblocks->[$x]->[1] . " color=red\n";
		print {$OUT} "link$linkNum $qidxchr " . $qblocks->[$x]->[0] . " " . $qblocks->[$x]->[1] . " color=red\n";
	}
	close $OUT;
}