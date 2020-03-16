#!/usr/bin/perl -w
use strict;

# Pass through valid-looking SAM records
# Suppress records with invalid CIGAR strings, sending warning to STDERR

while (<>) {
	if ($_ =~ /\@SQ/ || $_ =~ /\@PG/)
    {
		print $_;
		next;
	}
	my ($readName,$flag,$refName,$refOffset,$mapQuality,$cigar,$ignore1,$ignore2,$ignore3,$seq,$qual,$tag,@other) = split /\t/;
	unless (defined($seq))
    {
		warn "No sequence found in $_\n";
		next;
	}
	my $orig = $cigar;
	my $readbases = 0;
	while (length($cigar) > 0)
    {
		if ($cigar =~ /^([\d]+)([MISP\=X])(.*)/)
        {
			$readbases += int($1);
			$cigar = $3;
		}
        elsif ($cigar =~ /^([\d]+)([DNH])(.*)/)
        {
			# ignore
			$cigar = $3;
		}
        else
        {
			#warn "Failed to parse ($orig) in $_\n";
			print;
			last;
		}
	}
	my $seqlen = length($seq);
	if ($seqlen != $readbases && $orig ne "\*")
    {
		warn "$seqlen bases in $seq but $readbases bases parsed from CIGAR string $orig, $_";
	}
    else
    {
		print;
	}
}
