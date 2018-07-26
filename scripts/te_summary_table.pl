#!/usr/bin/perl

use strict; use warnings; use Getopt::Long; use Data::Dumper;

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fIn2,$fOut);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"d:s"=>\$fIn,
                "t:s"=>\$fIn2,
				) or &USAGE;
&USAGE unless ($fIn);

sub USAGE {#
	my $usage=<<"USAGE";
Description:	this program is used to generate look up table for embl id, flybase id and source id
Usage:
  Options:
  -d <file>  input dir, required
  -t <file>  te list, required
  -o <file>  output file, not required  
  -h         Help

USAGE
	print $usage;
	exit;
}

if ($fOut) {
    open (OUT, '>', $fOut);
} else {
    open(OUT, ">&STDOUT");
}

my @te_list=();
my @match_array=();
my @methods=qw(ngs_te_mapper temp popoolationte telocate retroseq relocate);
my @methods_used=();
my @results=();


# read line of TE file from standard input
open (TE, '<', $fIn2) or die $!;
while (my $line = <TE>) {
    chomp $line;
	push (@te_list, $line);
}
close TE;

my $join_te = join( '|', map quotemeta, @te_list );
my $re_all = qr/\b(?:$join_te)_\w+\b/;
my $re_ref = qr/\b(?:$join_te)_ref\w+\b/;
my $re_nonref = qr/\b(?:$join_te)_no\w+\b/;

# read files from input dir
my @files = glob "$fIn/*nonredundant.bed";
my $methods_len = @methods;
for (my $i=0; $i<$methods_len; $i++) {
    my @file_match = grep /$methods[$i]/, @files;
    if (@file_match) {
        push (@methods_used, $methods[$i]);
        my %hash_all = map {($_, 0)} @te_list;
        my %hash_ref = map {($_, 0)} @te_list;
        my %hash_nonref = map {($_, 0)} @te_list;
        my $file = $file_match[0];
        open (IN, '<', $file) or die $!;
        my $dummy=<IN>;
        while (<IN>) {
            my $string = $_;
            my @matches_all = $string =~ /($re_all)/g;
            if (@matches_all) {
                my $match_all = $matches_all[0];
                $match_all =~ s/_.*//g;
                $hash_all{$match_all} = $hash_all{$match_all} + 1;
            }
            my @matches_ref = $string =~ /($re_ref)/g;
            if (@matches_ref) {
                my $matches_ref = $matches_ref[0];
                $matches_ref =~ s/_.*//g;
                $hash_ref{$matches_ref} = $hash_ref{$matches_ref} + 1;
            }
            my @matches_nonref = $string =~ /($re_nonref)/g;
            if (@matches_nonref) {
                my $match_nonref = $matches_nonref[0];
                $match_nonref =~ s/_.*//g;
                $hash_nonref{$match_nonref} = $hash_nonref{$match_nonref} + 1;
            }
        }
        push (@results, \%hash_all);
        push (@results, \%hash_ref);
        push (@results, \%hash_nonref);
        close IN;
    }
}

my $suffix_all = '_all';
my $suffix_ref = '_ref';
my $suffix_nonref = '_nonref';

# my @methods_used=qw(ngs_te_mapper temp popoolationte telocate retroseq relocate);
my @methods_all = map "$_$suffix_all", @methods_used;
my @methods_ref = map "$_$suffix_ref", @methods_used;
my @methods_nonref = map "$_$suffix_nonref", @methods_used;

my @methods_out=();
for my $i ( 0 .. $#methods_used ) {
  push (@methods_out, $methods_all[$i], $methods_ref[$i], $methods_nonref[$i]);
}

my $methods_string=join(",",@methods_out);

# output
print OUT "TE Family,$methods_string\n"; 
my $te_len = @te_list;
my $results_len = @results;
for (my $i=0; $i<$te_len; $i++) {
    my $te = $te_list[$i];
    my @results_out=();
    for (my $j=0; $j<$results_len; $j++) {
        push (@results_out, "$results[$j]{$te}");
    }
    my $results_out_join=join(",",@results_out);
    print OUT "$te,$results_out_join\n";
}

close OUT;