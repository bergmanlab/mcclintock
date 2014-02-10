#!/usr/bin/perl

# Open the first (and only argument) - the fasta to split.
open(FASTA, $ARGV[0]) or die("Could not open  file.");
my $in_file = FALSE;

# Create list for RetroSeq input
my @seqnames_list;

# Loop through the file one line at a time
while(my $line = <FASTA>)
{	
	# If this is a sequence name...
	if($line =~ /^>/)
	{
		# Begin writing a new sequence file
		
		# If we have just written another sequence close that handle
		if($in_file)
		{
			close OUT;
		}

		my $seqname = substr $line, 1;
		chomp $seqname;	
		
		# Add sequence name to list for RetroSeq
		push (@seqname_list, $seqname);
		$filename = $seqname.".fa";
		
		# Create new file for output with 
		open(OUT, ">>$filename");
		$in_file = TRUE;
		
		# Write the name of the sequence in the new file
		print OUT ">$seqname\n";	
	}
	
	# If this is not sequence name then append to the current open handle
	elsif($in_file)
	{
		print OUT "$line";	
	}
	
	# If neither is true then something is broken...
	else
	{
		# Something went wrong, no sequences present.
		print "Woops, please supply one fasta as the only argument."
	}
}

# Print the element list required for RetroSeq
open(LIST, ">>elementlist");
foreach(@seqname_list)
{
	print LIST "$_\t$_.fa\n";
}
close LIST;
