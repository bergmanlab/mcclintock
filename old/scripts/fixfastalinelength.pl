#!/usr/bin/perl

# Usage: Enter none standard width fasta as first input, desired width as 
# second input and output file name as third input.
# Program reads each fasta sequence in to memory so make sure you have 
# enough.

open (FASTAIN, $ARGV[0]);
my $linelength = $ARGV[1];
open (FASTAOUT, ">$ARGV[2]");

# Set variable to indicate special first case.
my $first="true";
my $sequence;
while (my $line = <FASTAIN>)
{	
	# If this line is header
	if ($line=~m/^\>/)
	{
		# If this is not the first fasta entry
		if (!$first)
		{
			# Print the previous sequence
			# While there is content in the sequence split the sequence
			# in to linelength size pieces and print
			while (my $outputline = substr($sequence, 0, $linelength, ""))
			{
				print FASTAOUT "$outputline\n";
			}
		}

		# We have now passed the first header so set first to 0 (false)
		$first=0;
        	
		# Print the current header
		print FASTAOUT "$line";
		# Reset the sequence variable for this fasta entry
		my $sequence;
	}
	# If this is not header
	else 
	{
		# Remove the new line character
		chomp($line);
		# Add the current sequence line to the whole sequence 
		$sequence.=$line;
	}
}
# Once the file has all been read print the final sequence 
# While there is content in the sequence split the sequence
# in to linelength size pieces and print
while (my $outputline = substr($sequence, 0, $linelength, "")) 
{
	print FASTAOUT "$outputline\n";
}
