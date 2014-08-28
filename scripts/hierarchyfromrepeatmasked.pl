#! /bin/perl

# Creates a hierarchy file for McClintock from a GFF created by RepeatMasker and
# the fasta TE database used to run RepeatMasker.

# Open the GFF file output from a RepeatMasker run 
open (GFFIN, $ARGV[0]);
# Open the TE Fasta consensus sequence database used with RepeatMasker
open (FASTAIN, $ARGV[1]);
# Open a Gff to be the output GFF file with corrected IDs
open (GFFOUT, ">ID_$ARGV[0]");
# Open a file that will contain a ready made hierachy file for McClintock
open (HIERARCHY, ">$ARGV[2]");

# For each line in the fasta file
while (my $TE = <FASTAIN>)
{
	# Reset to the top of the GFF input
	seek (GFFIN, 0,0);
	# If this line is the header of a TE sequence
	if ($TE =~ m/^>/)
	{
		# Save the whole ID, split on the RepeatMasker dividing character
		my @teID = split (/[>#]/, $TE);
		my $family = $teID[1];
		my $family_count = 1;
		# For every GFF record in the input
		while (my $GFF_record = <GFFIN>)
		{
			# If the record is not header line
			if ($GFF_record !~ m/^#/)
			{
				my @GFF_fields = split (/\t/, $GFF_record);
				my @exact_family = split (/[:"]/, $GFF_fields[8]);
				
				# If the search family matches the family of the GFF record
				if ($exact_family[2] eq  $family)
				{
					# Print the GFF fields to the new GFF
					print GFFOUT join("\t", @GFF_fields[0 .. 7]), "\t";
					# Print the family name as the ID with a count appended to make it unique
					print GFFOUT "ID=$family","_","$family_count\n";
					# Print the TE ID and family name to the hierarchy file
					print HIERARCHY "$family","_","$family_count\t$teID[1]\#$teID[2]";
					$family_count++;
				} 
			}
		}
	}
}

