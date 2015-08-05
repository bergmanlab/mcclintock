#!/usr/bin/perl

use File::Basename;

# Script to calculate the numbers of concordant predictions in results from the mcclintock pipeline
# A prediction is counted as concordant if it is the same TE and is annotated within 100bp of the other method
# The output is a tab delimited file with the following columns:
# SampleID	Method1name	Method2name	Concordantnovelpredictions	Totalnovelpredictionshigh	Totalnovelpredictionslow	Concordantexistingpredictions	Totalexistingpredictionshigh	Totalexistingpredictionslow

my $window = @ARGV[0];
# For each directory supplied (these should be the directories in the results folder to be analysed)
foreach my $directory (@ARGV[1..$#ARGV])
{
	if ($directory =~ /reference/)
	{
		next;
	}
	# Open the relevant folder, this will be the results for one sample
	opendir my $sample, "$directory/results";
	
	$sampleID = basename($directory);
	
	# Read the results files (one for each different method)
	my @files = grep {/non/} readdir $sample;
	open (NEWOUTFILE, ">$directory/results/$sampleID.$window.new.concordance.tsv");	
	open (OLDOUTFILE, ">$directory/results/$sampleID.$window.old.concordance.tsv");	
	$count = 0;
	
	# Sort the files list to create more standardised output
	@files = sort @files;
	
	# Use indexes to sort into order used for rest of thesis.
	my @sortorder = qw/ 1 5 2 4 6 3/;	
	my @idx = sort { $sortorder[$a] <=> $sortorder[$b] } 0 .. $#sortorder;
	@files = @files[@idx];
	
	# For each results file
	foreach my $file1 (@files)
	{	
		# Make sure comparisons are only made once 
		for ($c = $count; $c < @files; $c ++)
		{ 
			$file2 = "$directory/results/@files[$c]";
			# Print the name of the sample being compared    
			print NEWOUTFILE "$sampleID \t";
			
			# Print the names of the methods being compared
			@method1 = split(/_/, $file1, 2);
			@method1name = split(/_non/, @method1[1], 2);
			@method2 = split(/_/, @files[$c], 2);
			@method2name = split(/_non/, @method2[1], 2);
			
			print NEWOUTFILE "@method1name[0] \t @method2name[0]\t";
			if (@method1name[0] ne "retroseq" && @method2name[0] ne "retroseq")
			{
				print OLDOUTFILE "$sampleID \t @method1name[0] \t @method2name[0] \t";
			}
			
			# Run bedtools window to detect overlapping predictions
			@result = `bedtools window -w $window -a $directory/results/$file1 -b $file2 2>&1`;
			
			$newcount = 0;
			$oldcount = 0;
			foreach my $line (@result)
			{
				my @fields = split "\t" , $line;
				my @te1 = split(/_r|_n/, @fields[3], 2);
				my @te2 = split(/_r|_n/, @fields[9], 2);
				# Only count predictions of the same TE as concordant
				if (@te1[0] eq @te2[0])
				{
					# Keep separate counts for new and old TE sites
					if (@fields[3] =~ /non-/ && @fields[9] =~ /non-/)
					{
						$newcount++;
					}
					if (@fields[3] =~ /_reference/ && @fields[9] =~ /_reference/)
					{
						$oldcount++;
					}
				}
			}
			
			# Print the count of concordant novel TE predictions
			print NEWOUTFILE "$newcount \t";
			
			# Print the total number of novel TE predictions for each method, high value first then low
			my $newpredictionsmethod1=`grep non- "$directory/results/$file1" | wc -l`;
			my $newpredictionsmethod2=`grep non- $file2 | wc -l`;
			
			chomp ($newpredictionsmethod1);
			chomp ($newpredictionsmethod2);
			
			if ($newpredictionsmethod1 >= $newpredictionsmethod2)
			{
				print NEWOUTFILE "$newpredictionsmethod1\t $newpredictionsmethod2\n";
			}
			else
			{
				print NEWOUTFILE "$newpredictionsmethod2\t $newpredictionsmethod1\n";
			}
				
			if (@method1name[0] ne "retroseq" && @method2name[0] ne "retroseq")
			{
				# Print the count of concordant exisitng TE predictions
				print OLDOUTFILE "$oldcount \t";
				
				# Print the total number of existing TE predictions for each method, high value first then low
				my $oldpredictionsmethod1=`grep _reference "$directory/results/$file1" | wc -l`;
				my $oldpredictionsmethod2=`grep _reference $file2 | wc -l`;
				
				chomp ($oldpredictionsmethod1);
				chomp ($oldpredictionsmethod2);
				
				if ($oldpredictionsmethod1 >= $oldpredictionsmethod2)
				{
					print OLDOUTFILE "$oldpredictionsmethod1\t $oldpredictionsmethod2\n";
				}
				else
				{
					print OLDOUTFILE "$oldpredictionsmethod2\t $oldpredictionsmethod1\n";
				}
			}
		}
        $count++;
	}
	closedir $sample;
	close (NEWOUTFILE);
	close (OLDOUTFILE);
}
exit;
