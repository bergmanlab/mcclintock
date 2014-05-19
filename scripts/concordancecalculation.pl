#!/usr/bin/perl

use File::Basename;

# Script to calculate the numbers of concordant predictions in results from the mcclintock pipeline
# A prediction is counted as concordant if it is the same TE and is annotated within 100bp of the other method
# The output is a tab delimited file with the following columns:
# SampleID	Method1name	Method2name	Concordantnovelpredictions	Totalnovelpredictionshigh	Totalnovelpredictionslow	Concordantexistingpredictions	Totalexistingpredictionshigh	Totalexistingpredictionslow

# For each directory supplied (these should be the directories in the results folder to be analysed)
foreach my $directory (@ARGV)
{
	# Open the relevant folder, this will be the results for one sample
	opendir my $sample, "$directory/results";
	
	$sampleID = basename($directory);
	
	# Read the results files (one for each different method), ignoring dotfiles
	my @files = grep {!/^\./} readdir $sample;
	for (my $i = 0; $i < @files; $i++)
	{	
		# Strip the decimal from any coordinates in the PoPoolationTE results file that have a non integer value	
		if (@files[$i] =~ /popoolationte/)
		{
			$old = @files[$i];
			open (TMP, ">$directory/results/@files[$i].nodec");
			print TMP "track name=\"$sampleID\_PoPoolationTE_nodec\" description=\"$sampleID\_PoPoolationTE_nodec\"\n";
			system ("awk '\$2=int(\$2),\$3=int(\$3)' OFS=\"\\t\" \"$directory/results/@files[$i]\" >> \"$directory/results/@files[$i].nodec\"");
			@files[$i] = "$old.nodec";
			close (TMP);
		}
	}
	
	open (NEWOUTFILE, ">$directory/results/$sampleID.new.concordance.tsv");	
	open (OLDOUTFILE, ">$directory/results/$sampleID.old.concordance.tsv");	
	$count = 0;
	
	# Sort the files list to create more standardised output
	@files = sort @files;	
	
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
			@method1name = split(/\./, @method1[1], 2);
			@method2 = split(/_/, @files[$c], 2);
			@method2name = split(/\./, @method2[1], 2);
			
			print NEWOUTFILE "@method1name[0] \t @method2name[0]\t";
			if (@method1name[0] ne "retroseq" && @method2name[0] ne "retroseq")
			{
				print OLDOUTFILE "$sampleID \t @method1name[0] \t @method2name[0] \t";
			}
			
			# Run bedtools window to detect overlapping predictions
			@result = `bedtools window -w 100 -a $directory/results/$file1 -b $file2 2>&1`;
			$newcount = 0;
			$oldcount = 0;
			foreach my $line (@result)
			{
				my @fields = split "\t" , $line;
				# Only count predictions of the same TE as concordant
				if (@fields[3] eq @fields[9])
				{
					# Keep separate counts for new and old TE sites
					if (@fields[3] =~ /new/)
					{
						$newcount++;
					}
					if (@fields[3] =~ /old/)
					{
						$oldcount++;
					}
				}
			}
			
			# Print the count of concordant novel TE predictions
			print NEWOUTFILE "$newcount \t";
			
			# Print the total number of novel TE predictions for each method, high value first then low
			my $newpredictionsmethod1=`grep new "$directory/results/$file1" | wc -l`;
			my $newpredictionsmethod2=`grep new $file2 | wc -l`;
			
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
				my $oldpredictionsmethod1=`grep old "$directory/results/$file1" | wc -l`;
				my $oldpredictionsmethod2=`grep old $file2 | wc -l`;
				
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