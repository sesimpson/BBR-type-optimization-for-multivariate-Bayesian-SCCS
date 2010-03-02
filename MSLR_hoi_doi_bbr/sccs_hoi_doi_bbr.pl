# to run, use names of drug table, condition table, desired output
#	file, e.g.
#
#	$ perl format_multivar.pl /Users/ses/Desktop sccs_drugs_per_period2.txt 
#										sccs_cases_per_period.txt
#										results.txt
#										sccs_multivar_Rin.txt
#										Rout.txt &

( $folder, $drugfile, $condfile, $outfile, $Rinfile, $Routfile, $os ) = @ARGV;

	$drugfile = $folder."/".$drugfile;
	$condfile = $folder."/".$condfile;

	$outfile = ">".$folder."/".$outfile;

	$Rscript = $folder."/sccs_hoi_doi_bbr.R";
	$Routfile = $folder."/".$Routfile;


$starttime = localtime;
print( "start, $starttime \n");		# START time
	
open(OUT, $outfile);
	print OUT "COND_ID\tDRUG_ID\tBETA_HAT\tSE_BETA_HAT\n";

	# get information from condition file
	open(COND, $condfile);
		$temp = <COND>;

		if ($temp =~ /\n/) { chomp $temp; }
		if ($temp =~ /\r/) { chop $temp; }
		@row = split("\t", $temp);

		undef %I; for($i=0; $i<@row; $i++) { $I{$row[$i]} = $i; }

		do { 
			$condline = <COND>;
			if ($condline) {
				if ($condline =~ /\n/) { chomp $condline; }
				if ($condline =~ /\r/) { chop $condline; }
				@condrow = split("\t", $condline);
	
				# load relevant information for a given drug/ae combo
				$cond = $condrow[$I{"CONDITION_CONCEPT_ID"}];
				
				# first time through the data
				if (!defined($lastcond)) { $lastcond = $cond; }

				$pid = $condrow[$I{"PERSON_ID"}];
				$start = $condrow[$I{"PERIOD_START_DATE"}];
				$end = $condrow[$I{"PERIOD_END_DATE"}];
				
				# we are still on the same condition
				if ($cond == $lastcond) {
					# will be null if there were no evts during that period
					$evt{$pid}{$start}{$end} = $condrow[$I{"PERIOD_CASE_COUNT"}];
				}
			}
			if (!$condline || ($cond != $lastcond)) {	# we just ended a block of conditions

				# if the drug info has already been recorded, remove the people who 
				#	did NOT have an occurrence of the condition ... ( maybe skip)

				# all the info for the last condition has been gathered, so we can
				# 	get the relevant drug information now
					open(DRUG, $drugfile);
						$temp = <DRUG>;

						if ($temp =~ /\n/) { chomp $temp; }
						if ($temp =~ /\r/) { chop $temp; }
						@row = split("\t", $temp);

						if (!defined(%D)) {
							 for($i=0; $i<@row; $i++) { $D{$row[$i]} = $i; }
						}

						$max_pd_drugs = 0;

						do {
							$drugline = <DRUG>;
		
							if ($drugline) {
								if ($drugline =~ /\n/) { chomp $drugline; }
								if ($drugline =~ /\r/) { chop $drugline; }
								@drugrow = split("\t", $drugline);
					
								$pid = $drugrow[$D{"PERSON_ID"}];
									
								# person has to have had at least one event of this condition
								if ($evt{$pid}) {

									$start = $drugrow[$D{"PERIOD_START_DATE"}];
									$end = $drugrow[$D{"PERIOD_END_DATE"}];
									$drug = $drugrow[$D{"DRUG_CONCEPT_ID"}];

									$pds{$pid}{$start}{$end} .= $drug.",";
										
									$pd_drugs = split(",", $pds{$pid}{$start}{$end});
					
									if ($pd_drugs > $max_pd_drugs) {
										$max_pd_drugs = $pd_drugs;		# largest num of drugs for any one period
									}
								}
							}
						} until !$drugline;
					close(DRUG);	# done reading in drugs for that condition
										

					# write out data file and 
					# 		combine periods with same pid, drugs and find lengths
					$condRinfile = $folder."/"."cond_".$lastcond."_".$Rinfile;
				
					open(R_IN, ">".$condRinfile);
						print R_IN "PID\tEVT\tOFFS";

						for ($i=1; $i <= $max_pd_drugs; $i++) {
							print R_IN "\tD".$i;
						}
						print R_IN "\n";

						# fix up the overlap periods + record
						foreach $pid (keys %pds) {
							
							# calculate for this pid
							foreach $start (keys %{$pds{$pid}}) {
								foreach $end (keys %{$pds{$pid}{$start}}) {

									$drugkey = join("\t", sort {$a <=> $b} split(",", $pds{$pid}{$start}{$end}));
	
									if (!defined($len{$drugkey})) { # add in extra "end" day
										$len{$drugkey} = 1;
									}

									$len{$drugkey} += $end - $start;	# add on the length of this period
									$evtnum{$drugkey} += $evt{$pid}{$start}{$end};
								}
							}
						
							# print out results for this pid
							foreach $drugkey (keys %len) {
								print R_IN $pid."\t";
								print R_IN $evtnum{$drugkey}."\t";
								print R_IN $len{$drugkey}."\t";
								print R_IN $drugkey."\n";
							}

							# clean up hashes for next pid
							undef %len;
							undef %evtnum;				
						}	
					close(R_IN);	

					
					# START: call out to R for analysis
						# $Rinfile = txt file with data to be read into R
						# $Rscript = name of script we want to call
						# $Routfile = name of (intermediate) file where R output is dumped
						#					so that we can get it back into perl

						$status = `R CMD BATCH --vanilla --slave --no-timing '--args $condRinfile' $Rscript $Routfile`;
				
						if ($os =~ /^windows$/) {	# windows... 
							@output = `type $Routfile`;	
						}
						else {						# unix/linux/mac
							@output = `cat $Routfile`;		
						}
		
						foreach $outline (@output) {
							print OUT "$lastcond\t$outline";
						}
					# END: call out to R for analysis

					# clean up hashes
						undef %evt;
						undef %pds;

					# record information from first row for new condition
					if ($condline) {
						$lastcond = $cond;

						$pid = $condrow[$I{"PERSON_ID"}];
						$start = $condrow[$I{"PERIOD_START_DATE"}];
						$end = $condrow[$I{"PERIOD_END_DATE"}];
		
						# will be null if there were no evts during that period
						$evt{$pid}{$start}{$end} = $condrow[$I{"PERIOD_CASE_COUNT"}];
					}
				#}
			}
		} until !$condline;
	close(COND);
close(OUT);

$endtime = localtime;
print( "end, $endtime \n");		# END time



