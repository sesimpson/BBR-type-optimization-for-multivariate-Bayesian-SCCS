# to run, use names of drug table, condition table, desired output
#	file, e.g.
#
#	$ perl format_multivar.pl sm_osim_sccs_drugs_per_period2.txt 
#										sm_osim_sccs_cases_per_period.txt
#										format_multivar_OUT_12drugs.txt
#

( $drugfile, $condfile, $outfile ) = @ARGV;

@user_cond = ( 373474 );
@user_drug = ( 
		332418, 721724, 722031, 725131, 731188, 
		740560, 781039, 789578, 903963, 904501, 
		904542, 913782, 918906, 929887, 940864, 
		941258, 956874, 961047, 974166, 975125, 
		1103006, 1103314, 1107830, 1112807, 1118084, 
		1124957, 1125315, 1137529, 1141018, 1154332, 
		1154343, 1163944, 1237049, 1307046, 1310149, 
		1322184, 1328165, 1341927, 1361711, 1501700, 
		1516800, 1518148, 1537655, 1545958, 1551099, 
		1557272, 1729720, 1734104, 1786621, 1797513, 
		19009540, 19035704, 19049105
 );

#@user_cond = ( 951 );
#@user_drug = ( 1 .. 50);
#@user_drug = ( 2406, 1687, 3846, 3209, 1749, 866, 2243, 2441, 3284, 4527, 2231, 688 );
#@user_drug = (	2406, 1687, 3846 );

@condlist{@user_cond} = (1) x @user_cond;
@druglist{@user_drug} = (1) x @user_drug;


# routine to do cho(m)ping
sub lineFix {
	$lineStr = shift @_;
	if ($lineStr =~ /\n/) { chomp $lineStr; }
	if ($lineStr =~ /\r/) { chop $lineStr; }
	@line = split("\t", $lineStr);
	return @line;
}


$starttime = localtime;
print( "start, $starttime \n");		# START time

$outfile = ">".$outfile;
#open(OUT, $outfile);

#	print OUT "PID\tEVT\tOFFS\tDRUGS\n";
	
	# get information from condition file
	open(COND, $condfile);
		$temp = <COND>;

		if ($temp =~ /\n/) { chomp $temp; }
		if ($temp =~ /\r/) { chop $temp; }
		@row = split("\t", $temp);

		undef %I; for($i=0; $i<@row; $i++) { $I{$row[$i]} = $i; }

		do { 
			$line = <COND>;
			if ($line) {
				if ($line =~ /\n/) { chomp $line; }
				if ($line =~ /\r/) { chop $line; }
				@row = split("\t", $line);
	
				# load relevant information for a given drug/ae combo
				$cond = $row[$I{"CONDITION_CONCEPT_ID"}];

				if ($condlist{$cond}) {
					$pid = $row[$I{"PERSON_ID"}];
					$start = $row[$I{"PERIOD_START_DATE"}];
					$end = $row[$I{"PERIOD_END_DATE"}];

					# will be null if there were no evts during that period
					$evt{$cond}{$pid}{$start}{$end} = $row[$I{"PERIOD_CASE_COUNT"}];
				}
			}
		} until !$line;
	close(COND);


	# get information from drug file
	open(DRUG, $drugfile);
		$temp = <DRUG>;

		if ($temp =~ /\n/) { chomp $temp; }
		if ($temp =~ /\r/) { chop $temp; }
		@row = split("\t", $temp);


		undef %I; for($i=0; $i<@row; $i++) { $I{$row[$i]} = $i; }

		$max_pd_drugs = 0;

		do {
			$line = <DRUG>;
		
			if ($line) {
				if ($line =~ /\n/) { chomp $line; }
				if ($line =~ /\r/) { chop $line; }
				@row = split("\t", $line);

				$pid = $row[$I{"PERSON_ID"}];

				$found = '';

				# only record if this person had at least one relevant event
				COND: foreach $cond (@user_cond) {
					if ($evt{$cond}{$pid}) { 
						$found = 1; 
						last COND;
					}
				}

				if ($found) {
					$start = $row[$I{"PERIOD_START_DATE"}];
					$end = $row[$I{"PERIOD_END_DATE"}];
					$drug = $row[$I{"DRUG_CONCEPT_ID"}];
					
					if (!defined($druglist{$drug})) { 
						$drug = 0; 
						$zero_drug{$pid} = 1;
					}
					$pds{$pid}{$start}{$end}{$drug} = 1;

					# record if person has at least one on-drug period
					if ($drug != 0) {
						$on_drug{$pid} = 1;
					}	

					$pd_drugs = keys %{$pds{$pid}{$start}{$end}};
					if ($zero_drug{$pid}) {
						$pd_drugs--;
					}

					if ($pd_drugs > $max_pd_drugs) {
						$max_pd_drugs = $pd_drugs;		# largest num of drugs for any one period
					}
				}
			}
		} until !$line;
	close(DRUG);


# merge periods with same pid, drugs and records number of evts during corresponding periods
	foreach $cond (@user_cond) {
		foreach $pid (keys %{$evt{$cond}}) {
			if ($on_drug{$pid}) {

				foreach $start (keys %{$pds{$pid}}) {
					foreach $end (keys %{$pds{$pid}{$start}}) {

						#$len_name = "length{".$pid."}";
						#$cond_name = "num_evts{".$cond."}{".$pid."}";

						@pd_drugs = sort {$a <=> $b} keys %{$pds{$pid}{$start}{$end}};
						
						if (($pd_drugs[0] == 0) && (@pd_drugs > 1)) {
							shift( @pd_drugs );				# remove the initial "zero" from the list
						}											# if other drugs are also present

						#foreach $drug (sort {$a <=> $b} keys %{$pds{$pid}{$start}{$end}}) {
							$pd_drugs_str = join("\t", @pd_drugs); 
							#$len_name = $len_name."{".$drug."}";
							#$cond_name = $cond_name."{".$drug."}";
						#}
					
						$length{$pid}{$pd_drugs_str} += $end-$start;
						$num_evts{$cond}{$pid}{$pd_drugs_str} += $evt{$cond}{$pid}{$start}{$end};

						#${$len_name} += $end-$start;									# add on to the length of period
						#${$cond_name} += $evt{$cond}{$pid}{$start}{$end};		# add onto num evts for period
					}
				}

			}
		}	
	}	


	# loop through and write out to text (rolled up format)

	open(OUT, $outfile);
		print OUT "PID\tEVT\tOFFS";

		for ($i=1; $i <= $max_pd_drugs; $i++) {
			print OUT "\tD".$i;
		}

		print OUT "\n";


		foreach $cond (@user_cond) {
			foreach $pid (keys %{$num_evts{$cond}}) {
				if ($on_drug{$pid}) {
						
					foreach $pd_drugs (sort {$a <=> $b} keys %{$length{$pid}}) {
						print OUT $pid."\t";
						print OUT $num_evts{$cond}{$pid}{$pd_drugs}."\t";
						print OUT $length{$pid}{$pd_drugs}."\t";
						print OUT $pd_drugs."\n";
					}
				}
			}
		}
	close(OUT);	




#	# loop through everything and write it out to text
#	foreach $cond (@user_cond) {
#		foreach $pid (keys %{$evt{$cond}}) {


#		# if this person does not have any on-drug periods, they will not 
#		# contribute to the likelihood and should be skipped
#		if ($on_drug{$pid}) {

#				foreach $start (keys %{$pds{$pid}}) {
#					foreach $end (keys %{$pds{$pid}{$start}}) {
#						print OUT $pid."\t";
	
#						$numevts = $evt{$cond}{$pid}{$start}{$end};
#						if (!$numevts) {
#							print OUT "0\t";
#						}					
#						else {
#							print OUT $numevts."\t";
#						}

#						$length = $end-$start;
#						print OUT $length;

#						$on_drug = '';
#						$pd_drugs = keys %{$pds{$pid}{$start}{$end}}; 
#						foreach $drug (keys %{$pds{$pid}{$start}{$end}}) {
#							if ($drug != 0) {
#								$on_drug = 1;
#								print OUT "\t".$drug;
#							}
#						}
					
#						if ($on_drug) {
#							foreach (1..($max_pd_drugs-$pd_drugs)) {
#								print OUT "\t";
#							}
#						}
#						else  {
#							print OUT "\t0";
							
#							foreach (1..($max_pd_drugs-1)) {
#								print OUT "\t";
#							}
#						}
						
#						print OUT "\n";
#					} # end END
#				} # end START

#			}	# end if (on-drug)
#		} # end PID
#	} # end COND

#close(OUT);


$endtime = localtime;
print( "end, $endtime \n");		# END time

