									PREVSTART: foreach $prevstart (keys %{$pds{$pid}}) {
										if ($pds{$pid}{$prevstart}{$start} && !$pds{$pid}{$start}{$end}) {
											
											$previous = $pds{$pid}{$prevstart}{$start};			

											# if previous period already includes the current drug
											if ( (($drug == 0) && ($previous !~ /^0,$/))
													|| ($previous =~ /(,$drug,|^$drug,.+)/) ) {

												# give extra day to previous period (it has more drugs)
												$start = $start + 1;
											} 
											elsif ( (($drug != 0) && ($previous =~ /^0,$/ )) ) {
								
												# give extra day to new period (it has more drugs)
												$pds{$pid}{$prevstart}{$start-1} = $previous;
												delete $pds{$pid}{$prevstart}{$start};		
											}
											elsif ( ($drug != 0) && ($previous != /(,$drug,|^$drug,.+)/) ) {

												# create new 1-day overlap period
												$pds{$pid}{$start}{$start} = $previous.$drug.",";

												# adjust previous period, if it was longer than one day
												if ($prevstart < $start) {
													$pds{$pids}{$prevstart}{$start-1} = $previous;
													delete $pds{$pid}{$prevstart}{$start};										
												}

												# update start day
												$start = $start + 1;
											}		

											last PREVSTART;
										}
									}




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


