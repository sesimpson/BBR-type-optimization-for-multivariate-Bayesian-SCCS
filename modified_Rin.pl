
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
							undef $lstdrugkey;
		
							foreach $start (sort {$a <=> $b} keys %{$pds{$pid}}) {
								foreach $end (sort {$a <=> $b} keys %{$pds{$pid}{$start}}) {

									$drugkey = $pds{$pid}{$start}{$end};
								
									if (!$lstdrugkey) {
										$lstdrugkey = $drugkey;
										$lststart = $start;
										$lstend = $end;
									}
									elsif ($lstdrugkey == $drugkey) {
										$lstend = $end;
									}
									else {	# $lstdrugkey != $drugkey

										@last = split(",", $lstdrugkey);
										@this = split(",", $drugkey);
	
										if (@last > @this) {
											$thisinlast = 1;
	
											THIS: foreach $this (@this) {
												if ($lstdrugkey !~ /(^$this|,$this),/) {
													$thisinlast = '';
													last THIS;
												}
											}		
										}
										else {
											$lastinthis = 1;
				
											LAST: foreach $last (@last) {
												if ($drugkey !~ /(^$last|,$last),/) {
													$lastinthis = '';
													last LAST;
												}
											}		
										}							
	
										if (($drugkey == "0,") || $thisinlast ) {
											$start++;
										}
										elsif (($lstdrugkey == "0,") || $lastinthis ) {
											$lstend--;
										}
										else {	# create new one day overlap period
											foreach (@last, @this) {
												if (!$combined{$_}) { $combined{$_} = 1; }
											}
											$overlap = join(",", sort {$a <=> $b} keys %combined);
	
											$adjpds{$lstend}{$start} = $overlap;
											$lstend--;
											$start++;
										}
	
										$adjpds{$lststart}{$lstend} = $lstdrugkey;
	
										$lstdrugkey = $drugkey;
										$lststart = $start;
										$lstend = $end;
									}	# END: $lstdrugkey != $drugkey

								}	# END: $end
							}	# END: $start
	
							# collect last period for this person
							$adjpds{$lststart}{$lstend} = $lstdrugkey;

							# sum up periods lengths by drugkey
							foreach $start (keys %adjpds) {
								foreach $end (keys %{$adjpds{$start}}) {
									$drugkey =  join("\t", sort {$a <=> $b} split(",", $adjpds{$start}{$end}));

									$len{$drugkey} += (1 + $end - $start);	# add on the length of this period
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
							undef %adjpds;
						}
					close(R_IN);	
