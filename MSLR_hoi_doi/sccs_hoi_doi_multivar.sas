* create SCCS_DATE_OF_CHANGE ;
	%macro create_sccs_date_of_change(inlib, outlib, drugs, exp_type);

		proc sql;
			CREATE TABLE &outlib..SCCS_DATE_OF_CHANGE AS
			SELECT PERSON_ID, observation_period_start_date AS CHANGE_DATE
			FROM &inlib..OBSERVATION_PERIOD
			UNION
				SELECT PERSON_ID, observation_period_end_date AS CHANGE_DATE
			FROM &inlib..OBSERVATION_PERIOD
			UNION
			SELECT PERSON_ID, DOI_START_DATE AS CHANGE_DATE
			FROM &inlib..DOI_ERA DE
			WHERE drug_exposure_type EQ "&exp_type"
			UNION
			SELECT PERSON_ID, DOI_END_DATE AS CHANGE_DATE
			FROM &inlib..DOI_ERA DE
			WHERE drug_exposure_type EQ "&exp_type"
			ORDER BY PERSON_ID, CHANGE_DATE ASC;
		quit;
		run;

	%mend create_sccs_date_of_change;



* create SCCS_PERIODS_ALL ;
	%macro create_sccs_periods_all(outlib);

		data &outlib..SCCS_PERIODS_ALL (keep=PERSON_ID PERIOD_START_DATE PERIOD_END_DATE);
			set &outlib..SCCS_DATE_OF_CHANGE;
			by PERSON_ID;
			
			retain this_start this_end CHANGE_DATE;
			
			if first.PERSON_ID then do;
				this_start = CHANGE_DATE;
			end;
	
			else do;
				this_end = CHANGE_DATE;
				PERIOD_START_DATE = this_start;
				PERIOD_END_DATE = this_end;
				this_start = this_end;
				output;
			end;
		run;

	%mend create_sccs_periods_all;



* create SCCS_PERIODS ;
	%macro create_sccs_periods(inlib, outlib);

		proc sql;
			CREATE TABLE &outlib..SCCS_PERIODS AS
			SELECT A.PERSON_ID,
				CASE 
					WHEN (A.PERIOD_START_DATE < O.observation_period_start_date) 
						THEN O.observation_period_start_date
					ELSE A.PERIOD_START_DATE
				END AS PERIOD_START_DATE,
				CASE
					WHEN (A.PERIOD_END_DATE > O.observation_period_end_date)
						THEN O.observation_period_end_date
					ELSE A.PERIOD_END_DATE
				END AS PERIOD_END_DATE 
			FROM &outlib..SCCS_PERIODS_ALL A
				INNER JOIN &inlib..OBSERVATION_PERIOD O
				ON A.PERSON_ID = O.PERSON_ID
			WHERE A.PERIOD_START_DATE < O.observation_period_end_date
				AND A.PERIOD_END_DATE > O.observation_period_start_date
				/*AND rx_data_availability IS NOT NULL*/
			ORDER BY PERSON_ID, PERIOD_START_DATE ASC;
		quit;
		run;

	%mend create_sccs_periods;




* create SCCS_DRUGS_PER_PERIOD ;
	%macro create_sccs_drugs_per_period(inlib, outlib, drugs, exp_type);

		proc sql;
			CREATE TABLE &outlib..SCCS_DRUGS_PER_PERIOD AS
			SELECT DISTINCT P.PERSON_ID,
				P.PERIOD_START_DATE,
				P.PERIOD_END_DATE, 
				DE.DOI_CONCEPT_ID

			FROM &outlib..SCCS_PERIODS P
				INNER JOIN
				&inlib..DOI_ERA DE
					ON P.PERSON_ID = DE.PERSON_ID
				INNER JOIN &outlib..&drugs D
					ON D.doi_concept_id = DE.doi_concept_id
					WHERE DE.DOI_START_DATE <= P.PERIOD_END_DATE
						AND DE.DOI_END_DATE > P.PERIOD_START_DATE
						AND DE.drug_exposure_type EQ "&exp_type"
			ORDER BY PERSON_ID, PERIOD_START_DATE ASC, DOI_CONCEPT_ID ASC;
		quit;
		run;	

	%mend create_sccs_drugs_per_period;



* create SCCS_DRUGS_PER_PERIOD2 ;
	%macro create_sccs_drugs_per_period2(outlib);

		proc sql;
			CREATE TABLE &outlib..SCCS_DRUGS_PER_PERIOD2 AS
			SELECT DISTINCT P.PERSON_ID,
				P.PERIOD_START_DATE,
				P.PERIOD_END_DATE,
				CASE 
					WHEN (R.DOI_CONCEPT_ID IS NULL) THEN 0
					ELSE R.DOI_CONCEPT_ID
				END AS DRUG_CONCEPT_ID
			FROM &outlib..SCCS_PERIODS P
				LEFT JOIN &outlib..SCCS_DRUGS_PER_PERIOD R
				ON P.PERSON_ID = R.PERSON_ID
					AND P.PERIOD_START_DATE = R.PERIOD_START_DATE
					AND P.PERIOD_END_DATE = R.PERIOD_END_DATE
			ORDER BY PERSON_ID, PERIOD_START_DATE ASC, DRUG_CONCEPT_ID ASC;
		quit;
		run;

	%mend create_sccs_drugs_per_period2;



* create SCCS_CASES_PER_PERIOD ;

	%macro create_sccs_cases_per_period(inlib, outlib, conditions, occur_type);

		proc sql;
			CREATE TABLE &outlib..SCCS_CASES_PER_PERIOD AS
			SELECT P.PERSON_ID,
				P.PERIOD_START_DATE,
				P.PERIOD_END_DATE,
				H.HOI_CONCEPT_ID,
				COUNT(H.HOI_ERA_ID) AS PERIOD_CASE_COUNT
					FROM &outlib..SCCS_PERIODS P
				LEFT JOIN
				&inlib..HOI_ERA H
					ON P.PERSON_ID = H.PERSON_ID
				INNER JOIN &outlib..&conditions C	
					ON C.hoi_concept_id = H.hoi_concept_id
			WHERE H.HOI_ERA_START_DATE > P.PERIOD_START_DATE
				AND H.HOI_ERA_END_DATE <= P.PERIOD_END_DATE	
				AND H.hoi_occurrence_type EQ "&occur_type"
			GROUP BY P.PERSON_ID,
				P.PERIOD_START_DATE,
				P.PERIOD_END_DATE,
				H.HOI_CONCEPT_ID
			ORDER BY HOI_CONCEPT_ID, PERSON_ID, PERIOD_START_DATE, PERIOD_END_DATE;
		quit;
		run;

	%mend create_sccs_cases_per_period;


* write out SCCS table to tab delimited text ;
	%macro to_text(lib, folder, tablename);

		proc export
			data = &lib..&tablename
			outfile = "&folder./&tablename..txt"
			dbms = tab replace;
		run;

	%mend to_text;



* read in table from tab delimited text ;
	%macro from_text(lib, folder, tablename);
		
		proc import
			datafile = "&folder./&tablename..txt"
			out = &lib..&tablename
			dbms = tab replace;
		run;

	%mend from_text;
	

	
* delete original table ;
	%macro del_table(lib, tablename);
	
		proc datasets lib=&lib nolist nodetails;
			delete &tablename;
		quit;
		
	%mend del_table;



* run without chunking ;	
%macro onetable_SCCS_multivar(infolder, inlib, outfolder, outlib, 
	conditions, occur_type, drugs, exp_type, del, syscall);

*	libname &inlib "&infolder";
	/* connect to oracle, infolder is path */
	libname &inlib oracle path="LSOMOP.CHNTVA1-DC2.CSCEHUB.COM"
		defer=no schema=mslr_cdm user=ssimps pass='Omop11'
		update_lock_type=row;

	libname &outlib "&outfolder";

	%if "&drugs" eq "all" %then %do;
		proc sql;
			CREATE TABLE &outlib..all_drugs AS
			SELECT DISTINCT doi_concept_id 
			FROM &inlib..doi_era;
		quit;
		run;
			
		%let drugs = all_drugs;
	%end;
	%else %do;
		%from_text(&outlib, &outfolder, &drugs);
		run;
	%end;
		
	%if "&conditions" EQ "all" %then %do;
		proc sql;
			CREATE TABLE &outlib..all_conditions AS
			SELECT DISTINCT hoi_concept_id 
			FROM &inlib..hoi_era;
		quit;
		run;
			
		%let conditions = all_conditions;
	%end;
	%else %do;
		%from_text(&outlib, &outfolder, &conditions);
		run;
	%end;


	%create_sccs_date_of_change(&inlib, &outlib, &drugs, &exp_type);
	run;

	%create_sccs_periods_all(&outlib);
	run;

	%create_sccs_periods(&inlib, &outlib);
	run;

	%create_sccs_drugs_per_period(&inlib, &outlib, &drugs, &exp_type);
	run;

	%create_sccs_drugs_per_period2(&outlib);
	run;

	/* remove the initial table */
	%del_table(&outlib, sccs_drugs_per_period);
	run;

	%create_sccs_cases_per_period(&inlib, &outlib, &conditions, &occur_type);
	run;

	/* export out drug and condition tables to text */
		%to_text(&outlib, &outfolder, sccs_cases_per_period);
		run;

		%to_text(&outlib, &outfolder, sccs_drugs_per_period2);
		run;


	%if "&del" EQ "deletetables" %then %do;
		%del_table(&outlib, sccs_date_of_change); run;
		%del_table(&outlib, sccs_periods_all); run; 
		%del_table(&outlib, sccs_periods); run; 
		%del_table(&outlib, sccs_drugs_per_period2); run;
		%del_table(&outlib, sccs_cases_per_period); run;
	%end;

	
	%if "&syscall" EQ "systemcall" %then %do;
		* call out to system to run perl script, which will perform the analyses ;

		x "nohup perl sccs_multivar_package.pl &outfolder sccs_drugs_per_period2.txt sccs_cases_per_period.txt sccs_multivar_RESULTS.txt sccs_multivar_Rin.txt Rout.txt &";
	%end;
	
	endsas;
	
%mend onetable_SCCS_multivar;


	
	
* sample invocation ;	
/*
	%onetable_SCCS_multivar(/omop1/dev/ssimpson/SMALL_OSIM, sm_osim, /omop1/dev/ssimpson/SMALL_OSIM, sm_osim,
		10_osim_conds, 65, 10_osim_drugs, 7, nodeletetables, systemcall);
	run;
*/
	
	%onetable_SCCS_multivar(lsomop.chnvta1-dc2.cscehub.com, MSLR_in, 
		/omop1/dev/ssimpson/MSLR_hoi_doi, MSLR_out,
		all, 65, all, 7, nodeletetables, nosystemcall);
	run;


