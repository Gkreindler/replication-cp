 * This do file: Master do file - cleans the SP data and calls all other do files
 * This version: 26 July, 2017
 * Authors: Ashwin MB
 
 *Things to do
 * 1. Remove duplicates values
 * 2. Call all the do files
 * 3. Merge with Recruit file

clear all
pause off
set more off

* Path based on global set in C:\ado\profile.do
	local path "$congestion_pricing_path"
	local pathcode "$cpcode_path"
	cd "`path'"


*** experimental sample
	use "data/treatment/treatment roster noPII.dta", clear
	keep uidp datestart
	rename uidp uidp_original

	* meet date and start date
	foreach va in datestart{
		cap drop temp
		gen temp = "2017-" + substr(`va',2,5)
		drop `va'
		gen `va' = date(temp, "YMD")
		format `va' %td
	}

	rename datestart datestart_experiment

	gen in_exp = 1
	tempfile exp_uidps
	save	`exp_uidps'

*** Which baseline survey to keep
	// import excel using "data/coded_cto/baseline/oneuidp_onesurvey_final.xlsx", clear firstrow
	// recode Keep (.=0)
	// rename Keep dum_keep
	// sort uidp_original time1 starttime
	// keep key time1 dum_keep
	// save "data/temp/baseline/oneuidp_onesurvey.dta", replace


*** Load
	use  "data/coded_cto/baseline_pooled.dta", clear


** Define the two sequences of the survey
	gen old_version = form <= 3
	assert subdate < mdy(4,1,2017) == old_version

** find out experimental respondents in the baseline data
	merge m:1 uidp_original using `exp_uidps'
	drop if _m==2
	drop _m

	gunique uidp_original
	gunique uidp_original if old_version == 0

	gunique uidp_original if in_exp == 1
	gunique uidp_original if in_exp == 1 & old_version == 0 

	duplicates report uidp_original if old_version == 0

**********************
*** CODE VARIABLES ***
**********************
		order am_group_intro_confirm, after(dest_flex)

	* departure times
		drop am_dep_time_? 
		gen am_dep_time_0 = am_dep_time_0_hh + am_dep_time_0_mm/60
		gen am_dep_time_1 = am_dep_time_1_hh + am_dep_time_1_mm/60
		gen am_dep_time_mean = (am_dep_time_0 + am_dep_time_1)/2
		drop am_dep_time_0_mm am_dep_time_0_hh am_dep_time_1_mm am_dep_time_1_hh am_dep_time_mean_mm am_dep_time_mean_hh
		order am_dep_time_0 am_dep_time_1 am_dep_time_mean, after(am_group_intro_confirm)

		drop pm_dep_time_? 
		gen pm_dep_time_0 = pm_dep_time_0_hh + pm_dep_time_0_mm/60
		gen pm_dep_time_1 = pm_dep_time_1_hh + pm_dep_time_1_mm/60
		gen pm_dep_time_mean = (pm_dep_time_0 + pm_dep_time_1)/2
		drop pm_dep_time_0_mm pm_dep_time_0_hh pm_dep_time_1_mm pm_dep_time_1_hh pm_dep_time_mean_mm pm_dep_time_mean_hh
		order pm_dep_time_0 pm_dep_time_1 pm_dep_time_mean, after(pm_group_intro_confirm)

	* duration
		drop am_mean_dur pm_mean_dur
		gen am_dur_mean = (am_dur_0 + am_dur_1)/2
		order am_dur_mean, after(am_dur_1)
		gen pm_dur_mean = (pm_dur_0 + pm_dur_1)/2
		order pm_dur_mean, after(pm_dur_1)

	* beliefs
		label var am_delta_t "AM time diff for belief question (15 or 30 minutes)"
		label var pm_delta_t "PM time diff for belief question (15 or 30 minutes)"

		* question asked in terms of "differnece to usual (as opposed to total time, as we did later)"
		gen am_beliefs_ask_delta = am_dep_delta_earlier != . & form == 2
		label var am_beliefs_ask_delta "AM question asked as difference"
		recode *_dep_earlier *_dep_later (99=.)
		drop am_dep_delta_earlier am_dep_delta_later
		rename am_dep_earlier 	am_beliefs_earlier
		rename am_dep_later 	am_beliefs_later
		rename am_dep_total_earlier am_beliefs_dur_earlier
		rename am_dep_total_later 	am_beliefs_dur_later
		order am_beliefs_ask_delta, after(am_beliefs_dur_later)

		gen pm_beliefs_ask_delta = pm_dep_delta_earlier != . & form == 2
		label var pm_beliefs_ask_delta "PM question asked as difference"
		recode *_dep_earlier *_dep_later (99=.)
		drop pm_dep_delta_earlier pm_dep_delta_later
		rename pm_dep_earlier 	pm_beliefs_earlier
		rename pm_dep_later 	pm_beliefs_later
		rename pm_dep_total_earlier pm_beliefs_dur_earlier
		rename pm_dep_total_later 	pm_beliefs_dur_later
		order pm_beliefs_ask_delta, after(pm_beliefs_dur_later)

*********************	
*** value of time ***
*********************
		/* 	VOT coding procedure: */
		order am_vot_101, before(am_vot_100)
		recode am_vot_101 (999=.)
		recode am_vot_100 am_vot_?? (99=.)

	* replacing Dont know (999/ 99) as missing values
		order pm_vot_101, before(pm_vot_100)
		recode pm_vot_101 (999=.)
		recode pm_vot_?? (99=.)

	* frequency of Sometimes
		foreach suf in 10 15 20 25 30 40 50 60 70 80 90 100{
			gen am_vot_`suf'_sometimes = am_vot_`suf' == 1
			gen am_vot_`suf'_charge    = am_vot_`suf' == 2
		}
		egen am_vot_any_sometimes = rowtotal(am_vot_*_sometimes)

		foreach suf in 10 15 20 25 30 40 50 60 70 80 90 100{
			gen pm_vot_`suf'_sometimes = pm_vot_`suf' == 1
			gen pm_vot_`suf'_charge    = pm_vot_`suf' == 2
		}
		egen pm_vot_any_sometimes = rowtotal(pm_vot_*_sometimes)

		foreach suf in 10 15 20 25 30 40 50 60 70 80 90 100{
			drop am_vot_`suf'_sometimes
			drop am_vot_`suf'_charge 
			drop pm_vot_`suf'_sometimes 
			drop pm_vot_`suf'_charge 
		}
		
	* earlier surveys didnt have am/pm_vot_101, replacing those values with 100
		replace am_vot_101 = 100 if (am_vot_100 == 2 | am_vot_100 == 1) & am_vot_101 == .
		replace pm_vot_101 = 100 if (pm_vot_100 == 2 | pm_vot_100 == 1) & pm_vot_101 == .

		gen highest_charged_am = 0 //if am_vot_* =="Charged"
		gen highest_charged_pm = 0 //if pm_vot_* =="Charged"
		gen lowest_free_am = 0 
		gen lowest_free_pm = 0

	* created a new variable which contains highest charged value
		foreach var in 100 90 80 70 60 50 40 30 25 20 15 10{
		  replace highest_charged_am = `var' if (am_vot_`var' == 2 & highest_charged_am == 0)
		  replace highest_charged_pm = `var' if (pm_vot_`var' == 2 & highest_charged_pm == 0)
		  replace lowest_free_am = `var' if am_vot_`var' == 0 
		  replace lowest_free_pm = `var' if pm_vot_`var' == 0
		}

		order lowest_free_am, after(highest_charged_am)
		order lowest_free_pm, after(highest_charged_pm)

		gen som_h_am=0 
		gen som_l_am=0
		gen som_h_pm=0 
		gen som_l_pm=0


	* created a new variable which contains highest and lowest "Sometimes" value
		foreach var in 100 90 80 70 60 50 40 30 25 20 15 10{
			replace som_h_am = `var' if (am_vot_`var' == 1 & som_h_am == 0)
			replace som_l_am = `var' if (am_vot_`var' == 1 & som_h_am != 0)
			replace som_h_pm = `var' if (pm_vot_`var' == 1 & som_h_pm == 0)
			replace som_l_pm = `var' if (pm_vot_`var' == 1 & som_h_pm != 0)
		}
	 
	* If VOT_100 is selected as sometimes, then replace it with the value they are ready to pay
		replace som_h_am = am_vot_101 if (am_vot_100 == 1 & am_vot_101 != . & am_vot_100 < am_vot_101)
		replace som_h_pm = pm_vot_101 if (pm_vot_100 == 1 & pm_vot_101 != . & pm_vot_100 < pm_vot_101)
	 
	* If "Charged" is higher than 100, then select that number
		gen highest_charged_vot101_am = 0
		replace highest_charged_vot101_am = am_vot_101 if (som_h_am ==0 & som_l_am==0)
	 
		gen highest_charged_vot101_pm = 0
		replace highest_charged_vot101_pm = pm_vot_101 if (som_h_pm ==0 & som_l_pm==0)
	 
	 * taking the average of "Sometimes"
		gen som_mean_am = 0
		replace som_mean_am = (som_h_am + som_l_am)/2
		gen som_mean_pm = 0
		replace som_mean_pm = (som_h_pm + som_l_pm)/2
		 
	 * taking the average of "charged and free"
		gen charge_free_mean_am = 0
		replace charge_free_mean_am = (highest_charged_am + lowest_free_am)/2 if (lowest_free_am!=0 & highest_charged_am!=0)
		gen charge_free_mean_pm = 0
		replace charge_free_mean_pm = (highest_charged_pm + lowest_free_pm)/2 if (lowest_free_pm!=0 & highest_charged_pm!=0)
		 
		gen VOT_am = 0
		gen VOT_pm = 0

	* VOT is average of "Sometimes" if it exists 
		replace VOT_am = som_mean_am if (som_h_am !=0 & som_l_am!=0)
		replace VOT_pm = som_mean_pm if (som_h_pm !=0 & som_l_pm!=0)
	 
	* Else, VOT is average of charged and free 
		replace VOT_am = charge_free_mean_am if (charge_free_mean_am!= . & VOT_am ==0)
		replace VOT_pm = charge_free_mean_pm if (charge_free_mean_pm!= . & VOT_pm ==0)

	* Else only "Charged"
		replace VOT_am = highest_charged_vot101_am if (VOT_am == 0 & highest_charged_vot101_am > VOT_am)
		replace VOT_pm = highest_charged_vot101_pm if (VOT_pm == 0 & highest_charged_vot101_pm > VOT_pm)
	 
	* VOT is 5 if they answer free for all options
		replace VOT_am = 5 if (am_vot_10 == 0 & am_vot_101 == .)
		replace VOT_pm = 5 if (pm_vot_10 == 0 & pm_vot_101 == .)

	* Removing uneccessary variables
		drop highest_charged* lowest_free* som_* charge_free*


	/* Removing some variables that are not used in the code */

	foreach var in 101 100 90 80 70 60 50 40 30 25 20 15 10 {
		drop am_vot_`var'
		drop pm_vot_`var'
	}

	drop   *_dt_noresp_check  *_label *_dt_intention_check
	drop dest_unique dest_flex
	order am_dt_p5, before(am_dt_p10)
	order pm_dt_p5, before(pm_dt_p10)


**************************
**** FURTHER CLEANING ****
**************************

/*** Removing and tagging absurd values ****/
						
/* am/pm_dep_time_0 value is absurd (too little like 0.5 , 0.333  but am/ pm_dep_time_1 seems proper between 7-11 am or 4-10 pm,
hence replacing dep_0 and dep_mean with dep_1 */

	#delimit ;
	replace pm_dep_time_mean = pm_dep_time_1 if inlist(key, "uuid:6a357e76-1461-4a5e-a653-64c148cd1459",
															"uuid:9f00bc3b-945d-40be-9e72-5fe609263fe0",
															"uuid:98357fd6-60f9-47f9-aefb-d1fb5a2d0f55");

															
	replace pm_dep_time_0 = pm_dep_time_1 if inlist(key, "uuid:6a357e76-1461-4a5e-a653-64c148cd1459", 
														 "uuid:9f00bc3b-945d-40be-9e72-5fe609263fe0",
														 "uuid:98357fd6-60f9-47f9-aefb-d1fb5a2d0f55");

	replace am_dep_time_mean = am_dep_time_1 if inlist(key, "uuid:ad0dcdc1-33bd-4046-94ca-fff395a5fa4a",
												            "uuid:7d62e346-37a7-4876-ad7d-d1414a87c0bc",
												            "uuid:f333f6a2-214a-4502-b7fe-2431a7593e63",
												            "uuid:98357fd6-60f9-47f9-aefb-d1fb5a2d0f55"); 

	replace am_dep_time_0 = am_dep_time_1 if inlist(key, "uuid:ad0dcdc1-33bd-4046-94ca-fff395a5fa4a",
												         "uuid:7d62e346-37a7-4876-ad7d-d1414a87c0bc",
												         "uuid:f333f6a2-214a-4502-b7fe-2431a7593e63",
												         "uuid:98357fd6-60f9-47f9-aefb-d1fb5a2d0f55");
	#delimit cr		

	/* Adding an indicator for absurd depature time for am below 4:00 and above 14:00 */
	gen am_absurd_dt =0
	replace am_absurd_dt =1 if (am_dep_time_0 < 4 | am_dep_time_0 > 14)
	replace am_absurd_dt =1 if ( am_dep_time_1 < 4 | am_dep_time_1 > 14)

	assert am_dep_time_0 <14 | am_dep_time_1 <14 if am_absurd_dt ==0
	assert am_dep_time_0 >4 | am_dep_time_1 >4 if am_absurd_dt ==0

	/* Adding an indicator for absurd depature time for am below 14:00 */
	gen pm_absurd_dt =0
	replace pm_absurd_dt =1 if (pm_dep_time_0 <14) 
	replace pm_absurd_dt =1 if (pm_dep_time_1 <14) 

	assert pm_dep_time_0 >14 | pm_dep_time_1 >14 if pm_absurd_dt ==0

	/* 95% confidence interval */

	sum am_dep_time_mean, d
	egen am_p5  = pctile(am_dep_time_mean), p(5) // finding 5% range for am
	egen am_p95 = pctile(am_dep_time_mean), p(95) // finding 95% range for am

	egen pm_p5  = pctile(pm_dep_time_mean), p(5) // finding 5% range for pm
	egen pm_p95 = pctile(pm_dep_time_mean), p(95) // finding 95% range for pm

	gen am_ci = am_dep_time_mean >am_p5 & am_dep_time_mean < am_p95 // applying 95% range for am
	gen pm_ci = pm_dep_time_mean >pm_p5 & pm_dep_time_mean < pm_p95 // applying 95% range for am

	/* removing values for beliefs_earlier and beliefs_later which are less than 5 minutes
	  some are negative values as well as some questions were asked as difference. It is challenging to ascertain which 
	  way the respondents answered. Hence, I have dropped all negative values. */
	  
	replace am_beliefs_dur_earlier 	= . if am_beliefs_dur_earlier 	< 5
	replace am_beliefs_dur_later 	= . if am_beliefs_dur_later 	< 5
	replace pm_beliefs_dur_earlier 	= . if pm_beliefs_dur_earlier 	< 5
	replace pm_beliefs_dur_later 	= . if pm_beliefs_dur_later 	< 5

	/* Finding the percentage difference from the mean duration in minutes*/
	gen am_early_test_perc = ((am_beliefs_dur_earlier 	- am_dur_mean)/am_dur_mean) * 100
	gen am_later_test_perc = ((am_beliefs_dur_later 	- am_dur_mean)/am_dur_mean) * 100
	gen pm_early_test_perc = ((pm_beliefs_dur_earlier 	- pm_dur_mean)/pm_dur_mean) * 100
	gen pm_later_test_perc = ((pm_beliefs_dur_later 	- pm_dur_mean)/pm_dur_mean) * 100

	/*Identifying the outliars with indicator which have change of less than -80% as lower bound. 
	-80% is randomly taken number (reasonably) as 1% cut off doesnt seem to apply for all the four below.
	eg.1% lower bound are am_early= -71%  am_later = -50% pm_early = - 73% pm_later = -60% */
	gen am_absurd_early_beliefs = am_early_test_perc < -80
	gen am_absurd_later_beliefs = am_later_test_perc < -80
	gen pm_absurd_early_beliefs = pm_early_test_perc < -80	
	gen pm_absurd_later_beliefs = pm_later_test_perc < -80

	/*Finding the difference of earlier/ later from the mean duration in minutes*/
	gen am_delta_early = (am_beliefs_dur_earlier - am_dur_mean)
	gen am_delta_later = (am_beliefs_dur_later - am_dur_mean)
	gen pm_delta_early = (pm_beliefs_dur_earlier - pm_dur_mean)
	gen pm_delta_later = (pm_beliefs_dur_later - pm_dur_mean)

	replace am_delta_early = . if am_absurd_early_beliefs ==1 // replacing absurd values with missing
	replace am_delta_later = . if am_absurd_later_beliefs ==1 
	replace pm_delta_early = . if pm_absurd_early_beliefs ==1 
	replace pm_delta_later = . if pm_absurd_later_beliefs ==1 

	/* Creating variable with beliefs_faster = 1 and beliefs_slower = 0 */
	/* Imp: beliefs_Sametime has been removed */

	/*
	gen am_early_faster_dum = am_delta_early <0 if am_delta_early !=. & am_delta_early !=0
	gen am_later_faster_dum = am_delta_later <0 if am_delta_later !=. & am_delta_later !=0
	gen pm_early_faster_dum = pm_delta_early <0 if pm_delta_early !=. & pm_delta_early !=0
	gen pm_later_faster_dum = pm_delta_later <0 if pm_delta_later !=. & pm_delta_later !=0
	*/

	/* Creating variable with beliefs_faster = 1 and beliefs_slower = 0 using the survey question*/
	/* Imp: beliefs_Sametime has been removed */
	tab am_beliefs_earlier
	gen am_early_faster_dum = am_beliefs_earlier 	< 0 if am_beliefs_earlier !=. & am_beliefs_earlier != 0
	gen am_later_faster_dum = am_beliefs_later 		< 0 if am_beliefs_earlier !=. & am_beliefs_earlier != 0
	gen pm_early_faster_dum = pm_beliefs_earlier 	< 0 if pm_beliefs_earlier !=. & pm_beliefs_earlier != 0
	gen pm_later_faster_dum = pm_beliefs_later 		< 0 if pm_beliefs_earlier !=. & pm_beliefs_earlier != 0


	/* Replacing DK/ NR values with missing values for DT question - would you leave earlier/sametime later/ sametime */
	replace am_dt_intention = . if am_dt_intention == 99
	replace pm_dt_intention = . if pm_dt_intention == 99

	/* Adding an indicator if its cheaper earlier and the respondent said he will leave later and vice versa */
	gen am_absurd_int=0
	replace am_absurd_int=1 if am_dt_cheaper == "earlier" & am_dt_intention ==1
	replace am_absurd_int=1 if am_dt_cheaper == "later"   & am_dt_intention == -1

	gen pm_absurd_int=0
	replace pm_absurd_int=1 if pm_dt_cheaper == "earlier" & pm_dt_intention ==1
	replace pm_absurd_int=1 if pm_dt_cheaper == "later"   & pm_dt_intention == -1

	/* Saving the cleaned file */
	// save "data/coded_cto/baseline/clean data_v1.dta", replace
	// save "data\analysis\analysis_v4_other\baseline_clean.dta"


**************************************
**** RESHAPE data to uidp x am/pm ****
**************************************
	/* Renaming variables to RESHAPE on am/ pm */
	rename (am_* pm_*) (*_am *_pm)

	#delimit ;
	reshape long  	dep_time_0  dep_time_1
				   	dep_time_mean dur_0 dur_1 dur_mean dur_c
	         		delta_t beliefs_earlier beliefs_later beliefs_c
	         		beliefs_dur_earlier beliefs_dur_later
					beliefs_ask_delta
					vot_base_extra vot_slower vot_dur vot_dur_slower vot_c vot_any_sometimes VOT
					dt_base_toll dt_c dt_cheaper dt_more_exp dt_p5 dt_p10 dt_p20 dt_intention 
					dt_delta dt_c_2 dt_delta_check dt_toll_5 has
					absurd_dt early_test_perc later_test_perc absurd_early_beliefs absurd_later_beliefs 
					delta_early delta_later early_faster_dum later_faster_dum ci absurd_int 
					, i(key) j(time1) string ;
	#delimit cr

	*** Stated behavior
	label var dep_time_0 "Earliest departure time"
	label var dep_time_1 "Latest departure time"
	label var dep_time_mean "Typical departure time"

	*** Stated beliefs
	label var dur_0 "Minimum trip duration (minutes)"
	label var dur_1 "Maximum trip duration (minutes)"
	label var dur_mean "Typical trip duration (minutes)"

	*** Conditions
	label var delta_t "Duration beliefs question: delta_t earlier/later"
	label var beliefs_ask_delta "Beliefs asked for _change_ (not levels). Only for form=2"

	*** Stated beliefs
	label var beliefs_earlier "Trip faster if leave earlier? (labeled)"
	label var beliefs_later   "Trip faster if leave later? (labeled)"
	label var early_faster_dum "Trip faster if leave earlier? (Dummy)"
	label var later_faster_dum "Trip faster if leave later? (labeled)"
	label var beliefs_dur_earlier "Trip duration if leave earlier (minutes)"
	label var beliefs_dur_later   "Trip duration if leave later (minutes)"
	label var delta_early "Earlier dep time duration minus typical duration"
	label var delta_later "Later dep time duration minus typical duration"

	*** Conditions
	label var vot_base_extra 	"VOT: extra minutes added to both options (0 or 20 minutes)"
	label var vot_slower 		"VOT: extra minutes on alternate route (10 or 15 minutes)"
	label var vot_dur 			"VOT: base route duration (dur_mean + vot_base_extra, min 30, max 120, rounded)"
	label var vot_dur_slower 	"VOT: alternate route duration (vot_dur + vot_slower)"
	label var vot_c 			"VOT: comments"
	label var vot_any_sometimes

	*** Stated behavior
	label var VOT "VOT: Indifference fare between (charged) base route and (free) alternate route"

	*** Conditions
	label var dt_base_toll 		"DT: charge for current departure time"
	label var dt_c 				"DT: comments"
	label var dt_cheaper 		"DT: cheaper earlier or later"
	label var dt_more_exp 	 	"DT: more exp earlier or later"
	label var dt_p5 			"DT: charge difference for 5 minutes"
	label var dt_p10 			"DT: charge difference for 10 minutes"
	label var dt_p20 			"DT: charge difference for 20 minutes"


	*** Stated behavior
	label var dt_intention 		"DT: stated departure time: earlier, same, later"
	label var dt_delta 			"DT: change in departure time (minutes)"
	label var dt_c_2 			"DT: comments"
	label var dt_delta_check 	"DT: verification question"
	label var dt_toll_5 		"DT: Toll for 5 minute change (aux)"

	*** Quality indicators
	label var absurd_dt 		"Departure Time outside usual range (4am-2pm for AM and 4pm-12am for PM)"
	label var early_test_perc 	"DT: early duration belief percentage difference from the mean duration in minutes"
	label var later_test_perc 	"DT: late  duration belief percentage difference from the mean duration in minutes"
	label var absurd_early_beliefs "Beliefs: huge change in travel time"
	label var absurd_later_beliefs "Beliefs: huge change in travel time"
	
	label var ci "Departure time within (5,95) of the distribution"
	label var absurd_int "DT: absurd reaction (travel earlier when it's more expensive, or vice versa)"

	gen     time = 1 if time1 == "_am"
	replace time = 2 if time1 == "_pm"
	assert inlist(time,1,2)

	label variable time "time label"
	label define time_label 1 "am" 2 "pm"
	label values  time time_label

*** Drop empty half-surveys
	assert (dep_time_0 == .) == (dep_time_1 == .) 
	assert (dep_time_0 == .) == (dur_0 == .)
	count if dep_time_0 == . & time == 1 
	count if dep_time_0 == . & time == 2
	drop if dep_time_0 == . /* 61 missing values in am and 67 missing values in pm */

	* ?
	assert has == 1
	drop has 

**************
*** SAMPLE ***
**************

*** Decide which surveys to keep (when the same person has multiple surveys)
	gen before_experiment = startdate <= datestart_experiment if in_exp == 1
	tab before_experiment // 40%

	/* 
	sample 1: earliest possible
	sample 2: earlist during late period
	sample 3: latest
	 */
	sort uidp_original time starttime
	by uidp_original time (starttime): gen baseline_sample_1 = _n==1
	by uidp_original time (starttime): gen baseline_sample_3 = _n==_N

	sort uidp_original time old_version starttime
	by uidp_original time old_version (starttime): gen baseline_sample_2 = _n==1 if old_version == 0
	recode baseline_sample_2 (.=0)

	// br uidp_original time old_version startdate baseline_sample_?

	label var baseline_sample_1 "Earliest survey (early wave included)"
	label var baseline_sample_2 "Earliest survey in late wave"
	label var baseline_sample_3 "Most recent survey"

	// sort uidp_original time starttime
	// by uidp_original time: gen baseline_first = _n==1
	// by uidp_original time: gen baseline_last  = _n==_N
	// tab baseline_first baseline_last

	// keep if old_version == 0
	// sort uidp_original time starttime
	// by uidp_original time: gen baseline_first = _n==1
	// by uidp_original time: gen baseline_last  = _n==_N
	// tab baseline_first baseline_last
	// tab before_experiment 
	// tab before_experiment if baseline_first == 1
	// tab before_experiment if baseline_last == 1

	/* Load dta with multiple observation per UIDP with tag to drop */
	// merge m:1 key time1 using "data/temp/baseline/oneuidp_onesurvey.dta" // merge on key time1 level
	// assert _m!=2
	// drop if _m == 1 
	// tab _m n_surveys, m
	// drop _m
	
	* check unique ID
	preserve
		keep if baseline_sample_1 == 1
		isid uidp_original time 
	restore
	preserve
		keep if baseline_sample_1 == 2
		isid uidp_original time 
	restore
	preserve
		keep if baseline_sample_1 == 3
		isid uidp_original time 
	restore

	mdesc

*** Save 
	saveold "data/coded_cto/baseline_coded.dta", replace 
