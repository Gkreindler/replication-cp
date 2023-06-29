
 * This do file: Codes chain trips at daily level (v0 - no HW information, no charges)
 * This version: 30 July, 2017
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

if "`1'" == ""{
	local chain_th "15"	
}
else{
	local chain_th "`1'"		
}

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"
	adopath ++ "code/ado/"


*** Visit H/W/W2
	use "data/coded_gps_dta/coded_trips_`chain_th'.dta", clear
	gen at_h  = oh==1  | dh==1		if !missing(oh, dh)
	gen at_w  = ow==1  | dw==1		if !missing(ow, dw)
	gen at_w2 = ow2==1 | dw2==1	if !missing(ow2, dw2)

	collapse (max) at_h at_w at_w2, by(uidp date)

	sum at_*

	tempfile visit_places
	save 	`visit_places'


********************
** Load Trip Data **
********************
	use "data/coded_gps_dta/coded_trips_`chain_th'.dta", clear

**************************
** DAILY OUTCOME CODING **
**************************
	isid uidp date chain
	sort uidp date chain
	
	** SUM OVER TRIPS	
	gen one = 1
	gen sample_trip_all = 1
	assert inlist(sample_trip_ok,0,1)
	assert am_sample * pm_sample * am_pk_sample * pm_pk_sample * am_pk1_sample * pm_pk1_sample != .

	forv i = 1/9{
		di "`i'"
		local prefix: 		word `i' of "nt_" 		"km_" 	"du_" 			"ra_"		"ch_"			"sm1_"		"s2_"		"sp1_"		"cha_"			
		local var_to_sum: 	word `i' of "any_trip" 	"plb" 	"dur_trips" 	"c_dt_rate"	"c_dt_charge" 	"c_dt_sm1"	"c_dt_s2"	"c_dt_sp1"	"c_area_charge"

		* no missing values!
		assert (`var_to_sum' != .) | ("`var_to_sum'" == "c_area_charge" & (c_area_charge != . | in_area_treatment==0))

		forv j=1/2{
			local prefix2: word `j' of "all" "ok"
			local trip_sample "sample_trip_`prefix2'"

			* no missing values!
			assert `trip_sample' != .

			forv k=1/3{
				di "i`i' j`j' k`k'"
				* the endpoints sample of the trip: all, home-work/2, and home-work/2 and work/2-work/2
				local prefix3: 	 word `k' of "" 	"hwp_" 		"hwpwp_"
				local trip_dest: word `k' of "one" 	"trip_hwp" 	"trip_hwpwp" 

				* no missing values!
				assert `trip_dest' != .

				gegen `prefix'`prefix2'_`prefix3' 	  = sum(`var_to_sum' * `trip_sample' * `trip_dest'), by(uidp date)

				* AM/PM samples (AM until 1pm, PM from 4pm)
				gegen `prefix'`prefix2'_`prefix3'am  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_sample), by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'pm  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_sample), by(uidp date)

				* full AM/PM samples (until/starting at 2pm)
				gegen `prefix'`prefix2'_`prefix3'aam  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_a_sample), by(uidp date) // area am/pm samples, by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'apm  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_a_sample), by(uidp date) // area am/pm samples, by(uidp date)

				* (narrow) peak samples: +/- 1.5h around the peak
				gegen `prefix'`prefix2'_`prefix3'cam = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_pk_sample * am_sample), by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'cpm = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_pk_sample * pm_sample), by(uidp date)

				* (wide) peak sample: +/- 2.5h around the peak
				gegen `prefix'`prefix2'_`prefix3'c1am = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_pk1_sample * am_sample), by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'c1pm = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_pk1_sample * pm_sample), by(uidp date)

				* (wide) PRE peak sample: -2.5h to 0 before the peak
				gegen `prefix'`prefix2'_`prefix3'pream = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_pk1_sample * am_pre_peak), by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'prepm = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_pk1_sample * pm_pre_peak), by(uidp date)

				* (wide) POST peak sample: 0 to 2.5h after the peak
				gegen `prefix'`prefix2'_`prefix3'postam = sum(`var_to_sum' * `trip_sample' * `trip_dest' * am_pk1_sample * (1 - am_pre_peak)), by(uidp date)
				gegen `prefix'`prefix2'_`prefix3'postpm = sum(`var_to_sum' * `trip_sample' * `trip_dest' * pm_pk1_sample * (1 - pm_pre_peak)), by(uidp date)

				if `i' == 1{
					* AM ramps
					gegen `prefix'`prefix2'_`prefix3'ramp_am_pre  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * ramp_am_pre), by(uidp date)
					gegen `prefix'`prefix2'_`prefix3'ramp_am_post = sum(`var_to_sum' * `trip_sample' * `trip_dest' * ramp_am_post), by(uidp date)

					gegen `prefix'`prefix2'_`prefix3'ramp_pm_pre  = sum(`var_to_sum' * `trip_sample' * `trip_dest' * ramp_pm_pre), by(uidp date)
					gegen `prefix'`prefix2'_`prefix3'ramp_pm_post = sum(`var_to_sum' * `trip_sample' * `trip_dest' * ramp_pm_post), by(uidp date)
				}
			}	
		}
	}
	drop one

* create new variables = at least one trip
	forv j=1/2{
		local prefix2: word `j' of "all" "ok"
		local trip_sample "sample_trip_`prefix2'"	

		forv k=1/3{
			di "j`j' k`k'"
			local prefix3: word `k' of "" "hwp_" "hwpwp_"
			// local trip_dest: word `k' of "one" 	"trip_hwp" 	"trip_hwpwp"

			assert nt_`prefix2'_`prefix3'   != .
			assert nt_`prefix2'_`prefix3'am != .
			assert nt_`prefix2'_`prefix3'pm != .
			assert nt_`prefix2'_`prefix3'aam != .
			assert nt_`prefix2'_`prefix3'apm != .
			assert nt_`prefix2'_`prefix3'cam != .
			assert nt_`prefix2'_`prefix3'cpm != .
			assert nt_`prefix2'_`prefix3'c1am != .
			assert nt_`prefix2'_`prefix3'c1pm != .

			assert nt_`prefix2'_`prefix3'pream != .
			assert nt_`prefix2'_`prefix3'prepm != .
			assert nt_`prefix2'_`prefix3'postam != .
			assert nt_`prefix2'_`prefix3'postpm != .

			gen nt0_`prefix2'_`prefix3' 	= nt_`prefix2'_`prefix3' > 0
			gen nt0_`prefix2'_`prefix3'am 	= nt_`prefix2'_`prefix3'am > 0
			gen nt0_`prefix2'_`prefix3'pm 	= nt_`prefix2'_`prefix3'pm > 0
			gen nt0_`prefix2'_`prefix3'aam 	= nt_`prefix2'_`prefix3'aam > 0 // area am/pm samples
			gen nt0_`prefix2'_`prefix3'apm 	= nt_`prefix2'_`prefix3'apm > 0
			gen nt0_`prefix2'_`prefix3'cam 	= nt_`prefix2'_`prefix3'cam > 0
			gen nt0_`prefix2'_`prefix3'cpm 	= nt_`prefix2'_`prefix3'cpm > 0
			gen nt0_`prefix2'_`prefix3'c1am = nt_`prefix2'_`prefix3'c1am > 0
			gen nt0_`prefix2'_`prefix3'c1pm = nt_`prefix2'_`prefix3'c1pm > 0

			gen nt0_`prefix2'_`prefix3'pream = nt_`prefix2'_`prefix3'pream > 0
			gen nt0_`prefix2'_`prefix3'prepm = nt_`prefix2'_`prefix3'prepm > 0
			gen nt0_`prefix2'_`prefix3'postam = nt_`prefix2'_`prefix3'postam > 0
			gen nt0_`prefix2'_`prefix3'postpm = nt_`prefix2'_`prefix3'postpm > 0

			assert cha_`prefix2'_`prefix3'   != .
			assert cha_`prefix2'_`prefix3'am != .
			assert cha_`prefix2'_`prefix3'pm != .
			assert cha_`prefix2'_`prefix3'aam != .
			assert cha_`prefix2'_`prefix3'apm != .
			assert cha_`prefix2'_`prefix3'cam != .
			assert cha_`prefix2'_`prefix3'cpm != .
			assert cha_`prefix2'_`prefix3'c1am != .
			assert cha_`prefix2'_`prefix3'c1pm != .

			assert cha_`prefix2'_`prefix3'pream != .
			assert cha_`prefix2'_`prefix3'prepm != .
			assert cha_`prefix2'_`prefix3'postam != .
			assert cha_`prefix2'_`prefix3'postpm != .

			gen cha0_`prefix2'_`prefix3' 	= cha_`prefix2'_`prefix3' > 0
			gen cha0_`prefix2'_`prefix3'am 	= cha_`prefix2'_`prefix3'am > 0
			gen cha0_`prefix2'_`prefix3'pm 	= cha_`prefix2'_`prefix3'pm > 0
			gen cha0_`prefix2'_`prefix3'aam = cha_`prefix2'_`prefix3'aam > 0 // area am/pm samples
			gen cha0_`prefix2'_`prefix3'apm = cha_`prefix2'_`prefix3'apm > 0
			gen cha0_`prefix2'_`prefix3'cam = cha_`prefix2'_`prefix3'cam > 0
			gen cha0_`prefix2'_`prefix3'cpm = cha_`prefix2'_`prefix3'cpm > 0
			gen cha0_`prefix2'_`prefix3'c1am = cha_`prefix2'_`prefix3'c1am > 0
			gen cha0_`prefix2'_`prefix3'c1pm = cha_`prefix2'_`prefix3'c1pm > 0

			gen cha0_`prefix2'_`prefix3'pream = cha_`prefix2'_`prefix3'pream > 0
			gen cha0_`prefix2'_`prefix3'prepm = cha_`prefix2'_`prefix3'prepm > 0
			gen cha0_`prefix2'_`prefix3'postam = cha_`prefix2'_`prefix3'postam > 0
			gen cha0_`prefix2'_`prefix3'postpm = cha_`prefix2'_`prefix3'postpm > 0

		}
	}

*** keep 
	keep 	uidp date dow ///
			qual_* data_quality report_outstation outstation_qual jump_dist tot_gap_hr ///
			nt0_* nt_* km_* du_* ra_* ch_* sm1_* s2_* sp1_* cha_* cha0_* ///
			date_exp_min date_exp_max cycle_type study_cycle study_day time_in_exp calendar_time date_wo_trips ///
			report_bal_current report_charge date_meeting in_area_treatment dropped

	duplicates drop 
	gisid date uidp

*** add 
	merge 1:1 uidp date using `visit_places'
	assert _m==3
	drop _m

****************************************************
**** Adding calculated charges at the DAY level **** - DT
****************************************************
	merge 1:1 uidp date using "data/coded_gps_dta/charges_days_dt_`chain_th'.dta"

	drop if date == mdy(4,26,2017) // already dropped from coded trips

	* these are all probably days with essentially no data (maybe one point) - drop
	count if _m==2
	assert r(N) == 127
	drop if _m==2  	

	* most non-matched from master are dropped (no subsequent data)
	count if dropped==1 & _m==1 
	assert r(N) <= 16080  // Master includes many dates after the last date with data (dropped == 1)

	count if dropped==0 & _m==1
	// assert r(N) == 106	
	assert r(N) == 0

	// list date data_quality dropped _m if uidp == "L0325191810"

	* all days with trips have charges - of course
	assert _m==3 if date_wo_trips == 0
	drop _m

****************************************************
**** Adding calculated charges at the DAY level **** - AREA
****************************************************

	merge 1:1 uidp date using "data/coded_gps_dta/charges_days_area_`chain_th'.dta"

	drop if date == mdy(4,26,2017) // already dropped from coded trips
	
	* these are all probably days with essentially no data (maybe one point) - drop
	count if _m==2
	assert r(N) == 81
	drop if _m==2  	

	* most non-matched from master are dropped (no subsequent data)
	count if dropped==1 & _m==1  & in_area_treatment == 1
	// assert r(N) == 457
	assert r(N) == 8012  // Master includes many dates after the last date with data (dropped == 1)

	count if dropped==0 & _m==1 & in_area_treatment == 1
	// assert r(N) == 83	
	assert r(N) == 0


	* all days with trips have charges - of course
	assert _m==3 if date_wo_trips == 0 & in_area_treatment==1
	drop _m

*** Checks
	mdesc

	* all  missing values for DT charges are from respondents who dropped out.
	* we don't need to replace these because anyway we'll restrict to available data
	count if cd_charge_am==. & dropped== 1
	assert r(N) <= 16080
	count if cd_charge_am==. & dropped== 0
	assert r(N) == 0

	mdesc cad_* if in_area_treatment == 1

	count if cad_charge_am==. & dropped==1 & in_area_treatment==1
	assert r(N) == 8012
	count if cad_charge_am==. & dropped==0 & in_area_treatment==1
	assert r(N) == 0

********************************
**** treatment and study dates
********************************

	merge m:1 uidp using "data/coded_gps_dta/treat_assign.dta"
	assert _m==3
	drop _m

****************
**** SAMPLE ****
****************

	* sample before experiment and during (drop spectial, outstation, and after)
	gen sample_analysis = inlist(study_cycle,0,1,2,3,4,5) & qual_gb == 1

	* sample except after experiment
	gen sample_attrition = inlist(study_cycle,0,1,2,3,4,5,99)

	* sample exclusing trial/special
	gen sample_study_proper = inlist(study_cycle,0,1,2,3,4,5)

*** TREATMENT VARs CODING
	gen dt_wcha = inlist(dt_treat,1,2)
	label var dt_wcha "Any DT treatment with charges (Low or High)"

	gen dt_noch = inlist(dt_treat,3,4)
	label var dt_noch "Any DT treatment without charges (Control or Info)"

	gen post = inrange(time_in_exp,1,25) //& study_cycle != 99
	assert post == inlist(study_cycle,1,2,3,4,5,99)  	// check that the definitions are equivalent
	assert post == study_cycle != 0 if sample_analysis == 1  // check that the definitions are equivalent

	gen dt_wcha_post 	= dt_wcha 	* post
	gen dt_noch_post 	= dt_noch 	* post

	gen dt_ctrl_post 	= dt_ctrl 	* post
	gen dt_info_post 	= dt_info 	* post
	gen dt_low_post 	= dt_low 	* post
	gen dt_high_post 	= dt_high 	* post

	*** AREA treatment
	assert a_treat!=.
	assert a_treat != 0 if a_late == 1
	assert a_treat != 0 if a_early == 1

	** Attention: remember to include study_cycle FE or restrict to inlist(study_cycle,1,4), otherwise not the correct control group
	gen 	area_treat = (a_late ==1 & study_cycle==4) | (a_early==1 & study_cycle==1)
	
	gen 	area_treat_early 	= area_treat * a_early
	gen 	area_treat_late  	= area_treat * a_late

	gen 	area_treat_high 	= area_treat * a_high
	gen 	area_treat_low 		= area_treat * a_low

	gen 	area_treat_long 	= area_treat * a_long
	gen 	area_treat_short 	= area_treat * a_short

	* 
	gen 	area_late_high = a_high * a_late
	gen 	area_late_low  = a_low  * a_late
	assert 	area_late_high + area_late_low == a_late

	* 
	gen 	area_late_short = a_short * a_late
	gen 	area_late_long  = a_long  * a_late
	assert 	area_late_short + area_late_long == a_late

*** TREATMENT SAMPLES
*** Define DT samples
	* full pre period
	* during experiment: 3 weeks based on DT late timing, including in control group
	gen 	sample_3w = 1 - post
	replace sample_3w = 1 if post == 1 & ((dt_late == 1 & inlist(study_cycle,2,3,4)) | (dt_late == 0 & inlist(study_cycle,1,2,3)))

	* full pre period
	* during experiment: 3 weeks in treatment (based on DT late timing), and full weeks 1-4 in control group
	gen 	sample_3w_fullc = 1 - post
	replace sample_3w_fullc = 1 if post == 1 & ((dt_late == 1 & inlist(study_cycle,2,3,4)) | (dt_late == 0 & inlist(study_cycle,1,2,3)))
	replace sample_3w_fullc = 1 if post == 1 & (dt_noch == 1 & inlist(study_cycle,1,2,3,4))
	// tab sample_3w sample_3w_fullc if post == 1, m

*** Define AREA samples
	* only for AREA treatment (a_treat != 0)
	gen sample_w14  = a_treat != 0 & inlist(study_cycle,1,4)
	gen sample_w1  	= a_treat != 0 & inlist(study_cycle,0,1)
	gen sample_w4  	= a_treat != 0 & inlist(study_cycle,0,4)
	gen sample_wall = a_treat != 0 & inlist(study_cycle,0,1,2,3,4)

********************
*** EXTRA CODING ***
********************

*** HETEROGENEITY ***
*** High versus low travellers
	* average daily pre-period KM travelled 
	* during the AM peak interval +/- 1 hour
	cap drop temp
	bys uidp: gegen temp = mean(km_ok_c1am) if qual_gb == 1 & study_cycle == 0
	bys uidp: gegen mean_km_ok_c1am = mean(temp)
	drop temp 

	* during the PM peak interval +/- 1 hour
	bys uidp: gegen temp = mean(km_ok_c1pm) if qual_gb == 1 & study_cycle == 0
	bys uidp: gegen mean_km_ok_c1pm = mean(temp)
	drop temp 

	* average daily pre-period KM travelled (total)
	gen mean_pre_km_ampm = mean_km_ok_c1am + mean_km_ok_c1pm

	* define dummy
	bys uidp:  gen o1=_n==1
	sum mean_pre_km_ampm if o1==1, d
	gen mprekm_high = mean_pre_km_ampm > r(p50)
	tab mprekm_high if o1 == 1

	* other testing	
	// scatter mean_km_ok_c1am mean_km_ok_c1pm if o1==1
	// reg     mean_km_ok_c1am mean_km_ok_c1pm if o1==1
	// reg mean_pre_ampm i.strat_cell_new, cl(uidp)
	// reg mean_pre_ampm strat3, cl(uidp)
	// tab mean_pre_ampm_low strat3 if o1==1, m // good overlap!
	// twoway (kdensity mean_pre_ampm if strat3==1 & o1==1) (kdensity mean_pre_ampm if strat3==0 & o1==1)

	* looks good - low and high are separate
		// twoway 	(kdensity mean_pre_ampm if strat3==1 & mean_pre_ampm_low==1 & o1==1) ///
		// 	(kdensity mean_pre_ampm if strat3==0 & mean_pre_ampm_low==0 & o1==1) ///
		// 	(kdensity mean_pre_ampm if strat3 != mean_pre_ampm_low & o1==1) ///
		// 	,legend(order(1 "low" 2 "high" 3 "disagree"))

*** HETEROGENEITY ***
*** Meeting after median
	sum date_meeting if o1==1,d
	di %td r(p50)
	gen recruit_late = date_meeting >= r(p50)
	tab recruit_late if o1==1


*** Format heterogeneity variables ***
	gen car_num = strat2
	gen mprekm_high_strat = 1 - strat3

	forv i=1/4{
		local hva: word `i' of car_num 	mprekm_high_strat 	mprekm_high 	recruit_late
		local suf: word `i' of "car" 	"km2" 				"km"			"rl"

		assert inlist(`hva',0,1)

		foreach va of varlist dt_info dt_low dt_high post study_cycle dt_info_post dt_low_post dt_high_post{
			gen `va'_`suf' = `va' * `hva'
		}
	}
	* s2 -> 1 = motor  		(originally 1 = car)
	* s3 -> 1 = high km  	(originally 1 = low)
	* km -> 1 = high km  	(originally 1 = low)
	* rl -> 1 = recruit late(originally 1 = early)




*** Define at any workplace:
	gen		at_w_all = at_w
	replace at_w_all = max(at_w, at_w2) if at_w2 != .
	tab at_w_all at_w, m

* check that all days are present
	sort uidp date
	by uidp: gegen date_max = max(date)
	by uidp: gegen date_min = min(date)
	by uidp: gegen n_dates = count(date)
	weekdays_bw date_min date_max, generate(n_weekdays_)
	assert inlist(n_weekdays_, n_dates, n_dates + 1)
	drop date_max date_min n_dates n_weekdays_

*** Final SAMPLE selection
	drop if date == mdy(4,26,2017)
	
*** Save 
	save "data/coded_gps_dta/coded_days_`chain_th'.dta", replace
