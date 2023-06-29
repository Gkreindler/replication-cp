
 * This do file: Code number of trips at day level, includes several variables counting the (smoothed) number of trips in bins relative to peak profile 
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


********************
** Load Trip Data **
********************
	use "data/coded_gps_dta/coded_trips_`chain_th'.dta", clear

	//10 //15//20//30
	local bw_minutes = 20


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

	* radius of Epanechnikov /\  <- 2xbw
	scalar bw = `bw_minutes' / 60
	scalar nbins = 36
	scalar bin_scale_factor = 3/nbins

	// every 5 minutes, +/- 3 hours around the peak
	forv tt = `=-nbins'/`=nbins'{

		local tt_name = `tt' + nbins

		local var_to_sum "any_trip"
		local prefix "ntr`tt_name'_"

	* Epanechnikov
		cap drop center
		gen center = `tt'*bin_scale_factor 
		cap drop kernel 
	 	gen kernel = max(0, 3/4 * (1 - ((t_rel - center) / bw)^2 ) ) / bw * 100
	 	// scatter kernel t_rel if t_rel < -2

	/* * trapezoidal	
		scalar t1 = `tt'*bin_scale_factor - 2 * ss
		scalar t2 = `tt'*bin_scale_factor -     ss
		scalar t3 = `tt'*bin_scale_factor
		scalar t4 = `tt'*bin_scale_factor +     ss
		scalar t5 = `tt'*bin_scale_factor + 2 * ss

		* kernel density   /--\
		cap drop kernel
		gen kernel = max(0, min(1, (t_rel - t1) / ss, (t5 - t_rel) / ss))
		replace kernel = . if t_rel == .
		// scatter kernel t_rel if t_rel < -3

		assert kernel == 1 if inrange(t_rel, t2, t4)
		assert inrange(kernel, 0, 1) if inrange(t_rel, t1, t2)
		assert inrange(kernel, 0, 1) if inrange(t_rel, t4, t5)
		assert inlist(kernel,0,.) if !inrange(t_rel,t1 - 0.001, t5 + 0.001) 
	// */

		* total sum:
		// cap drop kernel_sum
		// cap drop kernel_sum_am
		// cap drop kernel_sum_pm
		// egen kernel_sum = sum(kernel)
		// egen kernel_sum_am = sum(kernel * am_sample)
		// egen kernel_sum_pm = sum(kernel * pm_sample)

		* no missing values!
		assert (`var_to_sum' != .) | ("`var_to_sum'" == "c_area_charge" & (c_area_charge != . | in_area_treatment==0))

		forv j=2/2{
			local prefix2: word `j' of "all" "ok"
			local trip_sample "sample_trip_`prefix2'"

			* no missing values!
			assert `trip_sample' != .

			forv k=1/2{
				di ">>>>>>>>>>>>>>>>>>>>>>>>>   tt = `tt' j = `j' k = `k'"
				* the endpoints sample of the trip: all, home-work/2, and home-work/2 and work/2-work/2
				local prefix3: 	 word `k' of "" 	"hw_" 		"hwpwp_"
				local trip_dest: word `k' of "one" 	"trip_hwp" 	"trip_hwpwp" 

				* no missing values!
				assert `trip_dest' != .

				cap drop yvar
				gen yvar = `var_to_sum' * `trip_sample' * `trip_dest'
				by uidp date: egen `prefix'`prefix2'_`prefix3' = sum(yvar * kernel) 
				// replace 		   `prefix'`prefix2'_`prefix3' = `prefix'`prefix2'_`prefix3' / kernel_sum

				by uidp date: egen `prefix'`prefix2'_`prefix3'am = sum(yvar * kernel * am_sample) 
				// replace 		   `prefix'`prefix2'_`prefix3'am = `prefix'`prefix2'_`prefix3'am / kernel_sum_am

				by uidp date: egen `prefix'`prefix2'_`prefix3'pm = sum(yvar * kernel * pm_sample) 
				// replace 		   `prefix'`prefix2'_`prefix3'pm = `prefix'`prefix2'_`prefix3'pm / kernel_sum_pm

			} // k	
		} // j
	} // t_rel
	drop one

*** keep 
	keep 	uidp date dow ///
			qual_* data_quality report_outstation outstation_qual jump_dist tot_gap_hr ///
			ntr* ///
			date_exp_min date_exp_max cycle_type study_cycle study_day time_in_exp calendar_time date_wo_trips ///
			report_bal_current report_charge date_meeting in_area_treatment dropped

	duplicates drop 
	isid date uidp

********************************
**** treatment and study dates
********************************
	// local ver "v3"
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
* check that all days are present
	sort uidp date
	by uidp: egen date_max = max(date)
	by uidp: egen date_min = min(date)
	by uidp: egen n_dates = count(date)
	weekdays_bw date_min date_max, generate(n_weekdays_)
	assert inlist(n_weekdays_, n_dates, n_dates + 1)
	drop date_max date_min n_dates n_weekdays_

*** Final SAMPLE selection
	drop if date == mdy(4,26,2017)
	
*** Save 
	save "data/coded_gps_dta/coded_days_peak_`chain_th'_`bw_minutes'.dta", replace

