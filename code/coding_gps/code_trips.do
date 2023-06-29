
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

*** get peak times and date of meeting
	import delimited using "data/treatment/treatment roster noPII.csv", clear
	keep if meeting == "done"
	keep uidp dtamstart	dtpmstart datemeeting

	gen temp = clock(dtamstart, "hm")
	format %tc_hh:MM temp
	gen pk_am = hh(temp) + mm(temp) / 60 + 1.5 // this is the midpoint of the CC ramp

	gen temp2 = clock(dtpmstart, "hm")
	format %tc_hh:MM temp2
	gen pk_pm = hh(temp2) + mm(temp2) / 60 + 1.5 // this is the midpoint of the CC ramp

	* convert to date
	gen date_meeting = date("2017-" + substr(datemeeting,2,10), "YMD")
	format %td date_meeting

	keep uidp pk_am pk_pm date_meeting
	tempfile peak_times
	save 	`peak_times'

*** get area
	import delimited using "data/treatment/treatment roster noPII.csv", clear
	keep if meeting == "done"
	rename aradius area_radius
	assert (area_radius == . & atreat == "0 No treatment") | (area_radius != . & atreat != "0 No treatment")
	gen in_area_treatment = area_radius !=.
	keep uidp area_radius in_area_treatment
	tempfile area_radius
	save 	`area_radius'

******************
*** LOAD TRIPS ***
******************

	use "data/coded_gps_dta/chains_exp_`chain_th'.dta", clear

*** Merge DT-trip charges computed separately
	merge 1:1 uidp date chain using "data/coded_gps_dta/charges_trips_dt_`chain_th'.dta"
	assert _m!=2
	assert inrange(t_start, 0, 5) if _m==1 // only night time chains (charges not computed)
	drop _m

	* checks
	gen plb_rounded = floor(plb/100) / 10
	gen diff_plb = abs(plb_rounded - c_dt_distance)
	count if !(diff_plb == 0 | diff_plb == .)
	assert diff_plb < 0.1 if !(diff_plb == 0 | diff_plb == .)
	assert r(N) <= 1
	drop diff_plb c_dt_distance plb_rounded
	
	* checks
	gen dur_bg_mm_rounded = round(dur_bg_mm)
	gen diff_dur = abs(dur_bg_mm_rounded - c_dt_duration)
	assert inlist(diff_dur,0,1,.)
	drop diff_dur c_dt_duration dur_bg_mm_rounded

	* checks
	gen tstart_rounded = floor(t_start_bg_ / 60) / 60
	// br t_start* c_t* tstart_rounded
	gen diff_t = abs(tstart_rounded - c_dt_t_start)
	assert diff_t < 0.1 | diff_t == .
	drop diff_t tstart_rounded c_dt_t_start

	* generate ramp type dummies
	assert inlist(c_dt_ramp,.,-1,2,1,0)
	gen c_dt_sm1 = c_dt_rate * 100 / 24 * (c_dt_ramp == -1)  // minus 1 --- ascending ramp
	gen c_dt_s2  = c_dt_rate * 100 / 24 * (c_dt_ramp ==  2)  //		 --- peak
	gen c_dt_sp1 = c_dt_rate * 100 / 24 * (c_dt_ramp ==  1)  // plus 1	 --- descending ramp


*** add info on atreat and area radius
	merge m:1 uidp using `area_radius'
	assert _m==3
	drop _m

*** Merge AREA charges computed separately
	merge 1:1 uidp date chain using "data/coded_gps_dta/charges_trips_area_`chain_th'.dta"
	assert _m!=2
	assert inrange(t_start,0,5) | in_area_treatment==0 if _m==1
	drop _m
	
	// c_area_charge c_area_intersect c_area_min_dist c_area_min_edge c_area_time_intersect

* tag borderline crossings
	gen temp = c_area_min_dist / area_radius
	gen c_area_borderline = inrange(temp, 0.85, 1.15) if in_area_treatment == 1
	drop temp
	gen c_low_quality = c_area_min_edge > 1000
	drop in_area_treatment area_radius

	// lpoly c_area_borderline c_area_min_dist if c_area_min_dist < 2000, nosc ci
	// lpoly c_area_intersect c_area_min_dist if c_area_min_dist < 2000, nosc ci


************************
******* SAMPLE *********
************************
	* drop chains from outstation (based on trip data), weekends and nights (10pm-5am)
	keep if sample_drop_encoded > 3
	// drop if inlist(sample_drop_encoded,1,2,3) // drop weekends and night

	*** TESTING
	// tab date if uidp == "L0525103515_DA90"
	// sdfsd

************************

* merge with data quality to check number of days (good+bad) with missing trips (or with outstation)
	merge m:1 uidp date using "data/coded_gps_dta/quality_exp.dta" , ///
	keepusing(qual_* data_quality ///
			  date_exp_min date_exp_max cycle_type study_cycle study_day time_in_exp calendar_time ///
			  report_bal_current report_charge report_outstation outstation_qual jump_dist tot_gap_hr dropped)
	// dropped nodata_spell max_nodata_spell report_outstation 

	assert _m!=1 // all uidp - days appear in quality
	assert date > mdy(4,3,2017) if uidp == "L0314124709"
	assert date > mdy(4,3,2017) if uidp == "L0221092234" 
	* check no trip and outstation days for a particular respondent - looks ok!
	// tab date _m if uidp == "L0208100638"

	* if _m==2, either no data on that day (70%), or good/bad data but no trips
	gen date_wo_trips = _m==2 
	drop _m

	* NOTE: we have DROPPED outstation chains, so the outstation_qual = 1 (N=320) that are merged are discrepancies between chains and qual
	tab date_wo_trips data_quality, m
	tab date_wo_trips data_quality if sample_trip_ok==1, m
	gunique uidp date if date_wo_trips == 0 & qual_gb == 1 & outstation_qual == 0
	* --> 
	count 			 if date_wo_trips == 1 & qual_gb == 1 & outstation_qual == 0
	* --> 
	* ratio = 597 / (597 + 21154) = .02744701  <- includes night trips
	* ratio = 832 / (832 + 20919) = .03825111  < - excludes night trips (10pm-7am)
	* small number of - stay at home days

	* no days with trips and no_data, EXCEPT April 26
	count if qual_no_data == 1 & chain !=. & date!=mdy(4,26,2017)
	assert r(N) == 0

*** CODING
	* fix chain ID on days without trips
	assert chain == . if date_wo_trips == 1
	replace chain = 1 if date_wo_trips == 1

	* INCLUDE days without trips in the sample of good trips
	assert  sample_trip_ok == . if date_wo_trips == 1
	replace sample_trip_ok = 1 if date_wo_trips == 1
	assert  sample_trip_ok != .

*** Clean Vars
	replace plb = plb / 1000

	foreach va of varlist plb c_dt_charge c_dt_ramp dur_trips dur_bg_mm{
		rename `va' `va'_old
		winsor `va'_old, gen(`va') highonly p(0.01)
	}

*** Day-level coding
	isid uidp date chain

	* minor cleaning
	replace dow = dow(date) if dow==.
 	order dow, after(date)

 	* re-merge area treat 
 	merge m:1 uidp using `area_radius', keepusing(in_area_treatment)
	assert _m==3
	drop _m

	** Replace KM, duration, DT charges = 0 on days with >= halfday data but no trips (also for all locations)
	foreach va of varlist plb dur_bg_mm dur_trips c_dt_charge c_dt_rate c_dt_ramp trip_hw trip_hw2 trip_ww2 trip_hwp trip_hwpwp{
		assert `va' == . if date_wo_trips == 1
		replace `va' = 0 if date_wo_trips == 1
		assert `va' != .
	}


	** Replace area charges = 0 on days with >= halfday data but no trips (also for all locations)
	foreach va of varlist c_area_charge c_area_intersect{
		assert `va' == . if date_wo_trips == 1 & in_area_treatment == 1
		replace `va' = 0 if date_wo_trips == 1 & in_area_treatment == 1
		assert `va' != . if in_area_treatment == 1
	}

	** same for rate, ramp and the other indicators to zero on days without trips
	foreach va of varlist c_dt_sm1 c_dt_s2 c_dt_sp1{
		assert `va' == . if date_wo_trips == 1
		replace `va' = 0 if date_wo_trips == 1
		assert inrange(`va',0,100)
	}

	recode trip trips_dur_jump trips_plb_jump (.=0)

	* On days without trips, any_trip = 0. =1 otherwise
	gen any_trip = 1 - date_wo_trips
	assert inlist(any_trip,0,1)

	* normalize the DT rate to 0 to 100
	replace c_dt_rate = c_dt_rate / 24 * 100
	assert inrange(c_dt_rate, 0, 100) | c_dt_rate==.

	* normalize the AREA rate to 0 to 100
	replace c_area_charge = c_area_charge / 160 * 100
	assert inlist(c_area_charge, 0, 100) | c_area_charge==.

** CODING DT PEAK HOURS

	* add info for each person of the AM and PM peaks (whether they're in the DT treatment or not)
	merge m:1 uidp using `peak_times'
	assert _m!=1
	assert _m!=2
	drop _m

	* define the time intervals of interest
	gen am_sample = inrange(t_start,  7, 13)
	gen pm_sample = inrange(t_start, 16, 22)
	label var am_sample "(AM)  7:00-13:00"
	label var pm_sample "(PM) 16:00-22:00"

	gen am_a_sample = inrange(t_start,  7, 14)
	gen pm_a_sample = inrange(t_start, 14, 22) & t_start > 14
	label var am_a_sample "(AM)  7:00-14:00"
	label var pm_a_sample "(PM) 14:01-22:00"

	assert am_a_sample + pm_a_sample == 1 | t_start == .

	* peak (3 hours)
	gen am_pk_sample = abs(t_start - pk_am) < 1.5
	gen pm_pk_sample = abs(t_start - pk_pm) < 1.5
	label var am_pk_sample "+/- 1.5h around peak (AM)"
	label var pm_pk_sample "+/- 1.5h around peak (PM)"

	* peak plus minus half hour (4 hours)
	gen am_pk5_sample = abs(t_start - pk_am) < 2
	gen pm_pk5_sample = abs(t_start - pk_pm) < 2
	label var am_pk5_sample "+/- 2h around peak (AM)"
	label var pm_pk5_sample "+/- 2h around peak (PM)"

	* peak plus/minus 1 hour (5 hours)
	gen am_pk1_sample = abs(t_start - pk_am) < 2.5
	gen pm_pk1_sample = abs(t_start - pk_pm) < 2.5
	label var am_pk1_sample "+/- 2.5h around peak (AM)"
	label var pm_pk1_sample "+/- 2.5h around peak (PM)"

	* pre peak dummy
	gen am_pre_peak = t_start < pk_am & am_a_sample == 1
	gen pm_pre_peak = t_start < pk_pm & pm_a_sample == 1
	label var am_pre_peak "pre AM peak"
	label var pm_pre_peak "pre PM peak"

	assert am_pk1_sample <= am_sample
	assert pm_pk1_sample <= pm_sample

	gen 	t_rel = t_start - pk_am if t_start <= 14.5 
	replace t_rel = t_start - pk_pm if t_start  > 14.5 
	label var t_rel "Time relative to peak"

	gen t_absrel = abs(t_rel)
	label var t_absrel "Absolute time relative to peak"

	* ramp dummies
	gen ramp_am_pre  = inrange(t_rel, - 1.5, -0.5) & am_sample == 1
	gen ramp_am_post = inrange(t_rel,   0.5,  1.5) & am_sample == 1
	gen ramp_pm_pre  = inrange(t_rel, - 1.5, -0.5) & pm_sample == 1
	gen ramp_pm_post = inrange(t_rel,   0.5,  1.5) & pm_sample == 1

	label var ramp_am_pre	"Ramp am pre" 
	label var ramp_am_post	"Ramp am post" 
	label var ramp_pm_pre	"Ramp pm pre" 	
	label var ramp_pm_post	"Ramp pm post" 
	

	// twoway (kdensity t_rel if t_start <=14.5) (kdensity t_rel if t_start >14.5) ///
	// 		, legend(order(1 "am" 2 "pm")) xlabel(-5(0.5)5, grid)

*** Check missing values (see notes below)
	mdesc
	mdesc c_area* if in_area_treatment == 1
	assert c_area_charge!=. & c_area_intersect!=. if in_area_treatment==1

	/* Notes:
	- generally missing 12,178 obs for days without trips
	- try to add d_oh and others on days with location (but no trip) data using location location
	- same for inb_frac
	- segs/sids?

	- INFO: all sample_ have missing for date_wo_trips==1. (except sample_trip_ok)
	- INFO: DT ramp has missing values
	- INFO: Area charges missing except intersect and charge, 

	- INFO: quality: jump_pre_dist not defined for no_data
	 */

	 order uidp date chain trip date_wo_trips trip_hw trip_hw2 trip_ww2 trip_hwp trip_hwpwp d_dh d_dw d_dw2 d_oh d_ow d_ow2 oh ow ow2 dh dw dw2


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

*** TREATMENT VARs CODING
	gen dt_wcha = inlist(dt_treat,1,2)
	label var dt_wcha "Any DT treatment with charges (Low or High)"

	gen dt_noch = inlist(dt_treat,3,4)
	label var dt_noch "Any DT treatment without charges (Control or Info)"

	gen post = inrange(time_in_exp,1,25) //& study_cycle != 99

	* only 3 respondents for whom we skip this day. All on 4/26.
	count if post != inlist(study_cycle,1,2,3,4,5,99)
	assert r(N) == 3
	assert post == inlist(study_cycle,1,2,3,4,5,99) | date==mdy(4,26,2017) 	// check that the definitions are equivalent
	// assert post == study_cycle != 0 if sample_analysis == 1  // check that the definitions are equivalent

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

*** Final SAMPLE selection
	drop if date == mdy(4,26,2017)

********************
** Save Trip Data **
********************
	isid date uidp chain
	save "data/coded_gps_dta/coded_trips_`chain_th'.dta", replace
