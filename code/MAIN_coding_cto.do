
 * This do file: A guide document on the order to run coding and analysis for "Congestion Pricing" project
 * This version: 6 August, 2018
 * Authors: Anupriya Khemka, Gabriel Kreindler
 
clear all
set more off
pause on
// ssc install unique

* Path based on global set in C:/ado/profile.do
	local path "$congestion_pricing_path"
	di "`path'"
	cd "`path'"

******************************************
*** Survey and related data coding
******************************************

*** Coding Recruitment Surveys
	* Phone survey that collected app device ID for app installs
	do "$congestion_pricing_path/code/coding_cto/code_devid_v1.do"

	* Recruitment survey coding
	do "$congestion_pricing_path/code/coding_cto/code_recruitment.do"

	* UIDP-deviceid crosswalk
	do "$congestion_pricing_path/code/coding_cto\code_deviceid_uidp_crosswalk.do"	

*** Coding auxiliary data and merge with recruitment survey
	* Run the occupation, petrol pump and vehicle_prices do files 

   * collect manually cleaned occupation variables (3: detail, medium, coarse)
   * Coarse is a strict subset of indetail. indetail is not strictly included in medium (same indetail occ may appear in two medium occ)
   do "$congestion_pricing_path/code/coding_cto/code_occupation.do"
   * output: data/coded_cto/other/occupations_clean.dta

	* runs the petrol pump name do file that generates 1 new petrol pump variables
   do "$congestion_pricing_path/code/coding_cto/code_petrol_pump_codes.do"
	* output: data/coded_cto/other/pumps_clean.dta
   
   * runs the merge_recruit_olx_v8 
   do "$congestion_pricing_path/code/coding_cto/code_vehicle_values.do"
	* input: recruit_with_dict and olx_stats_brand/_model
	* output: data/vehicle_values/codeddata/vehicle_values_prices.dta (final recruitment with standardized vehicle brand and model and logprice)

	* NOT CHECKED IN DETAIL 
	do "$congestion_pricing_path/code/coding_cto/code_recruitment_other.do"	

*** Coding Stated Preferences Survey (called "baseline" for historical reasons)

	* cleaning Stated Preferences Survey
	do "$congestion_pricing_path/code/coding_cto/code_baseline.do"

	* coding variables and samples
	do "$congestion_pricing_path/code/coding_cto/code_baseline_part2.do"

******************************************
*** Google Maps data coding
******************************************

*** Google Maps (city-wide)
	do "$congestion_pricing_path/code/coding_maps/code_google_maps_citywide.do"
	* input: data/google_maps/gmaps_citywide_live.csv
	*    	 data\google_maps\gmaps_citywide_query.csv
	* output: data/google_maps/gmaps_citywide_live_coded.dta
	*		  data\google_maps\gmaps_citywide_query_btm.dta


*** Google Maps (home-work)
	do "$congestion_pricing_path/code/coding_maps/code_google_maps_hw.do"
	* input:	data/google_maps/gmaps_hw_query.csv
	*			data/google_maps/gmaps_hw_predicted.csv
	
	* output:	data\google_maps\gmaps_hw_predicted.dta
	*			data\google_maps\gmaps_hw_query.dta
