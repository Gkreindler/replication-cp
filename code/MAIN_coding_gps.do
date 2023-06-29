
 * This do file: MASTER for the coding of raw GPS data into trips
 * This version: August 11 2018
 * Authors: Gabriel Kreindler

clear all
pause on
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"


**********
* CODING Python
**********
 	* these scripts transform the raw GPS data in data on trips (with start and end)
 	* the raw data is not included in this replication package because

	*** Main script:
	* - coding_gps_py/tripcoding_main.py
	* 	  - input: 
	*			device-level raw gps data from data/raw_gps/coded_btr_...
	* 			crosswalk between survey and app IDs: data/coded_cto/crosswalk_uidp_deviceid.csv
	*	  - intermediate files in `data/coded_gps_byid/'
	*     - outputs: 
	*			everything in "data/coded_gps/..."
	*			data/coded_gps/sample_exp_chains_**.csv         // sample of trips for the experiment analysis
	*			data/coded_gps/sample_complete_chains_**.csv    // sample of trips for the road tech analysis
	*			data/coded_gps/quality_all.csv 					// user-day level information on GPS data quality

	*** Compute hypothetical charges for everyone in the experiment (also defined for control group). These depend on trips, type of trip, departure time, and whether the trip route intersects the congestion area.
	* - coding_gps_py/analysis_compute_shadow_charges.py
	*     - output: data/coded_gps/charges_daily_all_**.csv 	// charges at the DAY level
	*     - output: data/coded_gps/charges_trips_area_**.csv 	// area charges at the trip level
	*     - output: data/coded_gps/charges_trips_dt_**.csv 	 	// departure time charges at the trip level
	
	*** Compute number of trips by artery
	* - coding_gps_py/arteries/load_trips.py
	* 	- output: data/raw_other/road tech artery level/btm_volumes/...

	*** Note:
 	* 	The python script uses hard-coded paths. One option is to use symbolic links. For Windows users:
	*	1. Open Command Prompt in Administrator Mode: Search for **Command Prompt** > Right click on it and select **Run as administrator**.
	*	2. Enter the following command in Command Prompt (choose the one with your name):
	*		mklink /J C:\bang_cp_paper\ "C:\Users\Gabriel K\Dropbox (MIT)\0Re Projects\bang_cp_paper\"
	*	3. If successful, you should see the following message:
	*		`Junction created for C:\b_traffic_pilot2\ <<===>> <your path>`
	*	If you would like to remove the symbolic link, do not delete the "shortcut" as it would instead delete the original directory it was pointing to. Instead, use **Command Prompt** and run:
	*		`rmdir C:\bang_cp_paper\`


**********
* CODING STATA
**********

*** Process area candidates
	do "$congestion_pricing_path/code/coding_gps/load_area_candidates.do"
	* input:  data/treatment/all_area_candidates.csv
	* output: coding_gps/area_candidates.dta (uidp level, 258)


*** get reports from log files
	do "$congestion_pricing_path/code/coding_gps/study_dates.do"
	* (uidp-date level, only during the study)
	* outputs:
	* coded_gps_dta/study_dates_full.dta - inlcudes outstation, trial, special
	* coded_gps_dta/study_dates_wtrial.dta - inlcudes trial
	* coded_gps_dta/study_dates.dta - inlcudes only experiment

*** PREPARE TREATMENT sample and assignment
	do "$congestion_pricing_path/code/coding_gps/treatment_assignment.do"
	* (uidp level)
	* Output: coded_gps_dta/treat_assign.dta

*** load quality
	do "$congestion_pricing_path/code/coding_gps/load_quality_exp.do" 
	* (uidp-date level, includes pre period, drop dates after the experiment if the respondents dropped)
	* Output: coded_gps_dta/quality_exp.dta

	do "$congestion_pricing_path/code/coding_gps/load_quality_complete.do" 
	* (uidp-date level all users and days)
	* Output: coded_gps_dta/quality_complete.dta	


*** Load Trips (chains)

	* only experimental sample
	do "$congestion_pricing_path/code/coding_gps/load_chains_exp.do" 15  // Output: coded_gps_dta/chains_exp_15.dta
	do "$congestion_pricing_path/code/coding_gps/load_chains_exp.do" 30  // Output: 
	do "$congestion_pricing_path/code/coding_gps/load_chains_exp.do" 60  // Output: 

	* complete (all users, all days)
	do "$congestion_pricing_path/code/coding_gps/load_chains_complete.do" 15  // Output: coded_gps_dta/chains_complete_15.dta
	do "$congestion_pricing_path/code/coding_gps/load_chains_complete.do" 30  // Output: 
	do "$congestion_pricing_path/code/coding_gps/load_chains_complete.do" 60  // Output: 


*** Load Charges

	* DT charges at trip level
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_dt.do" 15   // Output: coded_gps_dta/charges_trips_dt_15.dta
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_dt.do" 30   // Output: 
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_dt.do" 60   // Output: 

	* AREA charges at trip level
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_area.do" 15   // Output: coded_gps_dta/charges_trips_area_15.dta
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_area.do" 30   // Output:
	do "$congestion_pricing_path/code/coding_gps/load_charge_trips_area.do" 60	 // Output:

	* DT charges at DAY level
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_dt.do"  15 	// Output: coded_gps_dta/charges_days_dt_`chain_th'.dta
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_dt.do"  30	// Output: 
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_dt.do"  60 	// Output:

	* AREA charges at DAY level
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_area.do" 15  // Output: coded_gps_dta/charges_days_area_`chain_th'.dta
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_area.do" 30  // Output: 
	do "$congestion_pricing_path/code/coding_gps/load_charge_days_area.do" 60  // Output: 


*** Code trips
	do "$congestion_pricing_path/code/coding_gps/code_trips.do" 15	// Output: coded_gps_dta/coded_trips_15.dta
	do "$congestion_pricing_path/code/coding_gps/code_trips.do" 30
	do "$congestion_pricing_path/code/coding_gps/code_trips.do" 60 // L0419164627
	* input: 
	*	treatment roster noPII.csv
	*	treat_assign.dta
	*	chains_**.dta
	*	charges_trips_**.dta
	*	charges_trips_area_**.dta
	*	quality.dta
	
	* output:
	* DESCRIPTION: 	trip level, includes days without trips (no data or no trips) with one observations,
	* 				where most vars are missing. Includes DT and AREA charges at trip level day sample 
	*				includes outstation, special and trial.
	*	coded_gps_dta/coded_trips_**.dta

	* For for entire sample
	do "$congestion_pricing_path/code/coding_gps/code_trips_complete.do" 15	// Output: coded_gps_dta/coded_trips_complete_15.dta
	do "$congestion_pricing_path/code/coding_gps/code_trips_complete.do" 30
	do "$congestion_pricing_path/code/coding_gps/code_trips_complete.do" 60

*** Code at DAY level
	do "$congestion_pricing_path/code/coding_gps/code_days.do" 15	// Output: coded_gps_dta/coded_days_15.dta
	do "$congestion_pricing_path/code/coding_gps/code_days.do" 30
	do "$congestion_pricing_path/code/coding_gps/code_days.do" 60
	* input: 
	*	coded_gps_dta/coded_trips.dta
	*	coded_gps_dta/charges_days_15.dta
	*	coded_gps_dta/charges_days_area_15.dta
	*	coded_gps_dta/treat_assign.dta
	
	* output:
	* DESCRIPTION: 	day level, includes days without trips (no data or no trips), includes DT and AREA 
	*				charges at DAY level, day sample includes outstation, special and trial.
	* 	coded_gps_dta/coded_days_**.dta

	do "$congestion_pricing_path/code/coding_gps/code_days_peak.do" 15 // Output: coded_gps_dta/coded_days_peak_15.dta
	do "$congestion_pricing_path/code/coding_gps/code_days_peak.do" 30
	do "$congestion_pricing_path/code/coding_gps/code_days_peak.do" 60
	* input: 
	*	coded_trips_`chain_th'.dta
	
	* output:
	* DESCRIPTION: 	day level, includes several variables counting the (smoothed) number of trips in bins relative to peak profile
	* 	coded_gps_dta/coded_days_peak_**.dta


*** Coding for route charges treatment analysis
	do "$congestion_pricing_path/code/coding_gps/code_area.do"
	* output:
	* 	- data/coded_gps_dta/area_coded.csv   (for reduced form analysis)
	* 	- data/coded_model/area_trip_vot.csv   (for GMM analysis)



*** Coding for road technology analysis (Figure 4, Table 3 and others)
	do "$congestion_pricing_path/code/coding_gps/code_road_tech_density_panel.do"
	* output:
	* - data/coded_road_tech/density_tm    // day minute
	* - data/coded_road_tech/density_th    // day hour
	* - data/coded_road_tech/density_th3   // day 20-minute
	* - data/coded_road_tech/density_t     // day
	

	do "$congestion_pricing_path/code/coding_gps/code_road_tech_density.do"
	* output:
	* - data/coded_road_tech/density_mdd  // minute level
	* - data/coded_road_tech/density_h    // hourly level 
	* - data/coded_road_tech/density_h3   // 20 minutes level

	do "$congestion_pricing_path/code/coding_gps/code_road_tech_speeds_volume.do"
	* - data/coded_road_tech/google_maps_th_level.dta
	* - data/coded_road_tech/google_maps_th3_level.dta
	* - data/coded_road_tech/google_maps_t_level.dta
	* - data/coded_road_tech/google_maps_h_level.dta
	* - data/coded_road_tech/google_maps_h_level.csv
	* - data/coded_road_tech/google_maps_h3_level.dta
	* - data/coded_road_tech/google_maps_h3_level.csv

	* - data/coded_road_tech/volumes_th_level   // day hour level
	* - data/coded_road_tech/volumes_th3_level  // day 20-minute level
	* - data/coded_road_tech/volumes_t_level 	// date level
	* - data/coded_road_tech/volumes_h_level 	// hour level
	* - data/coded_road_tech/volumes_h3_level 	// 20-minute level
