
 * This do file: preparing for demand estimation
 * Authors: Gabriel Kreindler
 
clear all
set more off
pause on

* Path based on global set in C:/ado/profile.do
	local path "$congestion_pricing_path"
	di "`path'"
	cd "`path'"

**************
*** CODING ***
**************

	*** computes the variance-covariance matrix of the estimated road technology
	do "code/analysis_simulations/code-road-tech-density-vcov.do"
	* output:
	* 	- data/coded_model/road_tech/vcov_density_linear.csv
	* 	- data/coded_model/road_tech/vcov_density_pow.csv


	*** code the route experiment results 
	do "code/coding_model/code_area_trip_vot.do"
	* output:
	*   - data/coded_model/area_trip_vot.csv

	*** code the heterogeneity route experiment results (note: I think these are no longer used in the final version of the paper)
	do "code/coding_model/code_struct_matlab_area_het.do"
	* output:
	*   - model/data_input/area_het_temp.dta

	*** code the departure time experiment results
	do "code/coding_model/code_struct_matlab_dt_het.do"
	* output:
	*   - model/data_input/dt_het.dta

	*** combine previous files and export to csv
	do "code/coding_model/code_struct.do"
	* output:
	*   - data/coded_model/coded_uidpn_a.csv    // area treatment participants
	*   - data/coded_model/coded_uidpn_na.csv   // non-area treatment participants

