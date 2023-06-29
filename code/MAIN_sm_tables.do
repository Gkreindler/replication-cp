
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
*** MAIN FIGURES
**************

	*** Table SM.I: Descriptive Statistics about Travel Behavior
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM1.do"

	*** Table SM.II: Experimental Design
		// by hand

	*** Table SM.III: Experimental Participant Sample Representativeness
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM3.do"

	*** Table SM.IV: Experimental Balance Checks
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM4 balance experiment.do"

	*** Table SM.V: GPS Data Quality at Daily Level (Attrition Check)
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM5 attrition.do"

	*** Table SM.VI: Impact of Departure Time Charges on Daily Total Hypothetical Rate: Commuting Trips
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM6 departure time trip level.do"

	*** Table SM.VII: Impact of Route Charges on Detour Route Usage
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM7 area.do"

	*** Table SM.VIII: Impact of Route Charge Sub-Treatments on Daily Outcomes
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Table SM8 area subtreats.do"

	*** Table SM.IX: Travel Demand Estimates: Additional Results
		/* estimation:
			run "code/analysis_estimation/estimation_robust*.jl" (6 files)
			
		assemble table:
			run "code/analysis_estimation/paper_table_robustness.jl"	
 */

	*** Table SM.X: Travel Demand Estimation: Discount Factor Robustness
/* 		estimation:
			run "code/analysis_estimation/estimation_delta_other.jl" 	  (columns 1-4)
			run "code/analysis_estimation/estimation_delta_d99.jl"   	  (column  4)
			run "code/analysis_estimation/estimation_delta_estimation.jl" (columns 5)
			
		assemble table:
			run "code/analysis_estimation/paper_table_delta.jl"
 */

	*** Table SM.XI: Dynamic Route Choice Model Identification
/* 		estimation:
			run "code/analysis_estimation/model_identification/model_sensitivity_full.jl"  -> also generates Figure SM.4 
			run "code/analysis_estimation/model_identification/model_sensitivity_vott.jl"
			run "code/analysis_estimation/model_identification/model_sensitivity_simple.jl"
			
		assemble table:
			run "code/analysis_estimation/model_identification/table_sensitivity_vott.jl"
 */

	*** Table SM.XII: Travel Demand Model Finite Sample Properties Check
		// see notes for SM Figure 5
