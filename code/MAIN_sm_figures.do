
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

	*** Figure SM.1. Example deadweight loss versus schedule cost (holding observed equilibrium fixed)
		// - available upon request

	*** Figure SM.2. Impact of Departure Time Charges on Departure Times (Commuting Trips)  
		// do "code/analysis_reduced_form/Figure 2 departure time.do"
		// note: already generated for Figure 2
 
	*** Figure SM.3. Travel Demand Model Fit
		// run  "code/analysis_estimation/paper_figure_model_fit.jl"

	*** Figure SM.4. Travel Demand Model: Understanding Identification
		// run "code/analysis_estimation/model_identification/model_sensitivity_full.jl"

	*** Figure SM.5: Travel Demand Model Numerical Identification Check and Finite Sample Properties estimation:
			// run "code/analysis_estimation/numerical_identification/simulate_data_numerical_check.jl"
			// run "code/analysis_estimation/numerical_identification/estimation_simulated_data.jl"
			
		// assemble figures and table:
			// run "code/analysis_estimation/numerical_identification/figure_sm_4.jl"
			// run "code/analysis_estimation/numerical_identification/table_sm_12_part1.jl" <- Julia
			do "${congestion_pricing_path}code/analysis_estimation/numerical_identification/paper_table_sm_12_part2.do" // <- Stata (to use quantile regression and output neatly to a latex table)

	*** Figure SM.6: Road Technology Estimation Robustness Checks
		*** panel A: 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM5A road tech percentiles.do"
		*** panel B: 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM5B road tech volume.do"
		*** panel C: 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM5C road tech recruitment.do"
		*** panel D: 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM5D road tech compare.do"
		*** panel E: 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM5E road tech pigou knight.do"

	*** Figure SM.7: Road Technology at the Daily Level 
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM6 road tech daily.do"

	*** Figure SM.8: Road Technology on Major Arteries
		do "${congestion_pricing_path}code/analysis_reduced_form/sm/Figure SM7 road tech artery level.do"

	*** Figure SM.9: Policy Counterfactual Additional Results
		// panel A: run "code/analysis_simulation/paper_smfig8A_charges.jl"
		// panel B: run "code/analysis_simulation/paper_figure5.jl"
		***panel C: 
		do "${congestion_pricing_path}code/analysis_simulations/paper_smfig8C_road_tech_powers.do" 
		// panel D: run "code/analysis_simulations/paper_smfig8D_factor.jl"

	*** Figure SM.10: Decomposing Gains and Losses in the Social Optimum
		// run "code/analysis_simulations/paper_smfig9_decomposition.jl" 

	*** Figure SM.11:
		*** panel A:    
		do "${congestion_pricing_path}code/analysis_simulations/paper_smfig10_panelA_road_tech_2routes.do"
		// panels B&C: run "code/analysis_simulations/paper_smfig10_2routes.jl" 
