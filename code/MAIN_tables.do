
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

	*** Table I: Impact of Departure Time Charges on Daily Outcomes
		do "${congestion_pricing_path}code/analysis_reduced_form/Table I departure time.do"


	*** Table II: Travel Demand Parameter Estimates
		/* estimation:
			run "code/analysis_estimation/estimation_main_benchmark.jl" (column 1)
			run "code/analysis_estimation/estimation_main_fixed_vott.jl" (columns 2 and 3)
			run "code/analysis_estimation/estimation_main_no_dt.jl" (column 4)

		assemble table:
			run "code/analysis_estimation/paper_table_II.jl"

		notes:
			- the estimation scripts can take a long time (especially for bootstrap)
			- the estimation scripts are set up to run in parallel 
			- shell script examples for running the code on a computer cluster are provided
 */

	*** Table III: Road Technology Supply Estimation: Travel Delay Linear in Traffic Density
		do "${congestion_pricing_path}code/analysis_reduced_form/Table III road tech.do"

	*** Table IV:  Policy Simulations: Commuter Welfare Gains from Optimal Peak-Hour Pricing
		/* simulations:
			run "code/analysis_simulations/runsims-main.jl"       // panel A
			run "code/analysis_simulations/runsims-main-boot.jl" 
			run "code/analysis_simulations/runsims-het.jl"        // panel B
			run "code/analysis_simulations/runsims-het-boot.jl"   
			run "code/analysis_simulations/runsims-varyparams.jl" // panel C
			run "code/analysis_simulations/runsims-power.jl"      // panel D
			run "code/analysis_simulations/runsims-2routes.jl"    // panel E

		assemble table:
			run "code/analysis_simulations/paper_table_IV.jl"  

		notes:
			- these rely on estimated parameters, so those scripts need to run first
			- also can take long (especially for bootstrap)
			- also can run in parallel and on computer cluster
 */

	*** Table V: Policy Simulations: Equilibrium with Extensive Margin Decision
		/* simulations:
			run "code/analysis_simulations/runsims-ext.jl"
		
		assemble table:
			run "code/analysis_simulations/paper_table_V.jl"
 */
