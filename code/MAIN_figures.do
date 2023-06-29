
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

	*** Figure 1: Average Travel Delay in the Study Region in Bangalore
		do "${congestion_pricing_path}code/analysis_reduced_form/Figure 1 google maps.do"

	*** Figure 2: Impact of Departure Time Charges on the Distribution of Departure Times
		run "${congestion_pricing_path}code/analysis_reduced_form/Figure 2 departure time.do"
		// manually change lines 27, 28 to generate all 4 combinations AM/PM and all/hw trips
		// note: this may take several minutes to run


	*** Figure 3: Impact of Route Charges on Average Detour Route Usage
		run "${congestion_pricing_path}code/analysis_reduced_form/Figure 3 area.do"


	*** Figure 4: Road Technology Supply Estimation: Travel Delay Linear in Traffic Density
		run "${congestion_pricing_path}code/analysis_reduced_form/Figure 4 road tech.do"


	*** Figure 5: Policy Simulation: Unpriced Nash Equilibrium and Social Optimum

	*** Run in Julia

/* 		simulations:
			run "code/analysis_simulations/runsims-main.jl"      
			run "code/analysis_simulations/runsims-main-boot.jl" 

		assemble figure:
			run "code/analysis_simulations/paper_figure5.jl" // also generates panel B in Figure SM9

 */	
