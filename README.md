# Replication instructions for "Peak-Hour Road Congestion Pricing: Experimental Evidence and Equilibrium Implications"
Gabriel Kreindler

June 2023

Questions? Email: gkreindler@g.harvard.edu


This folder contains intermediate data and code to replicate the results in the paper and the supplementary materials. Raw data with personally identifiable information (PII) is not included. The analysis uses Stata for the reduced form analysis and Julia for demand analysis and policy counterfactuals.

## Data sources

The paper uses three main data sources
1. Survey data (`data/raw_cto/`)
    - `recruit` has the recruitment survey, which also served as a baseline
    - `baseline` has a phone survey of beliefs and stated preferences

2. Raw GPS pings collected using the study smartphone app in (`data/raw_gps/`)
    - this data is not included here because of the concern that it may be used together with other data to identify specific individuals
    - intermediate coded data is available in (`data/coded_gps_dta/`)

3. Google Maps API data on travel times (`data/google_maps/`)
    - panel "real-time" data on a fixed set of routes in Bangalore (files with `gmaps_citywide_live`)
    - daily profiles between home and work for a sub-sample of respondents (files with `gmaps_hw`)
    - daily profiles across multiple routes between home and work (files with `all_unique_routes_by24h`)

Other data sources:
 - data from Akbar and Duranton (2017), digitized from Figrue 4 panel C. (`data/raw_other/Akbar and Duranton (2017)/`)
 - data on vehicle values compiled from posts on the online portal www.olx.in (`data/raw_other/vehicle_values/`)

 ### Pre-coded data
 
 Several datasets are provided in an intermediate form, primarily in order to avoid including data with PII.
 - survey data is provided in .dta format after coding from the raw .csv files, because this made it easier to remove PII (e.g. names and phone numbers)
 - coded trips data based on the GPS traces is provided in the folder `coded_gps_dta`. This data does not contain the precise routes.

Treatment assignment:
- "data/coded_trips_dta/treat_assign.dta" 
- "data/treatment/treatment roster noPII.csv" (and .dta)
- "data/treatment/all_area_candidates.csv"
- "data/treatment/time_in_exp_wtrial.xlsx"
- "data/treatment/prerandomized_assignment.dta"

Experiment logs:
- "data/treatment/respondent_stages noPII.csv"
- "data/treatment/full_daily_logs.csv"

## Software

### Stata
The reduced form analysis in this paper was done with Stata 17 and Stata 18 

#### Paths
Set the macro $congestion_pricing_path to point to the root of the replication package folder.
On Windows, add the following line to the file `C:\ado\profile.do` (this script is executed each time Stata starts. Create it if it does not exist.)
```
global congestion_pricing_path="C:/.../replication_package/" // fill in your path
```

#### Packages
The code relies on the following packages
- gtools
- ftools
- mdesc
- sdecode
- winsor
- winsor2 
- newey2
- ivreg2
- ranktest
- weakivtest
- avar
- grc1leg (net install grc1leg, from( http://www.stata.com/users/vwiggins/)
- keeporder

### Julia

#### Paths
Edit the file `code/paths.jl` and add the path pointing to the root of the replication package folder:
```
rootpath_base = "C:/.../replication_package/" # fill in your path
```

#### Versions and environments
The project uses environments that allows using the exact version of Julia and packages to replicate the results.

There are two environments (an environment has a Project.toml file and a Manifest.toml file), corresponding to a local machine and to a server. (The latter is recommended for the more time-intensive model estimation and policy simulations.)

Local environment:
```
	Julia version = "1.8.1"
	replication_package\code\analysis_estimation\cp_env\
```

Server environment:
```
	Julia version = "1.9.1 "
	replication_package\code\analysis_estimation\cp_env_server\
```


# Scripts: data coding

## Survey data coding
- run `MAIN_coding_cto.do` 

## GPS data and experimental coding
- run `MAIN_coding_gps.do` 

## Model coding
- run `MAIN_coding_model.do`

# Scripts: analysis (tables and figures)

For Stata, run the main file:
- run in Stata `MAIN_figures.do`
- run in Stata `MAIN_tables.do`
- run in Stata `MAIN_sm_figures.do`
- run in Stata `MAIN_sm_tables.do`

For Julia, run the main files:
- run in Julia `MAIN_estimation.jl` 
	- this takes long. it is better to run each sub-script independently, ideally on a computer cluster.
- run in Julia `MAIN_estimation_figures_tables.jl`
- run in Julia `MAIN_simulation.jl`
	- this takes long. it is better to run each sub-script independently, ideally on a computer cluster.
- run in Julia `MAIN_simulation_figures_tables.jl`


## Main Figures

```
	Figure 1: Average Travel Delay in the Study Region in Bangalore
		run "code/analysis_reduced_form/Figure 1 google maps.do"

	Figure 2: Impact of Departure Time Charges on the Distribution of Departure Times
		run "code/analysis_reduced_form/Figure 2 departure time.do"


	Figure 3: Impact of Route Charges on Average Detour Route Usage
		run "code/analysis_reduced_form/Figure 3 area.do"


	Figure 4: Road Technology Supply Estimation: Travel Delay Linear in Traffic Density
		run "code/analysis_reduced_form/Figure 4 road tech.do"


	Figure 5: Policy Simulation: Unpriced Nash Equilibrium and Social Optimum
		simulations:
			run "code/analysis_simulations/runsims-main.jl"      
			run "code/analysis_simulations/runsims-main-boot.jl" 

		assemble figure:
			run "code/analysis_simulations/paper_figure5.jl" // also generates panel B in Figure SM9
```


## Main Tables

```
	Table I: Impact of Departure Time Charges on Daily Outcomes
		run "code/analysis_reduced_form/Table 1 departure time.do"


	Table II: Travel Demand Parameter Estimates
		estimation:
			run "code/analysis_estimation/estimation_main_benchmark.jl" (column 1)
			run "code/analysis_estimation/estimation_main_fixed_vott.jl" (columns 2 and 3)
			run "code/analysis_estimation/estimation_main_no_dt.jl" (column 4)

		assemble table:
			run "code/analysis_estimation/paper_table_II.jl"

		notes:
			- the estimation scripts can take a long time (especially for bootstrap)
			- the estimation scripts are set up to run in parallel 
			- shell script examples for running the code on a computer cluster are provided


	Table III: Road Technology Supply Estimation: Travel Delay Linear in Traffic Density
		run "code/analysis_reduced_form/Table III road tech.do"

	Table IV:  Policy Simulations: Commuter Welfare Gains from Optimal Peak-Hour Pricing
		simulations:
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


	Table V: Policy Simulations: Equilibrium with Extensive Margin Decision
		simulations:
			run "code/analysis_simulations/runsims-ext.jl"
		
		assemble table:
			run "code/analysis_simulations/paper_table_V.jl"
```


## SM Figures

```
	Figure SM.1. Example deadweight loss versus schedule cost (holding observed equilibrium fixed)
		- available upon request

	Figure SM.2. Impact of Departure Time Charges on Departure Times (Commuting Trips)  
		run "code/analysis_reduced_form/Figure 2 departure time.do"
		note: run time around 10-20 minutes
 
	Figure SM.3. Travel Demand Model Fit
		run "code/analysis_estimation/paper_figure_model_fit.jl"

	Figure SM.4. Travel Demand Model: Understanding Identification
		run "code/analysis_estimation/model_identification/model_sensitivity_full.jl"

	Figure SM.5: Travel Demand Model Numerical Identification Check and Finite Sample Properties
		estimation:
			run "code/analysis_estimation/numerical_identification/simulate_data_numerical_check.jl"
			run "code/analysis_estimation/numerical_identification/estimation_simulated_data.jl"
			
		assemble figures and table:
			run "code/analysis_estimation/numerical_identification/figure_sm_4.jl"
			run "code/analysis_estimation/numerical_identification/table_sm_12_part1.jl" <- Julia
			run "code/analysis_estimation/numerical_identification/table_sm_12_part2.do" <- Stata (to use quantile regression and output neatly to a latex table)

	Figure SM.6: Road Technology Estimation Robustness Checks
		panel A: run "code/analysis_reduced_form/sm/Figure SM5A road tech percentiles.do"
		panel B: run "code/analysis_reduced_form/sm/Figure SM5B road tech volume.do"
		panel C: run "code/analysis_reduced_form/sm/Figure SM5C road tech recruitment.do"
		panel D: run "code/analysis_reduced_form/sm/Figure SM5D road tech compare.do"
		panel E: run "code/analysis_reduced_form/sm/Figure SM5E road tech pigou knight.do"

	Figure SM.7: Road Technology at the Daily Level 
		run "code/analysis_reduced_form/sm/Figure SM6 road tech daily.do"

	Figure SM.8: Road Technology on Major Arteries
		run "code/analysis_reduced_form/sm/Figure SM7 road tech artery level.do"

	Figure SM.9: Policy Counterfactual Additional Results
		panel A: run "code/analysis_simulation/paper_smfig8A_charges"
		panel B: run "code/analysis_simulation/paper_figure5.jl"
		panel C: run "code/analysis_simulations/paper_smfig8C_road_tech_powers.do" 
		panel D: run "code/analysis_simulations/paper_smfig8D_factor"

	Figure SM.10: Decomposing Gains and Losses in the Social Optimum
		run "code/analysis_simulations/paper_smfig9_decomposition.jl" 

	Figure SM.11:
		panel A:    run "code/analysis_simulations/paper_smfig10_panelA_road_tech_2routes.do"
		panels B&C: run "code/analysis_simulations/paper_smfig10_2routes.jl" 
```

## SM Tables

```
	Table SM.I: Descriptive Statistics about Travel Behavior
		run "code/analysis_reduced_form/sm/Table SM1.do"

	Table SM.II: Experimental Design
		no code provided

	Table SM.III: Experimental Participant Sample Representativeness
		run "code/analysis_reduced_form/sm/Table SM1.do"

	Table SM.IV: Experimental Balance Checks
		run "code/analysis_reduced_form/sm/Table SM4 balance experiment.do"

	Table SM.V: GPS Data Quality at Daily Level (Attrition Check)
		run "code/analysis_reduced_form/sm/Table SM5 attrition.do"

	Table SM.VI: Impact of Departure Time Charges on Daily Total Hypothetical Rate: Commuting Trips
		run "code/analysis_reduced_form/sm/Table SM6 departure time trip level.do"

	Table SM.VII: Impact of Route Charges on Detour Route Usage
		run "code/analysis_reduced_form/sm/Table SM7 area.do"

	Table SM.VIII: Impact of Route Charge Sub-Treatments on Daily Outcomes
		run "code/analysis_reduced_form/sm/Table SM8 area subtreats.do"

	Table SM.IX: Travel Demand Estimates: Additional Results
		estimation:
			run "code/analysis_estimation/estimation_robust*.jl" (6 files)
			
		assemble table:
			run "code/analysis_estimation/paper_table_robustness.jl"	


	Table SM.X: Travel Demand Estimation: Discount Factor Robustness
		estimation:
			run "code/analysis_estimation/estimation_delta_other.jl" 	  (columns 1-4)
			run "code/analysis_estimation/estimation_delta_d99.jl"   	  (column  4)
			run "code/analysis_estimation/estimation_delta_estimation.jl" (columns 5)
			
		assemble table:
			run "code/analysis_estimation/paper_table_delta.jl"


	Table SM.XI: Dynamic Route Choice Model Identification
		estimation:
			run "code/analysis_estimation/model_identification/model_sensitivity_full.jl"  -> also generates Figure SM.4 
			run "code/analysis_estimation/model_identification/model_sensitivity_vott.jl"
			run "code/analysis_estimation/model_identification/model_sensitivity_simple.jl"
			
		assemble table:
			run "code/analysis_estimation/model_identification/table_sensitivity_vott.jl"


	Table SM.XII: Travel Demand Model Finite Sample Properties Check
		see notes for SM Figure 5
```

## Numbers in the text

The exchange rates used in the paper are from 2017: 64.39 INR market rate and 20.65 INR PPP.

### section "Travel Demand Estimation Results"
```

(*) in the paragraph starting with "To explain what these numbers imply in terms of behavior,"

	run "code/analysis_estimation/calculations_in_text/text_dt_responsiveness.jl" // to get the response to charges

	run "code/analysis_simulations/paper_smfig9_decomposition.jl" // to show that commuters in the estimated (equilibrium) model arrive on average 2 minutes early before peak and 5 minutes late after the peak

(*) to get ratio and bootstrapped CIs in paragraph starting with 
		"There are few estimates of schedule costs in the transportation economics literature to compare with." 

	run "code/analysis_estimation/calculations_in_text/text_boot_ratio.jl"    // CI on ratio
	run "code/analysis_estimation/calculations_in_text/text_ratio_implied.do" // implies ratios

(*) to get new beta_E and beta_L estimates when reduced form departure time impact profile were attenuated by a factor of 2, in the paragraph starting with:
		"The model propagates reduced form moment uncertainty to parameter estimates in a relative transparent manner."

	run "code/analysis_estimation/calculations_in_text/text_half_DT_part1.jl" // estimation
	run "code/analysis_estimation/calculations_in_text/text_half_DT_part2.jl" // stats
```

### section on "Supply Estimation: Road Congestion Technology"
```
(*) sample size numbers
	run "code/coding_gps/code_road_tech_density_panel.do"

(*) motorcycle numbers
	run "code\analysis_reduced_form\other\text_moto_vs_car.do"



```

### section on "Policy Simulations"
```
(*) In the paragraph starting with "To find the social optimum, I use a nesting iterative procedure..." to get the correlation between partial- and equilibrium- externalities, run:
	"code/analysis_simulations/calculations_in_text/text_compare_msc_partial_vs_full_at_nash.jl"

(*) In the paragraph starting with "Optimal congestion charges lead to small but notable travel time reductions" to get the 10th and 90th percentiles, run:
	"code/analysis_simulations/paper_smfig9_decomposition.jl" (line 162)

(*) Paragraph starting with
		"Welfare gains from optimal pricing are very sensitive to the road technology, especially its slope at the peak."

	To get the elasticity and other numbers:

	run "code/analysis_simulations/paper_smfig8D_factor.jl"


(*) Paragraph starting with	"In panel D." To get the slopes at maximum density, run: 
	"code/analysis_simulations/paper_smfig8C_road_tech_powers.do"

(*) Paragraph starting with (and in Table V notes) "In Table V..." To get the min/max charges at the social optimum, run:
	"code/analysis_simulations/paper_table_V.jl"


```