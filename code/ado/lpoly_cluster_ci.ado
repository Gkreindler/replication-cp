* Description:  This ado file defines a program to generate clustered 
*		bootstrapped confidence bands for lpoly

* Date created: October 17, 2016.
* Created by: Gabriel Kreindler and Michael Fryar.

capture program drop lpoly_cluster_ci	// Drop old program

program define lpoly_cluster_ci		// Define new program
version 13				// Ensure backwards compatibility

/*
One dependent variable
[if] works as usual
cluster = a variable, required 
niter = bootstrap_iterations, optional--defaults to 100
bandwidth is fixed at 0.2 - should change for other applications
*/

syntax varlist(min=2 max=2) [if], cluster(varname) bw(real) [niter(integer 100)]

* Save data
tempfile original_data
qui save `original_data', replace

* One dependent variable	
tokenize `varlist'
local yvar "`1'"
local xvar "`2'"

* Keep subsample if provided
if("`if'" != ""){
	keep `if'
}

* Number of iterations
local bootstrap_iterations = `niter'
// local clustervar = "`cluster'"
if "`bw'" == ""{
	local bw = ""
}
else{
	local bw = "bw(`bw')"
}


* Prep grid
preserve
	keep `xvar'
	rename `xvar' `xvar'_lpoly
	qui duplicates drop
	sort `xvar'_lpoly
	gen n = _n

	tempfile atxvar
	qui save `atxvar'
restore

* Keep only essential variables!
keep `yvar' `cluster' `xvar'

forv iter = 1/`bootstrap_iterations'{
	preserve
		* Bootstrap sampling with replacement, at cluster level
		bsample, cluster(`cluster')

		gen n = _n
		qui merge 1:1 n using `atxvar'
		
		// Assert _m!=2
		qui count if _m==2
		if (r(N) > 0){
			noisily display "No data for `r(N)' x variable slots"
			drop if _m==2
		}
		drop _m

		* Run lpoly on current sample and save predicted at xvar
		qui lpoly `yvar' `xvar', `bw' ///
			gen(`yvar'_lpoly) at(`xvar'_lpoly) nograph
		
		noisily display "Iteration `iter' done."

		* keep only lpoly results
		keep `xvar'_lpoly `yvar'_lpoly

		* keep only lpoly results
		qui drop if `xvar'_lpoly == .
		gen iter = `iter'

		tempfile 	bsample`iter'
		qui save 	`bsample`iter''
	restore	
}

* Load results and compute confidence bands
	use `bsample1', clear

	forv iter=2/`bootstrap_iterations'{
		append using `bsample`iter''
	}

	/* Generate 2.5th, 50th and 97.5th percentiles (to compare with point 
	estimate) */
	bys `xvar'_lpoly: egen `yvar'_p025 = pctile(`yvar'_lpoly), p(2.5) 
	bys `xvar'_lpoly: egen `yvar'_median = median(`yvar'_lpoly) 
	bys `xvar'_lpoly: egen `yvar'_p975 = pctile(`yvar'_lpoly), p(97.5) 

	keep `xvar'_lpoly `yvar'_p* `yvar'_m*
	
	qui duplicates drop
	
	tempfile percentiles
	qui save `percentiles'

* Restore original data plus the confidence bands
	use `original_data', clear
	gen double `xvar'_lpoly = `xvar'
	
	qui merge m:1 `xvar'_lpoly using `percentiles'
	assert _m != 2
	drop _m

local closingmessage "Clustered bootstrapped confidence bands of `yvar' by"
local closingmessage "`closingmessage' `xvar' merged with original data."
noisily display "`closingmessage'"

end
