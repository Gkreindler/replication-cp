* Empirical Bayes shrinkage estimator based on Morris (1983)
* Adam Sacarny
* The latest version of this program can be found at http://www.sacarny.com/programs/

* ebayes fe se [vars] [if], [absorb() gen() bee() theta() var() uvar() rawvar() by()
*    simple tol() maxiter() ]

* fe - The fixed effects (or equivalently, the noisily estimated objects)
* se - Their standard errors (not variance!) (equivalently, the estimate of the
*      standard deviation of the noise in the objects)

* all arguments below are optional:

* vars - RHS vars in the WLS regression that is used to estimate the underlying mean
*        of the fe (i.e. theta=x*beta). This can be left blank in which case the
*        underlying mean is just a constant. Can be used with absorb() to allow the
*        underlying mean to equal a fixed effect at a higher level + unit level
*        covariates (e.g. theta = x*beta + market level FE)
* absorb() - allow the underlying mean to equal a fixed effect at the level of this
*            variable. can be used with or without listing any 'vars' (see above)
* gen() - variable to put empirical bayes adjusted fe
* bee() - variable to put attenuation factor ("B")
* theta() - variable to put underlying mean (theta)
*           this is the value toward which the fe was attenuated (e.g. x*beta)
* var() - variable to put underlying variance of the fe distribution
*         (conditional on the RHS vars)
* uvar() - variable to put unconditional underlying variance of the fe distribution
*         (underlying variance of the fe + variance of theta)
* rawvar() - variable to put the observed variance of the raw, noisily estimated fe
*            given as first argument to the program
*            (this variance is weighted using the same weights as var() and uvar())
* by() - perform empirical bayes procedure at this variable's level. the procedure will
*        be run separately for each value this variable takes on. e.g. by(year) would
*        run the procedure for each year, one at a time.
* simple - don't iterate the EB procedure. instead, use constant weights to estimate
*          the underlying mean, variance, and attenuation factor. equivalent to taking
*          the output from the first iteration of the estimation procedure.
* tol() - tolerance for convergence, taken off the estimate of the underlying
*         conditional variance (default 0.000001)
* maxiter() - maximum # of iterations (default 100)

program define ebayes, sortpreserve

	syntax varlist(numeric min=2) [if] , [ GENerate(name) ///
		Bee(name) Theta(name) Var(name) UVar(name) RAWVar(name) ///
		by(string) ABSorb(varname numeric) simple ///
		tol(real 0.000001) maxiter(integer 100) ]
	
	marksample touse
	
	gettoken fe se: varlist
	gettoken se coeff: se
	
	if ("`bee'"=="") {
		tempvar bee
	}
	if ("`theta'"=="") {
		tempvar theta
	}
	if ("`var'"=="") {
		tempvar var
	}
	if ("`uvar'"=="") {
		tempvar uvar
	}
	if ("`rawvar'"=="") {
		tempvar rawvar
	}
	
	qui gen `bee' = .
	qui gen `theta' = .
	qui gen `var' = .
	qui gen `uvar' = .
	qui gen `rawvar' = .
	
	tempvar w fevar tsum
	qui gen `w' = .
	qui gen `fevar' = `se'^2
	qui gen `tsum' = .
	
	if ("`by'"!="") {
		sort `by'
		local by = "by `by':"
	}
	
	if ("`absorb'"!="") {
		local absorb = "absorb(`absorb')"
	}
	
	`by' ebayes_by `fe' `coeff' if `touse', ///
		bee(`bee') theta(`theta') var(`var') uvar(`uvar') rawvar(`rawvar') ///
		w(`w') fevar(`fevar') tsum(`tsum') tol(`tol') maxiter(`maxiter') `absorb' `simple'
	
	if ("`generate'" != "") {
		qui gen `generate' = (1-`bee')*`fe' + `bee'*`theta'
	}

end

program ebayes_by, byable(recall)

	syntax varlist(numeric min=1) [if] [in], ///
		bee(varname) theta(varname) var(varname) uvar(varname) rawvar(varname) ///
		w(varname) fevar(varname) tsum(varname) tol(real) maxiter(integer) ///
		[ absorb(varname numeric) simple]

	marksample touse
	
	gettoken fe coeff: varlist
	
	local sigma_alpha
	
	* later we'll be estimating a beta
	* what regression command will we use for this?
	if ("`absorb'"=="") {
		local regcmd = "regress"
		local pstat = "xb"
	}
	else {
		local absorb = ", absorb(`absorb')"
		local regcmd = "areg"
		local pstat = "xbd"
	}

	****

	local sigma_alpha = 0
	local sigma_alpha_old = 0
	tempvar ttheta
	
	* iterate to find sigma_alpha and theta (x_ht*beta)
	local iter = 0
	* if we're NOT in simple mode...
	* -> we need at least 2 iterations
	* -> after the second iteration, stop if we've converged or hit MAXITER
	* if we ARE in simple mode
	* -> iterate once
	while ( ( "`simple'"=="" & ( (abs(`sigma_alpha' - `sigma_alpha_old') > `tol' & `iter'<`maxiter') | `iter'<2 ) ) | ///
			( "`simple'"=="simple" & `iter'<1) ) {

		* fix the new weights			
		if (`iter'==0) {
			* initial run: everyone gets the same weight
			qui replace `w' = 1 if `touse'
		}
		else {
			* ensuing runs: adjust weight for variance of fe estimate
			qui replace `w' = 1/(`sigma_alpha' + `fevar') if `touse'
		}
				
		* Run WLS to fix the new beta
		qui `regcmd' `fe' `coeff' [aw=`w'] if `touse' `absorb'
		
		* This is x_ht*beta
		capture drop `ttheta'
		qui predict `ttheta' if `touse', `pstat'
		
		* this is the variance of the error term
		* i.e. the weighted MSE of the regression
		local mse = e(rmse)^2
		
		// * this is the variance of x_ht*beta
		// local explvar = (e(tss) - e(rss))/e(df_r)
		// * this is wrong: the df adjustment is not right
		// * best to estimate var(xht*beta) by doing
		// * var(y) - var(e)
				
		* this is the variance of y (raw productivity)
		if ("`regcmd'"=="areg") {
			local yvar = e(tss)/(e(N)-1)
		}
		else {
			local yvar = (e(mss)+e(rss))/(e(N)-1)
		}
		* this is the weighted average of the measurement variance
		qui summ `fevar' [aw=`w'] if `touse'
		local avgfevar = r(mean)
		
		* this is the number of observations
		local N = e(N)
		* this is the degrees of freedom of the residual (N_obs - N_x - N_groups)
		local dof = e(df_r)
		
		* get the new sigma_alpha
		* if it's < 0, set it to 0 and force a break out of the loop
		if (`mse' - `avgfevar' < 0 ) {
			local sigma_alpha = 0
			local sigma_alpha_old = 0
			display "Iteration `iter': sigma_alpha^2 < 0 !"
		}
		else {
			local sigma_alpha_old = `sigma_alpha'
			local sigma_alpha = `mse' - `avgfevar'
			display "Iteration `iter': sigma_alpha^2: " string(`sigma_alpha',"%9.4g")
			
			local sigma_alpha_u = `yvar' - `avgfevar'
		}

		local iter = `iter' + 1
	}
	
	* indicator for failure to converge
	* condition: we're NOT in simple mode + we hit MAXITER + sigma_alpha isn't converged
	local failed = "`simple'"=="" & `iter'==`maxiter' & abs(`sigma_alpha' - `sigma_alpha_old') > `tol'

	* so long as we converged...
	if (!`failed' & `sigma_alpha' > 0) {
		qui replace `theta' = `ttheta' if `touse'
		qui replace `var' = `sigma_alpha' if `touse'
		qui replace `uvar' = `sigma_alpha_u' if `touse'
		qui replace `rawvar' = `yvar' if `touse'
		qui replace `bee' = ((`dof'-2)/(`dof'))*(`fevar'/(`fevar'+`sigma_alpha')) if `touse'
		display "Converged!"
		display "sigma_alpha^2: " string(`sigma_alpha',"%9.4g")
		display "sigma_alpha^2 [unconditional]: " string(`sigma_alpha_u',"%9.4g")
		display "raw var: " string(`yvar',"%9.4g")
		display "avg fe var: " string(`avgfevar',"%9.4g")
	}
	else {
		display "FAILED!"
		exit 198
	}
		
end
