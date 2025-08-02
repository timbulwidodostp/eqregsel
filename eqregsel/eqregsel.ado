* ado file for extremal quantile regression for selection models
* by D’Haultfoeuille, Maurel, Qiu and Zhang 
* This version: May 2019

quietly capture program drop eqregsel
program eqregsel,eclass 
	*version 12
	*version 12.0: set seed 123456789
	local boot_seed = 512
	if `c(stata_version)' >= 14{
		set rng kiss32
		set seed `boot_seed'
	}
	else{
		set seed `boot_seed'
	}
	syntax varlist(min=2 numeric) [if] [in] [, hom(integer 1) subs(integer 0) grid(integer 40) rep(integer 150) small ]
	marksample touse
	quietly count if `touse'
	if `r(N)' == 0 error 2000
	local N = r(N)
	
	* PARAMETERS
	tempfile _original_dt 
	qui save "`_original_dt'"
	qui keep if `touse'

	
	local d: word count `varlist'
	local num_indep = `d' - 1
	
	gettoken depvar varlist : varlist
	
	tempvar y_temp_  _id 
	gen `y_temp_' = -`depvar'
	gen `_id' = _n
	
	

	
	local G  `grid'
	local B  `rep'

	matrix l = (0.9,1.1)
	matrix l1 = (1,0.2)
	local J = colsof(l)
	local J1 = colsof(l1)
	
	if `G' <= 0 {
		di in red "grid() does not accept negative value."
		exit
	}
	
	if `B' <= 0 {
		di in red "rep() does not accept negative value."
		exit
	}
	
	if `subs' < 0{
		di in red "subs() does not accept negative value."
		exit
	}
	
	if `hom' < 0{
		di in red "hom() does not accept negative value."
		exit
	}
	else if `hom' >= `d'{
		di in red "hom() does not accept a value larger than the total number of independent variables."
		exit
	}
	

	if "`small'" != ""{
		tokenize `varlist'
		foreach var of numlist 1/`num_indep'{
			qui replace ``var'' = -``var''
			local tempvarlist `tempvarlist' ``var''
		}
	}
	sort `depvar' `varlist'

	* SPECIFY BOOTS
	if `subs' == 0{
		if `N'<=500                       local subs = 0.6*`N'
		else if (`N'>500)&(`N'<=1000)     local subs = 300 + 0.4*(`N'-500)
		else if (`N'>1000)&(`N'<=2000) 	  local subs = 500 + 0.2*(`N'-1000)
		else                              local subs = 700 + 0.2*log(2000)/log(`N')*(`N'-2000)
		local subs = floor(`subs')
	}
	
	local lower = min(80/`subs',0.1)
	local upper = 0.3
	local step = (`upper' - `lower')/`G'
	mat Phi = J(1,`d'-1,0)
	forvalues i = 1/`hom'{
		mat Phi[1,`i'] = 1
	}
	
	local phi Phi
	mat ss = J(1, colsof(`phi'),1) * (`phi')'
	local dbb = ss[1,1]

	mat par_bootb = J(`dbb',`B',0)
	mat chi_bootbb = J(`B',1,0)
	scalar disbb = 10000000000000000000     
	
	
	
	* SELECT OPTIMAL TAU
	forvalues gg = 1/`G'{
		local tau = `lower' + `step' * `gg'
		qui qreg_simplified12 `y_temp_' `varlist' if `touse',quantile(`tau') cformat(%10.0g)
		mat theta = e(b)'
		mat thetatemp = J(`d',`B',0)
		
		* BOOSTRAP THE FIRST STAGE ESTIMATOR TO COMPUTE OMEGA_0
		local tloop1 "$S_TIME"
		forvalues bb = 1/`B'{
			preserve
			qui bsample if `touse'
			qui qreg_simplified12 `y_temp_' `varlist',quantile(`tau') cformat(%10.0g)
			mat thetatemp[1,`bb'] = e(b)'
			restore	
		}
		
		* Adjust the order of coefficients of independent variables
		mat thetahat = (theta[`d',1]\theta[1..`d'-1,1])
		mat thetadagger = (thetatemp[`d',1..`B']\thetatemp[1..`d'-1,1..`B'])
		mat Sigma=(thetadagger-thetahat*J(1,`B',1))*(thetadagger-thetahat*J(1,`B',1))'/`B'
		
		* CALL MYFUN_HOM_NEW()
		* Instead of selecting optimal tau first, we compute beta each time
			
		mat pf = J(`d',`J1',0)
		forvalues j = 1/`J1'{
			local tau1 = `tau'*l1[1,`j']
			qui qreg_simplified12 `y_temp_' `varlist' if `touse',quantile(`tau1') cformat(%10.0g)
			mat pf[1,`j'] = e(b)'		
		} 
		
		local sigma Sigma
		local ll1 l1
		local ll l
		* The order of pf is adjusted in myfun_hom_new()
	    local PF  pf
		mata: myfun_hom_new( "`phi'","`sigma'", "`ll1'","`PF'", `d')
		mat beta = e(beta)
		local dis_b = e(dis_b)
		
			
		* Use subsample to compute the point estimator and J-test statistic
		forvalues bb = 1/`B'{
			preserve
			qui keep if `touse'
			qui sample `subs',count
			
			* CALL MYFUN_HOM_NEW()
			
			mat pf = J(`d',`J',0)
			forvalues j = 1/`J'{
				local tau1 = `tau'*l[1,`j']
				qui qreg_simplified12 `y_temp_' `varlist',quantile(`tau1') cformat(%10.0g)
				mat pf[1,`j'] = e(b)'		
			} 
			
			local PF  pf
			mata: myfun_hom_new("`phi'","`sigma'", "`ll'","`PF'", `d')
			mat par_bootb[1,`bb'] = e(beta)
			mat chi_bootbb[`bb',1] = e(dis_b)
			restore
		}
		
		local tloop2 "$S_TIME"
		
		* Selecting the optimal tau
		local Chi_bootbb  chi_bootbb
		local Par_bootb  par_bootb
		mata: IC("`Chi_bootbb'", "`Par_bootb'",`J',`dbb',`B',`subs',`N',`tau')
		local disbbtemp = e(disbbtemp)

		if(`disbbtemp' < disbb){
			scalar disbb = `disbbtemp'
			mat Sigmahat = Sigma
			mat beta_hom = beta
			local chibb = `dis_b'
			local tau0 = `tau'		
		}
		
		* Estimate computation time
		if(`gg'==1){
			local hdiff = real(substr("`tloop2'",1,2)) - real(substr("`tloop1'",1,2))
			local mdiff = real(substr("`tloop2'",4,2)) - real(substr("`tloop1'",4,2))
			local sdiff = real(substr("`tloop2'",7,2)) - real(substr("`tloop1'",7,2))
			
			local esttime = (`hdiff'*60 + `mdiff' + `sdiff' /60)* `G'
			local esttime = string(`esttime')
			di 
			di in ye "The estimation will take about " in wh "`esttime'" _c
			di in ye " minutes." 
		}
		
		* Printing progress bar if G >= 10
		if(`G'>=10 & `gg'==1){
			
			* Display the bar
			forvalues ss = 1/`G'{
				tempname per
				local per = (`ss'-1)/`G' 
				if(`per'==0 | `per'==0.2 | `per'==0.4 | `per'==0.6| `per'==0.8){ 
					di "|-" _c 
				}
				else if(`ss'==`G'){ 
					di "-|" 
				}
				else{ 
					di "--" _c 
				}
			}
			
			* Display progress
			forvalues ss = 1/`G'{
				tempname per
				local per = (`ss'-1)/`G' 
				if(`per'==0){   
					di "0" _c 
				}
				else if(`per'==0.2){  
					di "20" _c  
				}
				else if(`per'==0.4) {  
					di "40" _c   
				}
				else if(`per'==0.6) {  
					di "60" _c   
				}
				else if(`per'==0.8) {   
					di "80" _c  
				}
				else if(`ss'==`G'){  
					di "100"      
				}
				else{  
					di "  " _c  
				}
			}
			
			di ". " _c
		}
		else if(`G'>=10 & `gg'>1 & `gg'<=`G'){
			di ". " _c
		}
		else if(`G'>=10){
			di "."
		}
 	}
	
	
	ereturn post
	ereturn clear
		
	* BOOTSTRAP THE CONFIDENCE INTERVAL
	
	local sigmahat Sigmahat
	mat par_bootstraphomb = J(`dbb',`B',0)
	forvalues bb = 1/`B'{
		preserve
		qui bsample if `touse'	
			
		* CALL MYFUN_HOM_NEW()
			
		mat pf = J(`d',`J1',0)
		forvalues j = 1/`J1'{
			local tau1 = `tau0'*l1[1,`j']
			qui qreg_simplified12 `y_temp_' `varlist',quantile(`tau1') cformat(%10.0g)
			mat pf[1,`j'] = e(b)'		
		} 
			
		local PF  pf
		mata: myfun_hom_new("`phi'","`sigmahat'", "`ll1'","`PF'", `d')
		mat par_bootstraphomb[1,`bb'] = e(beta)
			
		restore
	}
	
	local Par_bootstraphomb par_bootstraphomb
	
	mata: stat("`Par_bootstraphomb'", `J1', `dbb',`chibb')
	local specificationtest e(specificationtest)
	mat std_b = e(std_b)
	mat V = diag(std_b)*diag(std_b)
	
	
	* RETURNS IN ECLASS
	
	ereturn scalar tau0 = `tau0'
	ereturn scalar specificationtest =`specificationtest'
	ereturn scalar subs = `subs'
	ereturn scalar homvar = `hom'
	
	
	* DISPLAY
	di ""
	di in gr "Number of observations = " %10.0g `N'
	di in gr "Optimal quantile index = " %10.0g e(tau0)
	//di in gr "Bootstrapped standard deviation = " %10.0g e(std_b)
	di in gr "J test(p-value) = " %10.0g e(specificationtest)
	di in gr "Subsampling size used in bootstrapping = " %10.0g e(subs)
	di in gr "Number of variables of interest = " %10.0g e(homvar)
	di ""
	* Display the results in a table

	mat Phi = `phi'
	tokenize `varlist'
	local i=1
	mat colnames beta_hom = "`depvar'"
	while "`1'" != "" {
		scalar index =  Phi[1,`i']
		if (index == 1){
			local names `names' `1'
		}
		mac shift
		local i=`i'+1
	}
	
	mat rownames V = `names'
	mat colnames V = `names'
	mat rownames beta_hom = `names'
	mat b = beta_hom'
	
	_coef_table , bmatrix(b) vmatrix(V)
	
	* RETURNS IN ECLASS
	* MATRIX
	ereturn matrix beta_hom = beta_hom
	ereturn matrix std_b = std_b
	ereturn matrix v = V

	use "`_original_dt'", clear
	
end

//-----------------------------------------------------------------------
//qreg_simplified12, Simplified version of qreg command
program qreg_simplified12, eclass byable(recall) sort prop(sw mi)
	version 6, missing
	local options "Level(cilevel)"
	if !replay() {
		local cmdline : copy local 0
		syntax varlist [aw fw] [if] [in] [, `options' /*
			*/ noConstant Quantile(real 0.5) /*
			*/ WLSiter(integer 1) noLOg * ]
		_get_diopts diopts options , `options'
		if "`constan'"!="" {
			di in red "nocons invalid"
			exit 198
		}
		if `quantil' >= 1 {
			local quant = `quantil'/100
		}
		else	local quant "`quantil'"
		if `quant' <= 0 | `quant' >= 1 {
			di in red "quantiles(`quantil') out of range"
			exit 198
		}
		if `wlsiter'<1 { error 198 }
		marksample touse
		qui count if `touse'
		if r(N)<2 { error cond(_n,2001,2000) }
		gettoken dep indep : varlist
		_rmcoll `indep' [`weight'`exp'] if `touse'
		local varlist `dep' `r(varlist)'
		tempvar r s2 p
		gen long `s2' = _n

		/* initial estimates via weighted least squares 	*/
		_qregwls `varlist' [`weight' `exp'] if `touse',	///
			iterate(`wlsiter') quant(`quant') r(`r') `log'

		quietly {
			if "`log'"=="" { local log "noisily" }
			else 	local log

			sort `r' `s2'
			drop `r'

			`log' _qreg `varlist' if `touse' [`weight'`exp'], /*
				*/ quant(`quant') `options'
			
end

//------------------------------------------------------------------------------
//Mata part: main function
mata:
// version 13
mata clear
mata set matastrict on

//------------------------------------------------------------------------------
//myfun_hom_new.m
void myfun_hom_new( string matrix Phi, string matrix sigma,string matrix ll, ///
					  string matrix PF, real scalar dd)
{

/*******************************************************************************\
*** Input                                                                       *
* myfun_hom_new computes the point estimator and the J-test statistic for a     *
  given quantile index tau.                                                     *
* tau: Quantile                                                                 *
* phi: User supplied vector to indicate which variable (in our application, it  *
  is "Black") is assumed to be homoskedastic. (The code allows for more than one*
  covariates to be homoskedastic.)                                              *
* X: Dependent variables                                                        *
* Y: Independent variable                                                       *
* l: Equations we will explore by minimum distance estimations are indexed by l.*
 In general, we will use quantile level tau, tau*l to estimate beta and delta   *
* Sigma: An estimator for \omega_0 in the paper                                 *
*** Output                                                                      *
* beta: The point estimator                                                     *
* dis_b: J-test statistic for this specification.                               *
\*******************************************************************************/

	real scalar  JJ, dbb, i, j
	real matrix  l, phi, Sigma, pf
	l = st_matrix(ll)
	JJ = cols(l)
	phi = st_matrix(Phi)
	dbb = sum(phi)
	Sigma = st_matrix(sigma)
	pf = st_matrix(PF)

    //Convert vector phi into a matrix which picks out the covariates that are 
	//homoskedastic.
	real scalar Count, counta
	real matrix phitemp
	phitemp = J(dbb,dd-1,0)
	Count = 1
	for(i=1;i<=dbb;i++){
		counta = 1
		for(j=Count;j<=dd-1;j++){
			if (phi[j] == 1  && counta == 1) {
				phitemp[i,j] = 1
				counta = counta + 1 
				Count = j + 1
			}
		}
	}
	
	
	//Compute matrix L, Gamma_2, Gamma_3 in the paper.
	real matrix l1, L, Gamma2, Gamma3
	l1 = l
	L = J(JJ,JJ,0)
	for(i = 1; i <= JJ; i++){
		for(j = 1; j <= JJ; j++){
			L[i,j] = min((l1[i],l1[j]))/sqrt(l1[i] * l1[j])
		}
	}
	
	Gamma2 = (J(dd-1,1,0),I(dd-1))
	Gamma3 = diag(1:/sqrt(l1))

	// The first step: extremal quantile regression
	real matrix tempy2
	tempy2 = J(dbb*JJ,1,0)
		
	for(j=1;j<=JJ;j++){
		// Be careful aboout the subscript of pf
		tempy2[dbb*(j-1)+1..dbb*j] = phitemp * pf[|1,j\dd-1,j|]
	}
	// The second step: minimum distance estimation 
	real matrix omega_0, W2, mom2,  GpG, beta
	real scalar dis_b
	omega_0 = Sigma
	GpG = Gamma3 # (phitemp * Gamma2)

	W2 = luinv(GpG * (L # omega_0) * GpG' )                                     //W2 is the optimal weighting matrix for homo beta

	beta = - luinv(J(JJ,1,I(dbb))' * W2 * J(JJ,1,I(dbb))) * J(JJ,1,I(dbb))' * W2 * tempy2
	mom2 = J(JJ*dbb, 1, 0)
	for(j = 1; j <= JJ; ++j){
		// Be careful aboout the subscript of pf
		mom2[dbb*(j-1)+1..dbb*j,1] = phitemp * pf[|1,j\dd-1,j|] + beta             //mom2 is the value of moments for homo beta
	}
	
	dis_b =  mom2' * W2 * mom2                                                  //here we compute the distance of homo beta evaluated at extremal quantile estimator of delta
	
	st_matrix("e(beta)",beta)
	st_numscalar("e(dis_b)",dis_b)
	
}

//------------------------------------------------------------------------------
//
void IC(string matrix Chi_bootbb, string matrix Par_bootb, real scalar J, ///
			real scalar dbb, real scalar B, real scalar boots, real scalar N, ///
			real scalar tau){
	real scalar median_b1, disbbbias, disbbvar, disbbtemp
	real matrix aux, aux1, tempchi, tempmean, chi_bootbb, par_bootb
	chi_bootbb = st_matrix(Chi_bootbb)
	par_bootb = st_matrix(Par_bootb)
	
	aux = rchi2(100000,1,(J-1)*dbb)
	aux1 = select(aux,aux:!=.)
	median_b1 = mm_median(aux1)
	//median_b1 = mm_median(rchi2(100000,1,(J-1)*dbb))
	//disbbbias=abs(mm_median(select(chi_bootbb,chi_bootbb:=.))*boots/N-median_b1)'
	tempchi = select(chi_bootbb,chi_bootbb:!=.)
	disbbbias=abs(mm_median(tempchi)*boots/N-median_b1)'
	tempmean = J(1,B,mean(par_bootb')')
	disbbvar = mean(colsum((par_bootb-tempmean):^2)')
	disbbtemp = disbbvar*boots/N + disbbbias/sqrt(tau*boots)
	st_numscalar("e(disbbtemp)",disbbtemp)
}

void stat(string matrix Par_bootstraphomb, real scalar J, real scalar dbb, ///
		real scalar chibb	){
	real scalar  specificationtest
	real matrix std_b,par_bootstraphomb
	par_bootstraphomb = st_matrix(Par_bootstraphomb)
	std_b = mm_colvar(par_bootstraphomb')'                                      //Compute the standard deviation.
	std_b = sqrt(std_b)
	specificationtest = 1 - chi2((J-1)*dbb,chibb)
	st_matrix("e(std_b)",std_b)
	st_numscalar("e(specificationtest)",specificationtest)
}

 
end
