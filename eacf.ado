*! version 1.0.0
cap program drop eacf
program define eacf, rclass
	version 16
	syntax varname(ts) [if] [in] [,ar(integer 7) ma(integer 13)]

	preserve
	marksample touse
	_ts tvar panelvar `if' `in', sort onepanel
	markout `touse' `tvar'
	quietly keep if `touse'
	
	tempvar zVar zVarMean
	// Demean time-series
	// DEMEAN is required indirectly to get the same result as that in the paper
	egen double `zVarMean' = mean(`varlist')
	gen double `zVar' = `varlist' - `zVarMean'

	local nrow = `ar' + `ma' + 1
	// NOT ar , because of the nature of using EQ.<2.7> in the paper
	// The reason is explained below

	// Calculate AR(k) coefficients Whose MA = 0
	matrix coefficients = J(`nrow',`nrow',0)
	forval ar_order = 1/`nrow' {
		quietly reg `zVar' L(1/`ar_order').`zVar', nocons
		matrix matret = e(b)
		matrix coefficients[`ar_order',1] = matret[1,1..`ar_order']
	}
	
	// Turn data to matrix, for MATA usage
	mkmat `zVar',matrix(z_matrix)
	restore
	// MATA calculation of EACF
	mata: eacf_calculate("z_matrix","coefficients",`ar',`ma')
	
	// Return key results
	return matrix seacf = seacf
	return matrix symbol = symbol
end

version 16
mata:

// ACF calculation of certain lag
// REFER to https://www.real-statistics.com/time-series-analysis/stochastic-processes/autocorrelation-function/#:~:text=Definition%201%3A%20The%20autocorrelation%20function%20%28ACF%29%20at%20lag,0%20is%20the%20variance%20of%20the%20stochastic%20process.
real eacf_acf(real matrix z, real scalar lag) {
    real len
	len = cols(z)
	return( sum(z[1..len-lag]:*z[1+lag..len])/sum(z:*z) )
}

// EACF Caluculator
void eacf_calculate(string scalar varname,string scalar coeffname,nar,nma) {
    real scalar len,nrow,ncol
	real matrix z,coefficients,zlag,seacf,symbol,n_minus_k,epsilon
	
    z = st_matrix(varname)'
	len = cols(z)
	coefficients = st_matrix(coeffname)
	nrow = rows(coefficients)
	coefficients = J(nrow,1,-1),coefficients
	ncol = cols(coefficients)
	
	// make lags for matrix acceleration
	zlag = J(nar+1,len,.)
	for (k=1;k<=nar+1;k++) {
		zlag[k,k..len] =  z[1..len-k+1]
	}
	
	seacf = J(nar + 1,nma + 1,0)
	// first row implies AR = 0
	symbol = J(nar + 1,nma + 1,.)
	n_minus_k = len :- range(0,nar,1)

	for (j =1;j<=nma+1;j++) {
	    // calculate eacf column by column
		// consistent to <symbol j> in the paper
		for (k=1;k<=nrow-1;k++) {
		    // consistent to <symbol k> in the paper
			// coefficient calculation REFER TO equation <2.7> in the paper
			coefficients[k,2..ncol] = coefficients[k+1,2..ncol] - coefficients[k,1..ncol-1] * coefficients[k+1,(k+1)+1] / coefficients[k,(k)+1]
		}
		seacf[1,j] = eacf_acf(z,j)
		// EACF of AR(0) equals ACF(z)
		for (k=1;k<=nar;k++) {
		    // calculate eacf row by row
			epsilon = coefficients[k,1..k+1] * zlag[1..k+1,k+1..len]
			// REFER TO equation <3.4> in the paper
			seacf[k+1,j] = eacf_acf(epsilon,j)
		}
		
		// ABS(eacf) exceeding 2 times of standard variance is marked 'x'
		// Variance is (n - k - j)^(-1), thus standard variance is (n - k - j)^(-0.5)
		symbol[1...,j] = abs(seacf[1...,j]) :> 2:/sqrt(n_minus_k :- j)

		nrow = nrow - 1
	    // the nature that coefficient calculation through EQ.<2.7> requires <k+1> information to calculate <k>'s coefficient
		// determines that each iteration(j) fails to calculate the last row's coefficient
		// so number of rows should minus 1 in each iteration(j)
	}
	
	// Printer
	//	Print Sample EACF
	printf("  AR\MA")
	for (j=1;j<=nma+1;j++) {
	    printf("%7.0f",j-1)
	}
	printf("\n")
	for (k=1;k<=nar+1;k++) {
	    printf("%7.0f",k-1)
		for (j=1;j<=nma+1;j++) {
		    printf("%7.2f",seacf[k,j])
		}
		printf("\n")
	}
	
	//	Print Sample EACF	
	printf("\n")
	printf(" AR\MA")
	for (j=1;j<=nma+1;j++) {
	    printf("%4.0f",j-1)
	}
	printf("\n")
	for (k=1;k<=nar+1;k++) {
	    printf("%6.0f",k-1)
		for (j=1;j<=nma+1;j++) {
		    if (symbol[k,j]) {
			    printf("%4s","x")
			} else {
			    printf("%4s","o")
			}
		}
		printf("\n")
	}
	
	// Return results to Main
	st_matrix("seacf",seacf)
	st_matrix("symbol",symbol)
}

end