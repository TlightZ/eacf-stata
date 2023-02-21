mata:
mata clear

// ACF calculation of certain lag
// REFER to https://www.real-statistics.com/time-series-analysis/stochastic-processes/autocorrelation-function/#:~:text=Definition%201%3A%20The%20autocorrelation%20function%20%28ACF%29%20at%20lag,0%20is%20the%20variance%20of%20the%20stochastic%20process.
real eacf_acf(real matrix z, real scalar lag) {
    real len
	len = cols(z)
	return( sum(z[1..len-lag]:*z[1+lag..len])/sum(z:*z) )
}

end

mata: mata mosave eacf_acf(),dir(PERSONAL) replace