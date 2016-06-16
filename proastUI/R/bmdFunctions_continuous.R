#' Main function for calculations with continuous response type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all 
#' @export
f.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.con")
	
	
	ans.all <- f.start.con(ans.all = ans.all, fitted = FALSE)
	
	
	for(main.ans.single in c(ans.all$main.ans, 13)){
		
		switch(as.character(main.ans.single), 
				
				"4" = { #Fit model
					
					ans.all <- f.mm4.con(ans.all, track = track)
					
				},
				
				"6" = { # Calculate CED
					
					if (!ans.all$fitted){
						
						stop("First fit the model")
						
					}
					
					if(ans.all$model.ans %in% c(1, 11)){
						
						stop("BMD not defined for null or full model")
						
					}
					
					ans.all <- f.mm6.con(ans.all, track = track)
					
				},
				
				"13" = { # Return results
					
#					main.ans <- 13
					if (track) 
						print("f.con:  END")
					return(ans.all)
					
				}) 
		
	}
	
}

#' Define initial parameter values for the chosen model, continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param fitted boolean TRUE if the model has already been fitted;
#'  default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.start.con <- function (ans.all, fitted = FALSE, track = FALSE) {
	
	if (track) 
		print("f.start.con")
	
	cc.upp <- 10000
	dd.upp <- 10
	
	ans.all$fitted <- fitted
	
	with(ans.all, {
				
				if (dtype == 5) {
					
					dtype <- 1
					
				}
				
				y.sort <- yNew[order(x)]
				x.sort <- sort(x)
				
				if (dtype %in% c(10, 15, 250, 260) & length(x) < 10) {
					# Continuous summary data
					
					xtmp <- x.sort
					ytmp <- y.sort
					
				} else {
					# Make 5 groups for linear model fit
					
					sub <- round(length(x.sort)/5)
					if (sub == 0) 
						sub <- 1
					xtmp <- numeric()
					ytmp <- numeric()
					e1 <- 0
					e2 <- 0
					for (k in 1:5) {
						e2 <- e2 + sub
						qq <- x.sort[e1:e2]
						xtmp[k] <- mean(qq, na.rm = TRUE)
						qq <- y.sort[e1:e2]
						if (dtype %in% c(25, 250)) 
							ytmp[k] <- mean(qq, na.rm = TRUE)
						if (dtype %in% c(26, 260)) 
							ytmp[k] <- mean(sqrt(qq, na.rm = TRUE))^2
						else if (!(dtype %in% c(25, 26, 250, 260))) 
							ytmp[k] <- exp(mean(log(qq[is.na(qq) == 
																	FALSE])))
						e1 <- e1 + sub
					}
				}
				if (length(ytmp[!is.na(ytmp)]) == 3) {
					ytmp[4:5] <- ytmp[3]
					xtmp[4:5] <- xtmp[3]
				}
				if (length(ytmp[!is.na(ytmp)]) == 4) {
					ytmp[5] <- ytmp[4]
					xtmp[5] <- xtmp[4]
				}
				top <- length(ytmp)
				intercept <- NA
				slope <- NA
				cc <- NA
				dd <- NA
				
				if (model.ans != 1) {
					
					lm.res <- f.start.lm.con(xtmp, ytmp, model.ans, dtype, track = track)
					
					intercept <- abs(lm.res$intercept)
					if (model.ans %in% c(2, 3, 7, 8, 17, 18)){
						
						slope <- lm.res$slope
						
					} else {
						
						slope <- abs(lm.res$slope)
						
					}
					
					increase <- lm.res$increase
					
				} else {
					
					increase <- 1
					
				}
				
				CES <- CES*increase
				
				switch(as.character(model.ans), 
						
						"1" = { #Null model
							
							intercept <- mean(ytmp)
							regr.par <- rep(intercept, nr.aa)
							
						}, 
						
						"11" = { #Full model
							
							if (dtype %in% c(10, 15, 250, 260)) {
								
								if (dtype %in% c(10, 15)) regr.par <- exp(mn.log)
								if (dtype == 250) regr.par <- (mn.log)
								if (dtype == 260) regr.par <- (mn.log)^2
								
							} else {
								
								regr.par <- numeric()
								nn <- numeric()
								if (nr.var == 1) { 
									
									for (jj in levels(factor(fct2))) for (ii in levels(factor(fct1))) {
											x.tmp <- x[fct1 == ii & fct2 == jj]
											y.tmp <- yNew[fct1 == ii & fct2 == jj]
											regr.tmp <- as.numeric(exp(tapply(log(y.tmp), 
																	x.tmp, mean)))
											regr.par <- c(regr.par, regr.tmp)
											nn.tmp <- tapply(y.tmp, x.tmp, length)
											nn <- c(nn, nn.tmp)
										}
									
								} else {
									
									regr.par <- numeric()
									for (jj in levels(factor(fct3))) {
										x.tmp <- x[fct3 == jj]
										y.tmp <- yNew[fct3 == jj]
										regr.tmp <- as.numeric(exp(tapply(log(y.tmp), 
																x.tmp, mean)))
										regr.par <- c(regr.par, regr.tmp)
										nn.tmp <- tapply(y.tmp, x.tmp, length)
										nn <- c(nn, nn.tmp)
									}
								}
								ans.all$nn <- nn
							}
							
							ans.all$lower <- regr.par
							ans.all$upper <- regr.par
							
						}, 
						
						"13" = { #E3 - CED
							
							if (fitted) {
								
								intercept <- par.start[nr.var + 1]
								slope <- par.start[nr.var + nr.aa + 1]
								
							}
							
							if (!fitted || is.na(slope) || slope == 0) {
								
								if (is.na(slope) || slope == 0) slope <- mean(xtmp)
								CED <- (1/slope) * log(CES + 1)
								CED <- abs(CED)
								if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
												min(x))/3
							}
							dd <- 1
							regr.par <- c(rep(intercept, nr.aa), rep(CED, nr.bb), 
									rep(dd, nr.dd))
							
						}, 
						
						"15" = { #E5 - CED
							
							if (fitted) {
								
								intercept <- par.start[nr.var + 1]
								slope <- par.start[nr.var + nr.aa + 1]
								
							}
							if (!fitted || is.na(slope) || slope == 0) {
								
								dd <- 1
								cc <- ytmp[top]/intercept
								if ((CES < 0) & (cc > 1 + CES)) {
									cc <- 1 + CES - (1 + CES)/100
								}
								if ((CES > 0) & (cc < 1 + CES)) {
									cc <- 1 + CES + CES/100
								}
								if (is.na(slope) || slope == 0) slope <- mean(xtmp)
								if (cc == 0) cc <- 0.1
								if (cc == 1) cc <- 0.9
								CED <- (-(1/slope) * log((-CES + 1 - cc)/(1 - 
																cc)))^(1/dd)
								CED <- abs(CED)
								if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
												min(x))/3
							}
							regr.par <- c(rep(intercept, nr.aa), rep(CED, nr.bb), 
									rep(cc, nr.cc), rep(dd, nr.dd))
						}, 
						
						"23" = { #H3 - CED
							if (CES > 0) slope <- -slope
							dd <- 1
							CED <- f.inv.con(model.ans = 18, c(intercept, slope, dd), track = track)
							regr.par <- c(rep(intercept, nr.aa), rep(CED, nr.bb), 
									rep(dd, nr.dd))
							y.expect <- f.expect.con(model.ans, xtmp, regr.par, 
									CES = CES, increase = increase)
							test <- abs(sum(sign(y.expect)))
							if (!any(is.na(test))) if (test != length(y.expect)) {
									slope <- 1.5 * slope
									regr.par <- c(rep(intercept, nr.aa), rep(slope, nr.bb), 
											rep(dd, nr.dd))
								}
							
						}, 
						
						"25" = { # H5 - CED
							
							cc <- ytmp[top]/intercept
							dd <- 1
							if (CES < 0 & cc > 1 + CES) {
								cc <- 1 + CES - (1 + CES)/100
							}
							if (CES > 0 & cc < 1 + CES) {
								cc <- 1 + CES + CES/100
							}
							if (is.na(slope) || slope == 0) slope <- mean(xtmp)
							if (cc == 0) cc <- 0.1
							if (cc == 1) cc <- 0.9
							CED <- f.inv.con(model.ans = 20, c(intercept, slope, cc, dd), 
									track = track)
							if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
											min(x))/3
							regr.par <- c(rep(intercept, nr.aa), rep(CED, nr.bb), 
									rep(cc, nr.cc), rep(dd, nr.dd))
						})
				
				ans.all$regr.par <- regr.par
				
				if (!(dtype %in% c(10, 15, 250, 260))) {
					
					expect.trans <- f.expect.con(model.ans, x, regr.par, 
							fct1 = fct1, fct2 = fct2, fct3 = fct3, CES = CES, 
							twice = twice, ttt = 0, y = yNew, increase = increase, 
							x.mn = x.mn, par.start = par.start)
					
					yy.trans <- yNew
					
					if (dtype == 26) {
						
						expect.trans <- sqrt(expect.trans)
						yy.trans <- sqrt(yNew)
						
					} else if (dtype %in% c(1, 5, 15)) {
						
						expect.trans <- log(expect.trans)
						yy.trans <- log(yNew)
						
					}
					
					if (!model.ans %in% c(1:10, 12:20, 22:25)) {
						
						ans.all$lower <- c(1e-06, ans.all$lower)
						ans.all$upper <- c(1, ans.all$upper)
						
					}
					
					if (any(is.na(expect.trans))) {
						
						ans.all$par.start <- c(rep(NA, nr.var), regr.par)
						warning("Adjust start values")
						ans.all$adjust.start <- TRUE
						
					}
					
					var.start <- var(yy.trans - expect.trans)
					
				} else {
					
					var.start <- mean(sd2.log)
					
				}
				
				ans.all$par.start <- c(rep(var.start, nr.var), regr.par)
				
				loglik.first <- -f.lik.con(ans.all$par.start, 
						x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
						sd2.log, nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
						ttt = 0, fct4 = fct4, fct5 = fct5, 
						cens.up = cens.up, par.tmp = NA, increase = increase, 
						x.mn = x.mn)
				
				if (is.na(loglik.first) | loglik.first < -1e+05) {
					
					warning("Adjust start values before fitting the model")
					ans.all$adjust.start <- TRUE
					
				}
				
				ans.all$loglik.first <- loglik.first
				
				
				ans.all$nr.var <- nr.var
				ans.all$npar <- length(ans.all$par.start)
				ans.all$CED <- NA
				ans.all$increase <- increase
				ans.all$CES <- CES
				ans.all$max.lev <- nr.aa * nr.bb
				
				ans.all <- f.constr.con(ans.all, track = track)
				
				if (track) 
					print("f.start.con:  END")
				
				return(ans.all)
				
			})
}



#' Calculate likelihood for continuous response
#' @param theta numeric vector, the initial regression parameter values
#' @param x numeric vector, the dose values
#' @param y numeric vector, the response values 
#' @param dtype integer, determines the type of response
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param fct3 numeric, value for parameter var
#' @param model.ans integer, determines the model that will be fitted 
#' @param mn.log numeric vector, transformation of the response values,
#' see f.execute()
#' @param sd2.log  numeric vector, transformation of the sd of the response
#' values, see f.execute()
#' @param nn numeric vector, the number of responses per dose level, for 
#' continuous summary data
#' @param Vdetlim numeric vector, values of detection limit
#' @param CES numeric, value for the CES
#' @param twice boolean, if some parameter values are equal, see f.execute()
#' @param ttt numeric, time variable 
#' @param fct4 numeric, value for parameter c
#' @param fct5 numeric, value for parameter d
#' @param cens.up numeric, value for right censoring
#' @param lb numeric vector, determines the lower bound for theta;
#' default value is -Inf
#' @param ub numeric vector, determines the upper bound for theta;
#' default value is Inf
#' @param par.tmp numeric vector, regression parameter values, see f.pars() 
#' @param increase boolean, whether the response values are increasing or 
#' decreasing for increasing dose values 
#' @param x.mn numeric value, the mean of dose values, see f.execute()
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric value, minus the sum of the scores (total likelihood)
#' @export
f.lik.con <- function (theta, x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
		sd2.log, nn, Vdetlim, CES, twice = TRUE, ttt = 0,  
		fct4 = 1, fct5 = 1, cens.up = NA, lb = -Inf, ub = Inf, par.tmp, 
		increase = increase, x.mn = NA, track = FALSE) {
	
	if (track) {
		
		print("f.lik.con:  begin")
		cat("Initial parameter values:", theta, "\n")
		
	} 
	
	
	if (any(is.na(theta))) {
		
		warning("Problem in f.lik.con: NAs in theta")
		return(NA)
		
	}
	
	if ((length(fct3) > 1) & (length(fct3) != length(x))) {
		
		stop("fct3 incorrect")
		
	}
	
	variance <- 0
	
	for (jj in 1:max(fct3)) variance <- variance + theta[jj] * (fct3 == jj)
	
	regr.par <- theta[(max(fct3) + 1):length(theta)]
	
	
	if (any(!is.finite(theta))) {
		
		theta[!is.finite(theta)] <- par.tmp[!is.finite(theta)]
		
	}
	
	if (sum(theta <= lb) > 0) {
		theta[theta < lb] <- 1.1 * par.tmp[theta < lb]
	}
	
	if (sum(theta >= ub) > 0) {
		theta[theta > ub] <- 0.9 * par.tmp[theta > ub]
	}
	
	expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
			fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, CES = CES, 
			twice = twice, ttt = ttt, y = y, increase = increase, 
			x.mn = x.mn)
	
	if (any(is.na(expect))) {
		
		warning("NAs in predicted response at parameter values", signif(theta, 4))
		
	}
	
	
	if (dtype %in% c(1, 5, 25, 26)) { 
		# Continuous response 
		
		if(dtype %in% c(1, 5)){
			
			expect <- log(expect)
			yTransformed <- 0
			yTransformed[y > 0] <- log(y)
			
			VdetlimTransformed <- log(Vdetlim)
			cens.upTransformed <- log(cens.up)
			
		} else if (dtype == 26) {
			
			expect <- sqrt(expect)
			yTransformed <- sqrt(y)
			
			VdetlimTransformed <- sqrt(Vdetlim)
			cens.upTransformed <- sqrt(cens.up)
			
		} else {
			
			yTransformed <-  y
			
			VdetlimTransformed <- Vdetlim
			cens.upTransformed <- cens.up
			
		}
		
		score1 <- (y > 0) * (-0.5 * log(2 * pi * variance) - 
					((yTransformed - expect)^2)/(2 * variance))
		score.detlim <- log(pnorm((VdetlimTransformed - expect)/sqrt(variance)))
		score.detlim[!is.finite(score.detlim)] <- 0
		score.censup <- log(1 - pnorm((cens.upTransformed - expect)/sqrt(variance)))
		score.censup[!is.finite(score.censup)] <- 0
		score2 <- (y == -1000) * score.detlim + (y == -2000) * 
				score.censup
		score <- score1 + score2
		
		
	} else if (dtype %in% c(10, 15, 250, 260)) {
		# Continuous summary response
		
		if(dtype %in% c(10, 15)){
			
			expect <- log(expect)
			
		} else if (dtype == 260){
			
			expect <- sqrt(expect)
			
		}
		
		if (model.ans != 11) {
			
			dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
			
		} else {
			
			dum <- (nn - 1) * sd2.log
			
		}
		
		score <- -(nn * log(sqrt(2 * pi * variance)) + dum/(2 * variance))
		
	} 
	
	if(track)
		print("f.lik.con END")
	
	return(-sum(score))
	
	
}



#' Calculate expected response values, for continuous response
#' @param model.ans integer, determines the type of model to be fitted
#' @param x numeric vector, the dose values
#' @param regr.par numeric vector, regression parameter values
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param fct3 numeric, value for parameter var
#' @param fct4 numeric, value for parameter c
#' @param fct5 numeric, value for parameter d
#' @param CES numeric, value for the CES
#' @param twice boolean, if some parameter values are equal, see f.execute()
#' @param ttt numeric, time variable
#' @param y numeric vector, the response values
#' @param increase boolean, whether the response values are increasing or 
#' decreasing for increasing dose values 
#' @param x.mn numeric value, the mean of dose values, see f.execute()
#' @param par.start numeric vector, regression parameter values (e.g. MLEs)
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, the expected response values under the estimated 
#' regression model
#' @export
f.expect.con <- function (model.ans, x, regr.par = 0, fct1 = 1, fct2 = 1, fct3 = 1, 
		fct4 = 1, fct5 = 1, CES = NA, twice, ttt = 0, y = 0, 
		increase, x.mn = NA, par.start = NA, track = FALSE) {
	
	if (track) 
		print("f.expect.con")
	
	nr.aa <- max(fct1)
	nr.bb <- max(fct2)
	nr.var <- max(fct3)
	nr.cc <- max(fct4)
	nr.dd <- max(fct5)
	
	
	if (model.ans != 11) {
		
		nrp <- length(regr.par)
		aa0 <- rep(0, length(x))
		aa.tmp <- regr.par[1:nr.aa]
		for (ii in (1:nr.aa)) aa0 <- aa0 + aa.tmp[ii] * (fct1 == 
						ii)
		bb0 <- rep(0, length(x))
		bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
		for (jj in (1:nr.bb)) bb0 <- bb0 + bb.tmp[jj] * (fct2 == 
						jj)
		par3 <- regr.par[nr.aa + nr.bb + 1]
		if (length(par3) == 0 || is.na(par3)) 
			par3 <- 0
		par4 <- regr.par[nr.aa + nr.bb + nr.cc + 1]
		cc0 <- par3
		dd0 <- par4
		if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 19, 
				20, 24, 25, 41, 42, 46)) {
			if (max(fct4) > 1) {
				cc0 <- rep(0, length(x))
				cc.tmp <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + 
									nr.bb + nr.cc)]
				for (kk in (1:nr.cc)) cc0 <- cc0 + cc.tmp[kk] * 
							(fct4 == kk)
			}
		}
		
		if (model.ans %in% c(3, 8, 13, 18, 23)) {
			dd0 <- par3
			if (max(fct5) > 1) {
				dd0 <- rep(0, length(x))
				dd.tmp <- regr.par[(nr.aa + nr.bb + 1):length(regr.par)]
				for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
							(fct5 == kk)
			}
		}
		
		if (model.ans %in% c(5, 6, 10, 15, 20, 25, 41, 42, 
				46)) {
			dd0 <- par4
			if (max(fct5) > 1) {
				dd0 <- rep(0, length(x))
				dd.tmp <- regr.par[(nr.aa + nr.bb + nr.cc + 
									1):length(regr.par)]
				for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
							(fct5 == kk)
			}
		}
	}
	
	switch(as.character(model.ans), 
			
			"1" = { # Null model
				
				y.expect <- aa0
				
			}, 
			
			"11" = { # Full model
				
				x.gr <- levels(factor(x))
				x.fact <- factor(x)
				y.tmp <- rep(NA, length(x))
				y.expect <- rep(0, length(x))
				if (nr.var > 1 & nr.aa == 1 & nr.bb == 1) {
					for (jj in (1:nr.var)) for (ii in (1:length(x.gr))) {
							y.tmp <- y[x.fact == x.gr[ii] & fct3 == jj]
							if (length(y.tmp) > 0) {
								y.mn <- exp(mean(log(y.tmp)))
								y.expect <- y.expect + y.mn * (x.fact == 
											x.gr[ii]) * (fct3 == jj)
							}
							y.mn <- y.tmp
						}
				} else {
					for (jj in (1:nr.aa)) for (kk in (1:nr.bb)) {
							for (ii in (1:length(x.gr))) {
								y.tmp <- y[x.fact == x.gr[ii] & fct1 == jj & 
												fct2 == kk]
								if (length(y.tmp) > 0) {
									y.mn <- exp(mean(log(y.tmp)))
									y.expect <- y.expect + y.mn * (x.fact == 
												x.gr[ii]) * (fct1 == jj) * (fct2 == kk)
								}
							}
						}
				}
			}, 
			
			"13" = { #E3 - CED
				y.expect <- aa0 * (CES + 1)^((x/bb0)^dd0)
				
			},
			
			"15" = { #E5 - CED
				y.expect <- aa0 * (cc0 - (cc0 - 1) * ((CES + 
									1 - cc0)/(1 - cc0))^((x/bb0)^(dd0)))
				
			}, 
			
			"23" = { #H3 - CED
				
				dum <- f.bb.con(model.ans, cc = NA, dd = dd0, CED = bb0, CES, track = track)
				
				if (increase == 1){
					
					y.expect <- aa0 * (1 - x^dd0/(-(-dum)^dd0 + x^dd0))
					
				} else {
					
					y.expect <- aa0 * (1 - x^dd0/((dum)^dd0 + x^dd0))
					
				}
			}, 
			
			"25" = { #H5 - CED
				
				dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, CES, track = track)
				y.expect <- aa0 * (1 + ((cc0 - 1) * x^dd0)/(dum^dd0 + 
								x^dd0))
				
			})
	
	if (track) 
		print("f.expect.con:  END")
	
	return(y.expect)
	
}

#' Determine value for parameter b for continuous response
#' @param model.ans integer, indicates the model that will be fitted;
#' one of \code{c(13, 15, 23, 25)}
#' @param cc numeric, value for parameter c
#' @param dd numeric, value for parameter d
#' @param CED numeric, value for the CED
#' @param CES numeric, value for the CES
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric, value for parameter b
#' @export
f.bb.con <- function (model.ans, cc, dd, CED, CES, track = FALSE) {
	
	if (track) 
		print("f.bb.con")
	
	switch(as.character(model.ans),
			
			"13" = {
				
				bb <- log(CES + 1)/CED^dd
				
			}, 
			
			"15" = {
				
				CES.tmp <- -abs(CES) * (cc < 1) + CES * (cc >= 1)
				dum <- (CES.tmp + 1 - cc)/(1 - cc)
				bb <- -log(dum)/(CED^dd)
				
			}, 
			
			"23" = {
				
				if (CES >= 0) {
					
					bb <- -CED/(CES/(1 + CES))^(1/dd)
					
				} else {
					
					bb <- CED/(-(CES/(1 + CES)))^(1/dd)
					
				}
				
			}, 
			
			"25" = {
				
				if (cc == 0) {
					
					bb <- -CED * ((1 + CES)/CES)^(1/dd)
					
				} else {
					
					CES <- abs(CES)
					if (cc >= 1) {
						
						bb <- CED/((CES/(cc - 1 - CES))^(1/dd))
						
					} else {
						
						bb <- CED/((-CES/(cc - 1 + CES))^(1/dd))
						
					}
					
				}
				
			})
	
	if (track) 
		print("f.bb.con:END")
	
	return(bb)
	
}

#' Determine value for the CED; continuous response
#' @param model.ans integer, indicates the model that will be fitted; 
#' one of \code{c(1, 11, 18, 20, 13, 15, 23, 25)}
#' @param params numeric vector, parameter values for the model parameters a,
#' b, c and d
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric, value for the CED 
#' @export
f.inv.con <- function (model.ans, params, track = FALSE) {
	
	
	if (track) 
		print("f.inv.con")
	
	aa <- params[1]
	bb <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	
	if(model.ans == 1){
		
		CED <- NA
		warning("For null model CED is not defined")
		
	} else if (model.ans == 11) {
		
		CED <- NA
		
#  } 
#  else if(model.ans == 18){
#    
#    CES.tmp <- CES * (par3 >= 1) - abs(CES) * (par3 < 1)
#    CED <- ((-sign(bb) * CES.tmp)/(1 + CES.tmp))^(1/par3) * (abs(bb))
#    
#  } else if(model.ans == 20){
#    
#    CES.tmp <- CES * (par3 >= 1) - abs(CES) * (par3 < 1)
#    
#    if (par3 == 0){
#      
#      CED <- ((-(bb^par4) * CES.tmp)/(1 + CES.tmp))^(1/par4)
#      
#    } else {
#      
#      CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))^(1/par4)
#      
#    }
#    
	} else {# Model 13, 15, 23, 25
		
		CED <- bb
		
	}
	
	if (is.na(CED) && bb < 1e-10) 
		CED <- 1e+10
	
	if (track) 
		print("f.inv.con : END")
	
	return(CED)
	
}

#' Determine values for the CED; continuous response
#' @param model.ans integer, indicates the model that will be fitted; 
#' one of \code{c(1, 11, 18, 20, 13, 15, 23, 25)}
#' @param regr.par.matr numeric matrix, each row contains a vector of parameter
#' values as needed by the function f.inv.con()
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, calculated CED values
#' @export
f.ced.con <- function (model.ans, regr.par.matr, track = FALSE) {
	
	if (track) 
		print("f.ced.con")
	
	CED <- apply(regr.par.matr, 1, function(par.tmp){
				
				f.inv.con(model.ans, par.tmp, track = track)
				
			})
	
	if (track) 
		print("f.ced.con: END")
	
	return(CED)
}





#' Fit linear model for continuous response
#' @param xtmp numeric vector, values for the independent variable
#' @param ytmp numeric vector, values for the dependent variable
#' @param model.ans integer, indicates which model will be fitted
#' @param dtype integer, indicates the type of response that is used
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return data frame with intercept and slope of the fitted model, the variable
#' increase is 1 if the slope is positive and -1 otherwise 
#' @export
f.start.lm.con <- function (xtmp, ytmp, model.ans, dtype, track = FALSE) {
	
	if (track) 
		print("f.start.lm.con")
	
	if (dtype %in% c(25, 250)) {
		
		yTransformed <- ytmp
		
	} else if (dtype %in% c(26, 260)) {
		
		yTransformed <- sqrt(ytmp)
		
	} else {
		
		yTransformed <- log(ytmp)
		
	}
	
	res.lm <- lm(yTransformed ~ xtmp)
	
	intercept <- res.lm[[1]][1]
	slope <- res.lm[[1]][2]
	
	if (dtype %in% c(25, 250)) {
		
		intercept <- intercept
		
	} else if (dtype %in% c(26, 260)) {
		
		intercept <- intercept^2
		
	} else {
		
		intercept <- exp(intercept)
		
	}
	
	increase <- 1 * (slope >= 0) - 1 * (slope < 0)
	
	if (!(model.ans %in% 2:15)) {
		
		yTransformed <- 1/ytmp
		res.lm <- lm(yTransformed ~ xtmp)
		intercept <- 1/res.lm[[1]][1]
		slope <- 1/res.lm[[1]][2]/intercept
		
	}
	
	if (track) 
		print("f.start.lm.con:  END")
	
	toReturn <- data.frame(intercept = intercept, slope = slope, increase = increase)
	
	return(toReturn)
	
}

#' Calculate residuals: difference of observed and expected response value
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.resid.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.resid.con")
	
	with(ans.all, {
				
				pred.value <- f.expect.con(model.ans, x, regr.par = regr.par, 
						fct1 = fct1, fct2 = fct2, fct3 = fct3, fct5 = fct5, 
						CES = CES, twice = twice, ttt = 0, y = yNew, 
						increase = increase, x.mn = x.mn)
				
				if (dtype %in% c(1, 5, 10)) {
					
					regr.resid <- log(yNew) - log(pred.value)
					
				} else if (dtype %in% c(25, 250)) { 
					
					regr.resid <- yNew - pred.value 
					
				} else if (dtype %in% c(26, 260)){
					
					regr.resid <- sqrt(yNew) - sqrt(pred.value)
					
				} 
				
				ans.all$regr.resid <- regr.resid
				ans.all$pred.value <- pred.value
				
				if (track) 
					print("f.resid.con:  END")
				
				return(ans.all)
				
			})
}



#' Fit the model for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm4.con <- function (ans.all, track = FALSE) {
	
	
	if (track) 
		print("f.mm4.con")
	
	ans.all <- with(ans.all, {
				
				ans.all$conf.int <- matrix(NA, ncol = 2)
				
				if (track) 
					cat("\n\n optimizing .... \n\n")
				
				ans.all <- f.nlminb(ans.all, track = track)
				f.hit.constr(ans.all, track = track)
				
				MLE <- ans.all$MLE
				ans.all$regr.par <- MLE[-(1:nr.var)]
				if (model.ans != 11){
					
					ans.all$regr.par.matr <- f.pars(ans.all, track = track)$regr.par.matr
					
				} 
				
				ans.all <- f.resid.con(ans.all, track = track)
				ans.all$boot <- FALSE
				
				if (cont) {
					
					if (track) 
						cat("\n\nNOTE: var is the variance of the residuals on (natural) log-scale\n")
					
				}
				
				f.check.cc(ans.all, track = track)
				
				covar.ans <- 1
				fit.res <- ans.all$fit.res
				varcov.matr <- NA
				corr.matr <- NA
				if (ans.all$converged) 
					if (length(MLE[MLE == 0]) != 0) {
						warning("The model has too many parameters, use a simpler model")
						covar.ans <- 1
					}
				if (covar.ans == 2) {
					fit.res$lower <- ans.all$lb
					fit.res$upper <- ans.all$ub
					fit.res$scale <- scale.dum
					varcov.matr <- vcov(fit.res)
					cat("\n\n variance-covariance matrix: \n")
					print(varcov.matr)
					if (varcov.matr[1] == -1000) {
						cat("\ntry fitting again, with scale adjusted\n")
						corr.matr <- 1000
					}
					else {
						tmp.matr <- sqrt(diag(varcov.matr))
						tmp.matr <- diag(1/tmp.matr)
						corr.matr <- tmp.matr %*% varcov.matr %*% tmp.matr
						cat("\n correlation matrix: \n")
						print(corr.matr)
						cat("\n maximum correlation:     ", round(max(abs(lower.tri(corr.matr) * 
																corr.matr)), 4), "\n")
					}
				}
				
				ans.all$par.start <- ans.all$MLE
				ans.all$list.logic <- FALSE
				if (model.ans %in% c(12:15, 22:25)) {
					
					ans.all$CED <- f.ced.con(model.ans, ans.all$regr.par.matr)
					
				}
				
				if (track) 
					print("f.mm4.con: END")
				
				return(ans.all)
				
			})        
	
	
	
}


#' Calculate the CED values for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm6.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.mm6.con")
	
	ans.all <- with(ans.all, {
				
				if (fitted || list.logic) {
					
					ans.all <- f.pars(ans.all, track = track)
					
				}
				
				CED <- f.ced.con(model.ans, regr.par.matr)
				CED.unscaled <- CED
				
				if (max(nr.aa, nr.bb) == 1) {
					
					CED.unscaled <- sf.x * CED[1]
					
				} else {
					
					CED.unscaled <- sf.x * CED
					
				}
				
				
				if (dtype %in% c(5, 15)) {
					
					warning("For calculating confidence intervals with litter effects 
									use the bootstrap method")
					
					ans.all$CED <- CED
					
					return(ans.all)
					
				}
				
				f.check.cc(ans.all, track = track)
				
				ans.all$group <- 0
				ans.all <- f.CI(ans.all, track = track)
				
				
				nr.CED <- max(fct2)
				
				if (!all(is.na(conf.int))) {
					
					CED.unique <- unique(CED)
					
					if (length(CED.unique) > 1) {
						
						CED.unscaled <- sf.x * CED.unique[nr.CED]
						
					} else {
						
						CED.unscaled <- sf.x * CED.unique
						
					}          
				}
				
				
				ans.all$CED <- CED
				ans.all$CED.origScale <- CED.unscaled
				ans.all$BMR <- NA
				
				if (track) 
					print("f.mm6.con:  END")
				
				return(ans.all)
				
			})
}

#' Define lower and upper bounds for the model parameters; continuous response 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.constr.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.constr.con")
	
	with(ans.all, {        
				
				# MV added this
				if(!is.null(ans.all$parameterConstraints)){
					
					lower <- ans.all$parameterConstraints[,"lowerBound"]
					upper <- ans.all$parameterConstraints[,"upperBound"]
					
					lower.var <- lower[1]
					lower.aa <- lower[2]
					lower.bb <- lower[3]
					upper.var <- upper[1]
					upper.aa <- upper[2]
					upper.bb <- upper[3]
					
					if (model.ans %in% c(15, 25)) {
						
						lower.cc <- lower[4]
						upper.cc <- upper[4]
						lower.dd <- lower[5]
						upper.dd <- upper[5]
						
					} else if (model.ans %in% c(13, 23)) {
						
						lower.dd <- lower[4]
						upper.dd <- upper[4]
						
					}
					
				} else {
					
					lower.var <- 1e-06
					upper.var <- 10
					
					if (dtype %in% c(25, 250)) {
						
						upper.var <- 1e+10
						
					} else if (dtype %in% c(26, 260)) {
						
						upper.var <- 1e+04
						
					}           
					
					if (model.ans %in% c(1, 13, 15, 23, 25)) {
						
						lower.aa <- min(yNew / 100)
						upper.aa <- max(yNew * 100)
						
						lower.bb <- 1e-06
						upper.bb <- Inf
						
						if (model.ans == 1) {
							
							lower.bb <- -Inf
							
						}
						
						if (increase == 1) {
							
							lower.cc <- 1.1
							upper.cc <- 1e+06
							
						} else {
							
							lower.cc <- 1e-06
							upper.cc <- 0.9
							
						}
						
						lower.dd <- 0.01
						upper.dd <- 10
						
						if (model.ans == 1) {
							
							lower <- lower.aa
							upper <- upper.aa
							
						} else if (model.ans %in% c(13, 23)) {
							
							lower <- c(lower.aa, lower.bb, lower.dd)
							upper <- c(upper.aa, upper.bb, upper.dd)
							
						} else {
							
							lower <- c(lower.aa, lower.bb, lower.cc, lower.dd)
							upper <- c(upper.aa, upper.bb, upper.cc, upper.dd)
							
						}
						
						lower <- c(lower.var, lower)
						upper <- c(upper.var, upper)
						
					} else if (model.ans == 11) {
						
						lower <- c(lower.var, regr.par)
						upper <- c(upper.var, regr.par)
						
					}
					
				}        
				
				
				# Repeat lower and upper bounds nr.<par> times
				if (model.ans == 1) {
					
					lb <- c(rep(lower.aa, nr.aa))
					ub <- c(rep(upper.aa, nr.aa))
					
				} else if (model.ans %in% c(13, 23)) {
					
					lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
							rep(lower.dd, nr.dd))
					ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
							rep(upper.dd, nr.dd))
					
				} else if (model.ans %in% c(15, 25)) {
					
					lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
							rep(lower.cc, nr.cc), rep(lower.dd, nr.dd))
					ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
							rep(upper.cc, nr.cc), rep(upper.dd, nr.dd))
					
				} else if (model.ans == 11) {
					
					lb <- regr.par
					ub <- regr.par
					
				}
				
				lb <- c(rep(lower.var, nr.var), lb)
				ub <- c(rep(upper.var, nr.var), ub)
				
				# TODO both lb and lower needed; ub and upper?
				ans.all$lb <- lb
				ans.all$ub <- ub
				ans.all$lower <- lower
				ans.all$upper <- upper
				
				if (track) 
					print("f.constr.con:  END")
				
				return(ans.all)
				
			})
}
