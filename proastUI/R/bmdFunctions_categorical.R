#' Main function for calculations with categorical response type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all 
#' @export
f.cat <- function (ans.all, track = FALSE) 
{
	if (track) 
		print("f.cat")
	
	# Note: contrary to f.cont, f.start.X not included at the beginning of the function
	# added, otherwise par.start is missing
	ans.all <- f.start.bin(ans.all) # user f.start.bin because model.type == 1

	for(main.ans.single in c(2, ans.all$main.ans, 15)){
		
		switch(as.character(main.ans.single), 
				
			# 4: fit regression
			'4' = {
					
				ans.all <- f.mm4.cat(ans.all, track = track)
		
			},
			
			# 6: Calculate BMD/CED
			'6' = {
					
				ans.all <- f.mm6.cat(ans.all, track = track)
				
			}, 
				
			# end session
				
			'15' = {
					
				return(ans.all)
				
			})
	
	}
	
}

#' Fit the model for a categorical response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm4.cat <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.mm4.cat")
	
	with(ans.all, {
				
		if (length(ans.all$text.par) == length(ans.all$text.par.0)) 
			ans.all$text.par <- ans.all$text.par.0
		
		# extract constraints -> to include in the UI
#		ans <- 1
#		while (ans <= npar) {
#			ans <- menu(c(paste(text.par, ":      ", lb, " ---- ", 
#				ub), "fit model"), title = "\n What parameter do you want to constraint? ")
#			if (ans <= npar) {
#				cat("\n Give infinite constraints as NA\n\n")
#				lb[ans] <- eval(parse(prompt = paste("Give lower bound for ", 
#					text.par[ans], ": ")))
#				ub[ans] <- eval(parse(prompt = paste("Give upper bound for ", 
#					text.par[ans], ": ")))
#			}
#		}
				
#		if ((model.type == 1) & (model.ans == 14)) 
#			cat("")
#		else for (ii in (1:npar)) {
#			if (lb[ii] == ub[ii]) 
#				text.par[ii] <- paste(text.par[ii], "(fixed)", sep = "")
#		}
				
#		cat("\n optimizing .... \n")
#		ans.all$lb <- lb
#		ans.all$ub <- ub
		ans.all <- f.nlminb(ans.all, track = track)
		f.hit.constr(ans.all)
		ans.all$CED.matr <- matrix(NA, ncol = 2, nrow = 1)
		loglik <- ans.all$loglik
		MLE <- ans.all$MLE
		
		if (model.type == 1 & model.ans == 14) {
			
			pi.full <- MLE
			if (dtype == 6) {
				alfa.start <- MLE[1:alfa.length]
				pi.full <- pi.full[-(1:alfa.length)]
			}
			
		}
		
		converged <- ans.all$converged
		aa <- MLE[1:max(fct1)]
		if ((min(x) == 0) & (model.type == 1) & (prod(aa) == 0)) {
			dum <- (y[x == 0] > 0)
			warning("ATTENTION: one of parameters a is estimated to be zero,",
				"which is impossible for nonzero observed response at dose zero.",
				"Therefore you should refit the model with positive lower bound for parameter a")
		}
		
		if (dtype == 6) 
			regr.par <- MLE[-(1:alfa.length)]
		else regr.par <- MLE

		if (nrp == 0) 
			regr.par <- numeric()
		
		ans.all$pi.full <- pi.full
		ans.all$regr.par <- regr.par
		ans.all$th.par <- th.par
		ans.all$sig.par <- sig.par
		ans.all$l.ty <- 1
		ans.all <- f.pars(ans.all)
		
		# remove all plots
		if (dtype %in% 2:3) {
			ans.all$combi <- FALSE
			ans.all$plot.type <- 1 # original: 5
			ans.all$categ.ans <- 0
		}

		if (!model.ans %in% 25) 
			if (any(MLE[1:nr.aa]) < 0) 
				warning("The parameter a is estimated to be negative",
					"Use maximum likelihood with constraints and refit the model")

		ans.all$odt <- odt
		ans.all$MLE <- MLE
		ans.all$alfa.start <- alfa.start
		ans.all$par.start <- MLE
		ans.all$categ.ans <- 1
		ans.all$fitted <- TRUE
		
		if (track) 
			print("f.mm4.cat:  END")
		
		return(ans.all)
		
	})
		
}

#' Calculate the CED values for a continuous response
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.mm6.cat <- function (ans.all, track = FALSE) 
{
	if (track) 
		print("f.mm6.cat")
	
	with(ans.all, {
				
		if (twice) 
			max.lev <- max(max(fct1), max(fct2))
		if (!twice) 
			max.lev <- max(fct1) * max(fct2)
		
		conf.int <- matrix(NA)
		
		if (fitted == FALSE) {
			
			cat("\nATTENTION: you did not yet fit the model!\n\n")
			regr.par <- regr.start
			th.par <- c(0, th.start)
			sig.par <- sig.start
		}
		
		profile.ans <- 2
		
		if (fitted & CED.model)	profile.ans <- 1
#					profile.ans <- menu(c("yes", "no"), title = "\nDo you want to calculate confidence intervals for BMD?")
		
		if (dtype == 3) {
			ces.ans <- 1
			ans.all$ces.ans <- 1
			ans.all$CES <- 0
		}
		
		# ask for CES value (done in the UI)
#			if (dtype != 3 && !CED.model) {
#				ces.ans <- menu(c("ED50", "Additional risk, i.e. P[BMD] - P[0]", 
#									"Extra risk, i.e. (P[BMD]-P[0])/(1-P[0])", "", 
#									""), title = "\nWhat type of Benchmark response do you want to consider?")
#				ans.all$ces.ans <- ces.ans
#				switch(ces.ans, CES <- 0, CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of additional risk > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of extra risk > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the BMR,\nin terms of difference in z-score > ")), 
#						CES <- eval(parse(prompt = "\nGive value for the percent change in risk > ")))
#			}
		
		ans.all$regr.par <- regr.par
		ans.all$th.par <- th.par
#		ans.all$CES <- CES
		CED.list <- f.ced.cat(ans.all)
		CED.matr <- CED.list$CED.matr
		gr.txt <- CED.list$gr.txt
		ans.all$gr.txt <- gr.txt
		ans.all$response.matr <- CED.list$response.matr
		pi.bg.matr <- CED.list$pi.bg.matr
			
		ans.all$CED.matr <- CED.matr
		nr.CED <- length(CED.matr[, 1])
		
		if (profile.ans == 1) {
			ans.all$trace <- T
			ans.all$trace.plt <- T
			ans.all$group <- 0

			if (model.type == 1 && model.ans == 25 && ces.ans == 1) 
				ans.all$group <- 1:nr.aa
			ans.all <- f.CI(ans.all)
			
			if (ans.all$update) {
				MLE <- ans.all$MLE
				switch(model.type, dum <- 1, dum <- 0)
				CED.list <- f.ced.cat(ans.all)
				CED.matr <- CED.list$CED.matr
				ans.all$regr.par <- MLE
				if (dtype == 6) 
					ans.all$regr.par <- ans.all$regr.par[-1]
				if (model.type == 2) {
					par.lst <- f.split.par(MLE, nrp, nth, dtype)
					ans.all$regr.par <- par.lst$regr.par
					ans.all$th.par <- c(0, par.lst$th.par)
					ans.all$sig.par <- par.lst$sig.par
				}
			}
			
			ans.all$CED.matr <- CED.matr
			ans.all$show <- f.show.cat(ans.all)
			conf.int <- ans.all$conf.int
			nr.CED <- length(conf.int[, 1])
			
			if (nr.CED > 1) 
				for (ii in 1:nr.CED) 
					cat("\nthe CED and the 90% confidence interval for group", 
						gr.txt[ii], "is: \n\n", signif(sf.x * CED.matr[ii, 1], 5), 
						"\n", signif(conf.int[ii, 1], 5), 
						"\n", signif(conf.int[ii, 2], 5), "\n")
			else cat("\nthe CED and the 90% confidence interval for category", 
					CES.cat, "is: \n\n", signif(sf.x * CED.matr[1, CES.cat], 5), 
					"\n", signif(conf.int[1], 5), 
					"\n", signif(conf.int[2], 5), "\n")
			
		}
		
		if (profile.ans == 2) {
			
			if (dtype == 3) {
				CED.categ <- cbind(round(t(sf.x * CED.matr), 3), scores.orig[2:length(scores.orig)])
				cat("\nCEDs (in original dose units) at (original) severity scores:\n")
				print(CED.categ)

			}else if (model.ans != 30) {
				if (nr.aa > 1 | nr.bb > 1) 
					for (ii in 1:nr.CED) 
						cat("\nthe CED for group", gr.txt[ii], "is: \n", 
							signif(sf.x * CED.matr[ii, 1], 5), "\n")
				else cat("\nthe CED is: \n\n", signif(sf.x * CED.matr[1, 1], 5), "\n")
			}
		}
		
		ans.scale <- 0
		
		# remove plot

		if (ces.ans > 1) {
			if (dtype %in% c(4, 6, 84)) {
				ans.all$CI.plt <- TRUE
				ans.all$combi.ans <- FALSE
			}
		}
		ans.all$pi.bg.matr <- pi.bg.matr
		ans.all$shift <- shift
		if (profile.ans == 1)	ans.all$conf.int <- conf.int
		ans.all$plot.type <- plot.type

		if (track)	print("f.mm6.cat:  END")
		
		return(ans.all)
		
	})
		
}

# 
#' Define initial parameter values for the chosen model, categorical response
#' 
#' From f.start.bin of the proast61.3 package, only: adjust == FALSE and tmp.quick == FALSE
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.start.bin <- function (ans.all, track = FALSE) 
{
	if (track) print("f.start.bin")
	
	ans.all$adjust <- FALSE
	ans.all$tmp.quick <- FALSE
	
	with(ans.all, {
				
		xx.tot <- x
		eps <- 1e-06
		
		if (dtype %in% c(4, 6, 84)) {
			
			xlmp <- x
			ylmp <- y
			nlmp <- nn
			fct1.lmp <- fct1
			fct2.lmp <- fct2
		}
		
		if (dtype == 2) {
			
			lmp.lst <- f.lump.cat(x, y, nn, fct1, fct2, dtype, twice, track = track)
			xlmp <- matrix(lmp.lst$xlmp, ncol = 1)
			ylmp <- lmp.lst$ylmp
			nlmp <- lmp.lst$nlmp
			if (length(fct1) > 1) 
				fct1.lmp <- lmp.lst$fct1.lmp
			if (length(fct2) > 1) 
				fct2.lmp <- lmp.lst$fct2.lmp
			
		}
		
		ylmp.corr <- ylmp + 0.01 * (ylmp == 0) - 0.01 * (ylmp == 1)
		
		top <- length(ylmp)
		
		switch(as.character(model.ans), 
				
			# 1: null model
			'1' = {
					
				if (max(fct2) > 1) {
					cat("\nATTENTION: \n\nThis model is not defined for two different b parameters!")
					ans.all$fct2 <- rep(1, length(x))
				} else {
					aa <- sum(ylmp)/length(ylmp)
					regr.start <- c(rep(aa, nr.aa))
					lb <- c(rep(eps, nr.aa))
					ub <- c(rep(1 - eps, nr.aa))
					scale.dum <- rep(1, nr.aa)
				}
			
			}, 
			
			# 14: full model
			'14' = {
				
				if (dtype == 2) {
					regr.start <- as.numeric(y)
					cat("\nNOTE:  Full model may not be applicable for dtype = 2\n\n")
				}
				
				if (nr.aa == 1 && nr.bb == 1) covar.tmp <- F else covar.tmp <- T
				
				if (max(covariate) > 1) {
					covar.tmp <- T
					twice <- T
					fct1 <- covariate
					fct2 <- covariate
					ans.all$fct1 <- covariate
					ans.all$fct2 <- covariate
				}
				
				if (dtype == 6 || dtype == 4) {
					if (!covar.tmp) {
						x.full <- tapply(x, x, mean)
						kk.tot <- as.numeric(tapply(kk, x, sum))
						nn.tot <- as.numeric(tapply(nn, x, sum))
						xx.tot <- as.numeric(tapply(x, x, sum))
						pi.full <- kk.tot/nn.tot
						fct1.full <- rep(1, length(x.full))
						fct2.full <- rep(1, length(x.full))
					}
					if (covar.tmp) {
						pi.full <- numeric()
						x.full <- numeric()
						xx.tot <- x
						fct1.full <- numeric()
						fct2.full <- numeric()
						if (twice) for (jj in levels(factor(fct2))) {
								x.tmp <- x[fct2 == jj]
								y.tmp <- y[fct2 == jj]
								xx.tmp <- as.numeric(tapply(x.tmp, x.tmp, mean))
								x.full <- c(x.full, xx.tmp)
								pi.tmp <- as.numeric(tapply(y.tmp, x.tmp, mean))
								pi.full <- c(pi.full, pi.tmp)
								fct1.tmp <- fct1[fct2 == jj]
								fct1.tmp <- as.numeric(tapply(fct1.tmp, x.tmp, mean))
								fct1.full <- c(fct1.full, fct1.tmp)
								fct2.tmp <- fct2[fct2 == jj]
								fct2.tmp <- as.numeric(tapply(fct2.tmp, x.tmp, mean))
								fct2.full <- c(fct2.full, fct2.tmp)
							}
						if (!twice) 
							for (jj in levels(factor(fct2))) 
								for (ii in levels(factor(fct1))) {
									x.tmp <- x[fct1 == ii & fct2 == jj]
									y.tmp <- y[fct1 == ii & fct2 == jj]
									xx.tmp <- as.numeric(tapply(x.tmp, x.tmp, mean))
									x.full <- c(x.full, xx.tmp)
									pi.tmp <- as.numeric(tapply(y.tmp, x.tmp, mean))
									pi.full <- c(pi.full, pi.tmp)
									fct1.tmp <- fct1[fct1 == ii & fct2 == jj]
									fct1.tmp <- as.numeric(tapply(fct1.tmp, x.tmp, mean))
									fct1.full <- c(fct1.full, fct1.tmp)
									fct2.tmp <- fct2[fct1 == ii & fct2 == jj]
									fct2.tmp <- as.numeric(tapply(fct2.tmp, x.tmp, mean))
									fct2.full <- c(fct2.full, fct2.tmp)
								}
					}
					regr.start <- pi.full
				}

				regr.start <- regr.start + eps * (regr.start == 0) - eps * (regr.start == 1)
				regr.start <- regr.start + eps * (regr.start == 0) - eps * (regr.start == 1)
				lb <- rep(eps, length(regr.start))
				ub <- rep(1 - eps, length(regr.start))
				scale.dum <- c(rep(1, length(regr.start)))
				
			}, 
			
			# 16: two-stage model
			'16' = {
				
				y.tmp <- logb(1 - ylmp.corr)
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				aa <- ylmp[1]/nlmp[1] + 1e-04
				bb <- -1/res.lm[[1]][2]
				bb <- abs(bb)
				cc <- 1
				BMD <- f.inv.bin(model.ans = 3, regr.par = c(aa, bb, cc), CES, ces.ans, track = track)
				if (ces.ans == 1) CES <- 0.5
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
				lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), eps)
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 1e+12)
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), 1)
				
			}, 
	
			# 18: log-logistic model
			'18' = {
				
				y.tmp <- logb(ylmp.corr/(1 - ylmp.corr))
				dose.tmp <- xlmp
				min.dose <- dose.tmp[dose.tmp != 0]
				dose.tmp[dose.tmp == 0] <- min.dose[2]/10
				dose.tmp <- logb(dose.tmp)
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				aa <- ylmp[1]/nlmp[1] + 1e-04
				cc <- -abs((res.lm[[1]][2]))
				bb <- exp(res.lm[[1]][1]/cc)
				cc <- -cc
				CES.tmp <- CES/(1 - aa)
				if (ces.ans == 1) CES.tmp <- 0.5
				BMD <- bb/exp(logb((1 - CES.tmp)/CES.tmp)/cc)
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
				lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
				
			}, 
			
			# 19: Weibull model
			'19' = {
				
				y.tmp <- logb(1 - ylmp.corr)
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = xlmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				aa <- ylmp[1]/nlmp[1] + 1e-04
				bb <- -1/res.lm[[1]][2]
				cc <- 1
				if (ces.ans == 1) CES <- 0.5
				BMD <- bb * (-logb(1 - CES))^(1/cc)
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
				lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
				
			}, 
	
			# 21: log-probit model
			'21' = {
				
				y.tmp <- qnorm(ylmp.corr)
				dose.tmp <- xlmp
				min.dose <- dose.tmp[dose.tmp != 0]
				dose.tmp[dose.tmp == 0] <- min.dose[2]/10
				dose.tmp <- logb(dose.tmp)
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				aa <- ylmp[1]/nlmp[1] + 0.01
				cc <- abs(res.lm[[1]][2])
				bb <- exp(-res.lm[[1]][1]/cc)
				if (ces.ans == 1) CES <- 0.5
				BMD <- bb * exp(qnorm(CES/(1 - aa))/cc)
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
				lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
				
			}, 

			# 24: Gamma model
			'24' = {
				
				BMD <- mean(x)
				aa <- 0.01
				cc <- 1
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb), cc)
				lb <- c(rep(eps, nr.aa), rep(eps, nr.bb), 0.01)
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb), 100)
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb), abs(1/cc))
				
			}, 
			
			# 25: probit model
			'25' = {
				
				y.tmp <- qnorm(ylmp.corr)
				dose.tmp <- xlmp
				min.dose <- dose.tmp[dose.tmp != 0]
				dose.tmp[dose.tmp == 0] <- min.dose[2]/10
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				bb <- res.lm[[1]][2]
				aa <- (-res.lm[[1]][1]/bb)
				if (ces.ans == 1) {
					BMD <- mean(xlmp)
					regr.start <- c(rep(BMD, nr.aa), rep(bb, nr.bb))
					lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
					ub <- c(rep(Inf, nr.aa), rep(Inf, nr.bb))
					scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
				}
				if (ces.ans %in% 2:3) {
					BMD <- qnorm(CES * (1 - pnorm(-aa * bb)) + pnorm(-aa * bb))
					BMD <- BMD/bb + aa
					nr.aa <- 1
					nr.bb <- 1
					regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
					lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
					ub <- c(rep(Inf, nr.aa), rep(Inf, nr.bb))
					scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
				}
				
			}, 
			
			# 26: logistic model
			'26' = {
				
				y.tmp <- logb(ylmp.corr/(1 - ylmp.corr))
				dose.tmp <- xlmp
				dt.fr <- data.frame(yyy = y.tmp, dose.tmp = dose.tmp)
				res.lm <- lm(yyy ~ dose.tmp, dt.fr)
				aa <- -res.lm[[1]][2]
				bb <- res.lm[[1]][2]
				if (ces.ans == 1) CES <- 0.5
				BMD <- -logb((1 - CES)/(CES + exp(aa))) - aa
				BMD <- BMD/bb
				regr.start <- c(rep(aa, nr.aa), rep(BMD, nr.bb))
				lb <- c(rep(-Inf, nr.aa), rep(eps, nr.bb))
				ub <- c(rep(1 - eps, nr.aa), rep(Inf, nr.bb))
				scale.dum <- c(rep(1, nr.aa), rep(1, nr.bb))
			}
			
		)
		
		par.start <- regr.start
		if (dtype == 6) {
			if (!(model.ans == 14 && model.type == 1)) 
				par.start <- c(alfa.start, regr.start)
			else par.start <- c(10, regr.start)
			lb <- c(1e-10, lb)
			ub <- c(Inf, ub)
		}
		
		nrp <- length(regr.start)
		
		# remove plot
#		if (plot.type > 0) {
#			ans.all.tmp <- ans.all
#			ans.all.tmp$regr.par <- regr.start
#			ans.all.tmp$x <- xlmp
#			ans.all.tmp$y <- ylmp
#			ans.all.tmp$nn <- nlmp
#			ans.all.tmp$heading <- "starting values"
#			ans.all.tmp$l.ty <- 2
#			ans.all.tmp$xy.lim[3] <- max(xlmp)
#			ans.all.tmp$show <- ""
#			ans.all.tmp$CED.matr <- NA
#			if (ans.all$ans.plt == 1) 
#				sep <- T
#			else sep <- F
#			aaa <- f.plot.all(ans.all.tmp, sep = sep)
#		}

		loglik.old <- -f.lik.cat(theta = par.start, x = x, y = y, kk = kk, 
			nn = nn, dtype = dtype, 
			fct1 = fct1, fct2 = fct2, 
			nrp = nrp, nth = 1, 
			nr.aa = nr.aa, nr.bb, 
			model.ans = model.ans, CES = CES, 
			ttt = ttt, twice = twice, cens = cens, alfa.length = alfa.length, 
			x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
			ces.ans = ces.ans, CES1 = CES1, CES2 = CES2, 
			nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot,
			track = track)
	
		cat("\nmodel: ", modelname, "\n")
		cat("\nLog-likelihood value at start values: ", signif(loglik.old, 5), "\n")
		cat("\nModel has not yet been fitted !\n")
		ans.all$loglik.old <- loglik.old

		ans.all$par.start <- par.start
		ans.all$regr.start <- regr.start
		ans.all$nrp <- nrp
		ans.all$npar <- length(ans.all$par.start)
		
		if (model.ans == 14 & dtype %in% c(4, 6)) {
			
			ans.all$x.full <- x.full
			ans.all$fct1.full <- fct1.full
			ans.all$fct2.full <- fct2.full
		
		}
		
		ans.all$lb <- lb
		ans.all$ub <- ub
		if (dtype == 6) 
			betabin <- 2
		else betabin <- 1
		if (model.ans %in% c(18:19, 24)) 
			if (constr != -Inf) 
				ans.all$lb[nr.aa + nr.bb + betabin] <- constr
		scale.dum <- abs(1/par.start)
		scale.dum[scale.dum < 0.001] <- 0.001
		scale.dum[scale.dum > 1000] <- 1000
		ans.all$scale.dum <- scale.dum
		ans.all$nn.tot <- nn.tot
		ans.all$kk.tot <- kk.tot
		ans.all$xx.tot <- xx.tot
		
		if (track)	print("f.start.bin:  END")
		
		return(ans.all)
		
	})
}

#' 
#' @param x vector independent variable
#' @param y vector response
#' @param nn vector number of observations
#' @param fct1 vector covariate on which parameter a is dependent
#' @param fct2 vector covariate on which parameter b is dependent
#' @param dtype integer response data type
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list
#' @export
f.lump.cat <- function (x, y, nn, fct1, fct2, dtype, twice, track = FALSE) {
	
	if (track)	print("f.lump.cat")
	
	x.ss <- tapply(x, x, length)
	if (mean(x.ss) > 8) 
		ngroups <- length(x.ss)
	else ngroups <- 5
	xlmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
	ylmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
	nlmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
	fct1.lmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
	fct2.lmp <- matrix(ncol = ngroups, nrow = max(fct1) * max(fct2))
	kk <- 0
	for (jj in 1:max(fct2)) {
		if (!twice) {
			for (ii in 1:max(fct1)) {
				kk <- kk + 1
				x.part <- x[fct1 == ii & fct2 == jj]
				y.part <- y[fct1 == ii & fct2 == jj]
				sub <- round(length(x.part)/ngroups)
				e1 <- 0
				e2 <- 0
				for (ee in 1:ngroups) {
					e2 <- e2 + sub
					xtmp <- x.part[e1:e2]
					ytmp <- y.part[e1:e2]
					xtmp <- xtmp[is.finite(xtmp)]
					ytmp <- ytmp[is.finite(ytmp)]
					xlmp[kk, ee] <- mean(xtmp)
					ylmp[kk, ee] <- mean(ytmp)
					e1 <- e1 + sub
				}
				nlmp[kk, ] <- rep(sub, ngroups)
				fct1.lmp[kk, ] <- rep(ii, ngroups)
				fct2.lmp[kk, ] <- rep(jj, ngroups)
			}
		}
		else {
			kk <- kk + 1
			x.part <- x[fct2 == jj]
			y.part <- y[fct2 == jj]
			sub <- round(length(x.part)/ngroups)
			e1 <- 0
			e2 <- 0
			for (ee in 1:ngroups) {
				e2 <- e2 + sub
				xtmp <- x.part[e1:e2]
				ytmp <- y.part[e1:e2]
				xtmp <- xtmp[is.finite(xtmp)]
				ytmp <- ytmp[is.finite(ytmp)]
				xlmp[kk, ee] <- mean(xtmp)
				ylmp[kk, ee] <- mean(ytmp)
				e1 <- e1 + sub
			}
			fct1.lmp[kk, ] <- rep(jj, ngroups)
			fct2.lmp[kk, ] <- rep(jj, ngroups)
			nlmp[kk, ] <- rep(sub, ngroups)
		}
	}
	xlmp <- matrix(t(xlmp), ncol = 1, nrow = ngroups * kk)
	ylmp <- matrix(t(ylmp), ncol = 1, nrow = ngroups * kk)
	nlmp <- matrix(t(nlmp), ncol = 1, nrow = ngroups * kk)
	fct1.lmp <- matrix(t(fct1.lmp), ncol = 1, nrow = ngroups * kk)
	fct2.lmp <- matrix(t(fct2.lmp), ncol = 1, nrow = ngroups * kk)
	ylmp <- ylmp[!is.na(xlmp)]
	nlmp <- nlmp[!is.na(xlmp)]
	fct1.lmp <- fct1.lmp[!is.na(xlmp)]
	fct2.lmp <- fct2.lmp[!is.na(xlmp)]
	xlmp <- xlmp[!is.na(xlmp)]
	lmp.lst <- list(xlmp = xlmp, ylmp = ylmp, nlmp = nlmp, fct1.lmp = fct1.lmp, 
			fct2.lmp = fct2.lmp, dtype = dtype)
	
	if (track)	print("f.lump.cat end: END")
	
	return(lmp.lst)
	
}

#' return the CED value for a set of parameters, CES and ces.ans
#' 
#' Originally the function can deal with model.ans in [1:13], 
#' but because currently this function is only used in the f.start.bin function
#' where model.ans == 16, in which model.ans is set to 3.
#' So only this part of the code is retained. 
#' @param model.ans integer, type of response model, only 3 available
#' @param regr.par vector with parameters values (a, b, c)
#' @param CES numeric, CES value
#' @param ces.ans integer, type of benchmark response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return CED value
#' @export
f.inv.bin <- function (model.ans = 3, regr.par, CES, ces.ans, track = FALSE) 
{
	if (track)	print("f.inv.bin")

	model.ans <- match(model.ans)
	
	aa <- regr.par[1]
	bb <- regr.par[2]
	cc <- regr.par[3]
	dd <- regr.par[4]
	
	if (ces.ans == 1){
		# only code for model.ans == 3
		qq <- log(2 - 2 * aa)
		CED <- sqrt((bb^2/cc) * (qq + 1/(4 * cc)))
		CED <- CED - 0.5 * bb/cc
	}
			
	if (ces.ans == 2){
		# only code for model.ans == 3
		if (bb == 0)
			cat("\n\nATTENTION: (one of the) parameter(s) bb  is zero!\n\n \nTherefore the CED cannot be calculated\n") 
		else if (cc == 0) CED <- -bb * logb(1 - CES/(1 - aa)) 
		else {
			dum1 <- (-0.5 * bb)/cc
			dum2 <- sqrt((0.25 * bb^2)/cc^2 - (bb^2/cc) * logb(1 - CES/(1 - aa)))
			CEDa <- dum1 + dum2
			CEDb <- dum1 - dum2
			CED <- max(CEDa, CEDb)
		}
	}

	if (ces.ans == 3){
		# only code for model.ans == 3	
		if (bb == 0)
			cat("\n\nATTENTION: (one of the) parameter(s) bb is zero!\n\nTherefore the CED cannot be calculated\n") 
		else if (cc < 1e-05) CED <- -bb * logb(1 - CES) 
		else {
			dum1 <- (-0.5 * bb)/cc
			dum2 <- sqrt((0.25 * bb^2)/cc^2 - (bb^2/cc) * logb(1 - CES))
			CEDa <- dum1 + dum2
			CEDb <- dum1 - dum2
			CED <- max(CEDa, CEDb)
		}
	}
	
	if (track)	print("f.inv.bin: END")
	
	return(CED)
	
}


#' Calculate likelihood for categorical response
#' @param theta numeric vector, the initial regression parameter values
#' @param x numeric vector, the dose values
#' @param y numeric vector, the response values 
#' @param kk 
#' @param nn 
#' @param dtype integer, determines the type of response
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param nrp 
#' @param nth 
#' @param nr.aa 
#' @param nr.bb 
#' @param model.ans 
#' @param CES numeric, value for the CES
#' @param CES.cat 
#' @param ttt numeric, time variable 
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param cens 
#' @param alfa.length 
#' @param x.full 
#' @param fct1.full 
#' @param fct2.full 
#' @param ces.ans index type of benchmark response
#' @param CES1 
#' @param CES2 
#' @param nn.tot 
#' @param kk.tot 
#' @param xx.tot  
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric value, minus the sum of the scores (total likelihood)
#' @export
f.lik.cat <- function (theta, x, y = NA, kk, nn, dtype, fct1 = 1, fct2 = 1, 
	nrp, nth = NA, nr.aa = 1, nr.bb = 1, model.ans, 
	CES = NA, ttt = 0, twice = NA, 
	cens = 0, alfa.length = 0, x.full = NA, 
	fct1.full = NA, fct2.full = NA, ces.ans = 1, 
	CES1 = CES1, CES2 = CES2, 
	nn.tot, kk.tot, xx.tot = xx.tot, track = FALSE){
	
	if (sum(is.na(theta)) > 0) 
		return(NA)
	
	if (model.ans != 25 & sum(theta == 0, na.rm = T) > 0) {
		theta[theta == 0] <- theta[theta == 0] + 0.1
	}
	
	if (sum(is.na(theta) > 0) > 0) {
		theta[is.na(theta)] <- 0.1
	}
	
	par.lst <- f.split.par(MLP = theta, nrp, nth, dtype, track = track)
	regr.par <- par.lst$regr.par
	th.par <- par.lst$th.par
	sig.par <- par.lst$sig.par
	
	if (dtype != 6) {
		
		if (prod(theta[1:nr.aa]) == 0) {
			kk <- kk[x != 0]
			nn <- nn[x != 0]
			fct1 <- fct1[x != 0]
			fct2 <- fct2[x != 0]
			x <- x[x != 0]
		}
		
		if (model.ans == 14 && nr.aa == 1 && nr.bb == 1) {
			kk <- kk.tot
			nn <- nn.tot
			x <- xx.tot
		}
		
		pipi <- f.expect.bin(model.ans = model.ans, x = x, 
			regr.par = theta, fct1 = fct1, fct2 = fct2, 
#			name = F, 
			kk = kk, nn = nn, dtype = dtype, 
			CES = CES, ttt = ttt, twice = twice, 
			x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
#			trace = FALSE, 
			ces.ans = ces.ans, CES1 = CES1, 
			CES2 = CES2, track = track)
	
		if (sum(!(is.na(pipi) == 0))) {
			return(1e+10)
		}
		
		lik <- numeric(0)
		
		if (dtype == 4 | dtype == 6 | dtype == 84) {
			
			if (sum(is.na(pipi)) > 0) {
				
				totlik <- -1e+12
				
			}else {
				
				if ((min(pipi) == 0) | (max(pipi) == 1)) {
					nth <- 1
					pi.tmp <- rep(pipi[1], nn[1])
					y <- c(rep(1, kk[1]), rep(0, nn[1] - kk[1]))
					for (ii in 2:length(x)) {
						pi.tmp <- c(pi.tmp, rep(pipi[ii], nn[ii]))
						y <- c(y, rep(1, kk[ii]), rep(0, nn[ii] - kk[ii]))
					}
					pipi <- matrix(pi.tmp, ncol = 1)
					dtype <- 400
					
				}else {
					
					if (max(pipi) > 1) {
						cat("\n\nATTENTION: pi > 1 \n")
						cat("current parameter values:\n")
						print(signif(theta, 3))
					}
					
					if (sum(cens) == 0) {
						loglik <- kk * logb(pipi) + (nn - kk) * logb(1 - pipi)
						totlik <- sum(loglik)
					}
					
					if (sum(cens) > 0) {
						lik <- dbinom(kk, nn, pipi) * (cens != 1) + 
							(1 - pbinom(kk, nn, pipi) + 
							dbinom(kk, nn, pipi)) * (cens == 1)
						totlik <- sum(logb(lik))
					}
					
				}
				
			}
			
		}
		
		if (dtype == 2) {
			nth <- 1
			pipi <- matrix(pipi, ncol = 1)
		}
		
		if (dtype == 2 | dtype == 3 | dtype == 400) {
			
			pi0 <- rep(1, length(y))
			for (k in 1:nth) pi0 <- pi0 - pipi[, k]
			pi0[pi0 == 0] <- 1e-12
			lik <- (y == 0) * pi0
			for (k in 1:nth) lik <- lik + (y == k) * pipi[, k]
			loglik <- logb(lik)
			totlik <- sum(loglik)
			
		}
		
		if (!is.na(totlik)) {
			if (totlik == -Inf) {
				totlik <- -2e+10
			}
		}
	}
	
	if (dtype == 6) {
		
		alfa <- theta[1:alfa.length]
		
		if (alfa.length > 1) {
			xtot <- as.numeric(tapply(x, x, mean))
			xtot <- round(xtot, 10)
			x.tmp <- round(x, 10)
			alfa.tmp <- rep(0, length(x.tmp))
			for (ii in 1:length(xtot)) 
				alfa.tmp <- alfa.tmp + alfa[ii] * (x.tmp == xtot[ii])
			alfa <- alfa.tmp
		}
		
		pipi <- f.expect.bin(model.ans = model.ans, x = x, regr.par = regr.par, 
			fct1 = fct1, fct2 = fct2, 
#			name = F, 
			kk = kk, nn = nn, dtype = dtype, 
			CES = CES, ttt = ttt, 
			twice = twice, x.full = x.full, 
			fct1.full = fct1.full, fct2.full = fct2.full, 
			#trace = trace, 
			ces.ans = ces.ans, CES1 = CES1, CES2 = CES2,
			track = track)
			
		pipi[pipi < 1e-04] <- 0.001
		pipi[pipi > 0.9999] <- 1 - 0.001
		lik <- numeric(0)
		beta.0 <- (alfa * (1 - pipi))/pipi
		log.beta.numerator <- lbeta(kk + alfa, nn - kk + beta.0)
		log.beta.denominator <- lbeta(alfa, beta.0)
		bin.coef <- 1
		loglik <- log(bin.coef) + log.beta.numerator - log.beta.denominator
		totlik <- sum(loglik)
		if (!is.finite(totlik)) 
			totlik <- NA
		
	}
	
	return(-totlik)
	
}


#' split the vector of parameters
#' @param MLE vector of parameters
#' @param nrp 
#' @param nth 
#' @param dtype integer, determines the type of response
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return 
#' @export
f.split.par <- function (MLE, nrp, nth, dtype, track = FALSE) 
{
	if (track)	print("f.split.par")
	
	if (dtype == 6) {
		MLE <- MLE[-1]
	}
	
	regr.par <- MLE[1:nrp]
	th.par <- MLE[(nrp + 1):(nrp + nth)]
	sig.par <- MLE[nrp + nth + 1]
	
	if (track)	print("f.split.par: END")
	
	return(list(regr.par = regr.par, th.par = th.par, sig.par = sig.par))
	
}

#' Calculate expected response values, for categorical response
#' @param model.ans integer, determines the type of model to be fitted
#' @param x numeric vector, the dose values
#' @param regr.par numeric vector, regression parameter values
#' @param fct1 numeric, value for parameter a
#' @param fct2 numeric, value for parameter b
#' @param kk 
#' @param nn 
#' @param dtype integer, determines the type of response
#' @param CES numeric, value for the CES
#' @param ttt numeric, time variable 
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @param x.full 
#' @param fct1.full 
#' @param fct2.full 
#' @param ces.ans 
#' @param CES1 
#' @param CES2 
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return numeric vector, the expected response values under the estimated 
#' regression model
#' @export
f.expect.bin <- function (model.ans = NA, x, regr.par = 0, fct1 = 1, fct2 = 1, 
	kk = -1, nn = -1, dtype = 6, CES = NA, ttt = 0, 
	twice = F, x.full = NA, fct1.full = NA, fct2.full = NA, #trace = F, 
	ces.ans = 3, CES1, CES2, track = FALSE) 
{
	if (track)	print("f.expect.bin")
	
	nr.aa <- max(fct1)
	nr.bb <- max(fct2)
	nrp <- length(regr.par)
	aa0 <- rep(0, length(x))
	bb0 <- rep(0, length(x))
	aa.tmp <- regr.par[1:nr.aa]
	bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
	
	for (ii in 1:nr.aa) 
		aa0 <- aa0 + aa.tmp[ii] * (fct1 == ii)
	
	for (jj in 1:nr.bb)
		bb0 <- bb0 + bb.tmp[jj] * (fct2 == jj)
	
	cc <- regr.par[nr.aa + nr.bb + 1]
	dd <- regr.par[nr.aa + nr.bb + 2]
	ee <- regr.par[nr.aa + nr.bb + 3]
	
	switch(as.character(model.ans), 
			
		# 1: null model
		'1' = pipi <- aa0, 
		
		# 14: full model
		'`4' = {
			
			if (dtype == 6) {
				
				pipi <- rep(0, length(x))
				
				if (twice) {
					for (ii in (1:max(fct1.full))) {
						regr.par.tmp <- regr.par[fct1.full == ii]
						x.full.tmp <- x.full[fct1.full == ii]
						for (kk in 1:length(x.full.tmp)) 
							pipi <- pipi + regr.par.tmp[kk] * (x == x.full.tmp[kk]) * (fct1 == ii)
					}
				}else{
					for (jj in (1:max(fct2.full))) 
						for (ii in (1:max(fct1.full))) {
							regr.par.tmp <- regr.par[fct1.full == ii & fct2.full == jj]
							x.full.tmp <- x.full[fct1.full == ii & fct2.full == jj]
							for (kk in 1:length(x.full.tmp))
								pipi <- pipi + regr.par.tmp[kk] * 
									(x == x.full.tmp[kk]) * (fct1 == ii) * (fct2 == jj)
						}
				}
			}
			
			if (dtype == 4)	pipi <- regr.par
			
		}, 
		
		# 16: two-stage model
		'16' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			
			if (ces.ans == -1) {
				if (cc == 0) {
					pipi <- aa0 + (1 - aa0) * (1 - exp((x * logb(1 - CES))/bb0))
				} else {
					dum0 <- ((-2 * bb0 * cc)/(1 - sqrt(1 - 4 * cc * logb(1 - CES))))
					pipi <- aa0 + (1 - aa0) * (1 - exp(-x/dum0 - cc * (x/dum0)^2))
				}
			} else {
				dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
					CED = bb0, CES = CES, ces.ans = ces.ans)
				pipi <- aa0 + (1 - aa0) * (1 - exp(-x/dum0 - cc * (x/dum0)^2))
			}
			
		}, 
		
		# 18: logistic model
		'18' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			if (ces.ans == -1)
				dum0 <- bb0 * exp(logb((1 - CES)/CES)/cc) 
			else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
						CED = bb0, CES = CES, ces.ans = ces.ans)
				
			pipi <- aa0 + (1 - aa0) * (1/(1 + exp(cc * (logb(dum0) - logb(x)))))
			
		}, 
		
		# 19: Weibull model
		'19' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			if (ces.ans == -1)
				dum0 <- bb0/(-logb(1 - CES))^(1/cc)
			else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
					CED = bb0, CES = CES, ces.ans = ces.ans)
		
			pipi <- aa0 + (1 - aa0) * (1 - exp(-(x/dum0)^cc))
		
		}, 
		
		# 21: log-probit model
		'21' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			if (ces.ans == -1)
				dum0 <- bb0 * exp(-qnorm(CES)/cc)
			else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
					CED = bb0, CES = CES, ces.ans = ces.ans)
		
			pipi <- aa0 + (1 - aa0) * pnorm((logb(x) - logb(dum0)) * cc)
		
		}, 
		
		# 24: Gammma model
		'24' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			if (ces.ans == -1)
				dum0 <- qgamma(CES, cc)/bb0
			else dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
					CED = bb0, CES = CES, ces.ans = ces.ans)
		
			pipi <- aa0 + (1 - aa0) * pgamma(dum0 * x, cc)
			
		}, 
		
		# 25: probit model
		'25' = {
			
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			
			if (ces.ans > 1 & ((nr.aa > 1) | (nr.bb > 1))) {
				
				pipi <- NA
				cat("\nATTENTION: covariates not implemented for probit model (f.expect.bin)!!\n")
				ans.all$fct1 <- 1
				ans.all$fct2 <- 1
				
			} else {
				
				aa <- aa0[1]
				bb <- bb0[1]
				if (ces.ans %in% 2:3) {
					dum <- f.bb.bin(model.ans, aa, cc = NA, dd = NA, 
						CED = bb, CES = CES, ces.ans = ces.ans)
					pipi <- pnorm(((x) - (aa)) * dum)
				}
				
				if (ces.ans == 1) {
					dum0 <- f.bb.bin(model.ans, aa, cc = cc, 
						dd = NA, CED = bb, CES = CES, ces.ans = ces.ans)
					pipi <- pnorm(((x) - (dum0)) * bb0)
				}
				
			}
			
		}, 
		
		# 26: logistic model
		{
			if (is.na(CES)) cat("\nValue for CES is not known !! \n")
			dum0 <- f.bb.bin(model.ans, aa0, cc = cc, dd = NA, 
				CED = bb0, CES = CES, ces.ans = ces.ans)
			pipi <- 1/(1 + exp(-aa0 - dum0 * x))
		}
		
	)
	
	return(pipi)
	
}