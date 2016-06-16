#' Wrapper function for plotting the results of proast
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.plot.all <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.plot.all")
	
	if (is.na(ans.all$xy.lim[1])) {
		
		x <- ans.all$x
		ans.all$xy.lim[1] <- min(x[x > 0])/5
		
	}
	
	with(ans.all, {
				
				if (cont) {
					
					ans.all$y.lim <- f.plot.con(ans.all, track = track)
					
					if (model.ans != 11){
						
						f.lines.con(ans.all, track = track)
						
					} 
					
					f.cedlines.con(ans.all, track = track)
					
				} else {
					
					if (dtype %in% c(2, 4, 6, 84)) {
						
						ans.all <- f.plot.frq(ans.all)
						if (model.ans != 14) 
							ans.all$gr.txt <- f.lines.frq(ans.all)
						
						f.cedlines.bin(ans.all)
						
						title(main = paste("\n\n\n", modelname), font.main = 1, 
								cex.main = 0.9)
						
					} else {
						
						ans.all.plt <- f.model.bb(ans.all)
						ans.all$y.lim <- f.lines.cat(ans.all.plt)
						
					}
					
				}
				
				if (track) 
					print("f.plot.all : END")
				
				return(ans.all)
				
			})
}


#' Plot the predicted response values (points) and estimated CIs 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return vector, range of the CI bounds on the predicted response values 
#' @export
f.plot.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.plot.con")
	
	with(ans.all, {
				
				if (is.null(x.leg)) 
					x.leg <- ""
				if (is.null(y.leg)) 
					y.leg <- ""
				
				nr.aa <- max(fct1)
				nr.bb <- max(fct2)
				nr.var <- max(fct3)
				
				twice <- FALSE
				
				if (sum(fct1 != fct2) == 0) 
					twice <- TRUE
				
				f.means <- function(xx, yNew, dtype, sd2.log, nn) {
					
					if (track) 
						print("f.means within f.plot.con")
					
					
					mean.y <- yNew
					
					if (dtype %in% c(1, 5, 15, 101)){
						
						mean.y <- exp(tapply(log(yNew), xx, mean))
						
					} else if (dtype == 26){
						
						mean.y <- tapply(sqrt(yNew), xx, mean)^2
						
					} else if (dtype == 25){
						
						mean.y <- tapply(yNew, xx, mean)
						
					} 
					
					
					mean.y.plt <- mean.y
					
					f.var <- function(dtype, xx, yNew, sd2.log, nn) {
						
						if (track) 
							print("f.var within f.plot.con")
						
						if (dtype %in% c(10, 15, 110, 250, 260)) {
							
							SS <- sum(sd2.log * (nn - 1))
							var.within <- SS/sum(nn - 1)
							df <- sum(nn - 1)
							
						} else {
							
							if (dtype %in% c(1, 5)) {
								
								yNew <- log(yNew)
								
							} else if (dtype == 26){
								
								yNew <- sqrt(yNew)
								
							} 
							
							mn <- tapply(yNew, xx, mean)
							nn <- tapply(yNew, xx, length)
							Vmn <- rep(mn[1], nn[1])
							
							
							if (dtype %in% c(1, 5)){
								
								for (ii in 2:length(mn)){
									Vmn <- c(Vmn, rep(mn[ii], nn[ii]))
								}
							} 
							
							resid <- yNew - Vmn
							df <- sum(nn - 1)
							var.within <- (sum(nn) - 1) * var(resid)/df
							
						}
						
						out.lst <- list(Vsem = sqrt(var.within/nn), df = df)
						return(out.lst)
						
					}
					
					var.out <- f.var(dtype, xx, yNew, sd2.log, nn)
					Vsem <- var.out$Vsem
					df <- var.out$df
					Vconf <- qt(0.975, df) * Vsem
					
					if (dtype %in% c(1, 5, 10, 15)) {
						
						conf.L <- mean.y/exp(Vconf)
						conf.U <- mean.y * exp(Vconf)
						
					} else if (dtype %in% c(25, 250)) {
						
						conf.L <- mean.y - Vconf
						conf.U <- mean.y + Vconf
						
					} else {
						
						conf.L <- (sqrt(mean.y) - Vconf)^2
						conf.U <- (sqrt(mean.y) + Vconf)^2
						
					}
					
					y.lim.tmp <- range(c(conf.L, conf.U))
					
					out.lst <- list(conf.L = conf.L, conf.U = conf.U, y.lim.tmp = y.lim.tmp)
					
					return(out.lst)
					
				}
				
				
				x.plt <- x
				y.plt <- yNew
				
				if (dtype %in% c(10, 110, 250, 260)) {
					
					plt.mns <- 3
					
				} else {
					
					plt.mns <- 1
					
				}
				
				mark <- c(1:2, 4:25, 33:500)
				
#        if (nr.aa == 1 && nr.bb == 1 && nr.var == 1) {
#          
#          if (max(sp.fact) > 1) {
#            
#            fct2 <- sp.fact
#            nr.bb <- max(fct2)
#            cat("\n\nno covariates, but subgroups are plotted distinctly\n\n")
#            
#          }
#          
#        }
				
				calculationsPlot <- f.means(x.plt, y.plt, dtype, sd2.log, nn)
				
				plot(x.plt, y.plt, main = heading, ylim = calculationsPlot$y.lim.tmp, 
						xlab = x.leg, ylab = y.leg, type = "n", col = color[1])
				
				
				if (nr.aa == 1 & nr.bb == 1 & nr.var > 1) {
					
					nr.bb <- nr.var
					fct2 <- fct3
					
				}
				
				if (nr.aa > 1 & nr.bb > 1 && twice) {
					
					nr.aa <- 1
					fct1 <- rep(1, length(x))
					
				}
				
				shift.tmp <- 0
				shift.abs <- 0
				bbb <- (max(x.plt, na.rm = TRUE) - min(x.plt, na.rm = TRUE))/70
				
				zz <- 0
				for (jj in 1:nr.bb) {
					
					for (ii in 1:nr.aa) {
						
						x.part <- x.plt[fct1 == ii & fct2 == jj]
						y.part <- y.plt[fct1 == ii & fct2 == jj]
						y.tmp <- yNew[fct1 == ii & fct2 == jj]
						
						if (length(y.part) > 0) {
							
							zz <- zz + 1
							if (zz > 500) 
								cat("\nATTENTION: number of subgroups too large for plotting\n\n")
							type.dum <- "p"
							lt.y <- 1
							x.part <- x.part + shift.tmp
							
							if (plt.mns < 3) {
								
								points(x.part, y.part, col = color[zz], pch = mark[zz], 
										cex = cex.1, type = "p", lty = lt.y)
								
							}
							
							mean.x <- NA
							
							if (dtype %in% c(1, 5, 15, 25, 26)) {
								
								mean.x <- as.numeric(tapply(x.part, x.part, mean))
								mean.y <- f.means(x.part, y.tmp, dtype)
								
							} else {
								
								mean.x <- x.part
								mean.y <- y.part
								
							}
							
							points(mean.x, mean.y, col = color[zz], pch = mark[zz], 
									cex = 1, type = type.dum, lty = lt.y)
							
							
							if (dtype %in% c(10, 250, 260)) {
								
								sd2.log.part <- sd2.log[fct1 == ii & fct2 == jj]
								nn.part <- nn[fct1 == ii & fct2 == jj]
								
							}
							
							out.lst <- f.means(x.part, y.tmp, dtype, 
									sd2.log = sd2.log.part, 
									nn = nn.part)
							
							conf.L.part <- out.lst$conf.L
							conf.U.part <- out.lst$conf.U
							
							for (qq in 1:length(conf.L.part)) {
								
								mean.x[qq] <- mean.x[qq]
								lines(rep(mean.x[qq], 2), c(conf.L.part[qq], 
												conf.U.part[qq]), col = color[zz])
								lines(c(mean.x[qq] - bbb, mean.x[qq] + 
														bbb), rep(conf.L.part[qq], 2), col = color[zz])
								lines(c(mean.x[qq] - bbb, mean.x[qq] + 
														bbb), rep(conf.U.part[qq], 2), col = color[zz])
								
							}
							
						}
						
						shift.tmp <- shift.tmp + shift.abs
						
					}
					
				}
				
				
				ans.all$y.lim <- range(y.plt)
				ans.all$y.lim.CI <- range(y.plt)
				
				if (track) 
					print("f.plot.con:  END")
				
				return(calculationsPlot$y.lim.tmp)
				
			})
}


#' Plot dotted lines for the determined CED value 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return NULL
#' @export
f.cedlines.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.cedlines.con")
	
	with(ans.all, {
				
				nr.subgr <- nrow(regr.par.matr)
				
				dum.contr <- xy.lim[1]
				low.y <- min(y.lim)
				
				for (jj in 1:nr.subgr) {
					
					CED.0 <- CED[jj]
					CES.0 <- CES
					ES.0 <- f.expect.con(model.ans, CED.0, regr.par.matr[jj,], 
							fct1 = 1, fct2 = 1, fct5 = 1, CES = CES.0, 
							increase = increase)
					
					ES.y <- rep(ES.0, 2)
					ES.x <- c(min(x), CED.0)
					CED.y <- c(low.y, ES.0)
					CED.x <- rep(CED.0, 2)
					
					lines(ES.x, ES.y, lty = 2, lwd = 1)
					lines(CED.x, CED.y, lty = 2, lwd = 1)
					
				}
				
				if (track) 
					print("f.cedlines.con: END")
				
				return(NULL)
				
			})
}


#' Plot the curve of the estimated model 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return NULL
#' @export
f.lines.con <- function (ans.all, track = FALSE) {
	
	if (track) 
		print("f.lines.con")
	
	with(ans.all, {
				
				f.lines.tmp <- function(ans.all.tmp) {
					
					with(ans.all.tmp, {
								
								nbins <- 1000
								nbins.0 <- nbins/10
								
								xline <- seq(from = min(x, na.rm = TRUE), 
										to = max(x, na.rm = TRUE), length = nbins)
								
								twice <- FALSE
								
								expect <- f.expect.con(model.ans, xline, regr.par, 
										fct1 = rep(1, length(xline)), fct2 = rep(1, 
												length(xline)), fct3 = 1, CES = CES, twice = twice, 
										increase = increase, x.mn = x.mn)
								
								lines(xline, expect, lty = 1)
								
								return(NULL)
								
							})
					
				}
				
				ans.all.tmp <- ans.all
				
				ans.all.tmp$dum.contr <- xy.lim[1]
				ans.all.tmp$min.x <- xy.lim[2]
				ans.all.tmp$max.x <- xy.lim[3]
				if ((nr.aa > 1 || nr.bb > 1) && (nr.cc > 1 || nr.dd > 1) && !twice) {
					
					warning("Combination of covariate on c or d and yet another covariate 
									\ndoes not allow plotting of all individual curves")
					
				}
				
				regr.par.matr <- matrix(regr.par, nrow = 1)
				
				count <- 0
				
				if (nr.aa == 1 && nr.bb == 1 && nr.var > 1){
					
					for (ii in 1:nr.var) {
						
						count <- count + 1
						ans.all.tmp$regr.par <- regr.par.matr[1, ]
						ans.all.tmp$col.tmp <- color[count]
						f.lines.tmp(ans.all.tmp)
						
					}
					
				} else {
					
					kk.max <- nrow(regr.par.matr)
					par.bb.old <- NA
					par.bb.new <- 0
					
					for (kk in 1:kk.max) {
						
						ans.all.tmp$regr.par <- regr.par.matr[kk, ]
						if (model.ans > 1) 
							par.bb.new <- ans.all.tmp$regr.par[2]
						
						if (!identical(par.bb.new, par.bb.old) || (model.ans == 1)) {
							count <- count + 1
							ans.all.tmp$col.tmp <- color[count]
						}
						
						f.lines.tmp(ans.all.tmp)
						par.bb.old <- par.bb.new
						
					}
				}
				
				if (track) 
					print("f.lines.con : END")
				
				return(NULL)
				
			})
}