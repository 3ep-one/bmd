#' Main function to perform proast analysis
#' @param odt list as returned by f.scan()
#' @param shinyInput list with values assigned in the shiny app
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return ans.all, list with all results that were obtained during the analysis
#' @export
f.proast <- function (odt = list(), shinyInput, track = FALSE) {
  
  if (track) 
    print("f.proast")
  
  ans.all <- f.ini(odt = odt, shinyInput = shinyInput, track = track)
  ans.all$user.nm <- deparse(substitute(odt))
  
  ans.all <- f.execute(ans.all, track = track)
  
  if (ans.all$cont){
    
    ans.all <- f.con(ans.all)
    
  } else {
    
    ans.all <- f.cat(ans.all)
    
  }
  
  if (track) 
    cat("\n\n\n    f.proast  END \n\n\n")
  
  return(ans.all)
  
}

#' Define default parameter values needed in the proast functions
#' @param odt list as returned by f.scan()
#' @param shinyInput list with values assigned in the shiny app
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return updated version of ans.all, list with all results that were obtained
#' during the analysis
#' @export
f.ini <- function (odt = NULL, shinyInput, track = FALSE) {
  
  if (track) 
    print("f.ini")
  
  
  ans.all <- shinyInput
  
  # TODO always 1?
  ans.all$fit.ans <- 1  # type of error in response
  
  ans.all$group <- 0
  
  
  # Likelihood optimization parameters
  ans.all$control.lst <- f.control(3)
  
  # lower/upper bounds for contrained parameters (only for categorical model)
  if(is.null(ans.all$lb))
    ans.all$lb <- NA
  if(is.null(ans.all$ub))
    ans.all$ub <- NA
  
  if(is.null(ans.all$sf.x)){
    
    ans.all$sf.x <- 1
    
  } 
  
  if(ans.all$cont){
    
    ans.all$model.name <- switch(as.character(ans.all$model.ans),
        
        "1" = "Null model: y = a",
        "11" = "Full model: y = group mean",
        "13" = "E3-CED: y = a*exp(bx^d)",
        "15" = "E5-CED: y = a * [c-(c-1)exp(-bx^d)]",
        "23" = "H3-CED: a * (1 - x^d/(b^d+x^d))",
        "25" = "H5-CED: a * (1 + (c-1)x^d/(b^d+x^d))"
    
    )
    
  } else {
    
    ans.all$model.name <- switch(as.character(ans.all$model.ans),
        
        "1" = "null model, y = a", 
        "14" = "full model: y = observed incidence",
        "16" = "two-stage in terms of BMD",
        "18" = "log-logistic in terms of BMD", 
        "19" = "Weibull in terms of BMD",
        "21" = "log-probit in terms of BMD",
        "24" = "gamma model in terms of BMD", 
        "25" = "probit model in terms of BMD",
        "26" = "logistic model in terms of BMD"
    ) 
          
  }
  
  
  # For the plots
  ans.all$xy.lim <- NA
  ans.all$color <- 1:100
  ans.all$heading <- ans.all$model.name
  
  
  # Save data information
  if (is.list(odt) && length(odt$nvar) == 1) {
    
    ans.all$varnames <- odt$varnames
    ans.all$nvar <- odt$nvar
    ans.all$odt <- odt$data
    ans.all$info <- odt$info
    
  }
  
  # Remove missing values & Check whether x and y are numeric vectors: in server.R
  ans.all$data.0 <- odt$data
  
  allFactors <- c(ans.all$fct1.no, ans.all$fct2.no, ans.all$fct3.no, 
      ans.all$fct4.no, ans.all$fct5.no)
  

  # from this part on, from execute function of proast 61.3
  # scale x
  ans.all$x <- as.numeric(ans.all$data.0[, ans.all$xans]) / ans.all$sf.x
  ans.all$x.leg <- ans.all$varnames[ans.all$xans]
  
  if (ans.all$sf.x != 1){
    
    ans.all$x.leg <- paste0(ans.all$x.leg, "/", ans.all$sf.x)
    
  } 
  
  # set y
  ans.all$y <- ans.all$data.0[, ans.all$yans]
  ans.all$y.leg <- ans.all$varnames[ans.all$yans]
  
  
  # Check for one level factors & provide warning
  if(length(allFactors) != 0){
    
    if (length(allFactors) == 1){
      
      oneLevel <- nlevels(ans.all$data.0[,c(ans.all$fct1.no, ans.all$fct2.no, 
                  ans.all$fct3.no, ans.all$fct4.no, ans.all$fct5.no)]) == 1
      
    } else {
      
      oneLevel <- apply(ans.all$data.0[,c(ans.all$fct1.no, ans.all$fct2.no, 
                  ans.all$fct3.no, ans.all$fct4.no, ans.all$fct5.no)], 2, 
          function(x) nlevels(x) == 1)
      
    }
    
    if(any(oneLevel)){
      
      parameters <- paste(c("a", "b", "var", "c", "d")[oneLevel], collapse = ",")
      warning("The factor you chose as covariate on parameter(s)", parameters, 
          "has only one level\n you might have selected a subgroup for this factor")
      
    }
    
  }
  
  ans.all$factor.name <- ""
  
  for(i in 1:5){
    
    if(is.null(ans.all[[paste0("fct", i, ".no")]])){
      
      ans.all[[paste0("fct", i, ".no")]] <- 0
      
      if(i == 3 & !ans.all$cont){
        
        ans.all[[paste0("fct", i)]] <- 1
        
      } else {
        
        ans.all[[paste0("fct", i)]] <- rep(1, nrow(ans.all$data.0))
        
      }
      
    } else {
      
      column <- ans.all[[paste0("fct", i, ".no")]]
      ans.all[[paste0("fct", i)]] <- as.numeric(factor(ans.all$data.0[, column]))
      ans.all[[paste0("fct", i, ".txt")]] <- levels(factor(ans.all$data.0[, column]))
      ans.all$factor.name <- paste0(ans.all$factor.name, " factor", i, ": ", ans.all$varnames[column])
      
    }
    
  } 
  
  if (track) 
    print("f.ini : END")
  
  return(ans.all)
  
}

#' Calculate values for the response variable given the response data type
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.execute <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.execute")
  
  with(ans.all, {
        
        yNew <- y
        data.0 <- data.0[order(x), ]
        
		# test if covariate is the same for parameters (a, b)|(a, d)|(a, c)
        twice <- FALSE
        if (fct1.no != 0) 
          twice <- fct1.no == fct2.no
        if (!twice & fct1.no > 1) 
          twice <- fct1.no == fct5.no
        if (!twice & fct1.no > 1) 
          twice <- fct1.no == fct4.no
        
        
        if (cont) {
          
          ans.all$low.y <- 0.98 * min(yNew)
          ans.all$upp.y <- 1.02 * max(yNew)
          
          nn <- 0
          
          
          if (dtype %in% c(10, 15, 250, 260)) {
            # Summary data
            
            sd <- data.0[, sans]
            nn <- data.0[, nans]
            
            if (sd.se == 2){
              
              sd <- sd * sqrt(nn)
              
            }               
            
            if (dtype %in% c(10, 15)) {
              # Log-scale
              
              y.mn <- y
              CV <- sd/y.mn
              mn.log <- log(y.mn/sqrt(1 + CV^2))
              yNew <- exp(mn.log)
              ans.all$mn.log <- mn.log
              ans.all$sd2.log <- log(CV^2 + 1)
              
            } else if (dtype == 250) {
              # Orig scale
              
              y.mn <- y
              sd2 <- sd^2
              ans.all$mn.log <- y.mn
              ans.all$sd2.log <- sd2
              ans.all$yNew <- y.mn
              
            }
          }
        }
        
        if (dtype %in% c(4, 6, 84)) {
          # Quantal data
          
          nn <- data.0[, nans]
          if (any(y > nn)) {
            
            stop("Number of responses is larger than sample size")
            
          }
          
          yNew <- y/nn
          
        }
        
        if (dtype %in% 2:3) {
          # Binary/ordinal
          
          nn <- 1
          y.original <- y
          scores.orig <- sort(unique(y))
          
          ymn <- tapply(y, x, mean)
          xmn <- tapply(x, x, mean)
          ymn <- as.numeric(ymn)
          
          if (length(xmn) > 1){
            if (var(xmn != 0) & var(ymn != 0)){ 
              if (cor(xmn, ymn) < 0) {
                
                warning("Proast assumes zero to represent normal, and higher scores abnormal
                        \nYour data seem to have the opposite direction and therefore will be reversed for analysis")
                
                scores.orig <- rev(scores.orig)
                
              }
            }
          } 
          
          score <- 0
          
          for (ii in seq_len(scores.orig)) {
            y[y.original == scores.orig[ii]] <- score
            score <- score + 1
          }
          
          scores.mtr <- cbind(scores.orig, levels(factor(y)))
          dimnames(scores.mtr) <- list(NULL, c("orig.scores", 
                  "temp.scores"))
          dum.ord <- sum(scores.mtr[, 1] == scores.mtr[, 2]) != 
              length(scores.mtr[, 1])
          if (dum.ord) {
            
            warning("The original scores have been transformed for analysis as follows", scores.mtr)
            
          }
        }
        
        # from f.adjust.saved
        ans.all$nr.aa <- max(fct1)
        ans.all$nr.bb <- max(fct2)
        if (cont) {
          
          ans.all$nr.var <- max(fct3)
          
        } else {
          
          ans.all$nr.var <- 0
          
        }
        ans.all$nr.cc <- max(fct4)
        ans.all$nr.dd <- max(fct5)
        
        if (!cont) {
          
          nr.gr.out <- f.nr.gr(ans.all$nr.aa, ans.all$nr.bb, twice)
          ans.all$nr.gr <- nr.gr.out[1]
          
        } else {
          
          ans.all$sd <- sd
          ans.all$x.mn <- mean(x)
          
        }
        
        ans.all$yNew <- yNew
        ans.all$nn <- nn
        
        if (dtype %in% c(4, 6, 84)) {
          
          ans.all$y <- yNew
          
        }
        
        ans.all$twice <- twice
        
        
        if (dtype == 3) {
          
          ans.all$scores.orig <- scores.orig
          ans.all$scores.mtr <- scores.mtr
          ans.all$nth <- max(y)
          
        }
        
        if (dtype == 4 && max(ans.all$fct3) > 1){
          
          nth <- max(as.numeric(fct3))
          
        } 
        
        if (dtype == 6){
          
          ans.all$alfa.length <- 1
          
        } 
        
        ans.all$dtype.0 <- dtype
        
        if (track) 
          print("f.execute:  END")
        
        return(ans.all)
        
        
      })
}







#' Minimization of the log-likelihood function 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.nlminb <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.nlminb")
  
  with(ans.all, {
        
        scale.dum <- abs(1/par.start)
        scale.dum[scale.dum == Inf] <- 1000
        if (dtype == 6) 
          scale.dum[1] <- 1
        
        count <- 0
        
        if (cont) {
          
          fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
			  # lower/upper for parameters
              lower = lb, upper = ub, control = control.lst, 
              x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
              fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
              mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
              Vdetlim = Vdetlim, CES = CES, twice = twice, 
              cens.up = cens.up, 
			  # lower/upper bound for theta (in objective function f.lik.con)
			  lb = lb, ub = ub, 
              par.tmp = par.start, increase = increase, x.mn = x.mn)
          
          loglik.check <- f.lik.con(fit.res$par, x, y, dtype, 
              fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
              nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
              ttt = ttt, fct4 = fct4, fct5 = fct5, 
              cens.up = cens.up, par.tmp = NA, increase = increase, 
              x.mn = x.mn)
          
          if (is.na(fit.res$obj)) 
            fit.res$obj <- 1e-12
          if (is.na(loglik.check)) 
            fit.res$obj <- 1e-12
          else if (fit.res$obj != loglik.check) 
            fit.res$obj <- 1e-12
          
          while (is.na(fit.res$par[1]) || (fit.res$obj == 0) || 
              abs(fit.res$obj) == 1e-12 || !is.finite(fit.res$obj)) {
            
            warning("Model is refitted with other scaling parameter")
            scale.dum <- rep(0.5, length(par.start))
            scale.dum <- 2 * scale.dum
            count <- count + 1
            print(count)
            
            fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum,
                lower = lb, upper = ub, control = control.lst, 
                x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                Vdetlim = Vdetlim, CES = CES, twice = TRUE, cens.up = cens.up, 
                lb = lb, ub = ub, par.tmp = fit.res$par, increase = increase, 
                x.mn = x.mn)
            
            if (is.na(fit.res$obj)) 
              fit.res$obj <- 1e-12
            if (count > 50) {
              fit.res$obj <- -1e+10
              fit.res$par <- 0
            }
            
          }
          
          ans.all$MLE <- fit.res$par
          ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
          
        } else {
          
          fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
			  # lower/upper bounds for parameters
              lower = lb, upper = ub, 
			  control = control.lst, 
              x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
              fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
              nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
              model.type = model.type, ans.nobg = 0, CES = CES, 
              CES.cat = CES.cat, ttt = ttt, twice = twice, 
              cens = cens, fct3 = fct3, fct4 = fct4, fct5 = fct5, 
              x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
              alfa.length = alfa.length, ces.ans = ces.ans, 
              decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
              nn.tot = nn.tot, kk.tot = kk.tot, xx.tot = xx.tot, 
              fct3.ref = fct3.ref)
          
          while (is.na(fit.res$par[1]) | (fit.res$obj == 0) | 
              (fit.res$obj == 1e-12) | is.na(fit.res$obj)) {
            
            warning("Model", model.ans, "refitted with other scaling parameter")
            
			# TODO: include in the UI?
            if (model.ans == 14 && model.type == 1 && dtype == 6) 
              par.start[1] <- eval(parse(prompt = paste("give start value for alfa", 
                          "  > ")))
            
            scale.dum <- rep(0.5, length(par.start))
            scale.dum <- 2 * scale.dum
            count <- count + 1
            
            fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                lower = lb, upper = ub, control = control.lst, 
                x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                model.type = model.type, ans.nobg = 0, CES = CES, 
                CES.cat = CES.cat, ttt = ttt, twice = twice, 
                cens = cens, fct3 = fct3, fct5 = fct5, x.full = x.full, 
                fct1.full = fct1.full, fct2.full = fct2.full, 
                alfa.length = alfa.length, ces.ans = ces.ans, 
                decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
                fct3.ref = fct3.ref)
            
            if (count > 50) {
              fit.res$obj <- -1e+10
              fit.res$par <- 0
            }
          }
          
          ans.all$MLE <- fit.res$par
          ans.all$regr.par <- ans.all$MLE
          
        }
        
        ans.all$loglik <- round(-fit.res$objective, 2)
        
        if (!is.finite(ans.all$loglik)) 
          ans.all$loglik <- -1e+12
        
        message <- fit.res$message
        ans.all$converged <- f.converged(message, fit.res$conv, track = track)
        
        if (dtype == 6 && (model.ans == 14 & model.type == 1)) {
          
          dum <- list()
          dum$loglik <- ans.all$loglik
          dum$MLE <- ans.all$MLE
          ans.all$full.model <- dum
          alfa.start <- ans.all$MLE[1]
          nn.dum <- table(x, x)
          nn.dum <- nn.dum[nn.dum != 0]
          if (length(nn.dum) > 15 && alfa.start > 20) 
            alfa.start <- 5
          ans.all$alfa.start <- alfa.start
          
        }
        
        ans.all$fitted <- TRUE
        
        if (track) 
          print("f.nlminb : END")
        
        return(ans.all)
        
      })
}


#' Define control parameters for nlminb(), which minimizes the log-likelihood function 
#' @param level integer, one of \code{1:3} determines which control parameters
#' are used; default value is 3
#' @return list, with control parameters 
#' @export
f.control <- function (level = 3) {
  
  lst <- list()
  switch(level, {
        
        lst$eval.max <- 50
        lst$iter.max <- 40
        lst$rel.tol <- 0.001
        lst$x.tol <- 0.015
        lst$step.min <- 0.00022
        
      }, {
        
        lst$eval.max <- 100
        lst$iter.max <- 75
        lst$rel.tol <- 1e-06
        lst$x.tol <- 0.00015
        lst$step.min <- 2.2e-07
        
      }, {
        
        lst$eval.max <- 1000
        lst$iter.max <- 750
        
      })
  
  return(lst)
  
}




#' Determine convergence type of the minimization algorithm
#' @param message object as returned from nlminb()$message
#' @param conv.out object as returned from nlminb()$conv
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return boolean, whether the optimization ended with convergence or not;
#' with attribute "message" which is the input object message
#' @export
f.converged <- function (message, conv.out, track = FALSE) {
  
  if (track) 
    print("f.converged")
  
  converged <- FALSE
  
  if (is.null(mode(message))) {
    
    warning("Convergence message is NULL")
    
  } else {
    
    converged <- 1 - conv.out
    
  }
  
  if (track) 
    print("f.converged:  END")
  
  attr(converged, "message") <- message
  
  return(converged)
  
}

#' Check whether any of the MLEs for the parameters hit the lower/upper constraint 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.hit.constr <- function (ans.all, track = FALSE) {
  
  if (track) 
    print("f.hit.constr")
  
  with(ans.all, {
        
        if (!(model.ans %in% c(11, 14))) {
          
          toReport <- c(MLE == lb && MLE != ub, MLE == ub && MLE != lb)
          
          if (any(toReport)) {
            
            warning("A parameter estimate was equal to the 
                    lower/upper constraint:", MLE[toReport])
            
          }
          
        }
        
        if (track) 
          print("f.hit.constr:   END")
        
        return(NULL)
        
      })
}


#' Construct regression parameter matrix 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.pars <- function (ans.all, track = FALSE) {
  if (track) 
    print("f.pars")
  
  with(ans.all, {
        
        nrp <- length(regr.par)
        regr.par.matr <- numeric()
        
        nr.aa <- max(fct1)
        nr.bb <- max(fct2)
        if (length(ans.all$nr.cc) == 0) 
          nr.cc <- 1
        if (length(ans.all$nr.dd) == 0) 
          nr.dd <- 1
        if (length(ans.all$nr.var) == 0) 
          nr.var <- 1
        nr.subgr <- max(nr.aa, nr.bb, nr.cc, nr.dd)
        
        if (nr.aa == 1) 
          fct1.txt <- rep("", nr.subgr)
        if (nr.bb == 1) 
          fct2.txt <- rep("", nr.subgr)
        if (nr.cc == 1) 
          fct4.txt <- rep("", nr.subgr)
        if (nr.dd == 1) 
          fct5.txt <- rep("", nr.subgr)
        if (identical(fct1, fct2)) 
          fct2.txt <- rep("", nr.subgr)
        if (identical(fct1, fct4)) 
          fct4.txt <- rep("", nr.subgr)
        if (identical(fct4, fct2)) 
          fct4.txt <- rep("", nr.subgr)
        
        if (all(c(nr.aa, nr.bb, nr.cc, nr.dd, nr.var) == 1)) {
          
          ans.all$regr.par.matr <- matrix(regr.par, nrow = 1)
          ans.all$nr.gr <- 1
          return(ans.all)
          
        }
        
        kk <- 1
        
        if (all(c(nr.aa, nr.bb, nr.cc, nr.dd) == 1)){
          
          regr.par.matr <- matrix(regr.par, nrow = 1)
          
        } else if (nr.dd == 1) {
          
          if ((dtype %in% c(1, 5)) && (model.ans == 11) && 
              (nr.var > 1)) 
            fct1 <- fct3
          
          par.tmp <- rep(NA, length(regr.par))
          dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
          gr.txt <- character(0)
          for (jj in 1:nr.bb) {
            f1 <- fct1[fct2 == jj]
            f1.lev <- levels(factor(f1))
            for (ii.index in f1.lev) {
              ii <- as.numeric(ii.index)
              f4 <- fct4[fct1 == ii & fct2 == jj]
              f4.lev <- levels(factor(f4))
              for (mm.index in f4.lev) {
                mm <- as.numeric(mm.index)
                par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                            jj], regr.par[nr.aa + nr.bb + mm], dd)
                par.tmp <- par.tmp[!is.na(par.tmp)]
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
                gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                    fct4.txt[mm], sep = "-")
                kk <- kk + 1
              }
            }
          }
          
          if (!cont && model.ans == 14 && model.type == 1) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          }
          else if (cont && model.ans == 11) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          }
          else {
            n.col <- nrp - nr.aa - nr.bb - nr.cc + 3
            regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
          }
        }
        if (nr.cc == 1 && nr.dd > 1) {
          if ((dtype %in% c(1, 5)) && (model.ans == 11) && 
              (nr.var > 1)) 
            fct1 <- fct3
          par.tmp <- rep(NA, length(regr.par))
          cc <- regr.par[nr.aa + nr.bb + 1]
          gr.txt <- character(0)
          for (jj in 1:nr.bb) {
            f1 <- fct1[fct2 == jj]
            f1.lev <- levels(factor(f1))
            for (ii.index in f1.lev) {
              ii <- as.numeric(ii.index)
              f5 <- fct5[fct1 == ii & fct2 == jj]
              f5.lev <- levels(factor(f5))
              for (mm.index in f5.lev) {
                mm <- as.numeric(mm.index)
                par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                            jj], cc, regr.par[nr.aa + nr.bb + 1 + mm])
                par.tmp <- par.tmp[!is.na(par.tmp)]
                regr.par.matr <- rbind(regr.par.matr, par.tmp)
                gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                    fct5.txt[mm], sep = "-")
                kk <- kk + 1
              }
            }
          }
          n.col <- nrp - nr.aa - nr.bb - nr.dd + 3
          if (model.ans == 14 && model.type == 1) {
            ans.all$regr.par.matr <- NA
            return(ans.all)
          }
          else regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                ncol = n.col)
        }
        if (identical(fct1, fct2) && identical(fct2, fct4)) 
          if (nr.aa > 1 && nr.bb > 1 && nr.cc > 1) {
            ii <- 0
            jj <- nr.aa
            if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 
                41, 42)) 
              dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
            else dd <- NULL
            kk <- nr.aa + nr.bb
            zz <- 1
            for (jj in 1:nr.cc) {
              if (nr.aa > 1) 
                ii <- ii + 1
              else ii <- 1
              if (nr.bb > 1) 
                jj <- jj + nr.aa
              else jj <- nr.aa + 1
              kk <- kk + 1
              par.tmp <- c(regr.par[ii], regr.par[jj], regr.par[kk], 
                  dd)
              regr.par.matr <- rbind(regr.par.matr, par.tmp)
              gr.txt[zz] <- paste(fct4.txt[zz])
              zz <- zz + 1
            }
          }
        if (identical(fct1, fct2) && identical(fct2, fct5)) 
          if (nr.aa > 1 && nr.bb > 1 && nr.dd > 1) {
            ii <- 0
            jj <- nr.aa
            if (model.ans %in% c(3, 8, 13, 18, 23, 25, 41, 
                42)) {
              cc <- NULL
              kk <- nr.aa + nr.bb
            }
            if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 
                41, 42)) {
              cc <- regr.par[nr.aa + nr.bb + 1]
              kk <- nr.aa + nr.bb + 1
            }
            zz <- 1
            for (jj in 1:nr.dd) {
              if (nr.aa > 1) 
                ii <- ii + 1
              else ii <- 1
              if (nr.bb > 1) 
                jj <- jj + nr.aa
              else jj <- nr.aa + 1
              kk <- kk + 1
              par.tmp <- c(regr.par[ii], regr.par[jj], cc, 
                  regr.par[kk])
              regr.par.matr <- rbind(regr.par.matr, par.tmp)
              gr.txt[zz] <- paste(fct5.txt[zz])
              zz <- zz + 1
            }
          }
        ans.all$regr.par.matr <- regr.par.matr
        ans.all$gr.txt <- gr.txt
        ans.all$nr.gr <- length(regr.par.matr[, 1])
        
        if (track) 
          print("f.pars: END ")
        
        return(ans.all)
        
      })
}





#' Check whether the parameter value for c is too close to CES
#' @param ans.all list, with all results that were obtained during the analysis
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return boolean, TRUE if the the parameter value for c is too close to CES; 
#' if so warning message is printed 
#' @export
f.check.cc <- function (ans.all, track = FALSE) {
  
  with(ans.all, {
        
        if (track) 
          print("f.check.cc")
        
        cc.OK <- TRUE
        
        if (model.ans %in% c(15, 25)) {
          
          cc <- MLE[nr.var + nr.aa + nr.bb + 1]
          
          if (increase == 1 && (cc < 0.02 + (1 + CES))) {
            
            warning("The value of parameter c is too close to CES. 
                    This indicates (nonrandom) errors in the data or too high value of CES")
            cc.OK <- FALSE
            
          } else if (increase == -1 && (cc > -0.02 + (1 - abs(CES)))) {
            
            warning("The value of parameter c is too close to CES. 
                    This indicates (nonrandom) errors in the data or too high value of CES")
            cc.OK <- FALSE
            
          }
          
        }       
        
        
        return(cc.OK)
        
      })
}


#' Wrapper function for calculating CI around parameter estimates, iterates 
#' f.profile.all()
#' @param ans.all list, with all results that were obtained during the analysis
#' @param nolog boolean; TRUE if calculations should be performed on the original
#' scale, FALSE if on the log-scale; default value is FALSE
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.CI <- function (ans.all, nolog = FALSE, track = FALSE) {
  
  if (track) 
    print("f.CI")
  
  with(ans.all, {
        
        if (length(ans.all$MLE) > 25) {
          
          stop("Number of parameters exceeds 25")
          
        }
        
        profile.out <- f.profile.all(ans.all, nolog, track = track)
        
        updated <- FALSE
        count <- 0
        loglik.tmp <- profile.out$loglik
        
        while (sum(profile.out$MLE.new) != 0) {
          
          count <- count + 1
          
          if (track)
            cat(" \n.... calculation of CI is re-started\n")
          
          ans.all$MLE <- profile.out$MLE.new
          
          if (cont) {
            
            ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
            
          }
          
          ans.all$loglik <- profile.out$loglik
          
          if (ans.all$loglik < loglik.tmp - 0.02) {
            
            warning("Log-likelihood fluctuates, CI is not calculated")
            profile.out$conf.int <- matrix(NA, ncol = 2)
            profile.out$MLE.new <- 0
            
          } else {
            
            loglik.tmp <- ans.all$loglik
            profile.out <- f.profile.all(ans.all, nolog, track = track)
            updated <- TRUE
            
            if (count > 50) {
              
              warning("No global optimum found")
              profile.out$conf.int <- matrix(NA, ncol = 2)
              profile.out$MLE.new <- 0
              
            }
          }
        }
        
        if (updated & !cont) {
          
          CED <- ans.all$MLE[(alfa.length + nr.aa) + 1:nr.bb]
          
          ans.all$CED <- CED
          ans.all$CED.matr <- matrix(CED, ncol = 1)
          ans.all$regr.par <- ans.all$MLE[1:nrp]
          
        }
        
        
        if (nolog) {
          
          conf.int <- profile.out$conf.int
          
        } else {
          
          conf.int <- 10^profile.out$conf.int
          
        }
        
        
        # Replace missing values, if not both limits are missing
        conf.int <- t(apply(conf.int, 1, function(row){
                  
                  isMissing <- is.na(row)
                  
                  if(sum(isMissing) == 1){
                    
                    row[isMissing] <- c(0, Inf)[isMissing]
                    
                  }
                  
                  return(row)
                  
                }))
        
        ans.all$conf.int <- conf.int * sf.x
        ans.all$update <- updated
        ans.all$profile <- profile.out$profile
        
        if (track) 
          print("f.CI END")
        
        return(ans.all)
        
      })
}


# TODO simplify code: lb and ub similar calculation twice
#' Profile likelihood method to calculate confidence intervals 
#' @param ans.all list, with all results that were obtained during the analysis
#' @param nolog boolean, whether the response is log-transformed or not;
#' default value is FALSE (ie log-transformed response)  
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list, updated version of ans.all
#' @export
f.profile.all <- function (ans.all, nolog = FALSE, track = FALSE) {
  
  if (track) 
    print("f.profile.all")
  
  with(ans.all, {
        
        max.runs <- 200
        count <- 0
        count.2 <- 0
        count.max <- 100
        crit <- 0.5 * qchisq(conf.lev, 1)
        jump.crit <- crit * 3
        small.step.low <- FALSE
        small.step.upp <- FALSE
        dist <- 1.01
        par.nr <- 0
        loglik.max <- loglik
        lb.orig <- lb
        ub.orig <- ub
        profile.out <- list()
        MLE.new <- 0
        tb <- "\t"
        
        sp.fact <- 1
        
        # Define number of groups
        if (cont) {
          
          sig2 <- mean(MLE[1:nr.var])
          
          if (group[1] == 0) {
            
            group <- nr.var + nr.aa + (1 : nr.bb)
            
          }
          
          if (nr.var > 1) 
            sp.fact <- fct3
          
          if (dtype == 5) {
            
            warning("Confidence intervals are calculated without accounting for \n
                    nested structure in the data (e.g. intralitter correlations) \n")
            
          }
          
        } else {
          
          if (dtype %in% c(4, 6)) 
            y <- kk/nn
          
          if (group[1] == 0) {
            
            group <- nr.aa + (1 : nr.bb)
            
            if (dtype == 6) {
              
              group <- group + 1
              
            }
            
          }
          
        }
        
        if (nr.bb > 1) 
          sp.fact <- fct2
        if (nr.aa > 1) 
          sp.fact <- fct1
        
        
        
        conf.int <- matrix(NA, nrow = length(group), ncol = 2)
        
        for (jj in group) {
          
          par.nr <- par.nr + 1
          
          if (track) {
            
            if (nr.bb == 1) 
              cat("\nCalculating C.I.......\n")
            else cat("\nCalculating C.I. for group", jj, "......\n")
            
          }
          
          
          # Calculating lower limit
          
          CED.upp.inf <- FALSE
          lb <- lb.orig
          ub <- ub.orig
          lb[jj] <- MLE[jj]
          
          start <- dist * MLE
          loglik.low <- loglik.max
          CED.low <- MLE[jj]
          
          loglik.low.old <- numeric()
          CED.low.old <- numeric()
          
          
          if (cont && dtype %in% c(1, 5, 10)){
            
            step.start <- 1.02 + 0.11 * sig2
            
          } else {
            
            step.start <- 1.08
            
          }
          
          step <- step.start
          stop <- FALSE
          run <- 0
          
          while (!stop) {
            
            run <- run + 1
            lb[jj] <- lb[jj]/step
            ub[jj] <- lb[jj]
            start[jj] <- lb[jj]
            
            step.tmp <- step
            
            if (lb[jj] < lb.orig[jj]) {
              
              step.tmp <- sqrt(step.tmp)
              lb[jj] <- lb.orig[jj] * step.tmp
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              step <- step.tmp
              
            }
            
            ans.all$par.start <- start
            ans.all$lb <- lb
            ans.all$ub <- ub
            
            nlminb.out <- f.nlminb(ans.all, track = track)
            MLE.current <- nlminb.out$MLE
            loglik.current <- nlminb.out$loglik
            
            if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
              
              cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
              start <- dist * nlminb.out$MLE
              loglik.low <- -Inf
              
            }
            
            if (loglik.current > loglik.max + 0.03) {
              
              ans.all$par.start <- MLE.current * dist
              cat("\nlocal optimum found....\n")
              
              profile.out <- list()
              profile.out$MLE.new <- nlminb.out$MLE
              profile.out$loglik <- loglik.current
              profile.out$conf.int <- conf.int
              return(profile.out)
              
            }
            
            if (length(loglik.low) <= 2 || run > max.runs) {
              
              loglik.low <- c(loglik.current, loglik.low)
              CED.low <- c(lb[jj], CED.low)
              
            } 
            
            if (length(loglik.low) > 2) {
              
              if (loglik.current > loglik.low[2]) {
                
                loglik.low <- c(loglik.current, loglik.low[-(1:2)])
                CED.low <- c(lb[jj], CED.low[-(1:2)])
                
              } else if (loglik.current > loglik.low[1]) {
                
                loglik.low <- c(loglik.current, loglik.low[-1])
                CED.low <- c(lb[jj], CED.low[-1])
                
              } else {
                
                loglik.low <- c(loglik.current, loglik.low)
                CED.low <- c(lb[jj], CED.low)
                
              }
              
            }
            
            if (loglik.current == -Inf) {
              
              loglik.low <- loglik.low[-1]
              CED.low <- CED.low[-1]
              
            }
            
            if (length(loglik.low) > 1) {
              
              if (abs(loglik.low[1] - loglik.low[2]) > 0.25 * crit) 
                step <- sqrt(step)
              if (abs(loglik.low[1] - loglik.low[2]) < 0.1 * crit) 
                step <- step^1.26
              
            }
            
            if (track) {
              cat("\n")
              cat(signif(MLE.current, 3), round(loglik.current, 
                      2), round(loglik.max - crit, 2), sep = tb)
              
            }
            
            if (step < 1.00001) {
              
              stop <- TRUE
              small.step.low <- TRUE
              
            }
            
            if (CED.low[1] < 1e-20) 
              stop <- TRUE
            
            if (length(loglik.low) > 30) 
              if (max(abs(diff(loglik.low[1:10]))) < 0.001) {
                stop <- TRUE
              }
            
            if (loglik.max - loglik.current > crit) {
              if (length(unique(loglik.low)) < 6) {
                if (track) 
                  cat("\n------- not enough points: ", length(loglik.low), 
                      " ------------------")
                
                CED.low <- c(CED.low, CED.low.old)
                loglik.low <- c(loglik.low, loglik.low.old)
                CED.low.old <- CED.low
                loglik.low.old <- loglik.low
                power <- length(CED.low)/15
                step <- step.start^power
                step.start <- step
                lb[jj] <- MLE[jj]
                start <- dist * MLE
                run <- 0
                count.2 <- count.2 + 1
                
                if (count.2 > 5) {
                  
                  count.2 <- -1
                  stop <- TRUE
                  loglik.low <- NA
                  
                }
                
              } else if ((loglik.max - loglik.current) > jump.crit) {
                
                lb[jj] <- lb[jj] * step
                step <- step.start^0.5
                step.start <- step
                count <- count + 1
                
                if (step.start < 1.01) {
                  
                  CED.low[1] <- lb[jj]
                  loglik.low[1] <- loglik.current
                  count <- count.max
                  
                }
                
                if (count == count.max) 
                  stop <- TRUE
                run <- 0
                
              } else {
                
                loglik.diff <- diff(loglik.low[order(CED.low)])
                if (sum(loglik.diff < -0.1) > length(CED.low) - 3) {
                  
                  step.start <- 1 + (step.start - 1)/2
                  step <- step.start
                  run <- 0
                  
                } else {
                  
                  stop <- TRUE
                  
                }
              }
            }
            
            if (run > max.runs) {
              
              warning("Maximum number of runs reached in establishing profile")
              if (abs(loglik.low[1] - loglik.low[2]) < 0.01) 
                stop <- TRUE
            }
            
            jump.low <- FALSE
            if (count == count.max) {
              
              warning("The jump in the log-likelihood could not be mitigated")
              jump.low <- TRUE
              
            }
          }
          
          # Calculating upper limit
          
          if (MLE.new[1] != 0) 
            MLE <- MLE.new
          lb <- lb.orig
          ub <- ub.orig
          lb[jj] <- MLE[jj]
          start <- dist * MLE
          loglik.upp <- loglik.max
          CED.upp <- MLE[jj]
          CED.upp.old <- numeric()
          loglik.upp.old <- numeric()
          
          if (cont && dtype %in% c(1, 5, 10)){
            
            step.start <- 1.02 + 0.11 * sig2
            
          } else {
            
            step.start <- 1.08
            
          }
          
          step <- step.start
          run <- 0
          stop <- FALSE
          
          count <- 0
          
          while (!stop) {
            
            run <- run + 1
            lb[jj] <- lb[jj] * step
            ub[jj] <- lb[jj]
            start[jj] <- lb[jj]
            
            
            step.tmp <- step
            if (lb[jj] < lb.orig[jj]) {
              
              step.tmp <- sqrt(step.tmp)
              lb[jj] <- lb.orig[jj] * step.tmp
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              step <- step.tmp
              
            }
            
            step.tmp <- step
            
            while (lb[jj] > ub.orig[jj]) {
              
              cat("f.profile.all: value of parameter exceeds upper constraint, attempting to avoid this\n")
              step.tmp <- sqrt(step.tmp)
              lb[jj] <- lb[jj]/step.tmp
              ub[jj] <- lb[jj]
              start[jj] <- lb[jj]
              step <- step.tmp
              
            }
            
            ans.all$par.start <- start
            ans.all$lb <- lb
            ans.all$ub <- ub
            
            nlminb.out <- f.nlminb(ans.all, track = track)
            MLE.current <- nlminb.out$MLE
            loglik.current <- nlminb.out$loglik
            
            if (!is.finite(loglik.current) || loglik.current <= -1e+10) {
              
              cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
              start <- dist * nlminb.out$MLE
              loglik.upp <- -Inf
              
            }
            
            if (loglik.current > loglik.max + 0.03) {
              
              ans.all$par.start <- MLE.current * dist
              cat("\nlocal optimum found....\n")
              
              profile.out <- list()
              profile.out$MLE.new <- nlminb.out$MLE
              profile.out$loglik <- loglik.current
              profile.out$conf.int <- conf.int
              return(profile.out)
              
            }
            
            if (length(loglik.upp) <= 2 || (run > max.runs)) {
              
              loglik.upp <- c(loglik.current, loglik.upp)
              CED.upp <- c(ub[jj], CED.upp)
              
            }
            
            if (length(loglik.upp) > 2) {
              
              if (loglik.current > loglik.upp[2]) {
                
                loglik.upp <- c(loglik.current, loglik.upp[-(1:2)])
                CED.upp <- c(ub[jj], CED.upp[-(1:2)])
                
              } else if (loglik.current > loglik.upp[1]) {
                
                loglik.upp <- c(loglik.current, loglik.upp[-1])
                CED.upp <- c(ub[jj], CED.upp[-1])
                
              } else if (loglik.current <= loglik.upp[1]) {
                
                loglik.upp <- c(loglik.current, loglik.upp)
                CED.upp <- c(ub[jj], CED.upp)
                
              }
            }
            
            if (loglik.current == -Inf) {
              
              loglik.upp <- loglik.upp[-1]
              CED.upp <- CED.upp[-1]
              
            }
            
            if (length(loglik.upp) > 1) {
              
              if (abs(loglik.upp[1] - loglik.upp[2]) > 0.25 * crit) 
                step <- sqrt(step)
              if (abs(loglik.upp[1] - loglik.upp[2]) < 0.1 * crit) 
                step <- step^1.26
            }
            
            if (track) {
              
              cat("\n")
              cat(signif(MLE.current, 3), round(loglik.current, 
                      2), round(loglik.max - crit, 2), sep = tb)
              
            }
            
            if (step < 1.00001) {
              
              stop <- TRUE
              small.step.upp <- TRUE
              
            }
            
            if (length(loglik.upp) > 30) 
              if (max(abs(diff(loglik.upp[1:10]))) < 0.001) 
                stop <- TRUE
            
            if (CED.upp[1] > 1e+20) 
              stop <- TRUE
            
            if (ub[jj] > 1e+50) {
              stop <- TRUE
            }
            
            if (loglik.max - loglik.current > crit) {
              if (length(unique(loglik.upp)) < 6) {
                if (track) 
                  cat("\n------- not enough points: ", length(loglik.upp), 
                      " ------------------")
                
                CED.upp <- c(CED.upp, CED.upp.old)
                loglik.upp <- c(loglik.upp, loglik.upp.old)
                CED.upp.old <- CED.upp
                loglik.upp.old <- loglik.upp
                power <- length(CED.upp)/15
                step <- step.start^power
                step.start <- step
                lb[jj] <- MLE[jj]
                start <- dist * MLE
                run <- 0
                count.2 <- count.2 + 1
                
                if (count.2 > 5) {
                  
                  count.2 <- -1
                  stop <- TRUE
                  loglik.upp <- NA
                  
                }
                
              } else if ((loglik.max - loglik.current) > jump.crit) {
                
                lb[jj] <- lb[jj]/step
                step <- step.start^0.5
                step.start <- step
                count <- count + 1
                
                if (step.start < 1.01) {
                  
                  CED.upp[1] <- lb[jj]
                  loglik.upp[1] <- loglik.current
                  count <- count.max
                  
                }
                
                if (count == count.max) 
                  stop <- TRUE
                
                run <- 0
                
              } else {
                
                loglik.diff <- diff(loglik.upp[order(CED.upp)])
                if (sum(loglik.diff < -0.1) > length(CED.upp) - 3) {
                  
                  step.start <- 1 + (step.start - 1)/2
                  step <- step.start
                  run <- 0
                  
                }
                else {
                  
                  stop <- TRUE
                  
                }
              }
            }
            
            if (run > max.runs) {
              
              warning("Maximum number of runs reached in establishing profile")
              if (abs(loglik.upp[1] - loglik.upp[2]) < 0.01) 
                stop <- TRUE
              
            }
            
          }
          
          jump.upp <- FALSE
          if (count == count.max) {
            
            warning("The jump in the log-likelihood could not be mitigated")
            jump.upp <- TRUE
            
          }
          
          
          # Modify lower bounds
          if (!nolog) 
            CED.low <- log10(CED.low)
          
          if (max(abs(diff(loglik.low))) > 0.01) {
            
            loglik.low <- loglik.low[order(CED.low)]
            CED.low <- sort(CED.low)
            spline.low <- spline(CED.low, loglik.low)
            
            if (min(spline.low$x) > min(CED.low)) {
              
              spline.low$x <- c(min(CED.low), spline.low$x)
              spline.low$y <- c(loglik.low[1], spline.low$y)
              
            }
            
            if (length(unique(spline.low$y)) == 1) {
              
              cat("\n\nonly one point in lower spline\n")
              ci.low <- list(y = Inf)
              conf.int[par.nr, 1] <- ci.low$y
              
            } else {
              
              ci.low <- approx(spline.low$y, spline.low$x, 
                  xout = loglik.max - crit)
              conf.int[par.nr, 1] <- ci.low$y
              
            }
            
            if (loglik.low[1] > loglik.max - crit) 
              conf.int[par.nr, 1] <- -Inf
            
            if (small.step.low) 
              conf.int[par.nr, 1] <- (CED.low[1])
            
          } else {
            
            spline.low <- list()
            spline.low$x <- CED.low
            spline.low$y <- rep(loglik.low[1], length(CED.low))
            ci.low <- list()
            ci.low$y <- min(CED.low)
            
          }
          
          # Modify upper bounds
          if (!nolog) 
            CED.upp <- log10(CED.upp)
          
          if (max(abs(diff(loglik.upp))) > 0.01) {
            
            loglik.upp <- loglik.upp[order(CED.upp)]
            CED.upp <- sort(CED.upp)
            spline.upp <- spline(CED.upp, loglik.upp)
            
            if (max(spline.upp$x) < max(CED.upp)) {
              
              spline.upp$x <- c(spline.upp$x, max(CED.upp))
              spline.upp$y <- c(spline.upp$y, loglik.upp[length(loglik.vec)])
              
            }
            
            if (length(unique(spline.upp$y)) == 1) {
              
              cat("\n\nonly one point in upper spline\n")
              ci.upp <- list(y = Inf)
              conf.int[par.nr, 2] <- ci.upp$y
              
            } else {
              
              ci.upp <- approx(spline.upp$y, spline.upp$x, 
                  xout = loglik.max - crit)
              conf.int[par.nr, 2] <- ci.upp$y
              
            }
            
            if (loglik.upp[length(loglik.upp)] > loglik.max - crit) 
              conf.int[par.nr, 2] <- Inf
            
            if (small.step.upp) 
              conf.int[par.nr, 2] <- log10(CED.upp[1])
            
            if (CED.upp.inf) 
              conf.int[par.nr, 2] <- Inf
            
          } else {
            
            spline.upp <- list()
            spline.upp$x <- CED.upp
            spline.upp$y <- rep(loglik.upp[1], length(CED.upp))
            ci.upp <- list()
            ci.upp$y <- max(CED.upp)
            
          }
          
          # Bind results
          CED.vec <- c(CED.low, CED.upp)
          loglik.vec <- c(loglik.low, loglik.upp)
          
          if (jump.low) 
            conf.int[par.nr, 1] <- CED.vec[1]
          if (jump.upp) 
            conf.int[par.nr, 2] <- CED.vec[length(CED.vec)]
          
          
        } # End for loop over groups
        
        
        if (max(abs(diff(loglik.low))) < 0.01) 
          conf.int[par.nr, 1] <- -Inf
        
        if (max(abs(diff(loglik.upp))) < 0.01) 
          conf.int[par.nr, 2] <- Inf
        
        if (count.2 == -1) {
          
          conf.int[par.nr, 1] <- NA
          conf.int[par.nr, 2] <- NA
          
        }
        
        profile.out <- list(conf.int = conf.int, MLE.new = MLE.new, 
            loglik = nlminb.out$loglik, profile = cbind(CED.vec, loglik.vec))
        
        if (track) 
          print("f.profile.all:  END")
        
        return(profile.out)
        
      })
}


#' Determine the number of group ?
#' @param nr.aa 
#' @param nr.bb 
#' @param twice logical, if TRUE two parameters are dependent of the same covariate
#' @return 
#' @export
f.nr.gr <- function (nr.aa, nr.bb, twice) 
{
#	ab.ans <- 1 # why used for?
	nr.gr <- 1
#	if (nr.aa == 1) 
#		ab.ans <- 1
	if (nr.aa > 1 && nr.bb == 1) {
		nr.gr <- nr.aa
	}
	if (nr.bb > 1 && nr.aa == 1) {
		nr.gr <- nr.bb
	}
	if (nr.aa > 1 && nr.bb > 1) {
		if (twice) 
			nr.gr <- nr.aa
		if (!twice) 
			nr.gr <- nr.aa * nr.bb
	}
	return(nr.gr)
}


