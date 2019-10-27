## stats: tests ==============================================================


#### TODO *************************
#### ******************************
#### ******TODO*TODO***************
#### ******xxxxxxxxx***************
#### ******************************

# original:

# https://github.com/nicebread/WRS
# Rand Wilcox,
# http://www.psychology.mcmaster.ca/bennett/boot09/rt2.pdf

#
#  Compute a 1-alpha confidence interval for the difference between
#  the trimmed means corresponding to two independent groups.
#  The bootstrap percentile t method is used.
#
#  The default amount of trimming is tr=.2
#  side=T indicates two-sided method using absolute value of the
#  test statistics within the bootstrap; otherwise the equal-tailed method
#  is used.
#
#  This function uses trimse.
#

# side<-as.logical(side)
# p.value<-NA
# yuenbt<-vector(mode="numeric",length=2)
# if (SEED)set.seed(2) # set seed of random number generator so that
# #             results can be duplicated.
# x<-x[!is.na(x)]  # Remove missing values in x
# y<-y[!is.na(y)]  # Remove missing values in y
# xcen<-x-mean(x,tr)
# ycen<-y-mean(y,tr)
# if (!side) {
#   if (pr)print("NOTE: p-value computed only when side=T")
# }
# test<-(mean(x,tr)-mean(y,tr))/sqrt(trimse(x,tr=tr)^2+trimse(y,tr=tr)^2)
# datax<-matrix(sample(xcen,size=length(x)*nboot,replace=TRUE),nrow=nboot)
# datay<-matrix(sample(ycen,size=length(y)*nboot,replace=TRUE),nrow=nboot)
# top<-apply(datax,1,mean,tr)-apply(datay,1,mean,tr)
# botx<-apply(datax,1,trimse,tr)
# boty<-apply(datay,1,trimse,tr)
# tval<-top/sqrt(botx^2+boty^2)
# if (plotit) {
#   if (op == 1)
#     akerd(tval)
#   if (op == 2)
#     rdplot(tval)
# }
# if (side)tval<-abs(tval)
# tval<-sort(tval)
# icrit<-floor((1-alpha)*nboot+.5)
# ibot<-floor(alpha*nboot/2+.5)
# itop<-floor((1-alpha/2)*nboot+.5)
# se<-sqrt((trimse(x,tr))^2+(trimse(y,tr))^2)
# yuenbt[1]<-mean(x,tr)-mean(y,tr)-tval[itop]*se
# yuenbt[2]<-mean(x,tr)-mean(y,tr)-tval[ibot]*se
# if (side) {
#   yuenbt[1]<-mean(x,tr)-mean(y,tr)-tval[icrit]*se
#   yuenbt[2]<-mean(x,tr)-mean(y,tr)+tval[icrit]*se
#   p.value<-(sum(abs(test)<=abs(tval)))/nboot
# }
# list(ci=yuenbt,test.stat=test,p.value=p.value,est.1=mean(x,tr),est.2=mean(y,tr),est.dif=mean(x,tr)-mean(y,tr),
#      n1=length(x),n2=length(y))



# getAnywhere(t.test.default)
#
# function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
#           mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
#           trim = 0, nboot = 599, na.rm = FALSE
#           ...)

.YuenTTestB <- function(x, y, trim = 0, conf.level = 0.95, nboot=599,
  alternative = c("two.sided", "less", "greater"), mu = 0, na.rm = FALSE) {


  TrimSE <- function(x, trim = 0, na.rm = FALSE) {

    #  Estimate the standard error of the gamma trimmed mean
    #  The default amount of trimming is trim = 0.2

    if (na.rm) x <- na.omit(x)

    winvar <- var(Winsorize(x, probs = c(trim, 1-trim)))

    trimse <- sqrt(winvar) / ((1 - 2 * trim) * sqrt(length(x)))
    trimse
  }


  alternative <- match.arg(alternative)
  method <- "Yuen Two Sample bootstrap t-test"
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  if (na.rm) x <- na.omit(x)
  if (na.rm) y <- na.omit(y)

  meanx <- mean(x, trim = trim)
  meany <- mean(y, trim = trim)

  tstat <- (meanx - meany ) / sqrt(TrimSE(x, trim = trim)^2 + TrimSE(y, trim = trim)^2)

  sampx <- matrix(sample(x - meanx, size=length(x) * nboot, replace=TRUE), nrow=nboot)
  sampy <- matrix(sample(y - meany, size=length(y) * nboot, replace=TRUE), nrow=nboot)

  top <- apply(sampx, 1, mean, trim) - apply(sampy, 1, mean, trim)
  botx <- apply(sampx, 1, TrimSE, trim)
  boty <- apply(sampy, 1, TrimSE, trim)
  tval <- top / sqrt(botx^2 + boty^2)


  alpha <- 1 - conf.level
  se <- sqrt((TrimSE(x, trim = trim))^2 + (TrimSE(y, trim = trim))^2)

  if (alternative == "two.sided") {
    tval <- abs(tval)
    icrit <- floor((1 - alpha) * nboot + .5)
    cint <- meanx - meany + c(-1, 1) * tval[icrit] * se
    pval <- (sum(abs(tstat) <= abs(tval))) / nboot

  } else {
    tval <- sort(tval)
    ibot <- floor(alpha/2 * nboot + .5)
    itop <- floor((1 - alpha/2) * nboot + .5)
    cint <- meanx - meany - tval[c(itop, ibot)] * se

  }

  names(tstat) <- "t"
  names(mu) <- "difference in means"
  estimate <- c(meanx, meany)
  names(estimate) <- c("mean of x", "mean of y")

  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)

}





#' Yuen t-Test For Trimmed Means
#' 
#' Performs one and two sample Yuen t-tests for trimmed means on vectors of
#' data.
#' 
#' 
#' @aliases YuenTTest YuenTTest.default YuenTTest.formula
#' @param x numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite
#' values will be omitted.
#' @param alternative is a character string, one of \code{"greater"},
#' \code{"less"}, or \code{"two.sided"}, or the initial letter of each,
#' indicating the specification of the alternative hypothesis. For one-sample
#' tests, \code{alternative} refers to the true median of the parent population
#' in relation to the hypothesized value of the mean.
#' @param paired a logical indicating whether you want a paired z-test.
#' @param mu a number specifying the hypothesized mean of the population.
#' @param conf.level confidence level for the interval computation.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of x before the mean is computed. Values of trim outside that range are
#' taken as the nearest endpoint.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return An object of class \code{htest} containing the following components:
#' \item{statistic}{the value of the t-statistic.} 
#' \item{parameter}{the degrees
#' of freedom for the t-statistic and the trim percentage used.}
#' \item{p.value}{the p-value for the test.} 
#' \item{conf.int}{a confidence
#' interval for the trimmed mean appropriate to the specified alternative
#' hypothesis.} 
#' \item{estimate}{the estimated trimmed mean or difference in
#' trimmed means depending on whether it was a one-sample test or a two-sample
#' test. } 
#' \item{null.value}{the specified hypothesized value of the trimmed
#' mean or trimmed mean difference depending on whether it was a one-sample
#' test or a two-sample test.} 
#' \item{alternative}{a character string describing
#' the alternative hypothesis.} 
#' \item{method}{a character string indicating
#' what type of test was performed.} 
#' \item{data.name}{a character string giving
#' the name(s) of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>, 
#' based on R-Core code of \code{\link{t.test}}
#' 
#' @seealso 
#' \code{\link{t.test}}, 
#' \code{\link{print.htest}}
#' 
#' @references Wilcox, R. R. (2005) Introduction to robust estimation and
#' hypothesis testing. \emph{Academic Press}.\cr
#' Yuen, K. K. (1974) The
#' two-sample trimmed t for unequal population variances. \emph{Biometrika},
#' 61, 165-170.
#' @keywords htest
#'
#'
#' @examples
#' 
#' x <- rnorm(25, 100, 5)
#' YuenTTest(x, mu=99)
#' 
#' # the classic interface
#' with(sleep, YuenTTest(extra[group == 1], extra[group == 2]))
#' 
#' # the formula interface
#' YuenTTest(extra ~ group, data = sleep)
#' 
#' 
#' # Stahel (2002), pp. 186, 196  
#' d.tyres <- data.frame(A=c(44.5,55,52.5,50.2,45.3,46.1,52.1,50.5,50.6,49.2),
#'                       B=c(44.9,54.8,55.6,55.2,55.6,47.7,53,49.1,52.3,50.7))
#' with(d.tyres, YuenTTest(A, B, paired=TRUE))
#' 
#' 
#' d.oxen <- data.frame(ext=c(2.7,2.7,1.1,3.0,1.9,3.0,3.8,3.8,0.3,1.9,1.9),
#'                      int=c(6.5,5.4,8.1,3.5,0.5,3.8,6.8,4.9,9.5,6.2,4.1))
#' with(d.oxen, YuenTTest(int, ext, paired=FALSE))

# @export
YuenTTest <- function(x, ...) {
  UseMethod("YuenTTest")
}

# @rdname YuenTTest
# @export
YuenTTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || 
      (length(formula) != 3L) || 
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- DoCall("YuenTTest", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L)
    names(y$estimate) <- paste("trimmed mean in group", levels(g))
  y
}

# @rdname YuenTTest
# @export
YuenTTest.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                               mu = 0, paired = FALSE, conf.level = 0.95, trim = 0.2, ...) {

  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired)
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }
    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]

  nx <- length(x)
  mx <- mean(x, trim = trim)
  vx <- var(Winsorize(x, probs = c(trim, 1-trim)))

  if (is.null(y) | paired) {
    if (nx < 2)
      stop("not enough 'x' observations")

    df <- nx - 2 * floor(trim * nx) - 1

    if (paired) {
      my <- mean(y, trim = trim)
      vy <- var(Winsorize(y, probs = c(trim, 1-trim)))
      covxy <- var(
        Winsorize(x, probs = c(trim, 1-trim)), 
        Winsorize(y, probs = c(trim, 1-trim))
      )
      stderr <- sqrt( (nx-1) * (vx + vy - 2 * covxy) / ((df + 1) * df) )
    } else {
      stderr <- sqrt(vx) / ((1 - 2 * trim) * sqrt(nx))
    }

    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")

    if (paired) {
      method <- "Yuen Paired t-test"
      tstat <- (mx - my - mu) / stderr
      estimate <- setNames(mx - my, "difference of trimmed means")

    } else {
      method <- "Yuen One Sample t-test"
      tstat <- (mx - mu)/stderr
      estimate <- setNames(mx, "trimmed mean of x")
    }

  }
  else {
    ny <- length(y)
    if (nx < 2)
      stop("not enough 'x' observations")
    if (ny < 2)
      stop("not enough 'y' observations")
    my <- mean(y, trim = trim)
    vy <- var(Winsorize(y, probs = c(trim, 1-trim)))
    method <- "Yuen Two Sample t-test"
    estimate <- c(mx, my)
    names(estimate) <- c("trimmed mean of x", "trimmed mean of y")

    dfx <- length(x) - 2 * floor(trim * length(x)) - 1
    dfy <- length(y) - 2 * floor(trim * length(y)) - 1

    stderrx <- (length(x) - 1) * vx / ((dfx + 1) * dfx)
    stderry <- (length(y) - 1) * vy / ((dfy + 1) * dfy)

    df <- (stderrx + stderry)^2 / (stderrx^2 / dfx + stderry^2 / dfy)

    stderr <- sqrt(stderrx + stderry)

    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu) / stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(trim) <- "trim"
  names(mu) <- if (paired || !is.null(y))
    "difference in trimmed means"
  else "trimmed mean"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = c(df, trim), p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}





#' Student's t-Test Based on Sample Statistics
#' 
#' Performs one and two sample t-tests based on user supplied summary
#' information instead of data as in \code{t.test()}.
#' 
#' \code{alternative = "greater"} is the alternative that \code{x} has a larger
#' mean than \code{y}.
#' 
#' If \code{paired} is \code{TRUE} then both \code{mx, sx} and \code{my, sy}
#' must be specified and \code{nx} must be equal to \code{ny}.  If
#' \code{var.equal} is \code{TRUE} then the pooled estimate of the variance is
#' used.  By default, if \code{var.equal} is \code{FALSE} then the variance is
#' estimated separately for both groups and the Welch modification to the
#' degrees of freedom is used.
#' 
#' If the input data are effectively constant (compared to the larger of the
#' two means) an error is generated.
#' 
#' @param mx a single number representing the sample mean of x.
#' @param my an optional single number representing the sample mean of y.
#' @param sx a single number representing the sample standard deviation of x.
#' @param sy an optional single number representing the sample standard
#' deviation of y.
#' @param nx a single number representing the sample size of x.
#' @param ny an optional single number representing the sample size of y.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}.  You can specify just the initial letter.
#' @param mu a number indicating the true value of the mean (or difference in
#' means if you are performing a two sample test).
#' @param paired a logical indicating whether you want a paired t-test.
#' @param var.equal a logical variable indicating whether to treat the two
#' variances as being equal. If \code{TRUE} then the pooled variance is used to
#' estimate the variance otherwise the Welch (or Satterthwaite) approximation
#' to the degrees of freedom is used.
#' @param conf.level confidence level of the interval.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \code{"htest"} containing the following
#' components: 
#' \item{statistic}{the value of the t-statistic.}
#' \item{parameter}{the degrees of freedom for the t-statistic.}
#' \item{p.value}{the p-value for the test.} 
#' \item{conf.int}{a confidence
#' interval for the mean appropriate to the specified alternative hypothesis.}
#' \item{estimate}{the estimated mean or difference in means depending on
#' whether it was a one-sample test or a two-sample test.}
#' \item{null.value}{the specified hypothesized value of the mean or mean
#' difference depending on whether it was a one-sample test or a two-sample
#' test.} 
#' \item{alternative}{a character string describing the alternative
#' hypothesis.} 
#' \item{method}{a character string indicating what type of t-test was 
#' performed.} 
#' \item{data.name}{a character string giving the name(s) of the data.}
#' 
#' @seealso 
#' \code{\link{t.test}}
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' ## Classical example: Student's sleep data
#' mx <- 0.75
#' my <- 2.33
#' sx <- 1.789010
#' sy <- 2.002249
#' nx <- ny <- 10
#' TTestA(mx=mx, my=my, sx=sx, sy=sy, nx=nx, ny=ny)
#' 
#' # compare to
#' with(sleep, t.test(extra[group == 1], extra[group == 2]))

# @export

TTestA <- function(mx, sx, nx, my=NULL, sy = NULL, ny=NULL,
                     alternative = c("two.sided", "less", "greater"),
          mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95,
          ...) {

  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")

  if (!is.null(my)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  } else {
    dname <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
  }

  vx <- sx^2

  if (is.null(my)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- if (paired)
      "Paired t-test"
    else "One Sample t-test"
    estimate <- setNames(mx, if (paired)
      "mean of the differences"
      else "mean of x")

  } else {
    # ny <- length(y)
    if (nx < 1 || (!var.equal && nx < 2))
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2))
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3)
      stop("not enough observations")
    # my <- mean(y)
    # vy <- var(y)
    vy <- sy^2
    method <- paste(if (!var.equal)
      "Welch", "Two Sample t-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if (var.equal) {
      df <- nx + ny - 2
      v <- 0
      if (nx > 1)
        v <- v + (nx - 1) * vx
      if (ny > 1)
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
    }
    else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx),
                                                abs(my)))
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- if (paired || !is.null(my))
    "difference in means"
  else "mean"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}





#' Sign Test
#' 
#' Performs one- and two-sample sign tests on vectors of data.
#' 
#' The formula interface is only applicable for the 2-sample test.
#' 
#' \code{SignTest} computes a \dQuote{Dependent-samples Sign-Test} if both
#' \code{x} and \code{y} are provided.  If only \code{x} is provided, the
#' \dQuote{One-sample Sign-Test} will be computed.
#' 
#' For the one-sample sign-test, the null hypothesis is that the median of the
#' population from which \code{x} is drawn is \code{mu}. For the two-sample
#' dependent case, the null hypothesis is that the median for the differences
#' of the populations from which \code{x} and \code{y} are drawn is \code{mu}.
#' The alternative hypothesis indicates the direction of divergence of the
#' population median for \code{x} from \code{mu} (i.e., \code{"greater"},
#' \code{"less"}, \code{"two.sided"}.)
#' 
#' The confidence levels are exact.
#' 
#' @aliases SignTest SignTest.default SignTest.formula
#' @param x numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite
#' values will be omitted.
#' @param mu a number specifying an optional parameter used to form the null
#' hypothesis. See Details.
#' @param alternative is a character string, one of \code{"greater"},
#' \code{"less"}, or \code{"two.sided"}, or the initial letter of each,
#' indicating the specification of the alternative hypothesis. For one-sample
#' tests, \code{alternative} refers to the true median of the parent population
#' in relation to the hypothesized value of the median.
#' @param conf.level confidence level for the returned confidence interval,
#' restricted to lie between zero and one.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ the S-statistic (the number of positive differences
#' between the data and the hypothesized median), with names attribute
#' \dQuote{S}.} 
#' \item{parameter}{ the total number of valid differences.}
#' \item{p.value}{ the p-value for the test.} 
#' \item{null.value}{is the value of
#' the median specified by the null hypothesis. This equals the input argument
#' \code{mu}. } 
#' \item{alternative}{a character string describing the alternative hypothesis.} 
#' \item{method}{ the type of test applied.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{conf.int}{ a confidence interval for the median.} 
#' \item{estimate}{ the sample median.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>
#' 
#' @seealso 
#' \code{\link{t.test}}, 
#' \code{\link{wilcox.test}},
#' \code{\link{ZTest}}, 
#' \code{\link{binom.test}}, 
#' \code{\link[BSDA]{SIGN.test}} in the package \pkg{BSDA} (reporting 
#' approximative confidence intervals).
#' 
#' @references 
#' Gibbons, J.D. and Chakraborti, S. (1992): \emph{Nonparametric
#' Statistical Inference}. Marcel Dekker Inc., New York.
#' 
#' Kitchens, L. J. (2003): \emph{Basic Statistics and Data Analysis}. Duxbury.
#' 
#' Conover, W. J. (1980): \emph{Practical Nonparametric Statistics, 2nd ed}.
#' Wiley, New York.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' x <- c(1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30)
#' y <- c(0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29)
#' 
#' SignTest(x, y)
#' wilcox.test(x, y, paired = TRUE)
#' 
#' 
#' d.light <- data.frame( 
#'   black = c(25.85,28.84,32.05,25.74,20.89,41.05,25.01,24.96,27.47),
#'   white <- c(18.23,20.84,22.96,19.68,19.5,24.98,16.61,16.07,24.59),
#'   d <- c(7.62,8,9.09,6.06,1.39,16.07,8.4,8.89,2.88)
#' )
#' 
#' d <- d.light$d
#' 
#' SignTest(x=d, mu = 4)
#' wilcox.test(x=d, mu = 4, conf.int = TRUE)
#' 
#' SignTest(x=d, mu = 4, alternative="less")
#' wilcox.test(x=d, mu = 4, conf.int = TRUE, alternative="less")
#' 
#' SignTest(x=d, mu = 4, alternative="greater")
#' wilcox.test(x=d, mu = 4, conf.int = TRUE, alternative="greater")
#' 
#' # test die interfaces
#' x <- runif (10)
#' y <- runif (10)
#' g <- rep(1:2, each=10) 
#' xx <- c(x, y)
#' 
#' SignTest(x ~ group, data=data.frame(x=xx, group=g ))
#' SignTest(xx ~ g)
#' SignTest(x, y)
#' 
#' SignTest(x - y)
#' 

# @export

SignTest <- function(x, ...) {
  UseMethod("SignTest")
}

# @rdname SignTest
# @export
SignTest.formula <- function(formula, data, subset, na.action, ...) {

  # this is designed just like wilcox.test.formula

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  names(DATA) <- c("x", "y")
  y <- DoCall("SignTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}

# test:
#  cbind( c(NA,sort(x)), 0:n, dbinom(0:n, size=n, prob=0.5),  pbinom(0:n, size=n, prob=0.5))

# @rdname SignTest
# @export
SignTest.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                             mu = 0, conf.level = 0.95, ...) {

  MedianCI_Binom <- function( x, conf.level = 0.95,
                              alternative = c("two.sided", "less", "greater"), na.rm = FALSE ) {
    # http://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
    # http://de.scribd.com/doc/75941305/Confidence-Interval-for-Median-Based-on-Sign-Test
    if (na.rm) x <- na.omit(x)
    n <- length(x)
    switch( match.arg(alternative)
            , "two.sided" = {
              k <- qbinom(p = (1 - conf.level) / 2, size=n, prob=0.5, lower.tail=TRUE)
              ci <- sort(x)[c(k, n - k + 1)]
              attr(ci, "conf.level") <- 1 - 2 * pbinom(k-1, size=n, prob=0.5)
            }
            , "greater" = {
              k <- qbinom(p = (1 - conf.level), size=n, prob=0.5, lower.tail=TRUE)
              ci <- c(sort(x)[k], Inf)
              attr(ci, "conf.level") <- 1 - pbinom(k-1, size=n, prob=0.5)
            }
            , "less" = {
              k <- qbinom(p = conf.level, size=n, prob=0.5, lower.tail=TRUE)
              ci <- c(-Inf, sort(x)[k])
              attr(ci, "conf.level") <- pbinom(k, size=n, prob=0.5)
            }
    )
    return(ci)
  }

  alternative <- match.arg(alternative)

  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")

  if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
    stop("'conf.level' must be a single number between 0 and 1")

  if (!is.numeric(x))
    stop("'x' must be numeric")

  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")

    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- x[OK]
    y <- y[OK]
    METHOD <- "Dependent-samples Sign-Test"
    x <- (x - y)

  } else {
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
    METHOD <- "One-sample Sign-Test"
  }

  d <- (x - mu)

  # Naive version:
  n.valid <- sum(d > 0) + sum(d < 0)
  if (n.valid > 0) {
    RVAL <- binom.test(x=sum(d > 0), n=n.valid, p=0.5, alternative = alternative, conf.level = conf.level )
  } else {
    RVAL <- binom.test(x=1, n=1)
  }

  RVAL$method <- METHOD
  RVAL$data.name <- DNAME
  names(mu) <- if (!is.null(y)) "median difference" else "median"

  names(RVAL$statistic) <- "S"
  RVAL$estimate <- median(d + mu, na.rm=TRUE)
  names(RVAL$parameter) <- "number of differences"
  mci <- MedianCI_Binom(d + mu, conf.level=conf.level, alternative=alternative, na.rm=TRUE)
  RVAL$conf.int <- mci
  attr(RVAL$conf.int, "conf.level") = round(attr(mci,"conf.level"), 3)

  names(RVAL$estimate) <- "median of the differences"
  RVAL$null.value <- mu
  class(RVAL) <- "htest"
  return(RVAL)

}





#' Z Test for Known Population Standard Deviation
#' 
#' Compute the test of hypothesis and compute confidence interval on the mean
#' of a population when the standard deviation of the population is known.
#' 
#' Most introductory statistical texts introduce inference by using the z-test
#' and z-based confidence intervals based on knowing the population standard
#' deviation. However statistical packages often do not include functions to do
#' z-tests since the t-test is usually more appropriate for real world
#' situations. This function is meant to be used during that short period of
#' learning when the student is learning about inference using z-procedures,
#' but has not learned the t-based procedures yet.  Once the student has
#' learned about the t-distribution the \code{t.test()} function should be used
#' instead of this one (but the syntax is very similar, so this function should
#' be an appropriate introductory step to learning \code{t.test()}).
#' 
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' @aliases ZTest ZTest.default ZTest.formula
#' @param x numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite
#' values will be omitted.
#' @param mu a number specifying the hypothesized mean of the population.
#' @param sd_pop a number specifying the known standard deviation of the
#' population.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}.  You can specify just the initial letter. \cr
#' For one-sample
#' tests, \code{alternative} refers to the true mean of the parent population
#' in relation to the hypothesized value of the mean.
#' @param paired a logical indicating whether you want a paired z-test.
#' @param conf.level confidence level for the interval computation.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} a factor with two levels giving the
#' corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class "\code{htest}" containing the following
#' components:
#' \item{statistic}{ the value of the z-statistic.} 
#' \item{p.value}{the p-value for the test} 
#' \item{conf.int}{a confidence interval for the mean
#' appropriate to the specified alternative hypothesis.}
#' \item{estimate}{the estimated mean or difference in means depending on 
#' whether it was a one-sample test or a two-sample test.} 
#' \item{null.value}{the specified
#' hypothesized value of the mean or mean difference depending on whether it
#' was a one-sample test or a two-sample test.} 
#' \item{alternative}{a character
#' string describing the alternative hypothesis.} 
#' \item{method}{ a character
#' string indicating what type of test was performed.} 
#' \item{data.name}{a character string giving the name(s) of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>, 
#' based on R-Core code of \code{\link{t.test}},\cr
#' documentation partly from Greg Snow <greg.snow@@imail.org>
#' 
#' @seealso 
#' \code{\link{t.test}}, 
#' \code{\link{print.htest}}
#' 
#' @references 
#' Stahel, W. (2002) \emph{Statistische Datenanalyse, 4th ed}, vieweg
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' x <- rnorm(25, 100, 5)
#' ZTest(x, mu=99, sd_pop=5)
#' 
#' # the classic interface
#' with(sleep, ZTest(extra[group == 1], extra[group == 2], sd_pop=2))
#' 
#' # the formula interface
#' ZTest(extra ~ group, data = sleep, sd_pop=2)
#' 
#' 
#' # Stahel (2002), pp. 186, 196
#' 
#' d.tyres <- data.frame(A=c(44.5,55,52.5,50.2,45.3,46.1,52.1,50.5,50.6,49.2),
#'                       B=c(44.9,54.8,55.6,55.2,55.6,47.7,53,49.1,52.3,50.7))
#' with(d.tyres, ZTest(A, B, sd_pop=3, paired=TRUE))
#' 
#' 
#' d.oxen <- data.frame(ext=c(2.7,2.7,1.1,3.0,1.9,3.0,3.8,3.8,0.3,1.9,1.9),
#'                      int=c(6.5,5.4,8.1,3.5,0.5,3.8,6.8,4.9,9.5,6.2,4.1))
#' with(d.oxen, ZTest(int, ext, sd_pop=1.8, paired=FALSE))
#' 
ZTest <- function(x, ...) {
  UseMethod("ZTest")
}

ZTest.formula <- function(formula, data, subset, na.action, ...)  {

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- DoCall("ZTest", c(DATA, list(...)))
  y$data.name <- DNAME
  if (length(y$estimate) == 2L)
    names(y$estimate) <- paste("mean in group", levels(g))
  y
}


ZTest.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
                           paired = FALSE, mu = 0, sd_pop, conf.level = 0.95,  ...)  {

  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                               conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (paired)
      xok <- yok <- complete.cases(x, y)
    else {
      yok <- !is.na(y)
      xok <- !is.na(x)
    }

    y <- y[yok]
  }
  else {
    dname <- deparse(substitute(x))
    if (paired)
      stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
  }
  x <- x[xok]

  if (paired) {
    x <- x - y
    y <- NULL
  }

  nx <- length(x)
  mx <- mean(x)
  # vx <- sd_pop^2

  if (is.null(y)) {
    if (nx < 2)
      stop("not enough 'x' observations")
    stderr <- sqrt(sd_pop^2/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx))
      stop("data are essentially constant")
    zstat <- (mx - mu)/stderr

    method <- method <- if (paired)
      "Paired z-test" else "One Sample z-test"
    estimate <- setNames(mx, if (paired)
      "mean of the differences"
      else "mean of x")
  }
  else {
    ny <- length(y)
    if (nx < 1)
      stop("not enough 'x' observations")
    if (ny < 1)
      stop("not enough 'y' observations")
    if (nx + ny < 3)
      stop("not enough observations")
    my <- mean(y)

    method <- paste("Two Sample z-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")

    stderr <- sqrt(sd_pop^2 * (1/nx + 1/ny))

    if (stderr < 10 * .Machine$double.eps * max(abs(mx),
                                                abs(my)))
      stop("data are essentially constant")
    zstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pnorm(zstat)
    cint <- c(-Inf, zstat + qnorm(conf.level))
  }
  else if (alternative == "greater") {
    pval <- pnorm(zstat, lower.tail = FALSE)
    cint <- c(zstat - qnorm(conf.level), Inf)
  }
  else {
    pval <- 2 * pnorm(-abs(zstat))
    alpha <- 1 - conf.level
    cint <- qnorm(1 - alpha/2)
    cint <- zstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(zstat) <- "z"
  names(mu) <- if (paired || !is.null(y))
    "difference in means"
  else "mean"
  names(sd_pop) <- "Std. Dev. Population"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = zstat, p.value = pval,
               parameter = sd_pop,
               conf.int = cint, estimate = estimate, null.value = mu,
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}







#' ChiSquare Test for One Variance and F Test to Compare Two Variances
#' 
#' Performs either a one sample ChiSquare test to compare the variance of a
#' vector with a given value or an F test to compare the variances of two
#' samples from normal populations.
#' 
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' The null hypothesis is that the ratio of the variances of the populations
#' from which \code{x} and \code{y} were drawn, or in the data to which the
#' linear models \code{x} and \code{y} were fitted, is equal to \code{ratio}.
#' 
#' @aliases VarTest VarTest.default VarTest.formula
#' @param x,y numeric vectors of data values.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}.  You can specify just the initial letter.
#' @param ratio the hypothesized ratio of the population variances of \code{x}
#' and \code{y}.
#' @param sigma.squared a number indicating the true value of the variance, if
#' one sample test is requested.
#' @param conf.level confidence level for the returned confidence interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} is a
#' numeric variable giving the data values and \code{rhs} a factor with two
#' levels giving the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return A list with class \code{"htest"} containing the following
#' components: 
#' \item{statistic}{the value of the F test statistic.}
#' \item{parameter}{the degrees of the freedom of the F distribution of the
#' test statistic.} 
#' \item{p.value}{the p-value of the test.} 
#' \item{conf.int}{a
#' confidence interval for the ratio of the population variances.}
#' \item{estimate}{the ratio of the sample variances of \code{x} and \code{y}.}
#' \item{null.value}{the ratio of population variances under the null.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.} 
#' \item{method}{the character string \code{"F test to compare two variances"}.} 
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> (One sample test)\cr
#' Two Sample test and help text from R-Core.
#' 
#' @seealso 
#' \code{\link{var.test}}, 
#' \code{\link{bartlett.test}} for testing homogeneity of variances in more than
#'  two samples from normal distributions;
#' \code{\link{ansari.test}} and 
#' \code{\link{mood.test}} for two rank based (nonparametric) two-sample tests 
#' for difference in scale.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' x <- rnorm(50, mean = 0, sd = 2)
#' 
#' # One sample test
#' VarTest(x, sigma.squared = 2.5)
#' 
#' # two samples
#' y <- rnorm(30, mean = 1, sd = 1)
#' VarTest(x, y)                  # Do x and y have the same variance?
#' VarTest(lm(x ~ 1), lm(y ~ 1))  # The same.
#' 
VarTest <- function(x, ...) {
  UseMethod("VarTest")
}


VarTest.default <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), ratio = 1,
                             sigma.squared = 1, conf.level = 0.95, ...) {

  if (is.null(y)) {
    # perform a one sample variance test

    alternative <- match.arg(alternative)
    if (!missing(sigma.squared) && (length(sigma.squared) != 1 || is.na(sigma.squared)))
      stop("'sigma.squared' must be a single number")

    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                                 conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")

    dname <- deparse(substitute(x))
    x <- na.omit(x)

    nx <- length(x)
    if (nx < 2)
      stop("not enough 'x' observations")
    df <- nx - 1
    vx <- var(x)

    xstat <- vx * df / sigma.squared

    method <- "One Sample Chi-Square test on variance"
    estimate <- vx

    if (alternative == "less") {
      pval <- pchisq(xstat, df)
      cint <- c(0, df * vx/qchisq((1 - conf.level), df))

    } else if (alternative == "greater") {
      pval <- pchisq(xstat, df, lower.tail = FALSE)
      cint <- c(df * vx/qchisq((1 - conf.level), df, lower.tail = FALSE), Inf)
    } else {
      pval <- 2 * min(pchisq(xstat, df), pchisq(xstat, df, lower.tail = FALSE))
      alpha <- 1 - conf.level
      cint <- qt(1 - alpha/2, df)
      cint <- xstat + c(-cint, cint)
      cint <- df * vx / c(qchisq((1 - conf.level)/2, df, lower.tail = FALSE),
                          qchisq((1 - conf.level)/2, df))
    }

    names(xstat) <- "X-squared"
    names(df) <- "df"
    names(sigma.squared) <- "variance"
    names(estimate) <- "variance of x"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = xstat, parameter = df, p.value = pval,
                 conf.int = cint, estimate = estimate, null.value = sigma.squared,
                 alternative = alternative, method = method, data.name = dname)
    class(rval) <- "htest"

    return(rval)

  } else {
    # perform a normal F-test
    var.test(x=x, y=y, ratio=ratio, alternative=alternative, conf.level=conf.level)
  }

}


VarTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  y <- do.call("VarTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}





# moved from Rcmdr 13 July 2004

# levene.test.default function slightly modified and generalized from Brian Ripley via R-help
# the original generic version was contributed by Derek Ogle
# last modified 28 January 2010 by J. Fox



#' Levene's Test for Homogeneity of Variance
#' 
#' Computes Levene's test for homogeneity of variance across groups.
#' 
#' 
#' @aliases LeveneTest LeveneTest.formula LeveneTest.lm LeveneTest.default
#' @param y response variable for the default method, or a \code{lm} or
#' \code{formula} object. If \code{y} is a linear-model object or a formula,
#' the variables on the right-hand-side of the model must all be factors and
#' must be completely crossed.
#' @param group factor defining groups.
#' @param center The name of a function to compute the center of each group;
#' \code{mean} gives the original Levene's test; the default, \code{median},
#' provides a more robust test (Brown-Forsythe-Test).
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param ... arguments to be passed down, e.g., \code{data} for the
#' \code{formula} and \code{lm} methods; can also be used to pass arguments to
#' the function given by \code{center} (e.g., \code{center=mean} and
#' \code{trim=0.1} specify the 10\% trimmed mean).
#' @return returns an object meant to be printed showing the results of the
#' test.
#' @note This function was previously published as leveneTest() in the
#' library(car) and has been integrated here without logical changes.
#' 
#' @author 
#' John Fox \email{jfox@@mcmaster.ca}; 
#' original generic version contributed by Derek Ogle\cr
#' adapted from a response posted by Brian Ripley to the r-help email list.
#' 
#' @seealso 
#' \code{\link{fligner.test}} for a rank-based (nonparametric)
#' \eqn{k}-sample test for homogeneity of variances; 
#' \code{\link{mood.test}} for another rank-based two-sample test for a 
#' difference in scale parameters;
#' \code{\link{var.test}} and 
#' \code{\link{bartlett.test}} for parametric tests for the homogeneity in
#' variance.
#' 
#' \code{\link[coin:ScaleTests]{ansari_test}} in package coin for exact and
#' approximate \emph{conditional} p-values for the Ansari-Bradley test, as well
#' as different methods for handling ties.
#' 
#' @references 
#' Fox, J. (2008) \emph{Applied Regression Analysis and Generalized
#' Linear Models}, Second Edition. Sage.
#' 
#' Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Regression},
#' Second Edition, Sage.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## example from ansari.test:
#' ## Hollander & Wolfe (1973, p. 86f):
#' ## Serum iron determination using Hyland control sera
#' ramsay <- c(111, 107, 100, 99, 102, 106, 109, 108, 104, 99,
#'             101, 96, 97, 102, 107, 113, 116, 113, 110, 98)
#' jung.parekh <- c(107, 108, 106, 98, 105, 103, 110, 105, 104,
#'             100, 96, 108, 103, 104, 114, 114, 113, 108, 106, 99)
#' 
#' LeveneTest( c(ramsay, jung.parekh),
#'   factor(c(rep("ramsay",length(ramsay)), rep("jung.parekh",length(jung.parekh)))))
#' 
#' LeveneTest( c(rnorm(10), rnorm(10, 0, 2)), factor(rep(c("A","B"),each=10)) )
#' 
#' \dontrun{
#' # original example from package car
#' 
#' with(Moore, LeveneTest(conformity, fcategory))
#' with(Moore, LeveneTest(conformity, interaction(fcategory, partner.status)))
#' 
#' LeveneTest(conformity ~ fcategory * partner.status, data = Moore)
#' LeveneTest(conformity ~ fcategory * partner.status, data = Moore, center = mean)
#' LeveneTest(conformity ~ fcategory * partner.status, data = Moore, center = mean, trim = 0.1)
#' 
#' LeveneTest(lm(conformity ~ fcategory*partner.status, data = Moore))
#' }
#' 
LeveneTest <- function(y, ...) {
  UseMethod("LeveneTest")
}

LeveneTest.default <- function(y, group, center=median, ...) { # original levene.test

  if (!is.numeric(y))
    stop(deparse(substitute(y)), " is not a numeric variable")

  if (!is.factor(group)) {
    warning(deparse(substitute(group)), " coerced to factor.")
    group <- as.factor(group)
  }

  valid <- complete.cases(y, group)
  meds <- tapply(y[valid], group[valid], center, ...)
  resp <- abs(y - meds[group])
  table <- anova(lm(resp ~ group))[, c(1, 4, 5)]
  rownames(table)[2] <- " "
  dots <- deparse(substitute(...))

  attr(table, "heading") <- paste("Levene's Test for Homogeneity of Variance (center = ",
                                  deparse(substitute(center)), if (!(dots == "NULL")) paste(":", dots),  ")", sep="")
  table
}


LeveneTest.formula <- function(formula, data, ...) {
  form <- formula
  mf <- if (missing(data)) model.frame(form) else model.frame(form, data)
  if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]]))))
    stop("Levene's test is not appropriate with quantitative explanatory variables.")
  y <- mf[,1]
  if (dim(mf)[2]==2) group <- mf[,2]
  else {
    if (length(grep("\\+ | \\| | \\^ | \\:",form))>0) stop("Model must be completely crossed formula only.")
    group <- interaction(mf[,2:dim(mf)[2]])
  }
  LeveneTest.default(y = y, group=group, ...)
}


# LeveneTest.formula <- function(formula, data, subset, na.action, ...) {
#
#   # replaced as the original did not support subsets
#
#   if (missing(formula) || (length(formula) != 3L))
#     stop("'formula' missing or incorrect")
#   m <- match.call(expand.dots = FALSE)
#   if (is.matrix(eval(m$data, parent.frame())))
#     m$data <- as.data.frame(data)
#   m[[1L]] <- quote(stats::model.frame)
#   mf <- eval(m, parent.frame())
#
#   if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]]))))
#     stop("Levene's test is not appropriate with quantitative explanatory variables.")
#
#   # from kruskal not to be applied here
#   # if (length(mf) > 2L)
#   #   stop("'formula' should be of the form response ~ group")
#
#   if (dim(mf)[2]==2)
#     group <- mf[, 2]
#   else {
#     if (length(grep("\\+ | \\| | \\^ | \\:", formula)) > 0)
#       stop("Model must be completely crossed formula only.")
#     group <- interaction(mf[, 2:dim(mf)[2]])
#   }
#
#   DNAME <- paste(names(mf), collapse = " by ")
#   names(mf) <- NULL
#   # y <- do.call("LeveneTest", as.list(mf))
#   LeveneTest.default(y=mf[, 1], group=group, ...)
#   # y$data.name <- DNAME
#   # y
# }



LeveneTest.lm <- function(y, ...) {
  LeveneTest.formula(formula(y), data=model.frame(y), ...)
}






#' Runs Test for Randomness
#' 
#' Performs a test whether the elements of \code{x} are serially independent -
#' say, whether they occur in a random order - by counting how many runs there
#' are above and below a threshold. If \code{y} is supplied a two sample
#' Wald-Wolfowitz-Test testing the equality of two distributions against
#' general alternatives will be computed.
#' 
#' The runs test for randomness is used to test the hypothesis that a series of
#' numbers is random. The 2-sample test is known as the Wald-Wolfowitz test.
#' \cr
#' 
#' For a categorical variable, the number of runs correspond to the number of
#' times the category changes, that is, where \eqn{x_{i}}{x_i} belongs to one
#' category and \eqn{x_{i+1}}{x_(i+1)} belongs to the other. The number of runs
#' is the number of sign changes plus one.\cr
#' 
#' For a numeric variable x containing more than two values, a run is a set of
#' sequential values that are either all above or below a specified cutpoint,
#' typically the median. This is not necessarily the best choice. If another
#' threshold should be used use a code like: \code{RunsTest(x > mean(x))}.
#' 
#' The exact distribution of runs and the p-value based on it are described in
#' the manual of SPSS "Exact tests"
#' \url{http://www.sussex.ac.uk/its/pdfs/SPSS_Exact_Tests_21.pdf}.
#' 
#' The normal approximation of the runs test is calculated with the expected
#' number of runs under the null \deqn{\mu_r=\frac{2 n_0 n_1}{n_0 + n_1} + 1}
#' and its variance \deqn{\sigma_r^2 = \frac{2 n_0 n_1 (2 n_0 n_1 - n_0 - n_1)
#' }{(n_0 + n_1)^2 \cdot (n_0 + n_1 - 1)}} as \deqn{\hat{z}=\frac{r - \mu_r +
#' c}{\sigma_r}} where \eqn{n_0, n_1} the number of values below/above the
#' threshold and \eqn{r} the number of runs.
#' 
#' Setting the continuity correction \code{correct = TRUE} will yield the
#' normal approximation as SAS (and SPSS if n < 50) does it, see
#' \url{http://support.sas.com/kb/33/092.html}. The c is set to \eqn{c = 0.5}
#' if \eqn{r < \frac{2 n_0 n_1}{n_0 + n_1} + 1} and to \eqn{c = -0.5} if \eqn{r
#' > \frac{2 n_0 n_1}{n_0 + n_1} + 1}.
#' 
#' \bold{Wald-Wolfowitz with Ties.} Ideally there should be no ties in the data
#' used for the Wald-Wolfowitz test. In practice there is no problem with ties
#' within a group, but if ties occur between members of the different groups
#' then there is no unique sequence of observations. For example the data sets
#' A: 10,14,17,19,34 and B: 12,13,17,19,22 can give four possible sequences,
#' with two possible values for r (7 or 9). The "solution" to this is to list
#' every possible combination, and calculate the test statistic for each one.
#' If all test statistics are significant at the chosen level, then one can
#' reject the null hypothesis. If only some are significant, then Siegel (1956)
#' suggests that the average of the P-values is taken. Help for finding all
#' permutations of ties can be found at:
#' \url{https://stackoverflow.com/questions/47565066/all-possible-permutations-in-factor-variable-when-ties-exist-in-r}
#' 
#' However this solutions seems quite coarse and in general, the test should
#' not be used if there are more than one or two ties. We have better tests to
#' distinguish between two samples!
#' 
#' @aliases RunsTest RunsTest.formula RunsTest.default
#' @param x a dichotomous vector of data values or a (non-empty) numeric vector
#' of data values.
#' @param y an optional (non-empty) numeric vector of data values.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"less"} or
#' \code{"greater"}.
#' @param exact a logical indicating whether an exact p-value should be
#' computed. By default exact values will be calculated for small vectors with
#' a total length <= 30 and the normal approximation for longer ones.
#' @param correct a logical indicating whether to apply continuity correction
#' when computing the test statistic. Default is \code{TRUE}. Ignored if
#' \code{exact} is set to \code{TRUE}. See details.
#' @param na.rm defines if \code{NA}s should be omitted. Default is
#' \code{FALSE}.
#' @param \dots further arguments to be passed to or from methods.
#' 
#' @return A list with the following components.  
#' \item{statistic}{z, the value of the standardized runs statistic, if not 
#' exact p-values are computed.}
#' \item{parameter}{the number of runs, the total number of zeros (m) and ones
#' (n)} 
#' \item{p.value}{the p-value for the test.} 
#' \item{data.name}{a character string giving the names of the data.} 
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>, 
#' exact p-values by Detlew Labes <detlewlabes@@gmx.de>
#' 
#' @seealso 
#' Run Length Encoding \code{\link{rle}}
#' 
#' @references 
#' 
#' Wackerly, D., Mendenhall, W. Scheaffer, R. L. (1986)
#' \emph{Mathematical Statistics with Applications}, 3rd Ed., Duxbury Press,
#' CA.
#' 
#' Wald, A. and Wolfowitz, J. (1940): On a test whether two samples are from
#' the same population, \emph{Ann. Math Statist}. 11, 147-162.
#' 
#' Siegel, S. (1956) \emph{Nonparametric Statistics for the Behavioural
#' Sciences}, McGraw-Hill Kogakusha, Tokyo.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' # x will be coerced to a dichotomous variable
#' x <- c("S","S", "T", "S", "T","T","T", "S", "T")
#' RunsTest(x)
#' 
#' 
#' x <- c(13, 3, 14, 14, 1, 14, 3, 8, 14, 17, 9, 14, 13, 2, 16, 1, 3, 12, 13, 14)
#' RunsTest(x)
#' # this will be treated as
#' RunsTest(x > median(x))
#' 
#' plot( (x < median(x)) - 0.5, type="s", ylim=c(-1,1) )
#' abline(h=0)
#' 
#' set.seed(123)
#' x <- sample(0:1, size=100, replace=TRUE)
#' RunsTest(x)
#' # As you would expect of values from a random number generator, the test fails to reject
#' # the null hypothesis that the data are random.
#' 
#' 
#' # SPSS example
#' x <- c(31,23,36,43,51,44,12,26,43,75,2,3,15,18,78,24,13,27,86,61,13,7,6,8)
#' RunsTest(x, exact=TRUE)       # exact probability
#' RunsTest(x, exact=FALSE)      # normal approximation
#' 
#' # SPSS example small dataset
#' x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)
#' RunsTest(x)
#' RunsTest(x, exact=FALSE)
#' 
#' # if y is not NULL, the Wald-Wolfowitz-Test will be performed
#' A <- c(35,44,39,50,48,29,60,75,49,66)
#' B <- c(17,23,13,24,33,21,18,16,32)
#' 
#' RunsTest(A, B, exact=TRUE)
#' RunsTest(A, B, exact=FALSE)
#' 
RunsTest <- function(x, ...) {
  UseMethod("RunsTest")
}

RunsTest.formula <- function(formula, data, subset, na.action, ...) {

  # this is a taken analogue to wilcox.test.formula

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  names(DATA) <- c("x", "y")
  y <- DoCall("RunsTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}


RunsTest.default <- function(x, y=NULL, alternative=c("two.sided", "less", "greater"), exact=NULL, correct=TRUE, na.rm = FALSE, ...) {

  # exact values:
  # http://www.reiter1.com/Glossar/Wald_Wolfowitz.htm

  # example:   x <- sample(c(0,1), size=20, r=TRUE)

  pruns <- function(r, n1, n2, alternative=c("two.sided", "less", "greater")) {

    # source: randomizeBE
    # author: D. Labes <detlewlabes at gmx.de>

    # function for calculating the denominator of the runs distribution
    .druns_nom <- function(r, n1, n2) {
      pp <- vector(mode="numeric",length=length(r))
      for (i in seq_along(r)) {
        if (2*r[i]%/%2==r[i]) {
          # even 2*k
          k <- r[i]/2
          pp[i] <- 2*choose(n1-1, k-1)*choose(n2-1, k-1)
        } else {
          # odd 2*k+1
          k <- (r[i]-1)/2
          pp[i] <- choose(n1-1,k-1) * choose(n2-1,k) +
            choose(n1-1,k)   * choose(n2-1,k-1)
        }
      }
      pp
    }

    alternative <- match.arg(alternative)

    n <- n1+n2

    if (r<=1) stop("Number of runs must be > 1")
    if (r>n) stop("Number of runs must be < (n1+n2")
    if (n1<1 | n2<1) return(0) #??? is not random!

    E <- 1 + 2*n1*n2/n

    denom <- choose(n,n1)
    # how long should we make the r vector?
    # in very unsymmetric cases only a few elements of
    # pp = density have values > 0 if rmax=n1+n2
    # number of runs possible: 2*m if n=m, 2*m+1 if m<n
    rmax <- ifelse(n1==n2, 2*n1, 2*min(n1,n2)+1)
    rv <- 2:rmax
    pp <- .druns_nom(rv, n1, n2)

    # pL is p(R<=r) -> left/lower tail
    pL <- sum(pp[rv<=r])/denom
    #pU is p(R>=r) -> right/upper tail
    pU <- 1 - sum(pp[rv<=(r-1)])/denom

    # Equn. 4.7 of the SPSS documentation
    p2 <- sum(pp[abs(rv-E)>=abs(r-E)])/denom

    # Next is the rule from:
    # Gibbons "Nonparametric Methods for Quantitative Analysis"
    # 0.5 is to avoid p>1 if both pL and pU are >0.5
    p2min <- 2*min(c(pL, pU, 0.5))

    # we are using the SPSS approach wich takes into account the
    # unsymmetric form of the distribution if n1 << n2

    return(
      switch( alternative
              , "less" = pL
              , "greater" = pU
              , "two.sided" = p2
      )
    )

  }



  if (!is.null(y)) {
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    # perform Wald-Wolfowitz-Test with 2 variables
    xy <- Sort(cbind(c(x,y), c(rep(0, length(x)), rep(1, length(y)))))[,2]

    TIES <- length(unique(u <- c(unique(x), unique(y)))) != length(u)
    rm(u)

    res <- RunsTest(x=xy, alternative=alternative, exact=exact, na.rm=na.rm)

    if (TIES)
      warning("cannot compute reliable p-values with inter-group ties between x and y")

    res$data.name <- dname
    res$method <- "Wald-Wolfowitz Runs Test "
    return(res)
  }

  alternative <- match.arg(alternative)
  dname <- deparse(substitute(x))

  # Strip NAs
  if (na.rm) x <- na.omit(x)

  # let's have a 0,1 vector if x is a numeric vector with more than 2 values
  if (is.numeric(x) & (length(unique(x))>2)) {
    est <- median(x, na.rm=TRUE)
    names(est) <- "median(x)"
    x <- ((x > est)*1)

  } else {
    est <- NULL
  }

  x <- factor(x)
  if ( nlevels(x) %nin% c(1,2) ) stop("Can only process dichotomous variables")
  x <- as.numeric(x) - 1

  # x <- sample(c(0,1), 100000000, replace=TRUE)
  # ### user  system elapsed
  #   9.39    6.76   16.30    system.time(rle(x))
  #   4.49    3.45    8.00    system.time(sum(diff(x) != 0) + 1)
  # so don't use rle...

  runs <- sum(diff(x) != 0) + 1
  m <- sum(x==0)
  n <- sum(x==1)

  if (is.null(exact)) { exact <- ((m +n) <= 30) }

  E <- 1 + 2*n*m / (n + m)
  s2 <- (2*n*m * (2*n*m - n - m)) / ((n + m)^2 * (n + m - 1))

  # this is the SPSS-Definition
  # http://publib.boulder.ibm.com/infocenter/spssstat/v20r0m0/index.jsp?topic=%2Fcom.ibm.spss.statistics.help%2Fidh_idd_npar_onesample_settings_tests_runs.htm
  # if ( n+m >= 50) {
  if (correct) {
    switch( as.character(cut(runs - E, breaks=c(-Inf, -0.5, 0.5, Inf), labels=c("a", "b", "c")))
            , "a" = statistic <- (runs - E + 0.5) / sqrt(s2)
            , "b" = statistic <- 0
            , "c" = statistic <- (runs - E - 0.5) / sqrt(s2)
    )
  } else {
    statistic <- (runs - E) / sqrt(s2)
  }

  switch( alternative
          , "less" = {
            p.value <- ifelse(exact, pruns(runs, m, n, alternative="less"), pnorm(statistic) )
            alternative <- "true number of runs is less than expected"
          }
          , "greater" = {
            p.value = ifelse(exact, pruns(runs, m, n, alternative="greater"), 1 - pnorm(statistic) )
            alternative <- "true number of runs is greater than expected"
          }
          , "two.sided" = {
            p.value = ifelse(exact, pruns(runs, m, n, alternative="two.sided"),
                             2 * min(pnorm(statistic), 1 - pnorm(statistic)) )
            alternative <- "true number of runs is not equal the expected number"
          }
  )

  method = "Runs Test for Randomness"
  names(statistic) <- "z"  # Standardized Runs Statistic

  # do not report statistic when exact p-value is calculated
  if (exact) statistic <- NULL

  structure(list(
    statistic = statistic,
    p.value = p.value,
    method = method,
    alternative = alternative,
    data.name = dname,
    estimate = est,
    parameter = c(runs=runs, m=m, n=n)),
    class = "htest")

}





#' Durbin-Watson Test
#' 
#' Performs the Durbin-Watson test for autocorrelation of disturbances.
#' 
#' The Durbin-Watson test has the null hypothesis that the autocorrelation of
#' the disturbances is 0. It is possible to test against the alternative that
#' it is greater than, not equal to, or less than 0, respectively. This can be
#' specified by the \code{alternative} argument.
#' 
#' Under the assumption of normally distributed disturbances, the null
#' distribution of the Durbin-Watson statistic is the distribution of a linear
#' combination of chi-squared variables. The p-value is computed using the
#' Fortran version of Applied Statistics Algorithm AS 153 by Farebrother (1980,
#' 1984). This algorithm is called "pan" or "gradsol". For large sample sizes
#' the algorithm might fail to compute the p value; in that case a warning is
#' printed and an approximate p value will be given; this p value is computed
#' using a normal approximation with mean and variance of the Durbin-Watson
#' test statistic.
#' 
#' Examples can not only be found on this page, but also on the help pages of
#' the data sets \code{\link[lmtest]{bondyield}},
#' \code{\link[lmtest]{currencysubstitution}},
#' \code{\link[lmtest]{growthofmoney}}, \code{\link[lmtest]{moneydemand}},
#' \code{\link[lmtest]{unemployment}}, \code{\link[lmtest]{wages}}.
#' 
#' For an overview on R and econometrics see Racine & Hyndman (2002).
#' 
#' @param formula a symbolic description for the model to be tested (or a
#' fitted \code{"lm"} object).
#' @param order.by Either a vector \code{z} or a formula with a single
#' explanatory variable like \code{~ z}. The observations in the model are
#' ordered by the size of \code{z}. If set to \code{NULL} (the default) the
#' observations are assumed to be ordered (e.g., a time series).
#' @param alternative a character string specifying the alternative hypothesis.
#' @param iterations an integer specifying the number of iterations when
#' calculating the p-value with the "pan" algorithm.
#' @param exact logical. If set to \code{FALSE} a normal approximation will be
#' used to compute the p value, if \code{TRUE} the "pan" algorithm is used. The
#' default is to use "pan" if the sample size is < 100.
#' @param tol tolerance. Eigenvalues computed have to be greater than
#' \code{tol} to be treated as non-zero.
#' @param data an optional data frame containing the variables in the model.
#' By default the variables are taken from the environment which
#' \code{DurbinWatsonTest} is called from.
#' 
#' @return An object of class \code{"htest"} containing: 
#' \item{statistic}{the test statistic.} 
#' \item{p.value}{the corresponding p-value.} 
#' \item{method}{a character string with the method used.} 
#' \item{data.name}{a character string with the data name.}
#' 
#' @note This function was previously published as \code{dwtest} in the
#' 
#' \pkg{lmtest} package and has been integrated here without logical changes.
#' 
#' @author 
#' Torsten Hothorn, Achim Zeileis, Richard W. Farebrother (pan.f),
#' Clint Cummins (pan.f), Giovanni Millo, David Mitchell
#' 
#' @seealso 
#' \code{\link{lm}}
#' 
#' @references
#' 
#' J. Durbin & G.S. Watson (1950), Testing for Serial Correlation in Least
#' Squares Regression I. \emph{Biometrika} \bold{37}, 409--428.
#' 
#' J. Durbin & G.S. Watson (1951), Testing for Serial Correlation in Least
#' Squares Regression II. \emph{Biometrika} \bold{38}, 159--178.
#' 
#' J. Durbin & G.S. Watson (1971), Testing for Serial Correlation in Least
#' Squares Regression III. \emph{Biometrika} \bold{58}, 1--19.
#' 
#' R.W. Farebrother (1980), Pan's Procedure for the Tail Probabilities of the
#' Durbin-Watson Statistic (Corr: 81V30 p189; AS R52: 84V33 p363- 366; AS R53:
#' 84V33 p366- 369). \emph{Applied Statistics} \bold{29}, 224--227.
#' 
#' R. W. Farebrother (1984), [AS R53] A Remark on Algorithms AS 106 (77V26
#' p92-98), AS 153 (80V29 p224-227) and AS 155: The Distribution of a Linear
#' Combination of \eqn{\chi^2} Random Variables (80V29 p323-333) \emph{Applied
#' Statistics} \bold{33}, 366--369.
#' 
#' W. Kr?mer & H. Sonnberger (1986), \emph{The Linear Regression Model under
#' Test}. Heidelberg: Physica.
#' 
#' J. Racine & R. Hyndman (2002), Using R To Teach Econometrics. \emph{Journal
#' of Applied Econometrics} \bold{17}, 175--189.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' 
#' ## generate two AR(1) error terms with parameter
#' ## rho = 0 (white noise) and rho = 0.9 respectively
#' err1 <- rnorm(100)
#' 
#' ## generate regressor and dependent variable
#' x <- rep(c(-1,1), 50)
#' y1 <- 1 + x + err1
#' 
#' ## perform Durbin-Watson test
#' DurbinWatsonTest(y1 ~ x)
#' 
#' err2 <- filter(err1, 0.9, method="recursive")
#' y2 <- 1 + x + err2
#' DurbinWatsonTest(y2 ~ x)
#' 
#' ## for a simple vector use:
#' e_t <- c(-32.33, -26.603, 2.215, -16.967, -1.148, -2.512, -1.967, 11.669,
#'          -0.513, 27.032, -4.422, 40.032, 23.577, 33.94, -2.787, -8.606,
#'           0.575, 6.848, -18.971, -29.063)
#' DurbinWatsonTest(e_t ~ 1)
#' 
DurbinWatsonTest <- function(formula, order.by = NULL, alternative = c("greater", "two.sided", "less"),
                             iterations = 15, exact = NULL, tol = 1e-10, data = list()) {

  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)

  if (!inherits(formula, "formula")) {
    if (!is.null(w <- weights(formula))) {
      if (!isTRUE(all.equal(as.vector(w), rep(1L, length(w)))))
        stop("weighted regressions are not supported")
    }
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }

  if (!is.null(order.by))
  {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  n <- nrow(X)
  if (is.null(exact)) exact <- (n < 100)
  k <- ncol(X)

  res <- lm.fit(X,y)$residuals
  dw <- sum(diff(res)^2)/sum(res^2)
  Q1 <- chol2inv(qr.R(qr(X)))
  if (n < 3) {
    warning("not enough observations for computing a p value, set to 1")
    pval <- 1
  } else {
    if (exact)
    {
      A <- diag(c(1,rep(2, n-2), 1))
      A[abs(row(A)-col(A))==1] <- -1
      MA <- diag(rep(1,n)) - X %*% Q1 %*% t(X)
      MA <- MA %*% A
      ev <- eigen(MA)$values[1:(n-k)]
      if (any(Im(ev)>tol)) warning("imaginary parts of eigenvalues discarded")
      ev <- Re(ev)
      ev <- ev[ev>tol]

      pdw <- function(dw) .Fortran("pan", as.double(c(dw,ev)), as.integer(length(ev)),
                                   as.double(0), as.integer(iterations), x=double(1),
                                   PACKAGE = "DescTools")$x
      pval <- switch(alternative,
                     "two.sided" = (2*min(pdw(dw), 1-pdw(dw))),
                     "less" = (1 - pdw(dw)),
                     "greater" = pdw(dw))

      if (is.na(pval) || ((pval > 1) | (pval < 0)))
      {
        warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
        exact <- FALSE
      }
    }
    if (!exact)
    {
      if (n < max(5, k)) {
        warning("not enough observations for computing an approximate p value, set to 1")
        pval <- 1
      } else {
        AX <- matrix(as.vector(filter(X, c(-1, 2, -1))), ncol = k)
        AX[1,] <- X[1,] - X[2,]
        AX[n,] <- X[n,] - X[(n-1),]
        XAXQ <- t(X) %*% AX %*% Q1
        P <- 2*(n-1) - sum(diag(XAXQ))
        Q <- 2*(3*n - 4) - 2* sum(diag(crossprod(AX) %*% Q1)) + sum(diag(XAXQ %*% XAXQ))
        dmean <- P/(n-k)
        dvar <- 2/((n-k)*(n-k+2)) * (Q - P*dmean)
        pval <- switch(alternative,
                       "two.sided" = (2*pnorm(abs(dw-dmean), sd=sqrt(dvar), lower.tail = FALSE)),
                       "less" = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE),
                       "greater" = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
      }
    }
  }

  alternative <- switch(alternative,
                        "two.sided" = "true autocorrelation is not 0",
                        "less" = "true autocorrelation is less than 0",
                        "greater" = "true autocorrelation is greater than 0")

  names(dw) <- "DW"
  RVAL <- list(statistic = dw, method = "Durbin-Watson test",
               alternative = alternative, p.value= pval, data.name=dname)
  class(RVAL) <- "htest"
  return(RVAL)
}






#' Von Neumann's Successive Difference Test 
#' 
#' A popular statistic to test for independence is the von Neumann ratio.
#' 
#' The VN test statistic is in the unbiased case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2 \cdot
#' n}{\sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2 \cdot (n-1)}
#' }{VN=\sum(x_i-x_{i+1})^2 / \sum(x_i-mean(x)^2 * n/n-1} It is known that
#' \eqn{(VN-\mu)/\sigma} is asymptotically standard normal, where
#' \eqn{\mu=\frac{2n}{n-1}}{\mu=2n/(n-1)} and \eqn{\sigma^2=4\cdot n^2
#' \frac{(n-2)}{(n+1)(n-1)^3}}{\sigma^2=[4*n^2 * (n-2)]/[(n+1)(n-1)^3]}.
#' 
#' The VN test statistic is in the original (biased) case
#' \deqn{VN=\frac{\sum_{i=1}^{n-1}(x_i-x_{i+1})^2}{\sum_{i=1}^{n}\left(x_i-\bar{x}\right)^2}}{VN=\sum(x_i-x_{i+1})^2
#' / \sum(x_i-mean(x)^2} The test statistic \eqn{(VN-2)/\sigma} is
#' asymptotically standard normal, where
#' \eqn{\sigma^2=\frac{4\cdot(n-2)}{(n+1)(n-1)}}{\sigma^2=[4*(n-2)]/[(n+1)(n-1)]}.
#' 
#' Missing values are silently removed.
#' 
#' @param x a numeric vector containing the observations
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param unbiased logical. In order for VN to be an unbiased estimate of the
#' true population value, the calculated value is multiplied by
#' \eqn{n/(n-1)}{n/(n-1)}. Default is TRUE. 
#' 
#' @return A list with class "htest" containing the components:
#' \item{statistic}{the value of the VN statistic and the normalized statistic
#' test.}
#' \item{parameter, n}{the size of the data, after the remotion of
#' consecutive duplicate values.}
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{a character string indicating the test
#' performed.}
#' \item{data.name}{a character string giving the name of the
#' data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>
#' 
#' @seealso 
#' \code{\link{BartelsRankTest}} 
#' 
#' @references von Neumann, J. (1941) Distribution of the ratio of the mean
#' square successive difference to the variance. \emph{Annals of Mathematical
#' Statistics} \bold{12}, 367-395.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' VonNeumannTest(d.pizza$temperature)
#' 
VonNeumannTest <- function(x, alternative = c("two.sided", "less", "greater"),
  unbiased=TRUE) {


  ## TODO: use incomplete beta for exact p-values
  ## ************************
  ## see: von Neumann Successive Difference 1941
  ##
  # n <- 50
  # vx <- 1
  #
  # mu2 <- (4 * (3*n - 4)/(n-1)^2) * vx^2
  #
  # q2 <- (3*n^4 - 10*n^3 -18*n^2 + 79*n - 60) / (8*n^3 - 50*n + 48)
  # q1 <- (4 - mu2 * (q2 + 1) * (q2 + 3)) / (4 - mu2 * (q2 + 1))
  # a2 <- 2 * (q1 - q2 - 2) / (q2 + 1)
  # cc <- a2 ^(q1 - q2 - 2) / beta(q1 - q2 -1, q2+1)
  #
  # c(q1, q2, a2, cc)
  #
  # pbeta(0.75, shape1 = q1 - q2 -1, shape2= q2+1)
  # pbeta(0.75, shape1 = q1 - q2 -1, shape2= q2+1)
  #
  # beta(q1 - q2 -1, q2+1)


  alternative <- match.arg(alternative)

  dname <- deparse(substitute(x))

  x <- x[!is.na(x)]

  d <- diff(x)
  n <- length(x)
  mx <- mean(x)

  if (unbiased) {

    # http://www.chegg.com/homework-help/detecting-autocorrelation-von-neumann-ratio-test-assuming-re-chapter-12-problem-4-solution-9780073375779-exc

    VN <- sum(d^2) / sum((x - mx)^2) * n/(n-1)
    Ex <- 2 * n/(n-1)
    Vx <- 4 * n^2 * (n-2) / ((n+1) * (n-1)^3)
    z <- (VN - Ex) / Vx

  } else {
    VN <- sum(d^2) / sum((x - mx)^2)
    z <- (1-(VN/2)) / sqrt((n-2)/(n^2 - 1))
  }


  if (alternative == "less") {
    pval <- pnorm(z)
  }
  else if (alternative == "greater") {
    pval <- pnorm(z, lower.tail = FALSE)
  }
  else {
    pval <- 2 * pnorm(-abs(z))
  }
  names(VN) <- "VN"
  method <- "Von Neumann Successive Difference Test"

  rval <- list(statistic = c(VN, z=z), p.value = pval,
               method = method,
               alternative = alternative, data.name = dname,
               z = z)

  class(rval) <- "htest"
  return(rval)

}







#' Bartels Rank Test of Randomness
#' 
#' Performs the Bartels rank test of randomness, which tests if a sample is
#' sampled randomly from an underlying population. Data must at least be
#' measured on an ordinal scale.
#' 
#' The RVN test statistic is
#' \deqn{RVN=\frac{\sum_{i=1}^{n-1}(R_i-R_{i+1})^2}{\sum_{i=1}^{n}\left(R_i-(n+1)/2\right)^2}}{RVN=\sum(R_i-R_{i+1})^2
#' / \sum(R_i-(n+1)/2)^2} where \eqn{R_i=rank(X_i), i=1,\dots,
#' n}{R_i=rank(X_i), i=1,...,n}. It is known that \eqn{(RVN-2)/\sigma} is
#' asymptotically standard normal, where
#' \eqn{\sigma^2=\frac{4(n-2)(5n^2-2n-9)}{5n(n+1)(n-1)^2}}{\sigma^2=[4(n-2)(5n^2-2n-9)]/[5n(n+1)(n-1)^2]}.
#' 
#' By using the alternative "\code{trend}" the null hypothesis of randomness is
#' tested against a trend. By using the alternative "\code{oscillation}" the
#' null hypothesis of randomness is tested against a systematic oscillation.
#' 
#' Missing values are silently removed.
#' 
#' Bartels test is a rank version of von Neumann's test.
#' 
#' @param x a numeric vector containing the observations
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "\code{two.sided}" (default), "\code{trend}" or
#' "\code{oscillation}".
#' @param method a character string specifying the method used to compute the
#' p-value. Must be one of \code{normal} (default), \code{beta} or \code{auto}.
#' 
#' @return A list with class "htest" containing the components:
#' \item{statistic}{the value of the normalized statistic test.}
#' \item{parameter, n}{the size of the data, after the remotion of consecutive
#' duplicate values.} 
#' \item{p.value}{the p-value of the test.}
#' \item{alternative}{a character string describing the alternative hypothesis.} 
#' \item{method}{a character string indicating the test performed.} 
#' \item{data.name}{a character string giving the name of the data.} 
#' \item{rvn}{the value of the RVN statistic (not show on screen).}
#' \item{nm}{the value of the NM statistic, the numerator of RVN (not show on
#' screen).} 
#' \item{mu}{the mean value of the RVN statistic (not show on screen).} 
#' \item{var}{the variance of the RVN statistic (not show on screen).}
#' 
#' @author 
#' Frederico Caeiro <fac@@fct.unl.pt>
#' 
#' @seealso
#' \code{\link[randtests]{rank.test}}, 
#' \code{\link{RunsTest}}
#' 
#' @references 
#' Bartels, R. (1982) The Rank Version of von Neumann's Ratio Test
#' for Randomness, \emph{Journal of the American Statistical Association},
#' \bold{77} (377), 40-46.
#' 
#' Gibbons, J.D. and Chakraborti, S. (2003) \emph{Nonparametric Statistical
#' Inference}, 4th ed. (pp. 97-98). URL:
#' \url{http://books.google.pt/books?id=dPhtioXwI9cC&lpg=PA97&ots=ZGaQCmuEUq}
#' 
#' von Neumann, J. (1941) Distribution of the ratio of the mean square
#' successive difference to the variance. \emph{Annals of Mathematical
#' Statistics} \bold{12}, 367-395.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' ## Example 5.1 in Gibbons and Chakraborti (2003), p.98.
#' ## Annual data on total number of tourists to the United States for 1970-1982.
#' 
#' years <- 1970:1982
#' tourists <- c(12362, 12739, 13057, 13955, 14123,  15698, 17523, 18610, 19842,
#'       20310, 22500, 23080, 21916)
#' plot(years, tourists, pch=20)
#' 
#' BartelsRankTest(tourists, alternative="trend", method="beta")
#' 
#' #  Bartels Ratio Test
#' #
#' # data:  tourists
#' # statistic = -3.6453, n = 13, p-value = 1.21e-08
#' # alternative hypothesis: trend
#' 
#' 
#' ## Example in Bartels (1982).
#' ## Changes in stock levels for 1968-1969 to 1977-1978 (in $A million), deflated by the
#' ## Australian gross domestic product (GDP) price index (base 1966-1967).
#' x <- c(528, 348, 264, -20, - 167, 575, 410, -4, 430, - 122)
#' 
#' BartelsRankTest(x, method="beta")
#' 
BartelsRankTest <- function(x, alternative = c("two.sided", "trend", "oscillation"),
                            method = c("normal", "beta", "auto")) {

  # Performs Bartels Ratio Test for Randomness.
  #
  # Args:
  #   x: a numeric vector containing the data.
  #   alternative hypothesis, must be one of "two.sided" (default), "left.sided" or "right.sided"
  #   pv.method: asymptotic aproximation method used to compute the p-value.
  #
  # Returns:
  #   statistic: the value of the RVN statistic test and the theoretical mean value and variance of the RVN statistic test.
  #   n: the sample size, after the remotion of consecutive duplicate values.
  #   p.value: the asymptotic p-value.
  #   method: a character string indicating the test performed.
  #   data.name: a character string giving the name of the data.
  #   alternative: a character string describing the alternative.

  pvalue <- match.arg(method, c("normal", "beta", "auto"))

  dname <- deparse(substitute(x))

  # Remove NAs
  x <- na.omit(x)
  stopifnot(is.numeric(x))
  n <- length(x)


  # if (alternative == "t") {alternative <- "two.sided"}
  # if (alternative == "l") {alternative <- "left.sided"}
  # if (alternative == "r") {alternative <- "right.sided"}
  # if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
  # {stop("must give a valid alternative")}

  alternative <- match.arg(alternative)

  if (n < 10) {stop("sample size must be greater than 9")}

  # unique
  rk <- rank(x)
  d <- diff(rk)
  d.rank <- sum(rk^2) - n * (mean(rk)^2)
  RVN <- sum(d^2) / d.rank
  mu <- 2
  vr <- (4*(n-2)*(5*n^2-2*n-9))/(5*n*(n+1)*(n-1)^2)

  # Computes the p-value
  if (pvalue == "auto") {
    pvalue <- ifelse(n <= 100, "beta", "normal")
  }

  if (pvalue == "beta") {
    btp <- (5*n*(n+1)*(n-1)^2)/(2*(n-2)*(5*n^2-2*n-9))-1/2
    pv0 <- pbeta(RVN/4, shape1=btp, shape2=btp)
  }
  if (pvalue=="normal") {
    pv0 <- pnorm((RVN - mu) / sqrt(vr))
  }

  if (alternative=="two.sided") {
    pv <- 2 * min(pv0, 1 - pv0)
    alternative <- "nonrandomness"
  }
  if (alternative == "trend") {
    pv <- pv0
    alternative <- "trend"
  }
  if (alternative == "oscillation") {
    pv <- 1 - pv0
    alternative <- "systematic oscillation"
  }

  test <- (RVN - mu) / sqrt(vr)
  rval <- list(statistic = c(RVN=RVN, z=test), nm=sum(d^2), rvn=RVN, mu=mu, var=vr, p.value = pv,
               method = "Bartels Ratio Test", data.name = dname, parameter=c(n=n), n=n, alternative=alternative)
  class(rval) <- "htest"
  return(rval)

}





#' Moses Test of Extreme Reactions
#' 
#' Perform Moses test of extreme reactions, which is a distribution-free
#' non-parametric test for the difference between two independent groups in the
#' extremity of scores (in both directions) that the groups contain.  Scores
#' from both groups are pooled and converted to ranks, and the test statistic
#' is the span of scores (the range plus 1) in one of the groups chosen
#' arbitrarily. An exact probability is computed for the span and then
#' recomputed after dropping a specified number of extreme scores from each end
#' of its range. The exact one-tailed probability is calculated.
#' 
#' For two independent samples from a continuous field, this tests whether
#' extreme values are equally likely in both populations or if they are more
#' likely to occur in the population from which the sample with the larger
#' range was drawn.
#' 
#' Note that the ranks are calculated in decreasing mode.
#' 
#' @aliases MosesTest MosesTest.default MosesTest.formula
#' @param x numeric vector of data values. \code{x} will be treated as control
#' group. Non-finite (e.g. infinite or missing) values will be omitted.
#' @param y numeric vector of data values. \code{y} will be treated as
#' experiment group. Non-finite (e.g. infinite or missing) values will be
#' omitted.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param extreme integer, defines the number of extreme values to be dropped
#' from the control group before calculating the span. Default (\code{NULL}) is
#' the integer part of \code{0.05 * length(x)} or \code{1}, whichever is
#' greater. If extreme is too large, it will be cut down to
#' \code{floor(length(x)-2)/2}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \dQuote{htest} containing the following
#' components: 
#' \item{statistic}{the value of the Moses Test statistic.}
#' \item{p.value}{the p-value for the test.} 
#' \item{method}{the character string \dQuote{Moses Test of Extreme Reactions}.} 
#' \item{data.name}{a character string giving the name(s) of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>
#' 
#' @seealso 
#' \code{\link{wilcox.test}}, 
#' \code{\link{ks.test}}
#' 
#' @references 
#' Moses, L.E. (1952) A Two-Sample Test, \emph{Psychometrika}, 17, 239-247.
#' \url{http://www-01.ibm.com/support/knowledgecenter/SSLVMB_20.0.0/com.ibm.spss.statistics.help/alg_npar_tests_moses.htm}
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#' 
#' MosesTest(x, y)
#' 
#' 
#' set.seed(1479)
#' x <- sample(1:20, 10, replace=TRUE)
#' y <- sample(5:25, 6, replace=TRUE)
#' 
#' MosesTest(x, y)
#' 
MosesTest <- function(x, ...) { 
  UseMethod("MosesTest")
}

# Extremreaktionen nach Moses: Nullhypothese: Die Spannweite der Werte ist
# in beiden Gruppen gleich gross. Die Werte beider Gruppen werden in eine gemeinsame
# Reihenfolge gebracht. Anschliessend werden ihnen Rangwerte zugewiesen.
# Eine der Gruppen (die Gruppe des Wertes, der in dem Dialogfeld
#                   Gruppen definieren als erstes angegeben ist) wird als Kontrollgruppe verwendet.
# Fuer diese Gruppe wird die Spannweite der Raenge als Differenz zwischen
# dem groessten und kleinsten Rangwert berechnet. Anhand dieser Spannweite errechnet
# sich die einseitige Signifikanz. Zusaetzlich wird der Test ein zweites
# Mal durchgefuehrt, wobei die Ausreisser der Gesamtstichprobe ausgeschlossen
# werden (insgesamt 5% der Faelle). Dabei werden sowohl die hoechsten als auch
# die niedrigsten Raenge entfernt. Das Testergebnis teilt die Anzahl der Faelle beider
# Gruppen, die Spannweiten und die entsprechenden einseitigen Signifikanzen
# fuer beide Tests (mit und ohne Ausreisser) mit. Fuer ein Beispiel siehe oben
# Abschnitt Moses-Test, S. 760.


MosesTest.formula <- function(formula, data, subset, na.action, ...) {

  # this is a taken analogue to wilcox.test.formula

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  names(DATA) <- c("x", "y")
  y <- DoCall("MosesTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}



MosesTest.default <- function(x, y, extreme = NULL, ...) {

  # example
  # x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
  # y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
  # MosesTest(y, x)

  if (is.null(extreme)) extreme <- pmax(floor(0.05 * length(x)), 1)
  h <- extreme
  if (2*h > length(x)-2) h <- floor((length(x)-2)/2)

  # no alternative for the moses.test
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

  nk <- length(x)
  ne <- length(y)
  # decreasing ranks following SPSS-calculation
  R1 <- rank(-c(x, y))[1:nk]
  R1 <- sort(R1)[(h+1):(length(R1)-h)]

  S <- ceiling(max(R1) - min(R1) + 1)

  tmp <- 0
  for( i in 0 : (S - nk + 2*h)) {
    tmp <- tmp + choose(i + nk - 2*h - 2, i) * choose(ne + 2*h + 1 - i, ne - i)
  }

  PVAL <- (tmp / choose(nk + ne, nk))

  structure(list(statistic = c(S = S),
                 p.value = PVAL,
                 method = "Moses Test of Extreme Reactions",
                 alternative = "extreme values are more likely in x than in y",
                 data.name = DNAME),
            class = "htest")

}





#' Siegel-Tukey Test For Equality In Variability
#' 
#' Non-parametric Siegel-Tukey test for equality in variability. The null
#' hypothesis is that the variability of x is equal between two groups. A
#' rejection of the null hypothesis indicates that variability differs between
#' the two groups. \code{SiegelTukeyRank} returns the ranks, calculated after
#' Siegel Tukey logic.
#' 
#' The Siegel-Tukey test has relatively low power and may, under certain
#' conditions, indicate significance due to differences in medians rather than
#' differences in variabilities (consider using the argument
#' \code{adjust.median}). Consider also using \code{\link{mood.test}} or
#' \code{\link{ansari.test}}.
#' 
#' @aliases SiegelTukeyTest SiegelTukeyRank SiegelTukeyTest.default
#' SiegelTukeyTest.formula
#' @param x,y numeric vector of data values. Non-finite (e.g. infinite or
#' missing) values will be omitted.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of x.
#' @param adjust.median Should between-group differences in medians be leveled
#' before performing the test? In certain cases, the Siegel-Tukey test is
#' susceptible to median differences and may indicate significant differences
#' in variability that, in reality, stem from differences in medians. Default
#' is \code{FALSE}.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param mu a number specifying an optional parameter used to form the null
#' hypothesis. See Details.
#' @param exact a logical indicating whether an exact p-value should be
#' computed. This is passed directly to \code{\link{wilcox.test}}.
#' @param correct a logical indicating whether to apply continuity correction
#' in the normal approximation for the p-value.
#' @param conf.int a logical indicating whether a confidence interval should be
#' computed.
#' @param conf.level confidence level of the interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param drop.median logical, defining whether the median of the combined
#' samples should be left out, ensuring that there's an even number of elements
#' (which is a requirement of the Siegel-Tukey test). Defaults to \code{TRUE}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ Siegel-Tukey test (Wilcoxon test on tie-adjusted
#' Siegel-Tukey ranks, after the median adjustment if specified).}
#' \item{p.value}{ the p-value for the test} 
#' \item{null.value}{is the value of
#' the median specified by the null hypothesis. This equals the input argument
#' \code{mu}. } 
#' \item{alternative}{a character string describing the
#' alternative hypothesis.} 
#' \item{method}{ the type of test applied}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @author 
#' Daniel Malter, Tal Galili <tal.galili@@gmail.com>, 
#' Andri Signorell <andri@@signorell.net>\cr
#' published on:
#' \url{http://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/}
#' 
#' @seealso 
#' \code{\link{mood.test}}, 
#' \code{\link{ansari.test}},
#' \code{\link{wilcox.test}}, 
#' \code{\link{LeveneTest}}
#' 
#' @references 
#' Siegel, S., Tukey, J. W. (1960): A nonparametric sum of ranks
#' procedure for relative spread in unpaired samples. \emph{Journal of the
#' American Statistical Association}.
#' 
#' Sheskin, D. J. (2004): \emph{Handbook of parametric and nonparametric
#' statistical procedures} 3rd edition. Chapman and Hall/CRC. Boca Raton, FL.
#' 
#' @keywords htest
#'
#'
#' @examples
#' 
#' # Duller, S. 183
#' x <- c(12, 13, 29, 30)
#' y <- c(15, 17, 18, 24, 25, 26)
#' SiegelTukeyTest(x, y)
#' SiegelTukeyTest(x, y, alternative="greater")
#' 
#' # Duller, S. 323
#' old <- c(870,930,935,1045,1050,1052,1055)
#' new <- c(932,970,980,1001,1009,1030,1032,1040,1046)
#' SiegelTukeyTest(old, new, alternative = "greater")
#' # compare to the recommended alternatives
#' mood.test(old, new, alternative="greater")
#' ansari.test(old, new, alternative="greater")
#' 
#' # Bortz, S. 250
#' x <- c(26.3,26.5,26.8,27.0,27.0,27.2,27.3,27.3,27.4,27.5,27.6,27.8,27.9)
#' id <- c(2,2,2,1,2,2,1,2,2,1,1,1,2)-1
#' SiegelTukeyTest(x ~ id)
#' 
#' 
#' # Sachs, Angewandte Statistik, 12. Auflage, 2007, S. 314
#' A <- c(10.1,7.3,12.6,2.4,6.1,8.5,8.8,9.4,10.1,9.8)
#' B <- c(15.3,3.6,16.5,2.9,3.3,4.2,4.9,7.3,11.7,13.1)
#' SiegelTukeyTest(A, B)
#' 
#' 
#' 
#' ### 1
#' x <- c(4,4,5,5,6,6)
#' y <- c(0,0,1,9,10,10)
#' SiegelTukeyTest(x, y)
#' 
#' ### 2
#' # example for a non equal number of cases:
#' x <- c(4,4,5,5,6,6)
#' y <- c(0,0,1,9,10)
#' SiegelTukeyTest(x, y)
#' 
#' ### 3
#' x <- c(33, 62, 84, 85, 88, 93, 97, 4, 16, 48, 51, 66, 98)
#' id <- c(0,0,0,0,0,0,0,1,1,1,1,1,1)
#' SiegelTukeyTest(x ~ id)
#' 
#' ### 4
#' x <- c(177,200,227,230,232,268,272,297,47,105,126,142,158,172,197,220,225,230,262,270)
#' id <- c(rep(0,8),rep(1,12))
#' SiegelTukeyTest(x ~ id, adjust.median=TRUE)
#' 
#' ### 5
#' x <- c(33,62,84,85,88,93,97)
#' y <- c(4,16,48,51,66,98)
#' SiegelTukeyTest(x, y)
#' 
#' ### 6
#' x <- c(0,0,1,4,4,5,5,6,6,9,10,10)
#' id <- c(0,0,0,1,1,1,1,1,1,0,0,0)
#' SiegelTukeyTest(x ~ id)
#' 
#' ### 7
#' x <- c(85,106,96, 105, 104, 108, 86)
#' id <- c(0,0,1,1,1,1,1)
#' SiegelTukeyTest(x ~ id)
#' 
SiegelTukeyTest <- function(x, ...)  {
  UseMethod("SiegelTukeyTest")
}

SiegelTukeyTest.formula <- function(formula, data, subset, na.action, ...)
{
  # this is a taken analogue to wilcox.test.formula

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  names(DATA) <- c("x", "y")
  y <- DoCall("SiegelTukeyTest", c(DATA, list(...)))
  y$data.name <- DNAME
  y

}



SiegelTukeyRank <- function(x, g, drop.median = TRUE) {

  # they do not drop the median in:
  # http://en.wikipedia.org/wiki/Siegel%E2%80%93Tukey_test
  # A <- c(33,62,84,85,88,93,97); B <- c(4,16,48,51,66,98)
  # this is wrong there, as the author did not leave the median out

  ord.x <- order(x, g)
  sort.x <- x[ord.x]
  sort.id <- g[ord.x]

  n <- length(x)
  if (drop.median) {
    if (n %% 2 > 0) {
      # gonna have to drop the (first) median value
      # as we sorted by the groupsize, this will be the one out of the bigger group (if existing)
      fm <- which( sort.x == median(sort.x))[1]
      sort.x <- sort.x[-fm]
      sort.id <- sort.id[-fm]
      n <- n-1
    }
  }

  base1 <- c(1, 4)
  iterator1 <- matrix(seq(from = 1, to = n, by = 4)) - 1
  rank1 <- apply(iterator1, 1, function(x) x + base1)

  iterator2 <- matrix(seq(from = 2, to = n, by = 4))
  base2 <- c(0, 1)
  rank2 <- apply(iterator2, 1, function(x) x + base2)

  if (length(rank1) == length(rank2)) {
    rank <- c(rank1[1:floor(n/2)], rev(rank2[1:ceiling(n/2)]))
  } else {
    rank <- c(rank1[1:ceiling(n/2)], rev(rank2[1:floor(n/2)]))
  }

  unique.ranks <- tapply(rank, sort.x, mean)
  unique.x <- as.numeric(as.character(names(unique.ranks)))

  ST.matrix <- merge(
    data.frame(sort.x, sort.id),          # this are the original values in x-order
    data.frame(unique.x, unique.ranks),   # this is the rank.matrix
    by.x = "sort.x", by.y = "unique.x")

  return(ST.matrix)
}


SiegelTukeyTest.default <- function(x, y, adjust.median = FALSE,
                                    alternative = c("two.sided","less","greater"), mu = 0,
                                    exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95, ...) {
  ###### published on:
  #   http://www.r-statistics.com/2010/02/siegel-tukey-a-non-parametric-test-for-equality-in-variability-r-code/
  #   Main author of the function:  Daniel Malter

  # Doku: http://www.crcnetbase.com/doi/abs/10.1201/9781420036268.ch14


  if (!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
    stop("'mu' must be a single number")

  if (conf.int) {
    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
      stop("'conf.level' must be a single number between 0 and 1")
  }

  if (!is.numeric(x))
    stop("'x' must be numeric")

  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
  }
  else {
    DNAME <- deparse(substitute(x))
    x <- x[is.finite(x)]
  }

  # adjusting median
  if (adjust.median) {
    x <- x - median(x)
    y <- y - median(y)
  }

  # the larger group comes first
  if ( length(x) > length(y) ) {
    xx <- c(x, y)
    id <- c(rep(0, length(x)), rep(1, length(y)))
  } else {
    xx <- c(y,x)
    id <- c(rep(0, length(y)), rep(1, length(x)))
  }

  strank <- SiegelTukeyRank(xx, g = id)
  ranks0 <- strank$unique.ranks[strank$sort.id == 0]
  ranks1 <- strank$unique.ranks[strank$sort.id == 1]

  RVAL <- wilcox.test(ranks0, ranks1, alternative = alternative,
                      mu = mu, paired = FALSE, exact = exact, correct = correct,
                      conf.int = conf.int, conf.level = conf.level)

  RVAL$statistic <- sum(ranks1)
  names(RVAL$statistic)  <- "ST"
  RVAL$data.name <- DNAME
  RVAL <- c(RVAL, list(stranks = strank, MeanRanks = c(mean(ranks0), mean(ranks1))))
  RVAL$method <- "Siegel-Tukey-test for equal variability"
  RVAL$null.value <- 1
  names(RVAL$null.value) <- "ratio of scales"
  class(RVAL) <- "htest"
  return(RVAL)

  if (suppressWarnings(wilcox.test(x,y)$p.value) < 0.05) warning("SiegelTukeyTest: wilcox.test(x, y) is significant! Consider setting adjust.median = TRUE." )

}






#' Exact Version of Jonckheere-Terpstra Test
#' 
#' Jonckheere-Terpstra test to test for ordered differences among classes.
#' 
#' JonckheereTerpstraTest is the exact (permutation) version of the
#' Jonckheere-Terpstra test.  It uses the statistic \deqn{\sum_{k<l} \sum_{ij}
#' I(X_{ik} < X_{jl}) + 0.5 I(X_{ik} = X_{jl}),} where \eqn{i, j} are
#' observations in groups \eqn{k} and \eqn{l} respectively.  The asymptotic
#' version is equivalent to \code{cor.test(x, g, method="k")}. The exact
#' calculation requires that there be no ties and that the sample size is less
#' than 100. When data are tied and sample size is at most 100 permutation
#' p-value is returned.\cr
#' 
#' If x is a list, its elements are taken as the samples to be compared, and
#' hence have to be numeric data vectors.  In this case, g is ignored, and one
#' can simply use JonckheereTerpstraTest(x) to perform the test.  If the
#' samples are not yet contained in a list, use JonckheereTerpstraTest(list(x,
#' ...)). \cr
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @aliases JonckheereTerpstraTest JonckheereTerpstraTest.default
#' JonckheereTerpstraTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of x. Ignored if x is a list.
#' @param alternative means are monotonic (\code{two.sided}),
#' \code{increasing}, or \code{decreasing}
#' @param nperm number of permutations for the reference distribution.  The
#' default is \code{NULL} in which case the permutation p-value is not
#' computed. It's recommended to set \code{nperm} to 1000 or higher if
#' permutation p-value is desired.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and rhs the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further argument to be passed to methods.
#' @note The function was previously published as \code{jonckheere.test()} in
#' the \pkg{clinfun} package and has been integrated here without logical
#' changes. Some argument checks and a formula interface were added.
#' 
#' @author 
#' Venkatraman E. Seshan <seshanv@@mskcc.org>, 
#' minor adaptations Andri Signorell
#' 
#' @references Jonckheere, A. R. (1954). A distribution-free k-sample test
#' again ordered alternatives. \emph{Biometrika} 41:133-145.
#' 
#' Terpstra, T. J. (1952). The asymptotic normality and consistency of
#' Kendall's test against trend, when ties are present in one ranking.
#' \emph{Indagationes Mathematicae} 14:327-333.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' set.seed(1234)
#' g <- ordered(rep(1:5, rep(10,5)))
#' x <- rnorm(50) + 0.3 * as.numeric(g)
#' 
#' JonckheereTerpstraTest(x, g)
#' 
#' x[1:2] <- mean(x[1:2]) # tied data
#' 
#' JonckheereTerpstraTest(x, g)
#' JonckheereTerpstraTest(x, g, nperm=5000)
#' 
#' # Duller, S. 222
#' coffee <- data.frame(
#'   time=c(
#'   447,396,383,410,
#'   438,521,468,391,504,472,
#'   513,543,506,489,407), 
#'   grp=Untable(c(4,6,5), type="ordered")[,1]
#' )  
#' 
#' # the formula interface:
#' JonckheereTerpstraTest(time ~ grp, data=coffee)
#' 
JonckheereTerpstraTest <- function(x, ...) {
  UseMethod("JonckheereTerpstraTest")
}

JonckheereTerpstraTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  y <- DoCall("JonckheereTerpstraTest", as.list(mf))
  y$data.name <- DNAME
  y
}

JonckheereTerpstraTest.default <- function(x, g, alternative = c("two.sided", "increasing", "decreasing"), nperm=NULL, ...) {

  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }
  n <- length(x)
  if (n < 2)
    stop("not enough observations")

  # start calculating

  jtpdf <- function(gsize) {
    ng <- length(gsize)
    cgsize <- rev(cumsum(rev(gsize)))
    mxsum <- sum(gsize[-ng]*cgsize[-1]) + 1
    zz <- .Fortran("jtpdf",
                   as.integer(mxsum),
                   pdf=double(mxsum),
                   as.integer(ng),
                   as.integer(cgsize),
                   double(mxsum),
                   double(mxsum))
    zz$pdf
  }

  jtperm.p <- function(x, ng, gsize, cgsize, alternative, nperm) {
    # this function computes the pdf using the convolution by Mark van de Wiel

    n <- length(x)
    pjtrsum <- rep(0, nperm)
    for (np in 1:nperm) {
      jtrsum <- 0
      for(i in 1L:(ng-1)) {
        na <- gsize[i]
        nb <- n-cgsize[i+1]
        # this jtrsum will be small if data are increasing and large if decreasing
        jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
      }
      pjtrsum[np] <- jtrsum
      # permute the data; this way the first value is the original statistic
      x <- sample(x)
    }
    # one-sided p-values
    # number of permuted values at least as small as original
    iPVAL <- sum(pjtrsum <= pjtrsum[1])/nperm
    # number of permuted values at least as large as original
    dPVAL <- sum(pjtrsum >= pjtrsum[1])/nperm
    # return p-value for the alternative of interest
    PVAL <- switch(alternative,
                   "two.sided" = 2*min(iPVAL, dPVAL, 1),
                   "increasing" = iPVAL,
                   "decreasing" = dPVAL)
    PVAL
  }


  # Alternative for the JT-Statistic
  # JT <- function(z) {
  #
  #   w <- function(x, y) {
  #     # verbatim from wilcox.test STATISTIC
  #     r <- rank(c(x, y))
  #     n.x <- as.double(length(x))
  #     n.y <- as.double(length(y))
  #     STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  #     STATISTIC
  #   }
  #
  #   k <- length(z)
  #   u <- 0
  #
  #   for(i in 2:k) {
  #     for(j in 1:(i-1))	{
  #       u <- u + w(z[[i]], z[[j]])
  #     } }
  #   u
  # }


  # see:
  #   library(NSM3)
  #   HOllander Wolfe pp 218
  #   piece <- c(40,35,38,43,44,41, 38,40,47,44,40,42, 48,40,45,43,46,44)
  #   grp <- factor(rep(c("ctrl","A","B"), each=6), ordered=TRUE, levels=c("ctrl","A","B"))
  #
  #   JonckheereTerpstraTest(piece, grp)
  #   pJCK(piece, grp)


  if (!is.numeric(x)) stop("data values should be numeric")
  if (!is.numeric(g) & !is.ordered(g)) stop("group should be numeric or ordered factor")
  alternative <- match.arg(alternative)
  METHOD <- "Jonckheere-Terpstra test"
  PERM <- !missing(nperm)
  n <- length(x)
  if (length(g) != n) stop("lengths of data values and group don't match")
  TIES <- length(unique(x)) != n
  gsize <- table(g)
  ng <- length(gsize)
  cgsize <- c(0,cumsum(gsize))
  x <- x[order(g)]
  jtrsum <- jtmean <- jtvar <- 0
  for(i in 1:(ng-1)) {
    na <- gsize[i]
    nb <- n-cgsize[i+1]
    jtrsum <- jtrsum + sum(rank(x[(cgsize[i]+1):n])[1:na]) - na*(na+1)/2
    jtmean <- jtmean + na*nb/2
    jtvar <- jtvar + na*nb*(na+nb+1)/12
  }
  # this jtrsum will be small if data are increasing and large if decreasing
  # to reverse this use 2*jtmean - jtrsum
  jtrsum <- 2*jtmean - jtrsum
  STATISTIC <- jtrsum
  names(STATISTIC) <- "JT"
  if (PERM) {
    PVAL <- jtperm.p(x, ng, gsize, cgsize, alternative, nperm)
  } else {
    if (n > 100 | TIES) {
      warning("Sample size > 100 or data with ties \n p-value based on normal approximation. Specify nperm for permutation p-value")
      zstat <- (STATISTIC-jtmean)/sqrt(jtvar)
      PVAL <- pnorm(zstat)
      PVAL <- switch(alternative,
                     "two.sided" = 2*min(PVAL, 1-PVAL, 1),
                     "increasing" = 1-PVAL,
                     "decreasing" = PVAL)
    } else {
      dPVAL <- sum(jtpdf(gsize)[1:(jtrsum+1)])
      iPVAL <- 1-sum(jtpdf(gsize)[1:(jtrsum)])
      PVAL <- switch(alternative,
                     "two.sided" = 2*min(iPVAL, dPVAL, 1),
                     "increasing" = iPVAL,
                     "decreasing" = dPVAL)
    }
  }

  RVAL <- list(statistic = STATISTIC,
               p.value = as.numeric(PVAL),
               alternative = alternative,
               method = METHOD,
               data.name = DNAME)
  class(RVAL) <- "htest"
  RVAL

}



# ***********************************
# Tests aus library(nortest)



#' Shapiro-Francia Test for Normality
#' 
#' Performs the Shapiro-Francia test for the composite hypothesis of normality.
#' 
#' The test statistic of the Shapiro-Francia test is simply the squared
#' correlation between the ordered sample values and the (approximated)
#' expected ordered quantiles from the standard normal distribution. The
#' p-value is computed from the formula given by Royston (1993).
#' 
#' @param x a numeric vector of data values, the number of which must be
#' between 5 and 5000. Missing values are allowed.
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the value of the Shapiro-Francia statistic.}
#' \item{p.value
#' }{the p-value for the test.} 
#' \item{method}{the character string \dQuote{Shapiro-Francia normality test}.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'  
#' @note The Shapiro-Francia test is known to perform well, see also the
#' comments by Royston (1993). The expected ordered quantiles from the standard
#' normal distribution are approximated by \code{qnorm(ppoints(x, a = 3/8))},
#' being slightly different from the approximation \code{qnorm(ppoints(x, a =
#' 1/2))} used for the normal quantile-quantile plot by \code{\link{qqnorm}}
#' for sample sizes greater than 10.
#' 
#' @author 
#' Juergen Gross <gross@@statistik.uni-dortmund.de>
#' 
#' @seealso 
#' \code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. 
#' \code{\link{AndersonDarlingTest}},
#' \code{\link{CramerVonMisesTest}}, 
#' \code{\link{LillieTest}},
#' \code{\link{PearsonTest}} for performing further tests for normality.
#' \code{\link{qqnorm}} for producing a normal quantile-quantile plot.
#' 
#' @references 
#' Royston, P. (1993): A pocket-calculator algorithm for the
#' Shapiro-Francia test for non-normality: an application to medicine.
#' \emph{Statistics in Medicine}, 12, 181--184.
#' 
#' Thode Jr., H.C. (2002): \emph{Testing for Normality}. Marcel Dekker, New
#' York. (2002, Sec. 2.3.2)
#' @keywords htest
#'
#'
#' @examples
#' 
#' ShapiroFranciaTest(rnorm(100, mean = 5, sd = 3))
#' ShapiroFranciaTest(runif (100, min = 2, max = 4))
#' 
ShapiroFranciaTest <- function(x) {

  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if ((n < 5 || n > 5000))
    stop("sample size must be between 5 and 5000")
  y <- qnorm(ppoints(n, a = 3/8))
  W <- cor(x, y)^2
  u <- log(n)
  v <- log(u)
  mu <- -1.2725 + 1.0521 * (v - u)
  sig <- 1.0308 - 0.26758 * (v + 2/u)
  z <- (log(1 - W) - mu)/sig
  pval <- pnorm(z, lower.tail = FALSE)
  RVAL <- list(statistic = c(W = W), p.value = pval, method = "Shapiro-Francia normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)

}




#' Pearson Chi-Square Test for Normality
#' 
#' Performs the Pearson chi-square test for the composite hypothesis of
#' normality.
#' 
#' The Pearson test statistic is \eqn{P=\sum (C_{i} - E_{i})^{2}/E_{i}}, where
#' \eqn{C_{i}} is the number of counted and \eqn{E_{i}} is the number of
#' expected observations (under the hypothesis) in class \eqn{i}. The classes
#' are build is such a way that they are equiprobable under the hypothesis of
#' normality. The p-value is computed from a chi-square distribution with
#' \code{n.classes}-3 degrees of freedom if \code{adjust} is \code{TRUE} and
#' from a chi-square distribution with \code{n.classes}-1 degrees of freedom
#' otherwise. In both cases this is not (!) the correct p-value, lying
#' somewhere between the two, see also Moore (1986).
#' 
#' @param x a numeric vector of data values. Missing values are allowed.
#' @param n.classes The number of classes. The default is due to Moore (1986).
#' @param adjust logical; if \code{TRUE} (default), the p-value is computed
#' from a chi-square distribution with \code{n.classes}-3 degrees of freedom,
#' otherwise from a chi-square distribution with \code{n.classes}-1 degrees of
#' freedom.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the value of the Pearson chi-square statistic.}
#' \item{p.value }{the p-value for the test.}
#' \item{method}{the character
#' string \dQuote{Pearson chi-square normality test}.}
#' \item{data.name}{a
#' character string giving the name(s) of the data.}
#' \item{n.classes}{the
#' number of classes used for the test.}
#' \item{df}{the degress of freedom of
#' the chi-square distribution used to compute the p-value.}
#' @note The Pearson chi-square test is usually not recommended for testing the
#' composite hypothesis of normality due to its inferior power properties
#' compared to other tests. It is common practice to compute the p-value from
#' the chi-square distribution with \code{n.classes} - 3 degrees of freedom, in
#' order to adjust for the additional estimation of two parameters. (For the
#' simple hypothesis of normality (mean and variance known) the test statistic
#' is asymptotically chi-square distributed with \code{n.classes} - 1 degrees
#' of freedom.) This is, however, not correct as long as the parameters are
#' estimated by \code{mean(x)} and \code{var(x)} (or \code{sd(x)}), as it is
#' usually done, see Moore (1986) for details. Since the true p-value is
#' somewhere between the two, it is suggested to run \code{PearsonTest} twice,
#' with \code{adjust = TRUE} (default) and with \code{adjust = FALSE}. It is
#' also suggested to slightly change the default number of classes, in order to
#' see the effect on the p-value. Eventually, it is suggested not to rely upon
#' the result of the test.
#' 
#' The function call \code{PearsonTest(x)} essentially produces the same result
#' as the S-PLUS function call \code{chisq.gof((x-mean(x))/sqrt(var(x)),
#' n.param.est=2)}.
#' 
#' @author 
#' Juergen Gross <gross@@statistik.uni-dortmund.de>
#' 
#' @seealso 
#' \code{\link{shapiro.test}} for performing the Shapiro-Wilk test for
#' normality. 
#' \code{\link{AndersonDarlingTest}},
#' \code{\link{CramerVonMisesTest}}, 
#' \code{\link{LillieTest}},
#' \code{\link{ShapiroFranciaTest}} for performing further tests for normality.
#' \code{\link{qqnorm}} for producing a normal quantile-quantile plot.
#' 
#' @references Moore, D.S., (1986) Tests of the chi-squared type. In:
#' D'Agostino, R.B. and Stephens, M.A., eds.: \emph{Goodness-of-Fit
#' Techniques}. Marcel Dekker, New York.
#' 
#' Thode Jr., H.C., (2002) \emph{Testing for Normality}. Marcel Dekker, New
#' York. Sec. 5.2
#' 
#' @keywords htest
#'
#' @examples
#' 
#' PearsonTest(rnorm(100, mean = 5, sd = 3))
#' PearsonTest(runif (100, min = 2, max = 4))
#' 
PearsonTest <- function(x, n.classes = ceiling(2 * (n^(2/5))), adjust = TRUE) {

  DNAME <- deparse(substitute(x))
  x <- x[complete.cases(x)]
  n <- length(x)
  if (adjust) {
    dfd <- 2
  }
  else {
    dfd <- 0
  }
  num <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
  count <- tabulate(num, n.classes)
  prob <- rep(1/n.classes, n.classes)
  xpec <- n * prob
  h <- ((count - xpec)^2)/xpec
  P <- sum(h)
  pvalue <- pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
  RVAL <- list(statistic = c(P = P), p.value = pvalue, method = "Pearson chi-square normality test",
               data.name = DNAME, n.classes = n.classes, df = n.classes -
                 1 - dfd)
  class(RVAL) <- "htest"
  return(RVAL)
}




#' Lilliefors (Kolmogorov-Smirnov) Test for Normality
#' 
#' Performs the Lilliefors (Kolmogorov-Smirnov) test for the composite
#' hypothesis of normality, see e.g. Thode (2002, Sec. 5.1.1).
#' 
#' The Lilliefors (Kolmogorov-Smirnov) test is an EDF omnibus test for the
#' composite hypothesis of normality. The test statistic is the maximal
#' absolute difference between empirical and hypothetical cumulative
#' distribution function. It may be computed as \eqn{D=\max\{D^{+}, D^{-}\}}
#' with \deqn{ D^{+} = \max_{i=1,\ldots, n}\{i/n - p_{(i)}\}, D^{-} =
#' \max_{i=1,\ldots, n}\{p_{(i)} - (i-1)/n\}, } where \eqn{p_{(i)} =
#' \Phi([x_{(i)} - \overline{x}]/s)}. Here, \eqn{\Phi} is the cumulative
#' distribution function of the standard normal distribution, and
#' \eqn{\overline{x}} and \eqn{s} are mean and standard deviation of the data
#' values. The p-value is computed from the Dallal-Wilkinson (1986) formula,
#' which is claimed to be only reliable when the p-value is smaller than 0.1.
#' If the Dallal-Wilkinson p-value turns out to be greater than 0.1, then the
#' p-value is computed from the distribution of the modified statistic \eqn{Z=D
#' (\sqrt{n}-0.01+0.85/\sqrt{n})}, see Stephens (1974), the actual p-value
#' formula being obtained by a simulation and approximation process.
#' 
#' @param x a numeric vector of data values, the number of which must be
#' greater than 4. Missing values are allowed.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the value of the Lilliefors (Kolomogorv-Smirnov)
#' statistic.}
#' \item{p.value }{the p-value for the test.}
#' \item{method}{the
#' character string \dQuote{Lilliefors (Kolmogorov-Smirnov) normality test}.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' @note The Lilliefors (Kolomorov-Smirnov) test is the most famous EDF omnibus
#' test for normality. Compared to the Anderson-Darling test and the Cramer-von
#' Mises test it is known to perform worse. Although the test statistic
#' obtained from \code{LillieTest(x)} is the same as that obtained from
#' \code{ks.test(x, "pnorm", mean(x), sd(x))}, it is not correct to use the
#' p-value from the latter for the composite hypothesis of normality (mean and
#' variance unknown), since the distribution of the test statistic is different
#' when the parameters are estimated.
#' 
#' The function call \code{LillieTest(x)} essentially produces the same result
#' as the S-PLUS function call \code{ks.gof(x)} with the distinction that the
#' p-value is not set to 0.5 when the Dallal-Wilkinson approximation yields a
#' p-value greater than 0.1. (Actually, the alternative p-value approximation
#' is provided for the complete range of test statistic values, but is only
#' used when the Dallal-Wilkinson approximation fails.)
#' 
#' @author 
#' Juergen Gross <gross@@statistik.uni-dortmund.de>
#' 
#' @seealso
#'  \code{\link{shapiro.test}} for performing the Shapiro-Wilk test for
#' normality. 
#' \code{\link{AndersonDarlingTest}},
#' \code{\link{CramerVonMisesTest}}, 
#' \code{\link{PearsonTest}},
#' \code{\link{ShapiroFranciaTest}} for performing further tests for normality.
#' \code{\link{qqnorm}} for producing a normal quantile-quantile plot.
#' 
#' @references 
#' Dallal, G.E. and Wilkinson, L. (1986) An analytic approximation
#' to the distribution of Lilliefors' test for normality. \emph{The American
#' Statistician}, 40, 294--296.
#' 
#' Stephens, M.A. (1974) EDF statistics for goodness of fit and some
#' comparisons. \emph{Journal of the American Statistical Association}, 69,
#' 730--737.
#' 
#' Thode Jr., H.C. (2002) \emph{Testing for Normality} Marcel Dekker, New York.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' LillieTest(rnorm(100, mean = 5, sd = 3))
#' LillieTest(runif (100, min = 2, max = 4))
#' 
LillieTest <- function(x) {

  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 5)
    stop("sample size must be greater than 4")
  p <- pnorm((x - mean(x))/sd(x))
  Dplus <- max(seq(1:n)/n - p)
  Dminus <- max(p - (seq(1:n) - 1)/n)
  K <- max(Dplus, Dminus)
  if (n <= 100) {
    Kd <- K
    nd <- n
  }
  else {
    Kd <- K * ((n/100)^0.49)
    nd <- 100
  }
  pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 *
                  Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598/sqrt(nd) +
                  1.67997/nd)
  if (pvalue > 0.1) {
    KK <- (sqrt(n) - 0.01 + 0.85/sqrt(n)) * K
    if (KK <= 0.302) {
      pvalue <- 1
    }
    else if (KK <= 0.5) {
      pvalue <- 2.76773 - 19.828315 * KK + 80.709644 *
        KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
    }
    else if (KK <= 0.9) {
      pvalue <- -4.901232 + 40.662806 * KK - 97.490286 *
        KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
    }
    else if (KK <= 1.31) {
      pvalue <- 6.198765 - 19.558097 * KK + 23.186922 *
        KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
    }
    else {
      pvalue <- 0
    }
  }
  RVAL <- list(statistic = c(D = K), p.value = pvalue, method = "Lilliefors (Kolmogorov-Smirnov) normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}




#' Cramer-von Mises Test for Normality
#' 
#' Performs the Cramer-von Mises test for the composite hypothesis of
#' normality, see e.g. Thode (2002, Sec. 5.1.3).
#' 
#' The Cramer-von Mises test is an EDF omnibus test for the composite
#' hypothesis of normality. The test statistic is \deqn{ W = \frac{1}{12 n} +
#' \sum_{i=1}^{n} \left (p_{(i)} - \frac{2i-1}{2n} \right), } where
#' \eqn{p_{(i)} = \Phi([x_{(i)} - \overline{x}]/s)}. Here, \eqn{\Phi} is the
#' cumulative distribution function of the standard normal distribution, and
#' \eqn{\overline{x}} and \eqn{s} are mean and standard deviation of the data
#' values. The p-value is computed from the modified statistic \eqn{Z=W (1.0 +
#' 0.5/n)} according to Table 4.9 in Stephens (1986).
#' 
#' @param x a numeric vector of data values, the number of which must be
#' greater than 7. Missing values are allowed.
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{the value of the Cramer-von Mises statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{the character string \dQuote{Cramer-von Mises normality test}.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' 
#' @author 
#' Juergen Gross <gross@@statistik.uni-dortmund.de>
#' 
#' @seealso 
#' \code{\link{shapiro.test}} for performing the Shapiro-Wilk test for normality. 
#' \code{\link{AndersonDarlingTest}}, 
#' \code{\link{LillieTest}},
#' \code{\link{PearsonTest}}, 
#' \code{\link{ShapiroFranciaTest}} for performing
#' further tests for normality. 
#' \code{\link{qqnorm}} for producing a normal quantile-quantile plot.
#' 
#' @references 
#' Stephens, M.A. (1986) Tests based on EDF statistics In:
#' D'Agostino, R.B. and Stephens, M.A., eds.: \emph{Goodness-of-Fit
#' Techniques}. Marcel Dekker, New York.
#' 
#' Thode Jr., H.C. (2002) \emph{Testing for Normality} Marcel Dekker, New York.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' CramerVonMisesTest(rnorm(100, mean = 5, sd = 3))
#' CramerVonMisesTest(runif(100, min = 2, max = 4))
#' 
CramerVonMisesTest <- function(x) {
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if (n < 8)
    stop("sample size must be greater than 7")
  p <- pnorm((x - mean(x))/sd(x))
  W <- (1/(12 * n) +
          sum(
            (p - (2 * seq(1:n) - 1)/(2 * n))^2
          ))
  WW <- (1 + 0.5/n) * W
  if (WW < 0.0275) {
    pval <- 1 - exp(-13.953 + 775.5 * WW - 12542.61 * WW^2)
  }
  else if (WW < 0.051) {
    pval <- 1 - exp(-5.903 + 179.546 * WW - 1515.29 * WW^2)
  }
  else if (WW < 0.092) {
    pval <- exp(0.886 - 31.62 * WW + 10.897 * WW^2)
  }
  else if (WW < 1.1) {
    pval <- exp(1.111 - 34.242 * WW + 12.832 * WW^2)
  }
  else {
    warning("p-value is smaller than 7.37e-10, cannot be computed more accurately")
    pval <- 7.37e-10
  }
  RVAL <- list(statistic = c(W = W), p.value = pval, method = "Cramer-von Mises normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}


#
# AndersonDarlingTest <- function(x) {
#
#     DNAME <- deparse(substitute(x))
#     x <- sort(x[complete.cases(x)])
#     n <- length(x)
#     if (n < 8)
#         stop("sample size must be greater than 7")
#     p <- pnorm((x - mean(x))/sd(x))
#     h <- (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
#     A <- -n - mean(h)
#     AA <- (1 + 0.75/n + 2.25/n^2) * A
#
#     if (AA < 0.2) {
#         pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
#     }
#     else if (AA < 0.34) {
#         pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
#     }
#     else if (AA < 0.6) {
#         pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
#     }
#     else if (AA < 160) {
#         pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
#     }
#     else {
#       pval <-0
#     }
#       RVAL <- list(statistic = c(A = A), p.value = pval, method = "Anderson-Darling normality test",
#         data.name = DNAME)
#     class(RVAL) <- "htest"
#     return(RVAL)
# }


##
## andarl.R
##
##  Anderson-Darling test and null distribution
##
## $Revision: 1.6 $ $Date: 2014/06/24 02:12:20 $
##



#' Anderson-Darling Test of Goodness-of-Fit
#' 
#' Performs the Anderson-Darling test of goodness-of-fit to a specified
#' continuous univariate probability distribution.
#' 
#' This command performs the Anderson-Darling test of goodness-of-fit to the
#' distribution specified by the argument \code{null}. It is assumed that the
#' values in \code{x} are independent and identically distributed random
#' values, with some cumulative distribution function \eqn{F}.  The null
#' hypothesis is that \eqn{F} is the function specified by the argument
#' \code{null}, while the alternative hypothesis is that \eqn{F} is some other
#' function.
#' 
#' @param x numeric vector of data values.
#' @param null a function, or a character string giving the name of a function,
#' to compute the cumulative distribution function for the null distribution.
#' @param \dots additional arguments for the cumulative distribution function.
#' @param nullname optional character string describing the null
#' distribution.\cr
#' The default is \code{"uniform distribution"}.
#' @return An object of class \code{"htest"} representing the result of the
#' hypothesis test.
#' 
#' @author 
#' Original C code by George Marsaglia and John Marsaglia
#' Interface by Adrian Baddeley.
#' @seealso \code{\link{shapiro.test}} and all other tests for normality.
#' @references Anderson, T.W. and Darling, D.A. (1952) Asymptotic theory of
#' certain 'goodness-of-fit' criteria based on stochastic processes.
#' \emph{Annals of Mathematical Statistics} \bold{23}, 193--212.
#' 
#' Anderson, T.W. and Darling, D.A. (1954) A test of goodness of fit.
#' \emph{Journal of the American Statistical Association} \bold{49}, 765--769.
#' 
#' Marsaglia, G. and Marsaglia, J. (2004) Evaluating the Anderson-Darling
#' Distribution.  \emph{Journal of Statistical Software} \bold{9} (2), 1--5.
#' February 2004.  \url{http://www.jstatsoft.org/v09/i02}
#' @keywords htest
#'
#' @examples
#' 
#' x <- rnorm(10, mean=2, sd=1)
#' AndersonDarlingTest(x, "pnorm", mean=2, sd=1)
#' 
AndersonDarlingTest <- function(x, null="punif", ..., nullname) {

  .recogniseCdf <- function(s="punif") {
    if (!is.character(s) || length(s) != 1) return(NULL)
    if (nchar(s) <= 1 || substr(s,1,1) != "p") return(NULL)
    root <- substr(s, 2, nchar(s))
    a <- switch(root,
                beta     = "beta",
                binom    = "binomial",
                birthday = "birthday coincidence",
                cauchy   = "Cauchy",
                chisq    = "chi-squared",
                exp      = "exponential",
                f        = "F",
                gamma    = "Gamma",
                geom     = "geometric",
                hyper    = "hypergeometric",
                lnorm    = "log-normal",
                logis    = "logistic",
                nbinom   = "negative binomial",
                norm     = "Normal",
                pois     = "Poisson",
                t        = "Student's t",
                tukey    = "Tukey (Studentized range)",
                unif     = "uniform",
                weibull  = "Weibull",
                NULL)
    if (!is.null(a))
      return(paste(a, "distribution"))
    b <- switch(root,
                AD     = "Anderson-Darling",
                CvM    = "Cramer-von Mises",
                wilcox = "Wilcoxon Rank Sum",
                NULL)
    if (!is.null(b))
      return(paste("null distribution of", b, "Test Statistic"))
    return(NULL)
  }


  xname <- deparse(substitute(x))
  nulltext <- deparse(substitute(null))
  if (is.character(null)) nulltext <- null
  if (missing(nullname) || is.null(nullname)) {
    reco <- .recogniseCdf(nulltext)
    nullname <- if (!is.null(reco)) reco else
      paste("distribution", sQuote(nulltext))
  }
  stopifnot(is.numeric(x))
  x <- as.vector(x)
  n <- length(x)
  F0 <- if (is.function(null)) null else
    if (is.character(null)) get(null, mode="function") else
      stop("Argument 'null' should be a function, or the name of a function")
  U <- F0(x, ...)
  if (any(U < 0 | U > 1))
    stop("null distribution function returned values outside [0,1]")
  U <- sort(U)
  k <- seq_len(n)
  ## call Marsaglia C code
  z <- .C("ADtestR",
          x = as.double(U),
          n = as.integer(n),
          adstat = as.double(numeric(1)),
          pvalue = as.double(numeric(1))
  )
  STATISTIC <- z$adstat
  names(STATISTIC) <- "An"
  PVAL <- z$pvalue
  METHOD <- c("Anderson-Darling test of goodness-of-fit",
              paste("Null hypothesis:", nullname))
  extras <- list(...)
  parnames <- intersect(names(extras), names(formals(F0)))
  if (length(parnames) > 0) {
    pars <- extras[parnames]
    pard <- character(0)
    for(i in seq_along(parnames))
      pard[i] <- paste(parnames[i], "=", paste(Format(pars[[i]], digits=DescToolsOptions("digits")), collapse=" "))
    pard <- paste("with",
                  ngettext(length(pard), "parameter", "parameters"),
                  "  ",
                  paste(pard, collapse=", "))
    METHOD <- c(METHOD, pard)
  }
  out <- list(statistic = STATISTIC,
              p.value = PVAL,
              method = METHOD,
              data.name = xname)
  class(out) <- "htest"
  return(out)
}

.pAD <- function(q, n=Inf, lower.tail=TRUE, fast=TRUE) {
  q <- as.numeric(q)
  p <- rep(NA_real_, length(q))
  if (any(ones <- is.infinite(q) & (q == Inf)))
    p[ones] <- 1
  if (any(zeroes <- (is.finite(q) & q <= 0) | (is.infinite(q) & (q == -Inf))))
    p[zeroes] <- 0
  ok <- is.finite(q) & (q > 0)
  nok <- sum(ok)
  if (nok > 0) {
    if (is.finite(n)) {
      z <- .C("ADprobN",
              a       = as.double(q[ok]),
              na      = as.integer(nok),
              nsample = as.integer(n),
              prob    = as.double(numeric(nok))
      )
      p[ok] <- z$prob
    } else if (fast) {
      ## fast version adinf()
      z <- .C("ADprobApproxInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok))
      )
      p[ok] <- z$prob
    } else {
      ## slow, accurate version ADinf()
      z <- .C("ADprobExactInf",
              a    = as.double(q[ok]),
              na   = as.integer(nok),
              prob = as.double(numeric(nok))
      )
      p[ok] <- z$prob
    }

  }
  if (!lower.tail)
    p <- 1 - p
  return(p)
}


# .qAD <- local({
#
#   f <- function(x, N, P, Fast) {
#     .pAD(x, N, fast=Fast) - P
#   }
#
#   .qAD <- function(p, n=Inf, lower.tail=TRUE, fast=TRUE) {
#     ## quantiles of null distribution of Anderson-Darling test statistic
#     stopifnot(all(p >= 0))
#     stopifnot(all(p <= 1))
#     if (!lower.tail) p <- 1-p
#     ans <- rep(NA_real_, length(p))
#     for(i in which(p >= 0 & p < 1))
#       ans[i] <- uniroot(f, c(0, 1), N=n, P=p[i], Fast=fast, extendInt="up")$root
#     return(ans)
#   }
#
#   .qAD
# })
#
#
#




# ***********************************
# Tests aus library(tseries)
#
# JarqueBeraTest <- function(x, robust=TRUE, na.rm=FALSE) {
#
#   # Author: Adrian Trapletti
#
#   if (NCOL(x) > 1)
#       stop("x is not a vector or univariate time series")
#
#   if (na.rm) x <- na.omit(x)
#
#   DNAME <- deparse(substitute(x))
#   n <- length(x)
#   m1 <- sum(x)/n
#   m2 <- sum((x-m1)^2)/n
#   m3 <- sum((x-m1)^3)/n
#   m4 <- sum((x-m1)^4)/n
#   b1 <- (m3/m2^(3/2))^2
#   b2 <- (m4/m2^2)
#   STATISTIC <- n * b1 / 6 + n * (b2 - 3)^2 / 24
#   names(STATISTIC) <- "X-squared"
#   PARAMETER <- 2
#   names(PARAMETER) <- "df"
#   PVAL <- 1 - pchisq(STATISTIC,df = 2)
#   METHOD <- "Jarque Bera Test"
#   structure(list(statistic = STATISTIC,
#                  parameter = PARAMETER,
#                  p.value = PVAL,
#                  method = METHOD,
#                  data.name = DNAME),
#             class = "htest")
# }
#
#





#' (Robust) Jarque Bera Test
#' 
#' This function performs the Jarque-Bera tests of normality either the robust
#' or the classical way.
#' 
#' The test is based on a joint statistic using skewness and kurtosis
#' coefficients. The robust Jarque-Bera (RJB) version of utilizes the robust
#' standard deviation (namely the mean absolute deviation from the median, as
#' provided e. g. by \code{\link{MeanAD}(x, FUN=median)}) to estimate sample
#' kurtosis and skewness. For more details see Gel and Gastwirth (2006). \cr
#'
#' Setting \code{robust} to \code{FALSE} will perform the original Jarque-Bera
#' test (see Jarque, C. and Bera, A (1980)).
#' 
#' @param x a numeric vector of data values.
#' @param robust defines, whether the robust version should be used.  Default
#' is \code{TRUE}.
#' @param method a character string out of \code{chisq} or \code{mc},
#' specifying how the critical values should be obtained. Default is
#' approximated by the chisq-distribution or empirically via Monte Carlo.
#' @param N number of Monte Carlo simulations for the empirical critical values
#' @param na.rm defines if \code{NAs} should be omitted. Default is
#' \code{FALSE}.
#' @return A list with class \code{htest} containing the following components:
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the
#' degrees of freedom.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{type of test was performed.}
#' \item{data.name}{a character
#' string giving the name of the data.}
#' 
#' @note 
#' This function is melted from the \code{jarque.bera.test} (in
#' \code{tseries} package) and the \code{rjb.test} from the package
#' \code{lawstat}.
#' 
#' @author 
#' W. Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @seealso Alternative tests for normality as 
#' \code{\link{shapiro.test}},
#' \code{\link{AndersonDarlingTest}}, 
#' \code{\link{CramerVonMisesTest}},
#' \code{\link{LillieTest}}, 
#' \code{\link{PearsonTest}},
#' \code{\link{ShapiroFranciaTest}}
#' 
#' \code{\link{qqnorm}}, 
#' \code{\link{qqline}} for producing a normal  quantile-quantile plot
#' 
#' @references 
#' Gastwirth, J. L.(1982) \emph{Statistical Properties of A Measure
#' of Tax Assessment Uniformity}, Journal of Statistical Planning and Inference
#' 6, 1-12.\cr
#' 
#' Gel, Y. R. and Gastwirth, J. L. (2008) \emph{A robust modification of the
#' Jarque-Bera test of normality}, Economics Letters 99, 30-32.\cr
#' 
#' Jarque, C. and Bera, A. (1980) \emph{Efficient tests for normality,
#' homoscedasticity and serial independence of regression residuals}, Economics
#' Letters 6, 255-259.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' x <- rnorm(100)    # null hypothesis
#' JarqueBeraTest(x)
#' 
#' x <- runif (100)    # alternative hypothesis
#' JarqueBeraTest(x, robust=TRUE)
#' 
JarqueBeraTest <- function(x, robust=TRUE, method=c("chisq", "mc"), N=0, na.rm=FALSE) {

  method <- match.arg(method)

  if (NCOL(x) > 1) { stop("x is not a vector or univariate time series") }
  if (na.rm) x <- na.omit(x)

  if ((method == "mc") & (N==0)) {
    stop("number of Monte Carlo simulations N should be provided for the empirical critical values")
  }

  DNAME <- deparse(substitute(x))

  ## Calculate the first 4 central moments
  n <- length(x)
  m1 <- sum(x)/n
  m2 <- sum((x - m1)^2)/n
  m3 <- sum((x - m1)^3)/n
  m4 <- sum((x - m1)^4)/n

  ## User can choose the Standard Jarque Bera Test or Robust Jarque Bera Test
  ## Robust Jarque Bera Test is default
  if (!robust) {
    b1 <- (m3/m2^(3/2))^2;
    b2 <- (m4/m2^2);
    statistic <- n * b1/6 + n * (b2 - 3)^2/24

  } else {
    J <- sqrt(pi/2) * mean(abs(x-median(x)))
    J2 <- J^2
    b1 <- (m3/(J2)^(3/2))^2
    b2 <- (m4/(J2)^2)
    vk<-64/n
    vs<-6/n
    ek<-3
    statistic <- b1/vs + (b2 - ek)^2/vk

  }

  if (method == "mc") {
    if (!robust) {
      ## computes empirical critical values for the JB statistic

      jb<-double(N)

      for (k in 1:N) {
        e <- rnorm(length(x), mean=0, sd = sqrt(1))
        m1 <- sum(e)/n
        m2 <- sum((e - m1)^2)/n
        m3 <- sum((e - m1)^3)/n
        m4 <- sum((e - m1)^4)/n
        b1 <- (m3/m2^(3/2))^2
        b2 <- (m4/m2^2)
        vk <- 24/n
        vs <- 6/n
        ek <- 3
        jb[k] <- b1/vs + (b2 - ek)^2/vk
      }

      y <- sort(jb)

      if (statistic >= max(y)) {
        p.value <- 0

      } else if (statistic<=min(y)) {
        p.value <- 1

      } else {
        bn <- which(y==min(y[I(y>=statistic)]))
        an <- which(y==max(y[I(y<statistic)]))
        a <- max(y[I(y<statistic)])
        b <- min(y[I(y>=statistic)])
        pa <- (an - 1) / (N - 1)
        pb <- (bn - 1) / (N - 1)
        alpha <- (statistic-a)/(b-a)
        p.value <- 1-alpha*pb-(1-alpha)*pa
      }

    } else {
      ## computes empirical critical values for the RJB statistic
      rjb <- double(N)

      for (k in 1:N) {
        e <- rnorm(length(x), mean=0, sd = sqrt(1))
        J <- sqrt(pi/2)*mean(abs(e-median(e)))
        J2 <- J^2
        m1 <- sum(e)/n
        m2 <- sum((e - m1)^2)/n
        m3 <- sum((e - m1)^3)/n
        m4 <- sum((e - m1)^4)/n
        b1 <- (m3/(J2)^(3/2))^2
        b2 <- (m4/(J2)^2)
        vk <- 64/n
        vs <- 6/n
        ek <- 3
        rjb[k] <- b1/vs + (b2 - ek)^2/vk
      }

      y <- sort(rjb)

      if (statistic >= max(y)) {
        p.value <- 0

      } else if (statistic <= min(y)) {
        p.value <- 1

      } else {
        bn <- which(y==min(y[I(y>=statistic)]))
        an <- which(y==max(y[I(y<statistic)]))
        a <- max(y[I(y<statistic)])
        b <- min(y[I(y>=statistic)])
        pa <- (an - 1) / (N - 1)
        pb <- (bn - 1) / (N - 1)
        alpha <- (statistic-a)/(b-a)
        p.value <- 1-alpha*pb-(1-alpha)*pa
      }
    }

  } else {
    p.value <- 1 - pchisq(statistic, df = 2)
  }

  METHOD <- ifelse(!robust, "Jarque Bera Test", "Robust Jarque Bera Test")
  STATISTIC=statistic
  names(STATISTIC) <- "X-squared"
  PARAMETER <- 2
  names(PARAMETER) <- "df"

  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = p.value, method = METHOD, data.name = DNAME),
            class = "htest")

}




# PageTest <- function(x) {
#
#   DNAME <- deparse(substitute(x))
#   x <- x[complete.cases(x),]
#
#   rnksum <- apply(apply(x, 1, rank),1, sum)
#   L <- sum(seq_along(rnksum) * rnksum)
#   nc <- ncol(x)
#   nr <- nrow(x)
#   mu <- nr * nc * (nc+1)^2 / 4
#   sig <- nr * nc^2 * (nc+1)^2*(nc-1) / 144
#   z <- (L - mu)/sqrt(sig)
#
#   pval <- pnorm(z, lower.tail = FALSE)
#   RVAL <- list(statistic = c(L = L), p.value = pval, method = "Page test for ordered alternatives",
#       data.name = DNAME)
#   class(RVAL) <- "htest"
#   return(RVAL)
#
# }



# PageTest<-function(x) {

# ### Alternative: package coin
# ### independence_test(scores ~ product | sitting, data = egg_data,
# ### scores = list(product = 1:10),
# ### ytrafo = yt)

# ### http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/pagesL


# if (missing(x))
# stop("Usage: PageTest(x)\n\twhere x is a matrix of ranks")

# dname <- deparse(substitute(x))

# dimx <- dim(x)

# ### This one only requires two dimensions
# page.crit3 <- array(
# c(28,41,54,66,79,91,104,116,128,141,153,165,178,190,202,215,227,239,251,
# NA,42,55,68,81,93,106,119,131,144,156,169,181,194,206,218,231,243,256,
# NA,NA,56,70,83,96,109,121,134,147,160,172,185,197,210,223,235,248,260),
# c(19,3)
# )

# ### the rest require three
# page.crit4plus <- array(
# c(58,84,111,137,163,189,214,240,266,292,317,
# 103,150,197,244,291,338,384,431,477,523,570,
# 166,244,321,397,474,550,625,701,777,852,928,
# 252,370,487,603,719,835,950,1065,1180,1295,1410,
# 362,532,701,869,1037,1204,1371,1537,1703,1868,2035,
# 500,736,971,1204,1436,1668,1900,2131,2361,2592,2822,
# 670,987,1301,1614,1927,2238,2549,2859,3169,3478,3788,
# 60,87,114,141,167,193,220,246,272,298,324,
# 106,155,204,251,299,346,393,441,487,534,581,
# 173,252,331,409,486,563,640,717,793,869,946,
# 261,382,501,620,737,855,972,1088,1205,1321,1437,
# 376,549,722,893,1063,1232,1401,1569,1736,1905,2072,
# 520,761,999,1236,1472,1706,1940,2174,2407,2639,2872,
# 696,1019,1339,1656,1972,2288,2602,2915,3228,3541,3852,
# NA,89,117,145,172,198,225,252,278,305,331,
# 109,160,210,259,307,355,403,451,499,546,593,
# 178,260,341,420,499,577,655,733,811,888,965,
# 269,394,516,637,757,876,994,1113,1230,1348,1465,
# 388,567,743,917,1090,1262,1433,1603,1773,1943,2112,
# 544,790,1032,1273,1512,1750,1987,2223,2459,2694,2929,
# 726,1056,1382,1704,2025,2344,2662,2980,3296,3612,3927),
# c(11,7,3)
# )

# mean.ranks <- apply(x, 2, mean)
# Lval <- NA
# p.table <- NA
# L <- sum(apply(x, 2, sum) * 1:dimx[2])

# if ((dimx[1] > 1 && dimx[1] < 13) && (dimx[2] > 3 && dimx[2] < 11))
# Lval <- page.crit4plus[dimx[1]-1,dimx[2]-3,]

# if ((dimx[1] > 1 && dimx[1] < 21) && dimx[2] == 3)
# Lval <- page.crit3[dimx[1]-1,]

# p.table <-
# ifelse(L > Lval[1],ifelse(L > Lval[2],ifelse(L > Lval[3],"<=.001","<=.01"),"<=.05"),"NS")
# #### print(Lval)

# ### if there was no tabled value, calculate the normal approximation
# if (length(Lval)<2) {
# munum <- dimx[1]*dimx[2]*(dimx[2]+1)*(dimx[2]+1)
# muL <- munum/4
# cat("muL =",muL,"\n")
# sigmaL <- (dimx[1]*dimx[2]*dimx[2]*(dimx[2]*dimx[2]-1)*(dimx[2]*dimx[2]-1))/
# (144*(dimx[2]-1))
# cat("sigmaL =",sigmaL,"\n")
# zL <- ((12*L-3*munum)/(dimx[2]*(dimx[2]-1)))*sqrt((dimx[2]-1)/dimx[1])
# pZ <- pnorm(zL,lower.tail=FALSE)
# } else {
# zL <- NA
# pZ <- NA
# }

# #### ptt <- list(ranks=x, mean.ranks=mean.ranks, L=L, p.table=p.table, Z=zL, pZ=pZ)
# #### class(ptt) <- "PageTest"
# #### return(ptt)

# if (is.na(p.table)) pval <- pZ else pval <- p.table

# RVAL <- list(statistic = c(L = L), p.value = pval, method = "Page test for ordered alternatives",
# data.name = dname)
# class(RVAL) <- "htest"
# return(RVAL)

# }

# print.PageTest<-function(x,...) {

# cat("\nPage test for ordered alternatives\n")
# cat("L =",x$L)

# if (is.na(x$p.table)) {
# plabel<-paste("Z =",x$Z,", p =",x$pZ,sep="",collapse="")
# cat(plabel,x$p.chisq,"\n\n")
# }
# else cat("  p(table) ",x$p.table,"\n\n")
# }




#' Exact Page Test for Ordered Alternatives
#' 
#' Performs a Page test for ordered alternatives using an exact algorithm by
#' Stefan Wellek (1989) with unreplicated blocked data.
#' 
#' \code{PageTest} can be used for analyzing unreplicated complete block
#' designs (i.e., there is exactly one observation in \code{y} for each
#' combination of levels of \code{groups} and \code{blocks}) where the
#' normality assumption may be violated.
#' 
#' The null hypothesis is that apart from an effect of \code{blocks}, the
#' location parameter of \code{y} is the same in each of the \code{groups}.\cr
#' The implemented alternative is, that the location parameter will be
#' monotonly greater along the groups, \cr
#' \eqn{H_{A}: \theta_{1} \le
#' \theta_{2} \le \theta_{3}} ... (where at least one inequality is strict).\cr
#'
#' If the other direction is required, the order of the groups has to be
#' reversed.  \cr
#'\cr
#' The Page test for ordered alternatives is slightly more
#' powerful than the Friedman analysis of variance by ranks.
#' 
#' If \code{y} is a matrix, \code{groups} and \code{blocks} are obtained from
#' the column and row indices, respectively.  \code{NA}'s are not allowed in
#' \code{groups} or \code{blocks}; if \code{y} contains \code{NA}'s,
#' corresponding blocks are removed.
#' 
#' For small values of k (methods) or N (data objects), \samp{PageTest} will
#' calculate the exact p-values.  For \samp{k, N > 15, Inf}, a normal
#' approximation is returned. Only one of these values will be returned.
#' 
#' @aliases PageTest PageTest.default PageTest.formula
#' @param y either a numeric vector of data values, or a data matrix.
#' @param groups a vector giving the group for the corresponding elements of
#' \code{y} if this is a vector; ignored if \code{y} is a matrix.  If not a
#' factor object, it is coerced to one.
#' @param blocks a vector giving the block for the corresponding elements of
#' \code{y} if this is a vector; ignored if \code{y} is a matrix.  If not a
#' factor object, it is coerced to one.
#' @param formula a formula of the form \code{a ~ b | c}, where \code{a},
#' \code{b} and \code{c} give the data values and corresponding groups and
#' blocks, respectively.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the L-statistic with names attribute
#' \dQuote{L}.}
#' \item{p.value}{the p-value of the test.}
#' \item{method}{the
#' character string \code{"Page test for ordered alternatives"}.}
#' \item{data.name}{a character string giving the names of the data.}
#' @note Special thanks to Prof. S. Wellek for porting old GAUSS code to R.
#' 
#' @author 
#' Stefan Wellek <stefan.wellek@@zi-mannheim.de> (exact p-values),
#' Andri Signorell <andri@@signorell.net> (interface) (strongly based on R-Core
#' code)
#' 
#' @seealso 
#' \code{\link{friedman.test}}
#' 
#' @references Page, E. (1963): Ordered hypotheses for multiple treatments: A
#' significance test for linear ranks. \emph{Journal of the American
#' Statistical Association}, 58, 216-230.
#' 
#' Siegel, S. & Castellan, N. J. Jr. (1988): \emph{Nonparametric statistics for
#' the behavioral sciences}. Boston, MA: McGraw-Hill.
#' 
#' Wellek, S. (1989): Computing exact p-values in Page's nonparametric test
#' against trend. \emph{Biometrie und Informatik in Medizin und Biologie 20},
#' 163-170
#' 
#' @keywords htest
#'
#' @examples
#' 
#'  # Craig's data from Siegel & Castellan, p 186
#'  soa.mat <- matrix(c(.797,.873,.888,.923,.942,.956,
#'   .794,.772,.908,.982,.946,.913,
#'   .838,.801,.853,.951,.883,.837,
#'   .815,.801,.747,.859,.887,.902), nrow=4, byrow=TRUE)
#'  PageTest(soa.mat)
#'  
#' 
#' # Duller, pg. 236 
#' pers <- matrix(c(
#' 1, 72, 72, 71.5, 69, 70, 69.5, 68, 68, 67, 68,
#' 2, 83, 81, 81, 82, 82.5, 81, 79, 80.5, 80, 81,
#' 3, 95, 92, 91.5, 89, 89, 90.5, 89, 89, 88, 88,
#' 4, 71, 72, 71, 70.5, 70, 71, 71, 70, 69.5, 69,
#' 5, 79, 79, 78.5, 77, 77.5, 78, 77.5, 76, 76.5, 76,
#' 6, 80, 78.5, 78, 77, 77.5, 77, 76, 76, 75.5, 75.5
#' ), nrow=6, byrow=TRUE) 
#' 
#' colnames(pers) <- c("person", paste("week",1:10))
#' 
#' # Alternative: week10 < week9 < week8 ... 
#' PageTest(pers[, 11:2])
#' 
#' 
#' # Sachs, pg. 464
#' 
#' pers <- matrix(c(
#'   3,2,1,4,
#'   4,2,3,1,
#'   4,1,2,3,
#'   4,2,3,1,
#'   3,2,1,4,
#'   4,1,2,3,
#'   4,3,2,1,
#'   3,1,2,4,
#'   3,1,4,2), 
#'   nrow=9, byrow=TRUE, dimnames=list(1:9, LETTERS[1:4]))  
#' 
#' # Alternative: B < C < D < A
#' PageTest(pers[, c("B","C","D","A")])
#' 
#' 
#' # long shape and formula interface
#' plng <- data.frame(expand.grid(1:9, c("B","C","D","A")), 
#'                    as.vector(pers[, c("B","C","D","A")]))
#' colnames(plng) <- c("block","group","x")
#' 
#' PageTest(plng$x, plng$group, plng$block)
#' 
#' PageTest(x ~ group | block, data = plng)
#' 
#' 
#' 
#' score <- matrix(c(
#'   3,4,6,9,
#'   4,3,7,8,
#'   3,4,4,6,
#'   5,6,8,9,
#'   4,4,9,9,
#'   6,7,11,10
#'   ), nrow=6, byrow=TRUE) 
#' 
#' PageTest(score)
#' 
PageTest <- function(y, ...) {
  UseMethod("PageTest")
}


PageTest.default <- function(y, groups, blocks, ...) {

  p.page <- function(k, n, L) {

    qvec <- .PageDF[k][[1]]
    f1 <- qvec

    for (i in 1:(n-1)) {
      erg <- convolve(f1, qvec, conj = TRUE, type = "open")
      f1 <- erg
    }
    p <- cumsum(erg)[n * k * (k+1) * (2*k+1)/6 + 1 - L]
    return(p)

  }


  DNAME <- deparse(substitute(y))
  if (is.matrix(y)) {
    groups <- factor(c(col(y)))
    blocks <- factor(c(row(y)))
  }
  else {
    if (any(is.na(groups)) || any(is.na(blocks)))
      stop("NA's are not allowed in 'groups' or 'blocks'")
    if (any(diff(c(length(y), length(groups), length(blocks))) !=
            0L))
      stop("'y', 'groups' and 'blocks' must have the same length")
    DNAME <- paste(DNAME, ", ", deparse(substitute(groups)),
                   " and ", deparse(substitute(blocks)), sep = "")
    if (any(table(groups, blocks) != 1))
      stop("not an unreplicated complete block design")
    groups <- factor(groups)
    blocks <- factor(blocks)
    o <- order(groups, blocks)
    y <- y[o]
    groups <- groups[o]
    blocks <- blocks[o]
  }
  k <- nlevels(groups)
  y <- matrix(unlist(split(y, blocks)), ncol = k, byrow = TRUE)
  y <- y[complete.cases(y), ]
  n <- nrow(y)


  rnksum <- apply(apply(y, 1, rank), 1, sum)
  L <- sum(seq_along(rnksum) * rnksum)
  nc <- ncol(y)
  nr <- nrow(y)

  if (nc < 16) {
    pval <- p.page(k=nc, n=nr, L=L)
  } else {
    mu <- nr * nc * (nc + 1)^2/4
    # sig <- nr * nc^2 * (nc + 1)^2 * (nc - 1)/144
    sigma <- nr * nc^2 * (nc+1) * (nc^2-1) / 144
    z <- (L - mu)/sqrt(sigma)
    pval <- pnorm(z, lower.tail = FALSE)

  }

  structure(list(statistic = c(L = L), p.value = pval, method = "Page test for ordered alternatives",
                 data.name = DNAME),
            class = "htest")
}


PageTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula))
    stop("formula missing")
  if ((length(formula) != 3L) || (length(formula[[3L]]) !=
                                  3L) || (formula[[3L]][[1L]] != as.name("|")) || (length(formula[[3L]][[2L]]) !=
                                                                                   1L) || (length(formula[[3L]][[3L]]) != 1L))
    stop("incorrect specification for 'formula'")
  formula[[3L]][[1L]] <- as.name("+")
  m <- match.call(expand.dots = FALSE)
  m$formula <- formula
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " and ")
  names(mf) <- NULL
  y <- DoCall("PageTest", as.list(mf))
  y$data.name <- DNAME
  y

}





#' Cochran's Q test
#' 
#' Perform the Cochran's Q test for unreplicated randomized block design
#' experiments with a binary response variable and paired data.
#' 
#' CochranQTest can be used for analyzing unreplicated complete block designs
#' (i.e., there is exactly one binary observation in y for each combination of
#' levels of groups and blocks) where the normality assumption may be violated.
#' 
#' The null hypothesis is that apart from an effect of blocks, the location
#' parameter of y is the same in each of the groups.
#' 
#' If y is a matrix, groups and blocks are obtained from the column and row
#' indices, respectively.  NA's are not allowed in groups or blocks; if y
#' contains NA's, corresponding blocks are removed.
#' 
#' Note that Cochran's Q Test is analogue to the Friedman test with 0, 1 coded
#' response. This is used here for a simple implementation.
#' 
#' @aliases CochranQTest CochranQTest.default CochranQTest.formula
#' @param y either a numeric vector of data values, or a data matrix.
#' @param groups a vector giving the group for the corresponding elements of y
#' if this is a vector; ignored if y is a matrix. If not a factor object, it is
#' coerced to one.
#' @param blocks a vector giving the block for the corresponding elements of y
#' if this is a vector; ignored if y is a matrix. If not a factor object, it is
#' coerced to one.
#' @param formula a formula of the form \code{a ~ b | c}, where a, b and c give
#' the data values and corresponding groups and blocks, respectively.
#' @param data an optional matrix or data frame (or similar: see model.frame)
#' containing the variables in the formula formula. By default the variables
#' are taken from environment(formula).
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to getOption("na.action").
#' @param \dots further arguments to be passed to or from methods.
#' @return
#' 
#' A list with class \code{htest} containing the following components:
#' 
#' \item{statistic}{the value of Cochran's chi-squared statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared
#' distribution of the test statistic.}
#' \item{p.value}{the p-value of the
#' test.}
#' \item{method}{the character string "Cochran's Q-Test".}
#' \item{data.name}{a character string giving the names of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>
#'
#' @examples
#' 
#' # example in: 
#' # http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf
#' # pp. S. 1824
#' 
#' # use expand.grid, xtabs and Untable to create the dataset
#' d.frm <- Untable(xtabs(c(6,2,2,6,16,4,4,6) ~ ., 
#'     expand.grid(rep(list(c("F","U")), times=3))), 
#'     colnames = LETTERS[1:3])
#' 
#' # rearrange to long shape    
#' d.long <- reshape(d.frm, varying=1:3, times=names(d.frm)[c(1:3)], 
#'                   v.names="resp", direction="long")
#' 
#' 
#' # after having done the hard work of data organisation, performing the test is a piece of cake....
#' CochranQTest(resp ~ time | id, data=d.long)
#' 
CochranQTest <- function(y, ...) {

  # Cochran's Q Test is analogue to the friedman.test with 0,1 coded response

  res <- friedman.test(y, ...)
  attr(res$statistic, "names") <- "Q"
  res$method <- "Cochran's Q test"
  return(res)
}

CochranQTest.default <- function(y, groups, blocks, ...) {
  res <- friedman.test(y, groups, blocks, ...)
  attr(res$statistic, "names") <- "Q"
  res$method <- "Cochran's Q test"
  return(res)
}

CochranQTest.formula <- function(formula, data, subset, na.action, ...) {
  res <- friedman.test(formula, data, subset, na.action, ...)
  attr(res$statistic, "names") <- "Q"
  res$method <- "Cochran's Q test"
  return(res)
}




#' Mantel-Haenszel Chi-Square Test
#' 
#' The Mantel-Haenszel chi-square statistic tests the alternative hypothesis
#' that there is a linear association between the row variable and the column
#' variable. Both variables must lie on an ordinal scale. 
#' 
#' The statistic is computed as \eqn{ Q_{MH} = (n-1) \cdot r^2}, where
#' \eqn{r^2} is the Pearson correlation between the row variable and the column
#' variable. The Mantel-Haenszel chi-square statistic use the scores specified
#' by srow and scol. Under the null hypothesis of no association, \eqn{Q_{MH}}
#' has an asymptotic chi-square distribution with one degree of freedom.
#' 
#' @param x a frequency table or a matrix.  
#' @param srow scores for the row variable, defaults to 1:nrow. 
#' @param scol scores for the colummn variable, defaults to 1:ncol. 
#' 
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value the Mantel-Haenszel chi-squared test
#' statistic.}
#' \item{parameter}{the degrees of freedom of the approximate
#' chi-squared distribution of the test statistic.}
#' \item{p.value}{the p-value
#' for the test.}
#' \item{method}{a character string indicating the type of test
#' performed.}
#' \item{data.name}{a character string giving the name(s) of the
#' data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso
#' \code{\link{chisq.test}}, for calculating correlation of a table:
#' \code{\link[boot]{corr}}
#' 
#' @references Agresti, A. (2002) \emph{Categorical Data Analysis}. John Wiley
#' & Sons, pp 86 ff. 
#' @keywords htest
#'
#' @examples
#' 
#' ## A r x c table  Agresti (2002, p. 57) Job Satisfaction
#' Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
#'               dimnames = list(income = c("< 15k", "15-25k", "25-40k", "> 40k"),
#'                               satisfaction = c("VeryD", "LittleD", "ModerateS", "VeryS"))
#'        )
#' 
#' MHChisqTest(Job, srow=c(7.5,20,32.5,60))
#' 
MHChisqTest <- function(x, srow=1:nrow(x), scol=1:ncol(x)) {

  # calculates Mantel-Haenszel Chisquare test

  # check for rxc 2-dim matrix
  p <- (d <- dim(x))[1L]
  if (!is.numeric(x) || length(d) != 2L)
    stop("'x' is not a rxc numeric matrix")

  DNAME <- deparse(substitute(x))

  STATISTIC <- (sum(x) - 1) * corr(d=CombPairs(srow, scol), as.vector(x))^2
  PARAMETER <- 1
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  METHOD <- "Mantel-Haenszel Chi-Square"

  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME), class = "htest")
}





#' G-Test for Count Data
#' 
#' \code{GTest} performs chi-squared contingency table tests and
#' goodness-of-fit tests.
#' 
#' The G-test is also called "Likelihood Ratio Test" and is asymptotically
#' equivalent to the Pearson ChiSquare-test but not usually used when analyzing
#' 2x2 tables. It is used in logistic regression and loglinear modeling which
#' involves contingency tables. The G-test is also reported in the standard
#' summary of \code{Desc} for tables.
#' 
#' If \code{x} is a matrix with one row or column, or if \code{x} is a vector
#' and \code{y} is not given, then a \emph{goodness-of-fit test} is performed
#' (\code{x} is treated as a one-dimensional contingency table).  The entries
#' of \code{x} must be non-negative integers.  In this case, the hypothesis
#' tested is whether the population probabilities equal those in \code{p}, or
#' are all equal if \code{p} is not given.
#' 
#' If \code{x} is a matrix with at least two rows and columns, it is taken as a
#' two-dimensional contingency table: the entries of \code{x} must be
#' non-negative integers.  Otherwise, \code{x} and \code{y} must be vectors or
#' factors of the same length; cases with missing values are removed, the
#' objects are coerced to factors, and the contingency table is computed from
#' these.  Then G-test is performed on the null hypothesis that the joint
#' distribution of the cell counts in a 2-dimensional contingency table is the
#' product of the row and column marginals.
#' 
#' % If \code{simulate.p.value} is \code{FALSE}, the p-value is computed % from
#' the asymptotic chi-squared distribution of the test statistic; % continuity
#' correction is only used in the 2-by-2 case (if \code{correct} % is
#' \code{TRUE}, the default).  Otherwise the p-value is computed for a % Monte
#' Carlo test (Hope, 1968) with \code{B} replicates.
#' 
#' % In the contingency table case simulation is done by random sampling % from
#' the set of all contingency tables with given marginals, and works % only if
#' the marginals are strictly positive.  (A C translation of the % algorithm of
#' Patefield (1981) is used.)  Continuity correction is % never used, and the
#' statistic is quoted without it.  Note that this is % not the usual sampling
#' situation assumed for the chi-squared test but % rather that for Fisher's
#' exact test.
#' 
#' % In the goodness-of-fit case simulation is done by random sampling from %
#' the discrete distribution specified by \code{p}, each sample being % of size
#' \code{n = sum(x)}.  This simulation is done in \R and may be % slow.  TOI
#' Yates' correction taken from Mike Camann's 2x2 G-test function.  GOF Yates'
#' correction as described in Zar (2000)
#' 
#' @param x a numeric vector or matrix. \code{x} and \code{y} can also both be
#' factors.
#' @param y a numeric vector; ignored if \code{x} is a matrix.  If \code{x} is
#' a factor, \code{y} should be a factor of the same length.
#' @param correct one out of \code{"none"} (default), \code{"williams"},
#' \code{"yates"} . See Details.
#' @param p a vector of probabilities of the same length of \code{x}.  An error
#' is given if any entry of \code{p} is negative.
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{statistic}{the value the chi-squared test statistic.}
#' \item{parameter}{the degrees of freedom of the approximate chi-squared
#' distribution of the test statistic, \code{NA} if the p-value is computed by
#' Monte Carlo simulation.}
#' \item{p.value}{the p-value for the test.}
#' \item{method}{a character string indicating the type of test performed, and
#' whether Monte Carlo simulation or continuity correction was used.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' \item{observed}{the observed counts.}
#' \item{expected}{the expected counts
#' under the null hypothesis.}
#' 
#' @author 
#' Pete Hurd <phurd@@ualberta.ca>
#' 
#' @seealso 
#' \code{\link{chisq.test}}.
#' 
#' 
#' @references 
#' Hope, A. C. A. (1968) A simplified Monte Carlo significance test
#' procedure. \emph{J. Roy, Statist. Soc. B} \bold{30}, 582--598.
#' 
#' Patefield, W. M. (1981) Algorithm AS159. An efficient method of generating
#' r x c tables with given row and column totals. \emph{Applied Statistics}
#' \bold{30}, 91--97.
#' 
#' Agresti, A. (2007) \emph{An Introduction to Categorical Data Analysis, 2nd
#' ed.}, New York: John Wiley & Sons.  Page 38.
#' 
#' Sokal, R. R., F. J. Rohlf (2012) \emph{Biometry: the principles and practice
#' of statistics in biological research}. 4th edition. W. H. Freeman and Co.:
#' New York. 937 pp.
#' 
#' @keywords htest distribution
#'
#' @examples
#' 
#' 
#' ## From Agresti(2007) p.39
#' M <- as.table(rbind(c(762, 327, 468), c(484,239,477)))
#' dimnames(M) <- list(gender=c("M","F"),
#'                     party=c("Democrat","Independent", "Republican"))
#' 
#' (Xsq <- GTest(M))   # Prints test summary
#' 
#' Xsq$observed        # observed counts (same as M)
#' Xsq$expected        # expected counts under the null
#' 
#' 
#' ## Testing for population probabilities
#' ## Case A. Tabulated data
#' x <- c(A = 20, B = 15, C = 25)
#' GTest(x)
#' GTest(as.table(x))             # the same
#' x <- c(89,37,30,28,2)
#' p <- c(40,20,20,15,5)
#' try(
#' GTest(x, p = p)                # gives an error
#' )
#' # works
#' p <- c(0.40,0.20,0.20,0.19,0.01)
#' # Expected count in category 5
#' # is 1.86 < 5 ==> chi square approx.
#' GTest(x, p = p)                #         maybe doubtful, but is ok!
#' 
#' ## Case B. Raw data
#' x <- trunc(5 * runif (100))
#' GTest(table(x))                # NOT 'GTest(x)'!
#' 
GTest <- function(x, y = NULL, correct=c("none", "williams", "yates"), p = rep(1/length(x), length(x))) {


  # Log-likelihood tests of independence & goodness of fit
  # Does Williams' and Yates' correction
  # does Monte Carlo simulation of p-values, via gtestsim.c
  #
  # G & q calculation from Sokal & Rohlf (1995) Biometry 3rd ed.
  # TOI Yates' correction taken from Mike Camann's 2x2 G-test fn.
  # GOF Yates' correction as described in Zar (2000)
  # more stuff taken from ctest's chisq.test()
  #
  # TODO:
  # 1) Beautify
  # 2) Add warnings for violations
  # 3) Make appropriate corrections happen by default
  #
  # V3.3 Pete Hurd Sept 29 2001. phurd@ualberta.ca


  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1)
      x <- as.vector(x)
  }
  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("x and y must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if ((nlevels(x) < 2) || (nlevels(y) < 2))
      stop("x and y must have at least 2 levels")
    x <- table(x, y)
  }
  if (any(x < 0) || any(is.na(x)))
    stop("all entries of x must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of x must be positive")

  correct <- match.arg(correct)

  #If x is matrix, do test of independence
  if (is.matrix(x)) {
    #Test of Independence
    nrows<-nrow(x)
    ncols<-ncol(x)
    if (correct=="yates") { # Do Yates' correction?
      if (dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
        stop("Yates' correction requires a 2 x 2 matrix")
      if ((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
      {
        #         x[1,1] <- x[1,1] - 0.5
        #         x[2,2] <- x[2,2] - 0.5
        #         x[1,2] <- x[1,2] + 0.5
        #         x[2,1] <- x[2,1] + 0.5
        #   this can be done quicker: 14.5.2015 AS
        x <- x + 0.5
        diag(x) <- diag(x) - 1

      } else {

        x <- x - 0.5
        diag(x) <- diag(x) + 1

        #         x[1,1] <- x[1,1] + 0.5
        #         x[2,2] <- x[2,2] + 0.5
        #         x[1,2] <- x[1,2] - 0.5
        #         x[2,1] <- x[2,1] - 0.5
      }
    }

    sr <- apply(x,1,sum)
    sc <- apply(x,2,sum)
    E <- outer(sr,sc, "*")/n
    # are we doing a monte-carlo?
    # no monte carlo GOF?
    #     if (simulate.p.value) {
    #       METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
    #       tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
    #                 as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
    #                 as.double(E), integer(nrows * ncols), double(n+1),
    #                 integer(ncols), results=double(B), PACKAGE= "ctest")
    #       g <- 0
    #       for (i in 1:nrows) {
    #         for (j in 1:ncols) {
    #           if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
    #         }
    #       }
    #       STATISTIC <- G <- 2 * g
    #       PARAMETER <- NA
    #       PVAL <- sum(tmp$results >= STATISTIC)/B
    #     }
    #     else {
    # no monte-carlo
    # calculate G
    g <- 0
    for (i in 1:nrows) {
      for (j in 1:ncols) {
        if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
      }
      # }
      q <- 1
      if (correct=="williams") { # Do Williams' correction
        row.tot <- col.tot <- 0
        for (i in 1:nrows) { row.tot <- row.tot + 1/(sum(x[i,])) }
        for (j in 1:ncols) { col.tot <- col.tot + 1/(sum(x[,j])) }
        q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
      }
      STATISTIC <- G <- 2 * g / q
      PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
      PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
      if (correct=="none")
        METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
      if (correct=="williams")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
      if (correct=="yates")
        METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
    }
  }
  else {
    # x is not a matrix, so we do Goodness of Fit
    METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
    if (length(x) == 1)
      stop("x must at least have 2 elements")
    if (length(x) != length(p))
      stop("x and p must have the same number of elements")
    E <- n * p

    if (correct=="yates") { # Do Yates' correction
      if (length(x)!=2)
        stop("Yates' correction requires 2 data values")
      if ( (x[1]-E[1]) > 0.25) {
        x[1] <- x[1]-0.5
        x[2] <- x[2]+0.5
      }
      else if ( (E[1]-x[1]) > 0.25) {
        x[1] <- x[1]+0.5
        x[2] <- x[2]-0.5
      }
    }
    names(E) <- names(x)
    g <- 0
    for (i in 1:length(x)) {
      if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
    }
    q <- 1
    if (correct=="williams") { # Do Williams' correction
      q <- 1+(length(x)+1)/(6*n)
    }
    STATISTIC <- G <- 2*g/q
    PARAMETER <- length(x) - 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  }
  names(STATISTIC) <- "G"
  names(PARAMETER) <- "X-squared df"
  names(PVAL) <- "p.value"
  structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
                 method=METHOD,data.name=DNAME, observed=x, expected=E),
            class="htest")
}








#' Stuart-Maxwell Marginal Homogeneity Test 
#' 
#' This function computes the marginal homogeneity test for a \eqn{k \times
#' k}{k x k} matrix of assignments of objects to k categories or an \eqn{n
#' \times 2 \times k}{n x 2} matrix of category scores for n data objects by
#' two raters. The statistic is distributed as chi-square with k-1 degrees of
#' freedom. \cr
#' It can be viewed as an extention of McNemar test to \eqn{k \times k}{k x k}
#' table
#' 
#' The null is that the probabilities of being classified into cells [i,j] and
#' [j,i] are the same.
#' 
#' If x is a matrix, it is taken as a two-dimensional contingency table, and
#' hence its entries should be nonnegative integers. Otherwise, both x and y
#' must be vectors or factors of the same length. Incomplete cases are removed,
#' vectors are coerced into factors, and the contingency table is computed from
#' these.
#' 
#' @param x either a 2-way contingency table in matrix form, or a factor
#' object. 
#' @param y a factor object; ignored if x is a matrix.
#' @return A list with class \code{"htest"} containing the following
#' components: 
#' \item{statistic}{the value of the test statistic.}
#' \item{parameter}{the degrees of freedom.} 
#' \item{p.value}{the p-value of the
#' test.} 
#' \item{method}{a character string indicating what type of test was
#' performed.} 
#' \item{data.name}{a character string giving the name of the
#' data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>, based on Code from Jim Lemon
#' 
#' @seealso 
#' \code{\link{mcnemar.test}}, 
#' \code{\link{chisq.test}},
#' \code{\link{MHChisqTest}}, 
#' \code{\link{BreslowDayTest}} 
#' 
#' @references 
#' Agresti, A. (2002) \emph{Categorical Data Analysis}. John Wiley & Sons, 
#' pp 86 ff.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' hyp <- as.table(matrix(c(20,3,0,10,30,5,5,15,40), nrow=3))
#' StuartMaxwellTest(hyp)
#' 
#' # Source: http://www.john-uebersax.com/stat/mcnemar.htm#stuart
#' mc <- as.table(matrix(c(
#'          732, 1524, 1575, 1577, 1602, 837, 1554, 1437, 
#'          1672, 1600, 841, 1363, 1385, 1484, 1524, 791), nrow=4))
#' 
#' StuartMaxwellTest(mc)
#' 
StuartMaxwellTest <- function(x, y = NULL) {

  # stuart.maxwell.mh computes the marginal homogeneity test for
  # a CxC matrix of assignments of objects to C categories or an
  # nx2 or 2xn matrix of category scores for n data objects by two
  # raters. The statistic is distributed as Chi-square with C-1
  # degrees of freedom.

  # The core code is form Jim Lemon, package concord
  # the intro is taken from mcnemar.test (core)

  if (is.matrix(x)) {
    r <- nrow(x)
    if ((r < 2) || (ncol(x) != r))
      stop("'x' must be square with at least two rows and columns")
    if (any(x < 0) || anyNA(x))
      stop("all entries of 'x' must be nonnegative and finite")
    DNAME <- deparse(substitute(x))
  }
  else {
    if (is.null(y))
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    r <- nlevels(x)
    if ((r < 2) || (nlevels(y) != r))
      stop("'x' and 'y' must have the same number of levels (minimum 2)")
    x <- table(x, y)
  }

  # get the marginals
  rowsums <- rowSums(x)
  colsums <- colSums(x)
  equalsums <- rowsums == colsums

  if (any(equalsums)) {
    # dump any categories with perfect agreement
    x <- x[!equalsums, !equalsums]
    # bail out if too many categories have disappeared
    if (dim(x)[1] < 2) stop("Too many equal marginals, cannot compute")
    # get new marginals
    rowsums <- rowSums(x)
    colsums <- colSums(x)
  }

  # use K-1 marginals
  Kminus1 <- length(rowsums) - 1
  smd <- (rowsums-colsums)[1:Kminus1]
  smS <- matrix(0, nrow=Kminus1, ncol=Kminus1)
  for(i in 1:Kminus1) {
    for(j in 1:Kminus1) {
      if (i == j) smS[i,j] <- rowsums[i] + colsums[j] - 2 * x[i,j]
      else smS[i,j] <- -(x[i,j] + x[j,i])
    }
  }

  STATISTIC <- t(smd) %*% solve(smS) %*% smd

  PARAMETER <- r - 1
  METHOD <- "Stuart-Maxwell test"

  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  names(STATISTIC) <- "chi-squared"
  names(PARAMETER) <- "df"
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER,
               p.value = PVAL, method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)

}






#' Breslow-Day Test for Homogeneity of the Odds Ratios 
#' 
#' Calculates the Breslow-Day test of homogeneity for a \eqn{2 \times 2 \times
#' k}{2 x 2 x k} table, in order to investigate if all \eqn{k} strata have the
#' same OR. If OR is not given, the Mantel-Haenszel estimate is used.
#' 
#' For the Breslow-Day test to be valid, the sample size should be relatively
#' large in each stratum, and at least 80\% of the expected cell counts should
#' be greater than 5. Note that this is a stricter sample size requirement than
#' the requirement for the Cochran-Mantel-Haenszel test for tables, in that
#' each stratum sample size (not just the overall sample size) must be
#' relatively large. Even when the Breslow-Day test is valid, it might not be
#' very powerful against certain alternatives, as discussed in Breslow and Day
#' (1980).
#' 
#' Alternatively, it might be better to cast the entire inference problem into
#' the setting of a logistic regression model. Here, the underlying question of
#' the Breslow-Day test can be answered by investigating whether an interaction
#' term with the strata variable is necessary (e.g. using a likelihood ratio
#' test using the \code{anova} function).
#' 
#' @param x a \eqn{2 \times 2 \times k}{2 x 2 x k} table. 
#' @param OR the odds ratio to be tested against. If left undefined (default)
#' the Mantel-Haenszel estimate will be used.
#' @param correct If TRUE, the Breslow-Day test with Tarone's adjustment is
#' computed, which subtracts an adjustment factor to make the resulting
#' statistic asymptotically chi-square. 
#' 
#' @author 
#' Michael Hoehle <hoehle@@math.su.se> 
#' 
#' @seealso 
#' \code{\link{WoolfTest}} 
#' 
#' @references 
#' Breslow, N. E., N. E. Day (1980) The Analysis of Case-Control
#' Studies \emph{Statistical Methods in Cancer Research: Vol. 1}. Lyon, France,
#' IARC Scientific Publications.
#' 
#' Tarone, R.E. (1985) On heterogeneity tests based on efficient scores,
#' \emph{Biometrika}, 72, pp. 91-95.
#' 
#' Jones, M. P., O'Gorman, T. W., Lemka, J. H., and Woolson, R. F. (1989) A
#' Monte Carlo Investigation of Homogeneity Tests of the Odds Ratio Under
#' Various Sample Size Configurations \emph{Biometrics}, 45, 171-181 \cr
#' 
#' Breslow, N. E. (1996) Statistics in Epidemiology: The Case-Control Study
#' \emph{Journal of the American Statistical Association}, 91, 14-26.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' migraine <- xtabs(freq ~ .,
#'             cbind(expand.grid(treatment=c("active", "placebo"),
#'                               response =c("better", "same"),
#'                               gender   =c("female", "male")),
#'                   freq=c(16, 5, 11, 20, 12, 7, 16, 19))
#'             )
#' 
#' # get rid of gender
#' tab <- xtabs(Freq ~ treatment + response, migraine)
#' Desc(tab)
#' 
#' # only the women
#' female <- migraine[,, 1]
#' Desc(female)
#' 
#' # .. and the men
#' male <- migraine[,, 2]
#' Desc(male)
#' 
#' BreslowDayTest(migraine)
#' BreslowDayTest(migraine, correct = TRUE)
#' 
#' 
#' salary <- array(
#'       c(38, 12, 102, 141, 12, 9, 136, 383),
#'       dim=c(2, 2, 2),
#'       dimnames=list(exposure=c("exposed", "not"),
#'                     disease =c("case", "control"),
#'                     salary  =c("<1000", ">=1000"))
#'                     )
#' 
#' # common odds ratio = 4.028269
#' BreslowDayTest(salary, OR = 4.02)
#' 
BreslowDayTest <- function(x, OR = NA, correct = FALSE) {

  # Function to perform the Breslow and Day (1980) test including the
  # corrected test by Tarone Uses the equations in Lachin (2000),
  # Biostatistical Methods, Wiley, p. 124-125.
  #
  # Programmed by Michael Hoehle <http://www.math.su.se/~hoehle>
  # Code taken originally from a Biostatistical Methods lecture
  # held at the Technical University of Munich in 2008.
  #
  # Params:
  #  x - a 2x2xK contingency table
  #  correct - if TRUE Tarones correction is returned
  #
  # Returns:
  #  a vector with three values
  #   statistic - Breslow and Day test statistic
  #   pval - p value evtl. based on the Tarone test statistic
  #               using a \chi^2(K-1) distribution
  #


  if (is.na(OR)) {
    #Find the common OR based on Mantel-Haenszel
    or.hat.mh <- mantelhaen.test(x)$estimate
  } else {
    or.hat.mh <- OR
  }

  #Number of strata
  K <- dim(x)[3]
  #Value of the Statistic
  X2.HBD <- 0
  #Value of aj, tildeaj and Var.aj
  a <- tildea <- Var.a <- numeric(K)

  for (j in 1:K) {
    #Find marginals of table j
    mj <- apply(x[,,j], MARGIN=1, sum)
    nj <- apply(x[,,j], MARGIN=2, sum)

    #Solve for tilde(a)_j
    coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),
              1-or.hat.mh)
    sols <- Re(polyroot(coef))
    #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
    tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
    #Observed value
    aj <- x[1,1,j]

    #Determine other expected cell entries
    tildebj <- mj[1] - tildeaj
    tildecj <- nj[1] - tildeaj
    tildedj <- mj[2] - tildecj

    #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
    Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)

    #Compute contribution
    X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)

    #Assign found value for later computations
    a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
  }

  # Compute Tarone corrected test
  # Add on 2015: The original equation from the 2008 lecture is incorrect
  # as pointed out by Jean-Francois Bouzereau.
  # X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
  X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.a) )

  DNAME <- deparse(substitute(x))

  STATISTIC <- if (correct) X2.HBDT else X2.HBD
  PARAMETER <- K - 1
  # Compute p-value based on the Tarone corrected test
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- if (correct) "Breslow-Day Test on Homogeneity of Odds Ratios (with Tarone correction)" else
    "Breslow-Day test on Homogeneity of Odds Ratios"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME
  ), class = "htest")

}




# BreslowDayTest <- function(x, OR = NA, correct = FALSE) {
#
#   # Function to perform the Breslow and Day (1980) test including
#   # the corrected test by Tarone
#   # Uses the equations in Lachin (2000) p. 124-125.
#   #
#   # Programmed by Michael Hoehle <http://www-m4.ma.tum.de/pers/hoehle>
#   # Note that the results of the Tarone corrected test do
#   # not correspond to the numbers in the Lachin book...
#   #
#   # Params:
#   #  x - a 2x2xK contingency table
#   #  correct - if TRUE Tarones correction is returned
#   #
#   # Returns:
#   #  a vector with three values
#   #   statistic - Breslow and Day test statistic
#   #   pval - p value evtl. based on the Tarone test statistic
#   #               using a \chi^2(K-1) distribution
#   #
#
#
#   if (is.na(OR)) {
#     #Find the common OR based on Mantel-Haenszel
#     or.hat.mh <- mantelhaen.test(x)$estimate
#   } else {
#     or.hat.mh <- OR
#   }
#
#   #Number of strata
#   K <- dim(x)[3]
#   #Value of the Statistic
#   X2.HBD <- 0
#   #Value of aj, tildeaj and Var.aj
#   a <- tildea <- Var.a <- numeric(K)
#
#   for (j in 1:K) {
#     #Find marginals of table j
#     mj <- apply(x[,,j], MARGIN=1, sum)
#     nj <- apply(x[,,j], MARGIN=2, sum)
#
#     #Solve for tilde(a)_j
#     coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),
#               1-or.hat.mh)
#     sols <- Re(polyroot(coef))
#     #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
#     tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
#     #Observed value
#     aj <- x[1,1,j]
#
#     #Determine other expected cell entries
#     tildebj <- mj[1] - tildeaj
#     tildecj <- nj[1] - tildeaj
#     tildedj <- mj[2] - tildecj
#
#     #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
#     Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)
#
#     #Compute contribution
#     X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)
#
#     #Assign found value for later computations
#     a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
#   }
#
#   # Compute Tarone corrected test
#   # This is incorrect as correctly pointed out by Jean-Francois Bouzereau..
#   # X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )
#   X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.a) )
#
#   DNAME <- deparse(substitute(x))
#
#   STATISTIC <- if (correct) X2.HBDT else X2.HBD
#   PARAMETER <- K - 1
#   # Compute p-value based on the Tarone corrected test
#   PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
#   METHOD <- if (correct) "Breslow-Day Test on Homogeneity of Odds Ratios (with Tarone correction)" else
#     "Breslow-Day test on Homogeneity of Odds Ratios"
#   names(STATISTIC) <- "X-squared"
#   names(PARAMETER) <- "df"
#   structure(list(statistic = STATISTIC, parameter = PARAMETER,
#                  p.value = PVAL, method = METHOD, data.name = DNAME
#   ), class = "htest")
#
# }



# the VCD package (available via CRAN) has a function called woolf_test()



#' Woolf Test For Homogeneity in 2x2xk Tables
#' 
#' Test for homogeneity on \eqn{2 \times 2 \times k}{2 x 2 x k} tables over
#' strata (i.e., whether the log odds ratios are the same in all strata).
#' 
#' 
#' @param x a \eqn{2 \times 2 \times k}{2 x 2 x k} table, where the last
#' dimension refers to the strata.
#' @return A list of class \code{"htest"} containing the following components:
#' \item{statistic}{the chi-squared test statistic.}
#' \item{parameter}{degrees
#' of freedom of the approximate chi-squared distribution of the test
#' statistic.}
#' \item{p.value}{\eqn{p}-value for the test.}
#' \item{method}{a
#' character string indicating the type of test performed.}
#' \item{data.name}{a
#' character string giving the name(s) of the data.}
#' \item{observed}{the
#' observed counts.}
#' \item{expected}{the expected counts under the null
#' hypothesis.}
#' @note This function was previously published as \code{woolf_test()} in the
#' \pkg{vcd} package and has been integrated here without logical changes.
#' 
#' @author 
#' David Meyer, Achim Zeileis, Kurt Hornik, Michael Friendly
#' 
#' @seealso 
#' \code{\link{mantelhaen.test}}, 
#' \code{\link{BreslowDayTest}}
#' 
#' @references
#'  Woolf, B. 1955: On estimating the relation between blood group
#' and disease. \emph{Ann. Human Genet.} (London) \bold{19}, 251-253.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' migraine <- xtabs(freq ~ .,
#'             cbind(expand.grid(treatment=c("active","placebo"),
#'                                response=c("better","same"),
#'                                gender=c("female","male")),
#'                   freq=c(16,5,11,20,12,7,16,19))
#'             )
#' 
#' WoolfTest(migraine)
#' 
WoolfTest <- function(x) {

  DNAME <- deparse(substitute(x))
  if (any(x == 0))
    x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1] * x[2,2]) / (x[1,2] * x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  o <- log(or)
  e <- weighted.mean(log(or), w)
  STATISTIC <- sum(w * (o - e)^2)
  PARAMETER <- k - 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Woolf Test on Homogeneity of Odds Ratios (no 3-Way assoc.)"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = o,
                 expected = e), class = "htest")

}




#' Lehmacher's Test for Marginal Homogenity
#' 
#' Performs Lehmacher's chi-squared test for marginal homogenity in a symmetric
#' two-dimensional contingency table.
#' 
#' The null is that the probabilities of being classified into cells [i,j] and
#' [j,i] are the same.
#' 
#' If x is a matrix, it is taken as a two-dimensional contingency table, and
#' hence its entries should be nonnegative integers. Otherwise, both x and y
#' must be vectors or factors of the same length. Incomplete cases are removed,
#' vectors are coerced into factors, and the contingency table is computed from
#' these.
#' 
#' @aliases LehmacherTest print.mtest
#' @param x either a two-dimensional contingency table in matrix form, or a
#' factor object.
#' @param y a factor object; ignored if x is a matrix.
#' @param digits a non-null value for digits specifies the minimum number of
#' significant digits to be printed in values. See details in
#' \code{\link{print.default}}.
#' @param \dots further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' @return A list with class \code{"mtest"} containing the following
#' components:
#' \item{statistic}{a vector with the value of the test
#' statistics.}
#' \item{parameter}{the degrees of freedom, which is always 1 in
#' LehmacherTest.}
#' \item{p.value}{a vector with the p-values of the single
#' tests.}
#' \item{p.value.corr}{a vector with the "hochberg" adjusted p-values
#' of the single tests. (See \code{\link{p.adjust}})}
#' \item{method}{a character
#' string indicating what type of test was performed.}
#' \item{data.name}{a
#' character string giving the name of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso 
#' \code{\link{mcnemar.test}} (resp. BowkerTest for a CxC-matrix),
#' \code{\link{StuartMaxwellTest}}, 
#' \code{\link{WoolfTest}}
#' 
#' @references 
#' Lehmacher, W. (1980) Simultaneous sign tests for marginal
#' homogeneity of square contingency tables \emph{Biometrical Journal}, Volume
#' 22, Issue 8, pages 795-798
#' 
#' @keywords htest
#'
#' @examples
#' 
#' x <- matrix(c(400,40,20,10, 
#'               50,300,60,20, 
#'               10,40,120,5, 
#'               5,90,50,80), nrow=4, byrow=TRUE)
#'               
#' LehmacherTest(x)
#' 
LehmacherTest <- function(x, y = NULL) {

  if (is.matrix(x)) {
    r <- nrow(x)
    if ((r < 2) || (ncol(x) != r))
      stop("'x' must be square with at least two rows and columns")
    if (any(x < 0) || anyNA(x))
      stop("all entries of 'x' must be nonnegative and finite")
    DNAME <- deparse(substitute(x))
  }
  else {
    if (is.null(y))
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    r <- nlevels(x)
    if ((r < 2) || (nlevels(y) != r))
      stop("'x' and 'y' must have the same number of levels (minimum 2)")
    x <- table(x, y)
  }

  rsum <- rowSums(x)
  csum <- colSums(x)

  STATISTIC <- (rsum-csum)^2 / (rsum + csum - 2*diag(x))
  PARAMETER <- 1
  PVAL <- 1 - pchisq(STATISTIC, PARAMETER)
  METHOD <- "Lehmacher-Test on Marginal Homogeneity"
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, p.value.corr = p.adjust(PVAL, "hochberg"),
                 method = METHOD, data.name = DNAME),
            class = "mtest")

}


print.mtest <- function(x, digits = 4L, ...) {

  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")

  out <- character()
  out <- cbind(format(round(x$statistic, 4)), format.pval(x$p.value, digits = digits),
               format.pval(x$p.value.corr, digits = digits),
               symnum(x$p.value.corr, corr = FALSE, na = FALSE,
                      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                      symbols = c("***", "**", "*", ".", " ")))
  colnames(out) <- c("X-squared", "pval", "pval adj", " ")
  rownames(out) <- if (is.null(rownames(x))) 1:length(x$statistic) else rownames(x)
  print.default(out, digits = 3, quote = FALSE, right = TRUE)

  cat("\n")
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  invisible(x)
}






#' Cochran-Armitage Test for Trend
#' 
#' Perform a Cochran Armitage test for trend in binomial proportions across the
#' levels of a single variable. This test is appropriate only when one variable
#' has two levels and the other variable is ordinal. The two-level variable
#' represents the response, and the other represents an explanatory variable
#' with ordered levels. The null hypothesis is the hypothesis of no trend,
#' which means that the binomial proportion is the same for all levels of the
#' explanatory variable.
#' 
#' 
#' @param x a frequency table or a matrix.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"increasing"} or
#' \code{"decreasing"}. You can specify just the initial letter. 
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ the z-statistic of the test.}
#' \item{parameter}{ the
#' dimension of the table.}
#' \item{p.value}{ the p-value for the test.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{the character string \dQuote{Cochran-Armitage
#' test for trend}.}
#' \item{data.name}{a character string giving the names of
#' the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> strongly based on code from
#' Eric Lecoutre <lecoutre@@stat.ucl.ac.be>\cr
#' \url{https://stat.ethz.ch/pipermail/r-help/2005-July/076371.html}
#' 
#' @seealso \code{\link{prop.trend.test}} 
#' @references 
#' Agresti, A. (2002) \emph{Categorical Data Analysis}. John Wiley & Sons
#' 
#' @keywords htest
#'
#' @examples
#' 
#' # http://www.lexjansen.com/pharmasug/2007/sp/sp05.pdf, pp. 4
#' dose <- matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2, 
#'   dimnames=list(resp=0:1, dose=0:3))
#' Desc(dose)
#' 
#' CochranArmitageTest(dose, "increasing")
#' CochranArmitageTest(dose)
#' CochranArmitageTest(dose, "decreasing")
#' 
#' 
#' # not exactly the same as in package coin:
#' # independence_test(tumor ~ dose, data = lungtumor, teststat = "quad")
#' lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
#'                         tumor = c(rep(c(0, 1), c(38, 2)),
#'                                   rep(c(0, 1), c(43, 7)),
#'                                   rep(c(0, 1), c(33, 15))))
#' tab <- table(lungtumor$dose, lungtumor$tumor)
#' CochranArmitageTest(tab)
#' 
#' # but similar to
#' prop.trend.test(tab[,1], apply(tab,1, sum))
#' 
CochranArmitageTest <- function(x, alternative = c("two.sided","increasing","decreasing")) {

  # based on:
  # http://tolstoy.newcastle.edu.au/R/help/05/07/9442.html
  DNAME <- deparse(substitute(x))

  if (!(any(dim(x)==2)))
    stop("Cochran-Armitage test for trend must be used with rx2-table", call.=FALSE)

  if (dim(x)[2]!=2) x <- t(x)

  nidot <- apply(x, 1, sum)
  n <- sum(nidot)

  # Ri <- scores(x, 1, "table")
  Ri <- 1:dim(x)[1]
  Rbar <- sum(nidot*Ri)/n

  s2 <- sum(nidot*(Ri-Rbar)^2)
  pdot1 <- sum(x[,1])/n
  z <- sum(x[,1]*(Ri-Rbar))/sqrt(pdot1*(1-pdot1)*s2)
  STATISTIC <- z

  alternative <- match.arg(alternative)

  PVAL <- switch(alternative,
                 two.sided = 2*pnorm(abs(z), lower.tail=FALSE),
                 increasing = pnorm(z),
                 decreasing = pnorm(z, lower.tail=FALSE) )

  PARAMETER <- dim(x)[1]
  names(STATISTIC) <- "Z"
  names(PARAMETER) <- "dim"

  METHOD <- "Cochran-Armitage test for trend"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, alternative = alternative,
                 p.value = PVAL, method = METHOD, data.name = DNAME
  ), class = "htest")

}





#' Barnard's Unconditional Test
#' 
#' Barnard's unconditional test for superiority applied to \eqn{2 \times
#' 2}{2x2} contingency tables using Score or Wald statistics for the difference
#' between two binomial proportions.
#' 
#' If \code{x} is a matrix, it is taken as a two-dimensional contingency table,
#' and hence its entries should be nonnegative integers.  Otherwise, both
#' \code{x} and \code{y} must be vectors of the same length.  Incomplete cases
#' are removed, the vectors are coerced into factor objects, and the
#' contingency table is computed from these.
#' 
#' For a 2x2 contingency table, such as \eqn{X=[n_1,n_2;n_3,n_4]}, the
#' normalized difference in proportions between the two categories, given in
#' each column, can be written with pooled variance (Score statistic) as
#' \deqn{T(X)=\frac{\hat{p}_2-\hat{p}_1}{\sqrt{\hat{p}(1-\hat{p})(\frac{1}{c_1}+\frac{1}{c_2})}},}
#' where \eqn{\hat{p}=(n_1+n_3)/(n_1+n_2+n_3+n_4)},
#' \eqn{\hat{p}_2=n_2/(n_2+n_4)}, \eqn{\hat{p}_1=n_1/(n_1+n_3)},
#' \eqn{c_1=n_1+n_3} and \eqn{c_2=n_2+n_4}. Alternatively, with unpooled
#' variance (Wald statistic), the difference in proportions can we written as
#' \deqn{T(X)=\frac{\hat{p}_2-\hat{p}_1}{\sqrt{\frac{\hat{p}_1(1-\hat{p}_1)}{c_1}+\frac{\hat{p}_2(1-\hat{p}_2)}{c_2}}}.}
#' The probability of observing \eqn{X} is
#' \deqn{P(X)=\frac{c_1!c_2!}{n_1!n_2!n_3!n_4!}p^{n_1+n_2}(1-p)^{n_3+n_4},}
#' where \eqn{p} is the unknown nuisance parameter.
#' 
#' Barnard's test considers all tables with category sizes \eqn{c_1} and
#' \eqn{c_2} for a given \eqn{p}. The p-value is the sum of probabilities of
#' the tables having a score in the rejection region, e.g. having significantly
#' large difference in proportions for a two-sided test. The p-value of the
#' test is the maximum p-value calculated over all \eqn{p} between 0 and 1.
#' 
#' @param x a numeric vector or a two-dimensional contingency table in matrix
#' form. \code{x} and \code{y} can also both be factors.
#' @param y a factor object; ignored if \code{x} is a matrix.
#' @param dp The resolution of the search space for the nuisance parameter
#' @param pooled Z statistic with pooled (Score) or unpooled (Wald) variance
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}.  You can specify just the initial letter.
#' @return A list with class \code{"htest"} containing the following
#' components:
#' \item{p.value}{the p-value of the test.}
#' \item{estimate}{an
#' estimate of the nuisance parameter where the p-value is maximized.}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{the character string \code{"Barnards
#' Unconditional 2x2-test"}.}
#' \item{data.name}{a character string giving the
#' names of the data.}
#' \item{statistic.table }{The contingency tables
#' considered in the analysis represented by 'n1' and 'n2', their scores, and
#' whether they are included in the one-sided (1), two-sided (2) tests, or not
#' included at all (0)}
#' \item{nuisance.matrix }{Nuisance parameters, \eqn{p},
#' and the corresponding p-values for both one- and two-sided tests}
#' 
#' @author 
#' Kamil Erguler, <k.erguler@@cyi.ac.cy>, 
#' Peter Calhoun <calhoun.peter@@gmail.com>, 
#' Rodrigo Duprat, 
#' Andri Signorell <andri@@signorell.net> (interface)
#' 
#' @seealso 
#' \code{\link{fisher.test}}
#' 
#' @references Barnard, G.A. (1945) A new test for 2x2 tables. \emph{Nature},
#' 156:177.
#' 
#' Barnard, G.A. (1947) Significance tests for 2x2 tables. \emph{Biometrika},
#' 34:123-138.
#' 
#' Suissa, S. and Shuster, J. J. (1985), Exact Unconditional Sample Sizes for
#' the 2x2 Binomial Trial, \emph{Journal of the Royal Statistical Society},
#' Ser. A, 148, 317-327.
#' 
#' Cardillo G. (2009) MyBarnard: a very compact routine for Barnard's exact
#' test on 2x2 matrix.
#' \url{http://ch.mathworks.com/matlabcentral/fileexchange/25760-mybarnard}
#' 
#' Galili T. (2010)
#' \url{http://www.r-statistics.com/2010/02/barnards-exact-test-a-powerful-alternative-for-fishers-exact-test-implemented-in-r/}
#' 
#' Lin C.Y., Yang M.C. (2009) Improved p-value tests for comparing two
#' independent binomial proportions. \emph{Communications in
#' Statistics-Simulation and Computation}, 38(1):78-91.
#' 
#' Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez, L.
#' Rodriguez-Cardozo N.A. Ramos-Delgado and R. Garcia-Sanchez. (2004).
#' Barnardextest:Barnard's Exact Probability Test. A MATLAB file. 
#' [WWW document]. \url{http://www.mathworks.com/}
#' 
#' @keywords nonparametric htest
#'
#' @examples
#' 
#' tab <- as.table(matrix(c(8, 14, 1, 3), nrow=2,
#'                 dimnames=list(treat=c("I","II"), out=c("I","II"))))
#' BarnardTest(tab)
#' 
#' # Plotting the search for the nuisance parameter for a one-sided test
#' bt <- BarnardTest(tab)
#' plot(bt$nuisance.matrix[, 1:2],
#'      t="l", xlab="nuisance parameter", ylab="p-value")
#' 
#' # Plotting the tables included in the p-value
#' ttab <- as.table(matrix(c(40, 14, 10, 30), nrow=2,
#'                  dimnames=list(treat=c("I","II"), out=c("I","II"))))
#' 
#' bt <- BarnardTest(ttab)
#' bts <- bt$statistic.table
#' plot(bts[, 1], bts[, 2],
#'      col=hsv(bts[, 4] / 4, 1, 1),
#'      t="p", xlab="n1", ylab="n2")
#' 
#' # Plotting the difference between pooled and unpooled tests
#' bts <- BarnardTest(ttab, pooled=TRUE)$statistic.table
#' btw <- BarnardTest(ttab, pooled=FALSE)$statistic.table
#' plot(bts[, 1], bts[, 2],
#'      col=c("black", "white")[1 + as.numeric(bts[, 4]==btw[, 4])],
#'      t="p", xlab="n1", ylab="n2")
#' 
BarnardTest <- function(x, y = NULL, alternative = c("two.sided", "less", "greater"), dp = 0.001, pooled = TRUE ) {

  if (is.matrix(x)) {
    r <- nrow(x)
    if ((r < 2) || (ncol(x) != r))
      stop("'x' must be square with at least two rows and columns")
    if (any(x < 0) || anyNA(x))
      stop("all entries of 'x' must be nonnegative and finite")
    DNAME <- deparse(substitute(x))
  } else {
    if (is.null(y))
      stop("if 'x' is not a matrix, 'y' must be given")
    if (length(x) != length(y))
      stop("'x' and 'y' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    OK <- complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    r <- nlevels(x)
    if ((r < 2) || (nlevels(y) != r))
      stop("'x' and 'y' must have the same number of levels (minimum 2)")
    x <- table(x, y)
  }

  nr <- nrow(x)
  nc <- ncol(x)

  # if ((nr == 2) && (nc == 2)) {
  #   alternative <- char.expand(alternative, c("two.sided", "less", "greater"))
  #   if (length(alternative) > 1L || is.na(alternative))
  #     stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
  # }
  alternative <- match.arg(alternative)

  method <- c("Wald", "Score")[1 + pooled]
  METHOD <- gettextf("Barnards Unconditional 2x2-test", method)

  vec.size <- 1.0 + 1.0 / dp
  mat.size <- 4.0 * prod(rowSums(x) + 1) # (n1 + n3 + 1) * (n2 + n4 + 1)


  if (pooled)
    ret1 <- .C( "ScoreS",
                as.integer(x[1]), as.integer(x[2]), as.integer(x[3]), as.integer(x[4]),
                as.numeric(dp),
                mat.size = as.integer(0),
                statistic.table = as.double(vector("double", mat.size)),
                statistic = as.double(0.0))
  else
    ret1 <- .C( "WaldS",
                as.integer(x[1]), as.integer(x[2]), as.integer(x[3]), as.integer(x[4]),
                as.numeric(dp),
                mat.size = as.integer(0),
                statistic.table = as.double(vector("double", mat.size)),
                statistic = as.double(0.0))


  xr <- seq(1, ret1$mat.size, 4) + 2

  ret1$statistic.table[xr + 1][
    (ret1$statistic <= 0 & ret1$statistic.table[xr] <= ret1$statistic) |
      (ret1$statistic >= 0 & ret1$statistic.table[xr] >= ret1$statistic)] <- 1

  ret1$statistic.table[xr + 1][
    (ret1$statistic <= 0 & ret1$statistic.table[xr] >= -ret1$statistic) |
      (ret1$statistic >= 0 & ret1$statistic.table[xr] <= -ret1$statistic)] <- 2

  ret2 <- .C("Barnard",
             as.integer(x[1]), as.integer(x[2]), as.integer(x[3]), as.integer(x[4]),
             as.numeric(dp),
             as.integer(ret1$mat.size),
             nuisance.vector.x = as.double(vector("double",vec.size)),
             nuisance.vector.y0 = as.double(vector("double",vec.size)),
             nuisance.vector.y1 = as.double(vector("double",vec.size)),
             statistic.table = as.double(ret1$statistic.table))

  np0 <- which.max(ret2$nuisance.vector.y0)
  np1 <- which.max(ret2$nuisance.vector.y1)

  nuisance.matrix <- matrix(cbind(ret2$nuisance.vector.x, ret2$nuisance.vector.y0, ret2$nuisance.vector.y1), ncol=3)
  statistic.table <- matrix(ret1$statistic.table, ncol=4, byrow=TRUE, dimnames=list(c(), c("n1", "n2", "statistic", "include.in.p.value")))


  STATISTIC <- ret1$statistic
  if (alternative == "two.sided") {
    PVAL <- ret2$nuisance.vector.y1[np1]
    ESTIMATE <- c(`Nuisance parameter` = ret2$nuisance.vector.x[np1])
  } else {
    PVAL <- ret2$nuisance.vector.y0[np0]
    ESTIMATE <- c(`Nuisance parameter` = ret2$nuisance.vector.x[np0])
  }

  names(STATISTIC) <- gettextf("%s statistic", method)
  RVAL <- list(statistic = STATISTIC, alternative = alternative, estimate = ESTIMATE,
               p.value = PVAL, method = METHOD, data.name = DNAME,
               statistic.table = statistic.table, nuisance.matrix = nuisance.matrix)

  class(RVAL) <- "htest"
  return(RVAL)
}






#' Breusch-Godfrey Test
#' 
#' \code{BreuschGodfreyTest} performs the Breusch-Godfrey test for higher-order
#' serial correlation.
#' 
#' Under \eqn{H_0} the test statistic is asymptotically Chi-squared with
#' degrees of freedom as given in \code{parameter}.  If \code{type} is set to
#' \code{"F"} the function returns a finite sample version of the test
#' statistic, employing an \eqn{F} distribution with degrees of freedom as
#' given in \code{parameter}.
#' 
#' By default, the starting values for the lagged residuals in the auxiliary
#' regression are chosen to be 0 (as in Godfrey 1978) but could also be set to
#' \code{NA} to omit them.
#' 
#' \code{BreuschGodfreyTest} also returns the coefficients and estimated
#' covariance matrix from the auxiliary regression that includes the lagged
#' residuals.  Hence, \code{CoefTest} (package: RegClassTools) can be used to
#' inspect the results. (Note, however, that standard theory does not always
#' apply to the standard errors and t-statistics in this regression.)
#' 
#' @param formula a symbolic description for the model to be tested (or a
#' fitted \code{"lm"} object).
#' @param order integer. maximal order of serial correlation to be tested.
#' @param order.by Either a vector \code{z} or a formula with a single
#' explanatory variable like \code{~ z}. The observations in the model are
#' ordered by the size of \code{z}. If set to \code{NULL} (the default) the
#' observations are assumed to be ordered (e.g., a time series).
#' @param type the type of test statistic to be returned. Either \code{"Chisq"}
#' for the Chi-squared test statistic or \code{"F"} for the F test statistic.
#' @param data an optional data frame containing the variables in the model. By
#' default the variables are taken from the environment which
#' \code{BreuschGodfreyTest} is called from.
#' @param fill starting values for the lagged residuals in the auxiliary
#' regression. By default \code{0} but can also be set to \code{NA}.
#' @return A list with class \code{"BreuschGodfreyTest"} inheriting from
#' \code{"htest"} containing the following components:
#' \item{statistic}{the
#' value of the test statistic.}
#' \item{p.value}{the p-value of the test.}
#' \item{parameter}{degrees of freedom.}
#' \item{method}{a character string
#' indicating what type of test was performed.}
#' \item{data.name}{a character
#' string giving the name(s) of the data.}
#' \item{coefficients}{coefficient
#' estimates from the auxiliary regression.}
#' \item{vcov}{corresponding
#' covariance matrix estimate.}
#' @note This function was previously published as \code{bgtest} in the
#' \pkg{lmtest} package and has been integrated here without logical changes.
#' 
#' @author 
#' David Mitchell <david.mitchell@@dotars.gov.au>,
#' Achim Zeileis
#' 
#' @seealso 
#' \code{\link[DescTools]{DurbinWatsonTest}}
#' 
#' @references 
#' Johnston, J. (1984): \emph{Econometric Methods}, Third Edition, McGraw Hill 
#' Inc.
#' 
#' Godfrey, L.G. (1978): `Testing Against General Autoregressive and Moving
#' Average Error Models when the Regressors Include Lagged Dependent
#' Variables', \emph{Econometrica}, 46, 1293-1302.
#' 
#' Breusch, T.S. (1979): `Testing for Autocorrelation in Dynamic Linear
#' Models', \emph{Australian Economic Papers}, 17, 334-355.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Generate a stationary and an AR(1) series
#' x <- rep(c(1, -1), 50)
#' 
#' y1 <- 1 + x + rnorm(100)
#' 
#' ## Perform Breusch-Godfrey test for first-order serial correlation:
#' BreuschGodfreyTest(y1 ~ x)
#' 
#' ## or for fourth-order serial correlation
#' BreuschGodfreyTest(y1 ~ x, order = 4)
#' 
#' ## Compare with Durbin-Watson test results:
#' DurbinWatsonTest(y1 ~ x)
#' 
#' y2 <- filter(y1, 0.5, method = "recursive")
#' BreuschGodfreyTest(y2 ~ x)
#' 
BreuschGodfreyTest <- function(formula, order = 1, order.by = NULL, type = c("Chisq", "F"),
                               data = list(), fill = 0) {

  # from lmtest

  dname <- paste(deparse(substitute(formula)))

  if (!inherits(formula, "formula")) {
    X <- if (is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if (is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }

  if (!is.null(order.by))
  {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  n <- nrow(X)
  k <- ncol(X)
  order <- 1:order
  m <- length(order)
  resi <- lm.fit(X,y)$residuals

  Z <- sapply(order, function(x) c(rep(fill, length.out = x), resi[1:(n-x)]))
  if (any(na <- !complete.cases(Z))) {
    X <- X[!na, , drop = FALSE]
    Z <- Z[!na, , drop = FALSE]
    y <- y[!na]
    resi <- resi[!na]
    n <- nrow(X)
  }
  auxfit <- lm.fit(cbind(X,Z), resi)

  cf <- auxfit$coefficients
  vc <- chol2inv(auxfit$qr$qr) * sum(auxfit$residuals^2) / auxfit$df.residual
  names(cf) <- colnames(vc) <- rownames(vc) <- c(colnames(X), paste("lag(resid)", order, sep = "_"))

  switch(match.arg(type),

         "Chisq" = {
           bg <- n * sum(auxfit$fitted^2)/sum(resi^2)
           p.val <- pchisq(bg, m, lower.tail = FALSE)
           df <- m
           names(df) <- "df"
         },

         "F" = {
           uresi <- auxfit$residuals
           bg <- ((sum(resi^2) - sum(uresi^2))/m) / (sum(uresi^2) / (n-k-m))
           df <- c(m, n-k-m)
           names(df) <- c("df1", "df2")
           p.val <- pf(bg, df1 = df[1], df2 = df[2], lower.tail = FALSE)
         })

  names(bg) <- "LM test"
  RVAL <- list(statistic = bg, parameter = df,
               method = paste("Breusch-Godfrey test for serial correlation of order up to", max(order)),
               p.value = p.val,
               data.name = dname,
               coefficients = cf,
               vcov = vc)
  class(RVAL) <- c("BreuschGodfreyTest", "htest")
  return(RVAL)
}


# vcov.BreuschGodfreyTest <- function(object, ...) object$vcov
# df.residual.BreuschGodfreyTest <- function(object, ...) if (length(df <- object$parameter) > 1L) df[2] else NULL







#' Effect Size Calculations for ANOVAs
#' 
#' Calculates eta-squared, partial eta-squared and generalized eta-squared
#' 
#' Calculates the eta-squared, partial eta-squared, and generalized eta-squared
#' measures of effect size that are commonly used in analysis of variance. The
#' input \code{x} should be the analysis of variance object itself. For
#' between-subjects designs, generalized eta-squared equals partial
#' eta-squared. The reported generalized eta-squared for repeated-measures
#' designs assumes that all factors are manipulated, i.e., that there are no
#' measured factors like gender (see references).
#' 
#' For unbalanced designs, the default in \code{EtaSq} is to compute Type II
#' sums of squares (\code{type=2}), in keeping with the \code{Anova} function
#' in the \code{car} package. It is possible to revert to the Type I SS values
#' (\code{type=1}) to be consistent with \code{anova}, but this rarely tests
#' hypotheses of interest. Type III SS values (\code{type=3}) can also be
#' computed. \code{EtaSq.aovlist} requires \code{type=1}.
#' 
#' @aliases EtaSq EtaSq.lm EtaSq.aovlist aovlDetails aovlErrorTerms
#' @param x An analysis of variance (\code{aov}, \code{aovlist}) object.
#' @param type What type of sum of squares to calculate? \code{EtaSq.aovlist}
#' requires \code{type=1}.
#' @param anova Should the full ANOVA table be printed out in addition to the
#' effect sizes?
#' @return If \code{anova=FALSE}, the output for \code{EtaSq.lm} is an M x 2
#' matrix, for \code{EtaSq.aovlist} it is an M x 3 matrix. Each of the M rows
#' corresponds to one of the terms in the ANOVA (e.g., main effect 1, main
#' effect 2, interaction, etc), and each of the columns corresponds to a
#' different measure of effect size. Column 1 contains the eta-squared values,
#' and column 2 contains partial eta-squared values. Column 3 contains the
#' generalized eta-squared values. If \code{anova=TRUE}, the output contains
#' additional columns containing the sums of squares, mean squares, degrees of
#' freedom, F-statistics and p-values. For \code{EtaSq.aovlist}, additional
#' columns contain the error sum of squares and error degrees of freedom
#' corresponding to an effect term.
#' 
#' @author 
#' Daniel Navarro <daniel.navarro@@adelaide.edu.au>,
#' Daniel Wollschlaeger <dwoll@@psychologie.uni-kiel.de>
#' 
#' @seealso 
#' \code{\link{aov}}, 
#' \code{\link{anova}}, 
#' \code{\link[car]{Anova}}
#' 
#' @references 
#' Bakeman, R. (2005). Recommended effect size statistics for
#' repeated measures designs. Behavior Research Methods 37(3), 379-384.
#' 
#' Olejnik, S. and Algina, J. (2003). Generalized Eta and Omega Squared
#' Statistics: Measures of Effect Size for Some Common Research Designs.
#' Psychological Methods 8(4), 434-447.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' #### Example 1: one-way ANOVA ####
#' 
#' outcome <- c(1.4,2.1,3.0,2.1,3.2,4.7,3.5,4.5,5.4)    # data
#' treatment1 <- factor(c(1,1,1,2,2,2,3,3,3))           # grouping variable
#' anova1 <- aov(outcome ~ treatment1)                  # run the ANOVA
#' summary(anova1)                                      # print the ANOVA table
#' EtaSq(anova1)                                        # effect size
#' 
#' #### Example 2: two-way ANOVA ####
#' 
#' treatment2 <- factor(c(1,2,3,1,2,3,1,2,3))       # second grouping variable
#' anova2 <- aov(outcome ~ treatment1 + treatment2) # run the ANOVA
#' summary(anova2)                                  # print the ANOVA table
#' EtaSq(anova2)                                    # effect size
#' 
#' #### Example 3: two-way ANOVA unbalanced cell sizes ####
#' #### data from Maxwell & Delaney, 2004              ####
#' #### Designing experiments and analyzing data       ####
#' 
#' dfMD <- data.frame(
#'   IV1 = factor(rep(1:3, c(3 + 5 + 7, 5 + 6 + 4, 5 + 4 + 6))),
#'   IV2 = factor(rep(rep(1:3, 3), c(3, 5, 7, 5, 6, 4, 5, 4, 6))),
#'   DV = c(
#'     c(41, 43, 50),
#'     c(51, 43, 53, 54, 46),
#'     c(45, 55, 56, 60, 58, 62, 62),
#'     c(56, 47, 45, 46, 49),
#'     c(58, 54, 49, 61, 52, 62),
#'     c(59, 55, 68, 63),
#'     c(43, 56, 48, 46, 47),
#'     c(59, 46, 58, 54),
#'     c(55, 69, 63, 56, 62, 67)
#'   )
#' )
#' 
#' # use contr.sum for correct sum of squares type 3
#' dfMD$IV1s <- C(dfMD$IV1, "contr.sum")
#' dfMD$IV2s <- C(dfMD$IV2, "contr.sum")
#' dfMD$IV1t <- C(dfMD$IV1, "contr.treatment")
#' dfMD$IV2t <- C(dfMD$IV2, "contr.treatment")
#' 
#' EtaSq(aov(DV ~ IV1s*IV2s, data=dfMD), type=3)
#' EtaSq(aov(DV ~ IV1t*IV2t, data=dfMD), type=1)
#' 
#' #### Example 4: two-way split-plot ANOVA -> EtaSq.aovlist ####
#' 
#' DV_t1 <- round(rnorm(3*10, -0.5, 1), 2)
#' DV_t2 <- round(rnorm(3*10,  0,   1), 2)
#' DV_t3 <- round(rnorm(3*10,  0.5, 1), 2)
#' dfSPF <- data.frame(id=factor(rep(1:(3*10), times=3)),
#'                     IVbtw=factor(rep(LETTERS[1:3], times=3*10)),
#' 					IVwth=factor(rep(1:3, each=3*10)),
#' 					DV=c(DV_t1, DV_t2, DV_t3))
#' spf <- aov(DV ~ IVbtw*IVwth + Error(id/IVwth), data=dfSPF)
#' EtaSq(spf, type=1, anova=TRUE)
#' 
EtaSq <- function(x, type = 2, anova = FALSE) {
  UseMethod("EtaSq")
}

EtaSq.lm <- function(x, type = 2, anova = FALSE) {

  # file:    etaSquared.R
  # author:  Dan Navarro
  # contact: daniel.navarro@adelaide.edu.au
  # changed: 13 November 2013
  # modified by Daniel Wollschlaeger 17.9.2014

  # etaSquared() calculates eta-squared and partial eta-squared for linear models
  # (usually ANOVAs). It takes an lm object as input and computes the effect size
  # for all terms in the model. By default uses Type II sums of squares to calculate
  # the effect size, but Types I and III are also possible. By default the output
  # only displays the effect size, but if requested it will also print out the full
  # ANOVA table.

  if (!is(anova, "logical") | length(anova) != 1) {
    stop("\"anova\" must be a single logical value")
  }
  if (!is(type, "numeric") | length(type) != 1) {
    stop("type must be equal to 1, 2 or 3")
  }
  if (type == 1) {
    ss <- anova(x)[, "Sum Sq", drop = FALSE]
    ss.res <- ss[dim(ss)[1], ]
    ss.tot <- sum(ss)
    ss <- ss[-dim(ss)[1], , drop = FALSE]
    ss <- as.matrix(ss)
  }
  else {
    if (type == 2) {
      ss.tot <- sum((x$model[, 1] - mean(x$model[, 1]))^2)
      ss.res <- sum((x$residuals)^2)
      terms <- attr(x$terms, "factors")[-1, , drop = FALSE]
      l <- attr(x$terms, "term.labels")
      ss <- matrix(NA, length(l), 1)
      rownames(ss) <- l
      for (i in seq_along(ss)) {
        vars.this.term <- which(terms[, i] != 0)
        dependent.terms <- which(apply(terms[vars.this.term, , drop = FALSE], 2, prod) > 0)
        m0 <- lm(x$terms[-dependent.terms], x$model)
        if (length(dependent.terms) > 1) {
          m1 <- lm(x$terms[-setdiff(dependent.terms, i)], x$model)
          ss[i] <- anova(m0, m1)$`Sum of Sq`[2]
        }
        else {
          ss[i] <- anova(m0, x)$`Sum of Sq`[2]
        }
      }
    }
    else {
      if (type == 3) {
        ## check if model was fitted with sum-to-zero contrasts
        ## necessary for valid SS type 3 (e.g., contr.sum, contr.helmert)
        IVs <- names(attr(model.matrix(x), "contrasts"))
        ## only relevant for more than one factor
        ## (and for unbalanced cell sizes and interactions, not tested here)
        if (length(IVs) > 1) {
          isSumToZero <- function(IV) {
            ## check if factor has directly associated contrasts
            if (!is.null(attr(x$model[, IV], "contrasts"))) {
              cm <- contrasts(x$model[, IV])
              all(colSums(cm) == 0)
            } else {
              ## check attributes from model matrix
              attr(model.matrix(x), "contrasts")[[IV]] %in% c("contr.sum", "contr.helmert")
            }
          }

          valid <- vapply(IVs, isSumToZero, logical(1))

          if (!all(valid)) {
            warning(c(ifelse(sum(!valid) > 1, "Factors ", "Factor "),
                      paste(IVs[!valid], collapse=", "),
                      ifelse(sum(!valid) > 1, " are", " is"),
                      " not associated with sum-to-zero contrasts",
                      " necessary for valid SS type III",
                      " when cell sizes are unbalanced",
                      " and interactions are present.",
                      " Consider re-fitting the model after setting",
                      " options(contrasts=c(\"contr.sum\", \"contr.poly\"))"))
          }
        }

        mod <- drop1(x, scope = x$terms)
        ss <- mod[-1, "Sum of Sq", drop = FALSE]
        ss.res <- mod[1, "RSS"]
        ss.tot <- sum((x$model[, 1] - mean(x$model[, 1]))^2)
        ss <- as.matrix(ss)
      }
      else {
        stop("type must be equal to 1, 2 or 3")
      }
    }
  }
  if (anova == FALSE) {
    eta2 <- ss/ss.tot
    eta2p <- ss/(ss + ss.res)
    E <- cbind(eta2, eta2p)
    rownames(E) <- rownames(ss)
    colnames(E) <- c("eta.sq", "eta.sq.part")
  }
  else {
    ss <- rbind(ss, ss.res)
    eta2 <- ss/ss.tot
    eta2p <- ss/(ss + ss.res)
    k <- length(ss)
    eta2p[k] <- NA
    df <- anova(x)[, "Df"]
    ms <- ss/df
    Fval <- ms/ms[k]
    p <- 1 - pf(Fval, df, rep.int(df[k], k))
    E <- cbind(eta2, eta2p, ss, df, ms, Fval, p)
    E[k, 6:7] <- NA
    colnames(E) <- c("eta.sq", "eta.sq.part", "SS", "df", "MS", "F", "p")
    rownames(E) <- rownames(ss)
    rownames(E)[k] <- "Residuals"
  }
  return(E)
}


EtaSq.aovlist <-  function(x, type = 2, anova = FALSE) {

  # author:  Daniel Wollschlaeger
  # contact: contact@dwoll.de
  # changed: 13 October 2014

  # EtaSq.aovlist() calculates partial eta-squared and generalized eta-squared
  # for aovlists

  if (!is(anova, "logical") | length(anova) != 1) {
    stop("\"anova\" must be a single logical value")
  }
  if (!is(type, "numeric") | length(type) != 1) {
    stop("type must be equal to 1, 2 or 3")
  }

  ## alternative: check design has balanced cell sizes
  if (type != 1) {
    stop("type must be equal to 1")
  }

  details <- aovlDetails(x)
  ss      <- details$Sum.Sq             # effect SS
  ss.res  <- sum(aovlErrorTerms(x)$SS)  # total error SS
  ss.tot  <- sum(ss) + sum(ss.res)

  # eta squared
  eta2 <- ss / ss.tot

  # partial eta squared
  # cf. Bakeman, R. (2005) Behavior Research Methods. 37(3), 379-384. Tables 1, 2
  eta2p <- ss / (ss + details$SSE)

  # generalized eta squared
  # if all factors are manipulated
  # cf. Bakeman, R. (2005) Behavior Research Methods. 37(3), 379-384. Tables 1, 2
  geta2 <- ss / (ss + sum(ss.res))

  if (anova == FALSE) {
    E <- cbind(eta2, eta2p, geta2)
    rownames(E) <- details$tt
    colnames(E) <- c("eta.sq", "eta.sq.part", "eta.sq.gen")
  } else {
    E <- data.frame(eta2=eta2,
                    eta2p=eta2p,
                    geta2=geta2,
                    ss=ss,
                    df=details$Df,
                    ms=details$Mean.Sq,
                    sse=details$SSE,
                    dfe=details$dfE,
                    Fval=details$F.value,
                    p=details$Pr..F.)
    colnames(E) <- c("eta.sq", "eta.sq.part", "eta.sq.gen", "SS", "df", "MS", "SSE", "dfE", "F", "p")
    rownames(E) <- details$tt
  }
  return(E)
}

# author:  Daniel Wollschlaeger
aovlDetails <- function(aovObj) {
  aovSum  <- summary(aovObj)
  etNames <- names(aovSum)  # error terms

  getOneRes <- function(tt, tab) {  # tab=anova table, tt = tested term
    ttIdx <- which(DescTools::StrTrim(rownames(tab)) == tt)
    list(df=tab[ttIdx,       "Df"],
         SS=tab[ttIdx,       "Sum Sq"],
         MS=tab[ttIdx,       "Mean Sq"],
         dfE=tab["Residuals", "Df"],
         SSE=tab["Residuals", "Sum Sq"],
         MSE=tab["Residuals", "Mean Sq"],
         F=tab[ttIdx, "F value"],
         p=tab[ttIdx, "Pr(>F)"])
  }

  getTermRes <- function(et) { # et = error term
    tab <- aovSum[[et]][[1]]
    at  <- DescTools::StrTrim(rownames(tab)) # all terms
    tt  <- at[-which(at == "Residuals")]     # tested terms only

    if (length(tt) > 0)
    {
      # error terms
      etRes <- list(df=tab["Residuals", "Df"],
                    SS=tab["Residuals", "Sum Sq"],
                    MS=tab["Residuals", "Mean Sq"])
      ttRes <- lapply(tt, getOneRes, tab=tab)
      ttRes <- setNames(ttRes, tt)
      ttIdx <- which(DescTools::StrTrim(rownames(tab)) %in% tt)
      return(data.frame(tt=tt, et=et,
                        tab[ttIdx, , drop=FALSE],
                        dfE=etRes$df, SSE=etRes$SS, MSE=etRes$MS,
                        stringsAsFactors=FALSE))
    } else {
      emptyDf <- data.frame(matrix(ncol=10, nrow=0))
      return(setNames(emptyDf, c("tt", "et", "Df", "Sum.Sq", "Mean.Sq", "F.value",
                                 "Pr..F.", "dfE", "SSE", "MSE")))
    }
  }

  detailsL  <- setNames(lapply(etNames, getTermRes), etNames)
  detailsDf <- do.call("rbind", detailsL)
  rownames(detailsDf) <- NULL
  return(detailsDf)
}

aovlErrorTerms <- function(aovObj) {
  aovSum  <- summary(aovObj)
  etNames <- names(aovSum)
  getSS <- function(x) {
    aovSum[[x]][[1]]["Residuals", "Sum Sq"]
  }

  getMS <- function(x) {
    aovSum[[x]][[1]]["Residuals", "Mean Sq"]
  }

  getDF <- function(x) {
    aovSum[[x]][[1]]["Residuals", "Df"]
  }

  SS <- vapply(etNames, getSS, numeric(1))
  MS <- vapply(etNames, getMS, numeric(1))
  DF <- vapply(etNames, getDF, numeric(1))
  return(list(SS=SS, MS=MS, DF=DF))
}





# TODO: examples are missing

#' Greenhouse-Geisser And Huynh-Feldt Epsilons 
#' 
#' Calculate Greenhouse-Geisser and Huynh-Feldt epsilons. 
#' 
#' @param S p x p covariance matrix
#' @param p dimension of observation vectors
#' @param g number of groups
#' @param n number of subjects
#' @return a numeric value
#' 
#' @author 
#' Hans Rudolf Roth <hroth@@retired.ethz.ch> 
#' 
#' @seealso 
#' \code{\link{aov}}
#' 
#' @references 
#' Vonesh, E.F., Chinchilli, V.M. (1997) \emph{Linear and Nonlinear
#' Models for the Analysis of Repeated Measurements} Marcel Dekker, New York,
#' p.84-86
#' 
#' Crowder, M.J., Hand, D.J. (1990) \emph{Analysis of Repeated Measures}.
#' Chapman & Hall, London, p.54-55 
#' @keywords models regression
#'
#' @examples
#' 
#' ## find!

Eps <- function(S, p, g, n) {

  ## Purpose: calculates the Greenhouse-Geisser and Huynh-Feldt epsilons
  ## ------------------------------------------------------------------~
  ## Arguments: S pxp covariance matrix
  ##            p dimension of observation vectors
  ##            g number of groups
  ##            n number of subjects

  ## Lit:    E.F. Vonesh + V.M. Chinchilli (1997), p.84-86
  ##         M.J. Crowder and D.J. Hand (1990), p.54-55

  ## Author: H.-R. Roth
  ## Date:   23.07.2002
  ## ------------------------------------------------------------------~

  # U is a matrix of (p-1) orthonormal contrasts
  U <- t(cbind(diag(p-1),0) - outer(1:(p-1), 1:p, "<") / ((p-1):1))
  a <- 1/sqrt(colSums(U^2))
  U <- U%*%diag(a)
  V <- t(U) %*% S %*% U
  e <- (sum(diag(V)))^2/sum(diag(V%*%V))/(p-1)

  GGepsilon <- e
  HFepsilon <- min(1, (n*(p-1)*e - 2) / ((p-1)* (n-g-(p-1)*e) ))
  t.output <- c(GGepsilon, HFepsilon)
  names(t.output)  <- c("G-G-epsilon", "H-F-epsilon")
  t.output

}





#' Power Calculations for ChiSquared Tests
#' 
#' Compute power of test or determine parameters to obtain target power (same
#' as \code{\link{power.anova.test}}).
#' 
#' Exactly one of the parameters \code{w}, \code{n}, \code{power} or
#' \code{sig.level} must be passed as NULL, and this parameter is determined
#' from the others. Note that the last one has non-NULL default, so \code{NULL}
#' must be explicitly passed, if you want to compute it.
#' 
#' @param n total number of observations.
#' @param w effect size.
#' @param df degree of freedom (depends on the chosen test.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @return Object of class "power.htest", a list of the arguments (including
#' the computed one) augmented with 'method' and 'note' elements.
#' @note \code{\link{uniroot}} is used to solve power equation for unknowns, so
#' you may see errors from it, notably about inability to bracket the root when
#' invalid arguments are given.
#' 
#' @author 
#' Stephane Champely <champely@@univ-lyon1.fr> \cr
#' but this is a mere copy of Peter Dalgaard's work on power.t.test
#' 
#' @seealso
#' \code{\link{power.t.test}}
#' 
#' @references 
#' Cohen, J. (1988) \emph{Statistical power analysis for the
#' behavioral sciences (2nd ed.)} Hillsdale, NJ: Lawrence Erlbaum.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Exercise 7.1 P. 249 from Cohen (1988) 
#' power.chisq.test(w=0.289, df=(4-1), n=100, sig.level=0.05)
#' 
#' ## Exercise 7.3 p. 251
#' power.chisq.test(w=0.346, df=(2-1)*(3-1), n=140, sig.level=0.01)
#' 
#' ## Exercise 7.8 p. 270
#' power.chisq.test(w=0.1, df=(5-1)*(6-1), power=0.80, sig.level=0.05)
#' 
power.chisq.test <- function(n = NULL, w = NULL, df = NULL, sig.level = 0.05, power = NULL) {

  if (sum(sapply(list(w, n, df, power, sig.level), is.null)) != 1)
    stop("exactly one of w, n, df, power or sig.level must be NULL")
  if (!is.null(w) && w < 0)
    stop("w must be positive")
  if (!is.null(n) && n < 1)
    stop("number of observations must be at least 1")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop(sQuote("sig.level"), " must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop(sQuote("power"), " must be numeric in [0, 1]")
  p.body <- quote({
    k <- qchisq(sig.level, df = df, lower = FALSE)
    pchisq(k, df = df, ncp = n * w^2, lower = FALSE)
  })
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(w))
    w <- uniroot(function(w) eval(p.body) - power, c(1e-10, 1e+05))$root
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(1 + 1e-10, 1e+05))$root
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) -
                           power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  METHOD <- "Chi squared power calculation"
  NOTE <- "n is the number of observations"
  structure(list(w = w, n = n, df = df, sig.level = sig.level,
                 power = power, method = METHOD, note = NOTE), class = "power.htest")
}






#' Pairwise Contrasts
#' 
#' Generate all pairwise contrasts for using in a post-hoc test, e.g.
#' ScheffeTest.
#' 
#' 
#' @param levs the levels to be used
#' 
#' @return A matrix with all possible pairwise contrasts, that can be built
#' with the given levels. 
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso 
#' \code{\link{ScheffeTest}} 
#' 
#' @keywords htest
#'
#' @examples
#' 
#' Contrasts(LETTERS[1:5])
#' 
#' #   B-A C-A D-A E-A C-B D-B E-B D-C E-C E-D
#' # A  -1  -1  -1  -1   0   0   0   0   0   0
#' # B   1   0   0   0  -1  -1  -1   0   0   0
#' # C   0   1   0   0   1   0   0  -1  -1   0
#' # D   0   0   1   0   0   1   0   1   0  -1
#' # E   0   0   0   1   0   0   1   0   1   1
#' 
#' 
Contrasts <- function(levs) {
  k = length(levs)
  M = data.frame(levs = levs)
  for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
      con = rep(0, k)
      con[i] = -1
      con[j] = 1
      nm = paste(levs[j], levs[i], sep = "-")
      M[[nm]] = con
    }
  }
  row.names(M) = levs

  return(M[-1])

}




#' Scheffe Test for Pairwise and Otherwise Comparisons
#' 
#' Scheffe's method applies to the set of estimates of all possible contrasts
#' among the factor level means, not just the pairwise differences considered
#' by Tukey's method.
#' 
#' 
#' @aliases ScheffeTest ScheffeTest.default ScheffeTest.aov
#' @param x either a fitted model object, usually an \code{\link{aov}} fit,
#' when g is left to \code{NULL} or a response variable to be evalutated by g
#' (which mustn't be \code{NULL} then). 
#' @param g the grouping variable.
#' @param which character vector listing terms in the fitted model for which
#' the intervals should be calculated. Defaults to all the terms.
#' @param contrasts a \eqn{r \times c}{r x c} matrix containing the contrasts
#' to be computed, while \code{r} is the number of factor levels and \code{c}
#' the number of contrasts. Each column must contain a full contrast ("sum")
#' adding up to 0. Note that the argument \code{which} must be defined, when
#' non default contrasts are used.  Default value of \code{contrasts} is
#' \code{NULL}. In this case all pairwise contrasts will be reported.
#' @param conf.level numeric value between zero and one giving the confidence
#' level to use.  If this is set to NA, just a matrix with the p-values will be
#' returned. 
#' @param \dots further arguments, currently not used.
#' 
#' @return A list of classes \code{c("PostHocTest")}, with one component for
#' each term requested in \code{which}. Each component is a matrix with columns
#' \code{diff} giving the difference in the observed means, \code{lwr.ci}
#' giving the lower end point of the interval, \code{upr.ci} giving the upper
#' end point and \code{pval} giving the p-value after adjustment for the
#' multiple comparisons.
#' 
#' There are print and plot methods for class \code{"PostHocTest"}. The plot
#' method does not accept \code{xlab}, \code{ylab} or \code{main} arguments and
#' creates its own values for each plot.
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso 
#' \code{\link{pairwise.t.test}},
#' \code{\link{TukeyHSD}} 
#' 
#' @references 
#' Robert O. Kuehl, Steel R. (2000) \emph{Design of experiments}. Duxbury
#' 
#' Steel R.G.D., Torrie J.H., Dickey, D.A. (1997) \emph{Principles and
#' Procedures of Statistics, A Biometrical Approach}. McGraw-Hill
#' 
#' @keywords htest
#'
#' @examples
#' 
#' fm1 <- aov(breaks ~ wool + tension, data = warpbreaks)
#' 
#' ScheffeTest(x=fm1)
#' ScheffeTest(x=fm1, which="tension")
#' 
#' TukeyHSD(fm1)
#' 
#' # some special contrasts
#' y <- c(7,33,26,27,21,6,14,19,6,11,11,18,14,18,19,14,9,12,6,
#'        24,7,10,1,10,42,25,8,28,30,22,17,32,28,6,1,15,9,15,
#'        2,37,13,18,23,1,3,4,6,2)
#' group <- factor(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,
#'        3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6))
#' 
#' r.aov <- aov(y ~ group)
#' 
#' ScheffeTest(r.aov, contrasts=matrix( c(1,-0.5,-0.5,0,0,0,
#'                                        0,0,0,1,-0.5,-0.5), ncol=2))
#' 
#' # just p-values:
#' ScheffeTest(r.aov, conf.level=NA)
#' 
ScheffeTest <- function(x, ...) {
  UseMethod("ScheffeTest")
}

ScheffeTest.default <- function(x, g = NULL, which = NULL, contrasts = NULL, conf.level = 0.95, ...) {
  ScheffeTest(x=aov(x~g), which=which, contrasts=contrasts, conf.level=conf.level, ...)
}



ScheffeTest.aov <- function(x, which=NULL, contrasts = NULL, conf.level=0.95, ...) {

  mm <- model.tables(x, "means")
  if (is.null(mm$n))
    stop("no factors in the fitted model")
  tabs <- mm$tables[-1L]

  if (is.null(which)) which <- seq_along(tabs)

  tabs <- tabs[which]
  nn <- mm$n[names(tabs)]
  nn_na <- is.na(nn)
  if (all(nn_na))
    stop("'which' specified no factors")
  if (any(nn_na)) {
    warning("'which' specified some non-factors which will be dropped")
    tabs <- tabs[!nn_na]
    nn <- nn[!nn_na]
  }
  out <- setNames(vector("list", length(tabs)), names(tabs))
  MSE <- sum(x$residuals^2)/x$df.residual

  autoContr <- is.null(contrasts)
  if (!is.null(contrasts)) {
    contrasts <- data.frame(contrasts)
  }

  # nm <- "tension"
  for (nm in names(tabs)) {
    tab <- tabs[[nm]]
    means <- as.vector(tab)

    nms <- if (length(d <- dim(tab)) > 1L) {
      dn <- dimnames(tab)
      apply(do.call("expand.grid", dn), 1L, paste, collapse = ":")
    } else names(tab)

    n <- nn[[nm]]
    if (length(n) < length(means))
      n <- rep.int(n, length(means))

    if (autoContr) contrasts <- Contrasts(nms)

    psi <- apply(contrasts * means, 2, sum)
    sscoeff <- apply(contrasts * contrasts / n, 2, sum)
    mspsi <- (psi * psi) / sscoeff

    # Korrektur von Daniel Wollschlaeger 9.9.2014:
    #     psi <- contrasts %*% means
    #     sscoeff <- contrasts * contrasts %*% (1/n)

    dferr <- x$df.residual
    dfgrp <- length(x$residuals) - dferr - 1

    pval <- pf(psi^2/(MSE*sscoeff*dfgrp),
               df1=dfgrp, df2=dferr, lower.tail=FALSE)

    critvalue <- dfgrp * qf(1-conf.level, dfgrp, dferr, lower.tail=FALSE)

    lwr <- psi - sqrt(critvalue) * sqrt(MSE * sscoeff)
    upr <- psi + sqrt(critvalue) * sqrt(MSE * sscoeff)

    out[[nm]] <- cbind(diff=psi, lwr, upr, pval)
    colnames(out[[nm]]) <- c("diff","lwr.ci","upr.ci","pval")

    if (!autoContr) {
      # define contrasts rownames
      rownames(out[[nm]]) <-  apply(contrasts, 2, function(x)
        gettextf("%s-%s", paste(nms[x>0], collapse=","),
                 paste(nms[x<0], collapse=",")) )
      if (is.na(conf.level)) out[[nm]] <- out[[nm]][,-c(2:3)]
    }

    if (autoContr & is.na(conf.level)) {
      out[[nm]] <- matrix(NA, nrow=length(means), ncol=length(means))
      out[[nm]][lower.tri(out[[nm]], diag = FALSE)] <- pval
      dimnames(out[[nm]]) <- list(nms, nms)
      out[[nm]] <- out[[nm]][-1, -ncol(out[[nm]])]
    }

  }

  class(out) <- c("PostHocTest")
  attr(out, "orig.call") <- x$call
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- FALSE
  attr(out, "method") <- "Scheffe Test"
  attr(out, "method.str") <- gettextf("\n  Posthoc multiple comparisons of means : %s \n", attr(out, "method"))


  return(out)

}




#' Post-Hoc Tests
#' 
#' A convenience wrapper for computing post-hoc test after having calculated an
#' ANOVA.
#' 
#' The function is designed to consolidate a couple of post-hoc tests with the
#' same interface for input and output.
#' 
#' Choosing Tests\cr
#' Different Post Hoc tests use different methods to control
#' FW and PE. Some tests are very conservative. Conservative tests go to great
#' lengths to prevent the user from committing a Type I error.  They use more
#' stringent criterion for determining significance. Many of these tests become
#' more and more stringent as the number of groups increases (directly limiting
#' the FW and PE error rate). Although these tests buy you protection against
#' Type I error, it comes at a cost. As the tests become more stringent, you
#' loose Power (1-B).  More Liberal tests, buy you Power but the cost is an
#' increased chance of Type I error.  There is no set rule for determining
#' which test to use, but different researchers have offered some guidelines
#' for choosing. Mostly it is an issue of pragmatics and whether the number of
#' comparisons exceeds K-1.
#' 
#' Fisher's LSD\cr
#' The Fisher LSD (Least Significant Different) sets Alpha
#' Level per comparison. Alpha = .05 for every comparison. df = df error (i.e.
#' df within). This test is the most liberal of all Post Hoc tests. The
#' critical t for significance is unaffected by the number of groups. This test
#' is appropriate when you have 3 means to compare. In general the alpha is
#' held at .05 because of the criterion that you can't look at LSD's unless the
#' Anova is significant. This test is generally not considered appropriate if
#' you have more than 3 means unless there is reason to believe that there is
#' no more than one true Null Hypothesis hidden in the means.
#' 
#' Dunn's (Bonferroni)\cr
#' Dunn's t-test is sometimes referred to as the
#' Bonferroni t because it used the Bonferroni PE correction procedure in
#' determining the critical value for significance. In general, this test
#' should be used when the number of comparisons you are making exceeds the
#' number of degrees of freedom you have between groups (e.g. K-1). This test
#' sets alpha per experiment; Alpha = (.05)/c for every comparison. df = df
#' error (c = number of comparisons (K(K-1))/2) This test is extremely
#' conservative and rapidly reduces power as the number of comparisons being
#' made increase.
#' 
#' Newman-Keuls\cr
#' Newman-Keuls is a step down procedure that is not as
#' conservative as Dunn's t test. First, the means of the groups are ordered
#' (ascending or descending) and then the largest and smallest means are tested
#' for significant differences. If those means are different, then test
#' smallest with next largest, until you reach a test that is not significant.
#' Once you reach that point then you can only test differences between means
#' that exceed the difference between the means that were found to be
#' non-significant. Newman-Keuls is perhaps one of the most common Post Hoc
#' test, but it is a rather controversial test. The major problem with this
#' test is that when there is more than one true Null Hypothesis in a set of
#' means it will overestimate they FW error rate. In general we would use this
#' when the number of comparisons we are making is larger than K-1 and we don't
#' want to be as conservative as the Dunn's test is.
#' 
#' Tukey's HSD\cr
#' Tukey HSD (Honestly Significant Difference) is essentially
#' like the Newman-Keul, but the tests between each mean are compared to the
#' critical value that is set for the test of the means that are furthest apart
#' (rmax e.g. if there are 5 means we use the critical value determined for the
#' test of X1 and X5). This Method corrects for the problem found in the
#' Newman-Keuls where the FW is inflated when there is more than one True Null
#' Hypothesis in a set of means. It buys protection against Type I error, but
#' again at the cost of Power. It tends to be the most common test and
#' preferred test because it is very conservative with respect to Type I error
#' when the Null hypothesis is true. In general, HSD is preferred when you will
#' make all the possible comparisons between a large set of means (Six or more
#' means).
#' 
#' Scheffe\cr
#' The Scheffe Test is designed to protect against a Type I error
#' when all possible complex and simple comparisons are made. That is we are
#' not just looking the possible combinations of comparisons between pairs of
#' means. We are also looking at the possible combinations of comparisons
#' between groups of means. Thus Scheffe is the most conservative of all tests.
#' Because this test does give us the capacity to look at complex comparisons,
#' it essentially uses the same statistic as the Linear Contrasts tests.
#' However, Scheffe uses a different critical value (or at least it makes an
#' adjustment to the critical value of F). This test has less power than the
#' HSD when you are making Pairwise (simple) comparisons, but it has more power
#' than HSD when you are making Complex comparisons. In general, only use this
#' when you want to make many Post Hoc complex comparisons (e.g. more than
#' K-1).
#' 
#' Tables\cr
#' For tables pairwise chi-square test can be performed, either
#' without correction or with correction for multiple testing following the
#' logic in \code{\link{p.adjust}}.
#' 
#' @aliases PostHocTest PostHocTest.aov PostHocTest.table PostHocTest.matrix
#' print.PostHocTest plot.PostHocTest
#' @param x an aov object. 
#' @param method one of \code{"hsd"}, \code{"bonf"}, \code{"lsd"},
#' \code{"scheffe"}, \code{"newmankeuls"}, defining the method for the pairwise
#' comparisons.\cr
#' For the post hoc test of tables the methods of
#' \code{\link{p.adjust}} can be supplied. See the detail there. 
#' @param which a character vector listing terms in the fitted model for which
#' the intervals should be calculated. Defaults to all the terms. 
#' @param conf.level a numeric value between zero and one giving the
#' family-wise confidence level to use.  If this is set to NA, just a matrix
#' with the p-values will be returned. 
#' @param ordered a logical value indicating if the levels of the factor should
#' be ordered according to increasing average in the sample before taking
#' differences. If ordered is \code{TRUE} then the calculated differences in
#' the means will all be positive. The significant differences will be those
#' for which the lower end point is positive. \cr
#' This argument will be ignored
#' if method is not either \code{hsd} or \code{newmankeuls}.
#' @param digits controls the number of fixed digits to print.
#' @param \dots further arguments, not used so far. 
#' @return an object of type "PostHocTest", which will either be \cr
#' A) a list of data.frames containing the mean difference, lower ci, upper ci 
#' and the
#' p-value, if a conf.level was defined (something else than NA) or \cr 
#' B) a list of matrices with the p-values, if conf.level has been set to NA.
#'  
#' @author 
#' Andri Signorell <andri@@signorell.net> 
#' 
#' @seealso 
#' \code{\link{TukeyHSD}}, 
#' \code{\link{aov}},
#' \code{\link{pairwise.t.test}}, 
#' \code{\link{ScheffeTest}}
#' 
#' @keywords htest
#'
#' @examples
#' 
#' PostHocTest(aov(breaks ~ tension, data = warpbreaks), method = "lsd")
#' PostHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd")
#' PostHocTest(aov(breaks ~ tension, data = warpbreaks), method = "scheffe")
#' 
#' r.aov <- aov(breaks ~ tension, data = warpbreaks)
#' 
#' # compare p-values:
#' round(cbind(
#'     lsd= PostHocTest(r.aov, method="lsd")$tension[,"pval"]
#'   , bonf=PostHocTest(r.aov, method="bonf")$tension[,"pval"]
#' ), 4)
#' 
#' # only p-values by setting conf.level to NA
#' PostHocTest(aov(breaks ~ tension, data = warpbreaks), method = "hsd",
#'             conf.level=NA)
#' 
PostHocTest <- function(x, ...) {
  UseMethod("PostHocTest")
}


PostHocTest.aov <- function(x, which = NULL,
                             method=c("hsd","bonferroni","lsd","scheffe","newmankeuls","duncan"),
                             conf.level = 0.95, ordered = FALSE, ...) {

  method <- match.arg(method)

  if (method=="scheffe") {
    out <- ScheffeTest(x=x, which=which, conf.level=conf.level, ...)

  } else {

    mm <- model.tables(x, "means")
    if (is.null(mm$n))
      stop("no factors in the fitted model")
    tabs <- mm$tables[-1L]

    if (is.null(which)) which <- seq_along(tabs)
    tabs <- tabs[which]

    nn <- mm$n[names(tabs)]
    nn_na <- is.na(nn)
    if (all(nn_na))
      stop("'which' specified no factors")
    if (any(nn_na)) {
      warning("'which' specified some non-factors which will be dropped")
      tabs <- tabs[!nn_na]
      nn <- nn[!nn_na]
    }
    out <- setNames(vector("list", length(tabs)), names(tabs))
    MSE <- sum(x$residuals^2)/x$df.residual
    for (nm in names(tabs)) {
      tab <- tabs[[nm]]
      means <- as.vector(tab)
      nms <- if (length(d <- dim(tab)) > 1L) {
        dn <- dimnames(tab)
        apply(do.call("expand.grid", dn), 1L, paste, collapse = ":")
      }
      else names(tab)
      n <- nn[[nm]]
      if (length(n) < length(means))
        n <- rep.int(n, length(means))

      # this will be ignored for bonferroni, lsd
      if (method %in% c("hsd", "newmankeuls", "duncan") & as.logical(ordered)) {
        ord <- order(means)
        means <- means[ord]
        n <- n[ord]
        if (!is.null(nms))
          nms <- nms[ord]
      }

      center <- outer(means, means, "-")
      keep <- lower.tri(center)
      center <- center[keep]

      switch(method
             ,"bonferroni" = {
               width <-  qt(1 - (1 - conf.level)/(length(means) * (length(means) - 1)), x$df.residual) *
                 sqrt(MSE * outer(1/n, 1/n, "+"))[keep]
               est <- center/sqrt(MSE * outer(1/n, 1/n, "+")[keep])

               pvals <- pmin(2 * pt(abs(est), df = x$df.residual, lower.tail = FALSE)
                             * ((length(means)^2 - length(means))/2), 1)
               method.str <- "Bonferroni"

             }
             ,"lsd" = {
               width <-  qt(1 - (1 - conf.level)/2, x$df.residual) *
                 sqrt(MSE * outer(1/n, 1/n, "+"))[keep]
               est <- center/sqrt(MSE * outer(1/n, 1/n, "+")[keep])
               pvals <- 2 * pt(abs(est), df = x$df.residual, lower.tail = FALSE)
               method.str <- "Fisher LSD"
             }
             ,"hsd" = {
               width <- qtukey(conf.level, length(means), x$df.residual) *
                 sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep]
               est <- center/(sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep])
               pvals <- ptukey(abs(est), length(means), x$df.residual,
                               lower.tail = FALSE)
               method.str <- "Tukey HSD"

             }
             ,"newmankeuls" ={
               nmean <- (abs(outer(rank(means), rank(means), "-")) + 1)[keep]

               width <- qtukey(conf.level, nmean, x$df.residual) *
                 sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep]

               est <- center/(sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep])

               pvals <- ptukey(abs(est), nmean, x$df.residual, lower.tail = FALSE)
               method.str <- "Newman-Keuls"

             }
             ,"duncan" = {
               # same as newmankeuls, but with bonferroni corrected alpha
               nmean <- (abs(outer(rank(means), rank(means), "-")) + 1)[keep]

               width <- qtukey(conf.level^(nmean-1), nmean, x$df.residual) *
                 sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep]

               est <- center/(sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep])
               pvals <- 1-(1-ptukey(abs(est), nmean, x$df.residual,
                                    lower.tail = FALSE))^(1/(nmean - 1))

               method.str <- "Duncan's new multiple range test"

             }
             ,"dunnett" = {
               method.str <- "Dunnett"
             }
             ,"scottknott" = {
               method.str <- "Scott Knott"
             }
             ,"waller" = {
               method.str <- "Waller"
             }
             ,"gabriel" = {
               method.str <- "Gabriel"
             }
      )

      if (!is.na(conf.level)) {
        dnames <- list(NULL, c("diff", "lwr.ci", "upr.ci", "pval"))
        if (!is.null(nms))
          dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
        out[[nm]] <- array(c(center, center - width,
                             center + width, pvals), c(length(width), 4L), dnames)
      } else {
        out[[nm]] <- matrix(NA, nrow=length(means), ncol=length(means))
        out[[nm]][lower.tri(out[[nm]], diag = FALSE)] <- pvals
        dimnames(out[[nm]]) <- list(nms, nms)
        out[[nm]] <- out[[nm]][-1, -ncol(out[[nm]])]

      }
    }

    class(out) <- c("PostHocTest")
    attr(out, "orig.call") <- x$call
    attr(out, "conf.level") <- conf.level
    attr(out, "ordered") <- ordered
    attr(out, "method") <- method.str
    attr(out, "method.str") <- gettextf("\n  Posthoc multiple comparisons of means : %s \n", attr(out, "method"))

  }

  return(out)

}


PostHocTest.matrix <- function(x, method = c("none","fdr","BH","BY","bonferroni","holm","hochberg","hommel"),
                               conf.level = 0.95, ...) {

  # http://support.sas.com/resources/papers/proceedings14/1544-2014.pdf

  # no conf.level supported so far
  conf.level  <- NA

  method <- match.arg(method)

  #  out <- setNames(vector("list", length(tabs)), names(tabs))

  pvals <- DescTools::PairApply(t(as.matrix(x)), FUN = function(y1, y2) chisq.test(cbind(y1,y2))$p.value, symmetric=TRUE)
  pvals[upper.tri(pvals, diag=TRUE)] <- NA

  if (method != "none")
    pvals[] <- p.adjust(pvals, method=method)

  #  pvals[] <- format.pval(pvals, digits = 2, na.form = "-")
  pvals <- pvals[-1, -ncol(pvals)]
  out <- list()
  out[[deparse(substitute(x))]] <- pvals

  class(out) <- c("PostHocTest")
  attr(out, "orig.call") <- "table"
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- FALSE
  attr(out, "method") <- method
  attr(out, "method.str") <- gettextf("\n  Posthoc multiple comparisons on chi-square test : %s \n", attr(out, "method"))

  return(out)

}


PostHocTest.table <- function(x, method = c("none","fdr","BH","BY","bonferroni","holm","hochberg","hommel"),
                              conf.level = 0.95, ...) {
  class(x) <- "matrix"
  PostHocTest(x, method=method, conf.level=conf.level, ...)
}




print.PostHocTest <- function(x, digits = getOption("digits", 3), ...) {

  cat(attr(x, "method.str"))
  if (!is.na(attr(x, "conf.level")))
    cat("    ", format(100 * attr(x, "conf.level"), 2), "% family-wise confidence level\n",
        sep = "")
  if (attr(x, "ordered"))
    cat("    factor levels have been ordered\n")
  if (!is.language(attr(x, "orig.call")) && !is.null(attr(x, "orig.call")))
    cat("\nFit: ", deparse(attr(x, "orig.call"), 500L), "\n\n", sep = "")
  else
    cat("\n")
  xx <- unclass(x)

  attr(xx, "orig.call") <- attr(xx, "conf.level") <-
    attr(xx, "ordered") <-  attr(xx, "method.str") <-  attr(xx, "method") <- NULL

  xx["data.name"] <- NULL

  if (!is.na(attr(x, "conf.level"))) {
    xx <- lapply(xx, as.data.frame)
    for(nm in names(xx)) {
      xx[[nm]]$" " <- Format(xx[[nm]]$"pval", fmt="*")
      xx[[nm]]$"pval" <- format.pval(xx[[nm]]$"pval", digits=2, nsmall=4)
    }

    print.default(xx, digits=digits, ...)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  } else {
    for(nm in names(xx)) {
      xx[[nm]][] <- format.pval(xx[[nm]], 2, na.form = "-")
    }
    #     attributes(pp) <- attributes(x$p.value)
    print(xx, digits=digits, quote = FALSE, ...)
  }
  cat("\n")

  invisible(x)
}



plot.PostHocTest <- function(x, ...) {
  # original:   stats:::plot.TukeyHSD(x, ...)

  # don't need that here..
  x$data.name <- NULL

  for (i in seq_along(x)) {
    xi <- x[[i]][, -4L, drop = FALSE]
    yvals <- nrow(xi):1L
    dev.hold()
    on.exit(dev.flush())
    plot(c(xi[, "lwr.ci"], xi[, "upr.ci"]), rep.int(yvals, 2L),
         type = "n", axes = FALSE, xlab = "", ylab = "", main = NULL,
         ...)
    axis(1, ...)
    axis(2, at = nrow(xi):1, labels = dimnames(xi)[[1L]],
         srt = 0, ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0.5, ...)
    segments(xi[, "lwr.ci"], yvals, xi[, "upr.ci"], yvals, ...)
    segments(as.vector(xi), rep.int(yvals - 0.1, 3L), as.vector(xi),
             rep.int(yvals + 0.1, 3L), ...)
    title(main = paste0(format(100 * attr(x, "conf.level"),
                               digits = 2L), "% family-wise confidence level\n"),
          xlab = paste("Differences in mean levels of", names(x)[i]))
    box()
    dev.flush()
    on.exit()
  }

}





#' Dunn's Test of Multiple Comparisons
#' 
#' Performs Dunn's test of multiple comparisons using rank sums.
#' 
#' \code{DunnTest} performs the post hoc pairwise multiple comparisons
#' procedure appropriate to follow the rejection of a Kruskal-Wallis test. The
#' Kruskal-Wallis test, being a non-parametric analog of the one-way ANOVA, is
#' an omnibus test of the null hypothesis that none of k groups stochastically
#' dominate one another. Dunn's test is constructed in part by summing jointly
#' ranked data. The rank sum test, itself a non-parametric analog of the
#' unpaired t-test, is possibly intuitive, but inappropriate as a post hoc
#' pairwise test, because (1) it fails to retain the dependent ranking that
#' produced the Kruskal-Wallis test statistic, and (2) it does not incorporate
#' the pooled variance estimate implied by the null hypothesis of the
#' Kruskal-Wallis test.
#' 
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors.  In this case, \code{g} is
#' ignored, and one can simply use \code{DunnTest(x)} to perform the test.  If
#' the samples are not yet contained in a list, use \code{DunnTest(list(x,
#' ...))}.
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @aliases DunnTest DunnTest.default DunnTest.formula print.DunnTest
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param method the method for adjusting p-values for multiple comparisons.
#' The function is calling \code{\link{p.adjust}} and this parameter is
#' directly passed through.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param out.list logical, indicating if the results should be printed in list
#' mode or as a square matrix. Default is list (TRUE).
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param digits controls the number of fixed digits to print.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \code{"DunnTest"} containing the following
#' components:
#' \item{res}{an array containing the mean rank differencens and
#' the according p-values}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>, the interface is based on R-Core code
#' 
#' @seealso 
#' \code{\link{kruskal.test}}, 
#' \code{\link{wilcox.test}},
#' \code{\link{p.adjust}}
#' 
#' @references 
#' Dunn, O. J. (1961) Multiple comparisons among means
#' \emph{Journal of the American Statistical Association}, 56(293):52-64.
#' 
#' Dunn, O. J. (1964) Multiple comparisons using rank sums
#' \emph{Technometrics}, 6(3):241-252.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' DunnTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' # do the kruskal.test first
#' kruskal.test(x, g)
#' 
#' # ...and the pairwise test afterwards
#' DunnTest(x, g)
#' 
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' DunnTest(Ozone ~ Month, data = airquality)
#' 
DunnTest <- function(x, ...) {
  UseMethod("DunnTest")
}


DunnTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) ||
      (length(formula) != 3L) || 
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  if (length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  y <- DoCall("DunnTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME
  y
}




DunnTest.default <- function(x, g, method = c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"),
                              alternative = c("two.sided", "less", "greater"),
                              out.list = TRUE, ...) {

  alternative <- match.arg(alternative)

  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }
  N <- length(x)
  if (N < 2)
    stop("not enough observations")

  method <- match.arg(method)

  nms <- levels(g)

  n <- tapply(g, g, length)
  rnk <- rank(x)
  mrnk <- tapply(rnk, g, mean)

  tau <- table(rnk[AllDuplicated(rnk)])
  tiesadj <- sum(tau^3 - tau) / (12*(N-1))
  mrnkdiff <- outer(mrnk, mrnk, "-")

  z <- mrnkdiff / sqrt( ((N*(N+1)/12) - tiesadj) * outer(1/n, 1/n, "+"))

  # extension for alternative in v. 0.99.16:
  if (alternative == "less") {
    pvals <- pnorm(abs(z))
  }
  else if (alternative == "greater") {
    pvals <- pnorm(abs(z), lower.tail=FALSE)
  }
  else {
    pvals <- 2 * pnorm(abs(z), lower.tail=FALSE)
  }


  keep <- lower.tri(pvals)
  pvals <- pvals[keep]
  m <- sum(keep)

  out <- list()

  pvals <- p.adjust(pvals, method=method)
  method.str <- method

  if (out.list) {
    dnames <- list(NULL, c("mean rank diff", "pval"))
    if (!is.null(nms))
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    out[[1]] <- array(c(mrnkdiff[keep], pvals), c(length(mrnkdiff[keep]), 2L), dnames)

  } else {
    out[[1]] <- matrix(NA, nrow=length(nms), ncol=length(nms))
    out[[1]][lower.tri(out[[1]], diag = FALSE)] <- pvals
    dimnames(out[[1]]) <- list(nms, nms)
    out[[1]] <- out[[1]][-1, -ncol(out[[1]])]

  }

  class(out) <- c("DunnTest")
  attr(out, "main") <- gettextf("Dunn's test of multiple comparisons using rank sums : %s ", method.str)
  attr(out, "method") <- method.str
  attr(out, "out.list") <- out.list

  return(out)

}




print.DunnTest <- function(x, digits = getOption("digits", 3), ...) {

  cat("\n", attr(x, "main"), "\n\n")
  xx <- unclass(x)

  if (attr(x, "out.list")==TRUE) {
    xx <- data.frame(x[1])
    xx$" " <- Format(xx$"pval", fmt="*")
    xx$"pval" <- format.pval(xx$"pval", digits=2, nsmall=4)

    print.data.frame(xx, digits=digits, ...)
    cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  } else {
    xx[[1]][] <- format.pval(xx[[1]], 2, na.form = "-")
    #     attributes(pp) <- attributes(x$p.value)
    print(xx[[1]], digits=digits, quote = FALSE, ...)
  }
  cat("\n")

  invisible(x)
}





#' Conover's Test of Multiple Comparisons
#' 
#' Perform Conover's test of multiple comparisons using rank sums as post hoc
#' test following a significant \code{\link{kruskal.test}}.
#' 
#' \code{ConoverTest} performs the post hoc pairwise multiple comparisons
#' procedure appropriate to follow the rejection of a Kruskal-Wallis test.
#' Conover's test is more powerful than Dunn's post hoc multiple comparisons
#' test (\code{\link{DunnTest}}). The interpretation of stochastic dominance
#' requires an assumption that the CDF of one group does not cross the CDF of
#' the other. \cr
#' ConoverTest makes m = k(k-1)/2 multiple pairwise comparisons
#' based on the Conover-Iman t-test-statistic for the rank-sum differences:
#' \deqn{\left | \bar{R}_{i}-\bar{R}_{j} \right | > t_{1-\alpha/2, n-k} \cdot
#' \sqrt{ s^2 \cdot \left [ \frac{n-1-\hat{H}^*}{n-k} \right ] \cdot \left [
#' \frac{1}{n_i} + \frac{1}{n_j} \right ] } } with the (tie corrected)
#' statistic of the Kruskal Wallis test \deqn{\hat{H}^* = \frac{\frac{12}{n
#' \cdot (n+1)} \cdot \sum_{i=1}^{k}\frac{R_{i}^2}{n_i} - 3\cdot(n+1) }
#' {1-\frac{\sum_{i=1}^{r} \left ( t_i^3-t_i \right )}{n^3-n}} } and the
#' \eqn{s^2} being \deqn{s^2 = \frac{1}{n-1} \cdot \left [ \sum{R_i^2} - n
#' \cdot \frac{(n+1)^2}{4} \right ]}
#' 
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors.  In this case, \code{g} is
#' ignored, and one can simply use \code{ConoverTest(x)} to perform the test.
#' If the samples are not yet contained in a list, use
#' \code{ConoverTest(list(x, ...))}.
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @aliases ConoverTest ConoverTest.default ConoverTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param method the method for adjusting p-values for multiple comparisons.
#' The function is calling \code{\link{p.adjust}} and this parameter is
#' directly passed through.
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of \code{"two.sided"} (default), \code{"greater"} or
#' \code{"less"}. You can specify just the initial letter.
#' @param out.list logical, indicating if the results should be printed in list
#' mode or as a square matrix. Default is list (TRUE).
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class \code{"DunnTest"} containing the following
#' components:
#' \item{res}{an array containing the mean rank differencens and
#' the according p-values}
#' 
#' @author
#' Andri Signorell <andri@@signorell.net>, the interface is based on R-Core code
#' 
#' @seealso 
#' \code{\link{DunnTest}},
#' \code{\link{NemenyiTest}},
#' \code{\link{kruskal.test}}, 
#' \code{\link{wilcox.test}},
#' \code{\link{p.adjust}}
#' 
#' @references 
#' Conover W. J., Iman R. L. (1979) On multiple-comparisons
#' procedures, \emph{Tech. Rep.} LA-7677-MS, Los Alamos Scientific Laboratory.
#' 
#' Conover, W. J. (1999) Practical Nonparametric Statistics \emph{Wiley},
#' Hoboken, NJ. 3rd edition.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' ConoverTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' # do the kruskal.test first
#' kruskal.test(x, g)
#' 
#' # ...and the pairwise test afterwards
#' ConoverTest(x, g)
#' 
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' ConoverTest(Ozone ~ Month, data = airquality)
#' 
ConoverTest <- function(x, ...) {
  UseMethod("ConoverTest")
}

ConoverTest.default <- function(x, g,
  method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  alternative = c("two.sided", "less", "greater"), out.list = TRUE, ...) {

  alternative <- match.arg(alternative)

  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  } else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }

  N <- length(x)
  if (N < 2)
    stop("not enough observations")
  method <- match.arg(method)
  nms <- levels(g)
  n <- tapply(g, g, length)
  rnk <- rank(x)
  mrnk <- tapply(rnk, g, mean)
  tau <- table(rnk[AllDuplicated(rnk)])
  tiesadj <- 1-sum(tau^3 - tau)/(N^3 - N)
  mrnkdiff <- outer(mrnk, mrnk, "-")

  # Kruskal-Wallis H statistic
  H <- (12 / (N * (N + 1))) * sum(tapply(rnk, g, sum)^2 / n) - 3 * (N + 1)
  if (tiesadj == 1) {
    s2 <- N * (N + 1) / 12
  } else {
    s2 <-   ( 1 / (N - 1)) * (sum(rnk^2) - (N * (((N + 1)^2) / 4)))
  }

  tval <- mrnkdiff/sqrt(s2 * ((N - 1 - H/tiesadj) / (N - k)) * outer(1/n, 1/n, "+"))

  if (alternative == "less") {
    pvals <- pt(abs(tval), df=N - k)

  } else if (alternative == "greater") {
    pvals <- pt(abs(tval), df=N - k, lower.tail = FALSE)

  } else {
    pvals <- 2 * pt(abs(tval), df=N - k, lower.tail = FALSE)

  }

  keep <- lower.tri(pvals)
  pvals <- pvals[keep]
  m <- sum(keep)
  out <- list()
  pvals <- p.adjust(pvals, method = method)
  method.str <- method
  if (out.list) {
    dnames <- list(NULL, c("mean rank diff", "pval"))
    if (!is.null(nms))
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    out[[1]] <- array(c(mrnkdiff[keep], pvals),
      c(length(mrnkdiff[keep]), 2L), dnames)
  } else {
    out[[1]] <- matrix(NA, nrow = length(nms), ncol = length(nms))
    out[[1]][lower.tri(out[[1]], diag = FALSE)] <- pvals
    dimnames(out[[1]]) <- list(nms, nms)
    out[[1]] <- out[[1]][-1, -ncol(out[[1]])]
  }

  class(out) <- c("ConoverTest", "DunnTest")
  attr(out, "main") <- gettextf("Conover's test of multiple comparisons : %s ",
                                method.str)
  attr(out, "method") <- method.str
  attr(out, "out.list") <- out.list

  return(out)

}




ConoverTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) ||
      (length(formula) != 3L) ||
      (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  if (length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")

  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  y <- DoCall("ConoverTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME

  y

}




# Test  NemenyiTest
#
# d.frm <- data.frame(x=c(28,30,33,35,38,41, 36,39,40,43,45,50, 44,45,47,49,53,54),
#                     g=c(rep(LETTERS[1:3], each=6)))
# NemenyiTest(x~g, d.frm)
#
# library(coin)
# library(multcomp)
# nem <- oneway_test(x ~ g, data=d.frm,
#                    ytrafo = function(data) trafo(data, numeric_trafo=rank),
#                    xtrafo = function(data) trafo(data, factor_trafo = function(x)
#                      model.matrix(~x - 1) %*% t(contrMat(table(x), "Tukey"))),
#                    teststat="max")
# nem
# pvalue(nem, method="single-step")



# TODO!! Tell when to use this test. References needed!
# 
#' Nemenyi Test 
#' 
#' Performs Nemenyi's test of multiple comparisons.
#' 
#' Nemenyi proposed a
#' test based on rank sums and the application of the family-wise error method
#' to control Type I error inflation, if multiple comparisons are done. The
#' Tukey and Kramer approach uses mean rank sums and can be employed for
#' equally as well as unequally sized samples without ties.
#' 
#' @aliases NemenyiTest NemenyiTest.default NemenyiTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param dist the distribution used for the test. Can be \code{tukey}
#' (default) or \code{chisq}.
#' @param out.list logical, defining if the output should be organized in
#' listform.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{htest}, containing the following components:
#' \item{statistic}{ Nemenyi test}
#' \item{p.value}{ the p-value for the test}
#' \item{null.value}{is the value of the median specified by the null
#' hypothesis. This equals the input argument \code{mu}. }
#' \item{alternative}{a
#' character string describing the alternative hypothesis.}
#' \item{method}{ the
#' type of test applied}
#' \item{data.name}{a character string giving the names
#' of the data.}
#' 
#' @author 
#' Andri Signorell <andri@@signorell.net>
#' 
#' @seealso 
#' \code{\link{DunnTest}}, 
#' \code{\link{ConoverTest}} 
#' 
#' @references Nemenyi, P. B. (1963) \emph{Distribution-Free Multiple
#' Comparisons} New York, State University of New York, Downstate Medical
#' Center
#' 
#' Hollander, M., Wolfe, D.A. (1999) \emph{Nonparametric Statistical Methods}
#' New York, Wiley, pp. 787
#' 
#' Friedman, M. (1937) The use of ranks to avoid the assumption of normality
#' implicit in the analysis of variance \emph{Journal of the American
#' Statistical Association}, 32:675-701
#' 
#' Friedman, M. (1940) A comparison of alternative tests of significance for
#' the problem of m rankings \emph{Annals of Mathematical Statistics}, 11:86-92
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' 
#' NemenyiTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' NemenyiTest(x, g)
#' 
#' ## Formula interface.
#' boxplot(Ozone ~ Month, data = airquality)
#' NemenyiTest(Ozone ~ Month, data = airquality)
#' 
#' # Hedderich & Sachs, 2012, p. 555
#' d.frm <- data.frame(x=c(28,30,33,35,38,41, 36,39,40,43,45,50, 44,45,47,49,53,54),
#'                     g=c(rep(LETTERS[1:3], each=6)))
#' 
#' NemenyiTest(x~g, d.frm)
#' 
NemenyiTest <- function(x, ...) {
  UseMethod("NemenyiTest")
}

NemenyiTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  if (length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  y <- DoCall("NemenyiTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME
  y
}




NemenyiTest.default <- function(x, g,
                                 dist = c("tukey", "chisq"), out.list = TRUE, ...) {

  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  }
  else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }
  N <- length(x)
  if (N < 2)
    stop("not enough observations")

  dist <- match.arg(dist, c("tukey", "chisq"))

  nms <- levels(g)

  n <- tapply(g, g, length)
  rnk <- rank(x)
  mrnk <- tapply(rnk, g, mean)

  tau <- table(rnk[AllDuplicated(rnk)])
  tiesadj <- min(1, 1 - sum(tau^3 - tau) / (N^3 - N))
  mrnkdiff <- outer(mrnk, mrnk, "-")

  if (dist == "chisq") {
    chi <- mrnkdiff^2 / ((N*(N+1)/12) * outer(1/n, 1/n, "+"))
    pvals <- pchisq(tiesadj * chi, df=k-1, lower.tail=FALSE)
  } else {
    z <- abs(mrnkdiff) / sqrt( (N*(N+1)/12) * outer(1/n, 1/n, "+"))
    pvals <- ptukey(z * sqrt(2), nmeans=k, df=Inf, lower.tail=FALSE)
  }


  keep <- lower.tri(pvals)
  pvals <- pvals[keep]
  m <- sum(keep)

  out <- list()

  # no p.adjustment in this test
  # pvals <- p.adjust(pvals, method=method)
  method.str <- "none" #method

  if (out.list) {
    dnames <- list(NULL, c("mean rank diff", "pval"))
    if (!is.null(nms))
      dnames[[1L]] <- outer(nms, nms, paste, sep = "-")[keep]
    out[[1]] <- array(c(mrnkdiff[keep], pvals), c(length(mrnkdiff[keep]), 2L), dnames)

  } else {
    out[[1]] <- matrix(NA, nrow=length(nms), ncol=length(nms))
    out[[1]][lower.tri(out[[1]], diag = FALSE)] <- pvals
    dimnames(out[[1]]) <- list(nms, nms)
    out[[1]] <- out[[1]][-1, -ncol(out[[1]])]

  }

  class(out) <- c("DunnTest")
  attr(out, "main") <- gettextf("Nemenyi's test of multiple comparisons for independent samples (%s) ", dist)
  attr(out, "method") <- method.str
  attr(out, "out.list") <- out.list

  return(out)

}






#' Dunnett's Test for Comparing Several Treatments With a Control
#' 
#' Performs Dunnett's test for comparing several treatments with a control.
#' 
#' \code{DunnettTest} does the post hoc pairwise multiple comparisons
#' procedure.
#' 
#' If \code{x} is a list, its elements are taken as the samples to be compared,
#' and hence have to be numeric data vectors.  In this case, \code{g} is
#' ignored, and one can simply use \code{DunnettTest(x)} to perform the test.
#' If the samples are not yet contained in a list, use
#' \code{DunnettTest(list(x, ...))}.
#' 
#' Otherwise, \code{x} must be a numeric data vector, and \code{g} must be a
#' vector or factor object of the same length as \code{x} giving the group for
#' the corresponding elements of \code{x}.
#' 
#' @aliases DunnettTest DunnettTest.default DunnettTest.formula
#' @param x a numeric vector of data values, or a list of numeric data vectors.
#' @param g a vector or factor object giving the group for the corresponding
#' elements of \code{x}.  Ignored if \code{x} is a list.
#' @param control the level of the control group against which the others
#' should be tested. If there are multiple levels the calculation will be
#' performed for every one.
#' @param conf.level confidence level of the interval.
#' @param formula a formula of the form \code{lhs ~ rhs} where \code{lhs} gives
#' the data values and \code{rhs} the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list of class \code{c("PostHocTest")}, containing one matrix named
#' after the control with columns \code{diff} giving the difference in the
#' observed means, \code{lwr.ci} giving the lower end point of the interval,
#' \code{upr.ci} giving the upper end point and \code{pval} giving the p-value
#' after adjustment for the multiple comparisons.
#' 
#' There are print and plot methods for class \code{"PostHocTest"}. The plot
#' method does not accept \code{xlab}, \code{ylab} or \code{main} arguments and
#' creates its own values for each plot.
#' 
#' @author Andri Signorell <andri@@signorell.net>, the interface is based on
#' R-Core code
#' 
#' @seealso 
#' \code{\link{PostHocTest}}
#' 
#' @references 
#' Dunnett C. W. (1955) A multiple comparison procedure for
#' comparing several treatments with a control, \emph{Journal of the American
#' Statistical Association}, 50:1096-1121.
#' 
#' @keywords htest
#'
#' @examples
#' 
#' ## Hollander & Wolfe (1973), 116.
#' ## Mucociliary efficiency from the rate of removal of dust in normal
#' ##  subjects, subjects with obstructive airway disease, and subjects
#' ##  with asbestosis.
#' x <- c(2.9, 3.0, 2.5, 2.6, 3.2) # normal subjects
#' y <- c(3.8, 2.7, 4.0, 2.4)      # with obstructive airway disease
#' z <- c(2.8, 3.4, 3.7, 2.2, 2.0) # with asbestosis
#' 
#' DunnettTest(list(x, y, z))
#' 
#' ## Equivalently,
#' x <- c(x, y, z)
#' g <- factor(rep(1:3, c(5, 4, 5)),
#'             labels = c("Normal subjects",
#'                        "Subjects with obstructive airway disease",
#'                        "Subjects with asbestosis"))
#' 
#' DunnettTest(x, g)
#' 
#' ## Formula interface
#' boxplot(Ozone ~ Month, data = airquality)
#' DunnettTest(Ozone ~ Month, data = airquality)
#' 
#' DunnettTest(Ozone ~ Month, data = airquality, control="8", conf.level=0.9)
#' 
DunnettTest <- function(x, ...) {
  UseMethod("DunnettTest")
}


DunnettTest.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  if (length(mf) > 2L)
    stop("'formula' should be of the form response ~ group")
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  y <- DoCall("DunnettTest", c(as.list(mf), list(...)))
  y$data.name <- DNAME
  y
}



DunnettTest.default <- function(x, g, control = NULL
                                 , conf.level = 0.95, ...) {

  if (is.list(x)) {
    if (length(x) < 2L)
      stop("'x' must be a list with at least 2 elements")
    DNAME <- deparse(substitute(x))
    x <- lapply(x, function(u) u <- u[complete.cases(u)])
    k <- length(x)
    l <- sapply(x, "length")
    if (any(l == 0))
      stop("all groups must contain data")
    g <- factor(rep(1:k, l))
    x <- unlist(x)
  } else {
    if (length(x) != length(g))
      stop("'x' and 'g' must have the same length")
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    OK <- complete.cases(x, g)
    x <- x[OK]
    g <- g[OK]
    if (!all(is.finite(g)))
      stop("all group levels must be finite")
    g <- factor(g)
    k <- nlevels(g)
    if (k < 2)
      stop("all observations are in the same group")
  }
  N <- length(x)
  if (N < 2)
    stop("not enough observations")

  # just organisational stuff so far, got a fine x and g now

  if (is.null(control)) control <- levels(g)[1]

  ctrls <- control
  out <- list()

  for(ii in seq_along(ctrls)) {

    control <- ctrls[ii]

    ni <- tapply(x, g, length)

    means <- tapply(x, g, mean)
    meandiffs <- means[names(means) != control] - means[control]

    fittedn <- ni[names(ni) != control]
    controln <- ni[control]

    s <- sqrt( sum(tapply(x, g, function(x) sum((x - mean(x))^2) )) /
                 (N - k))

    Dj <- meandiffs / (s * sqrt((1/fittedn) + (1/controln)))
    Rij <- sqrt(fittedn/(fittedn + controln))

    R <- outer(Rij, Rij, "*")
    diag(R) <- 1

    set.seed(5)  # for getting consistent results every run
    qvt <- mvtnorm::qmvt((1 - (1 - conf.level)/2), df = N - k, sigma = R, tail = "lower.tail")$quantile

    lower <- meandiffs - s * sqrt((1/fittedn) + (1/controln)) * qvt
    upper <- meandiffs + s * sqrt((1/fittedn) + (1/controln)) * qvt

    pval <- c()
    for (i in 1:(k-1)) {
      pval[i] <- 1 - mvtnorm::pmvt(-abs(Dj[i]), abs(Dj[i]), corr=R, delta=rep(0, k-1), df=N - k)[1]
    }

    out[[ii]] <- cbind(diff=meandiffs, lower, upper, pval)
    dimnames(out[[ii]]) <- list(paste(names(meandiffs), control, sep="-"), c("diff", "lwr.ci", "upr.ci","pval"))
  }

  names(out) <- ctrls

  class(out) <- c("PostHocTest")
  #  attr(out, "orig.call") <- NA
  attr(out, "conf.level") <- conf.level
  attr(out, "ordered") <- FALSE
  attr(out, "method") <- ""
  attr(out, "method.str") <- gettextf("\n  Dunnett's test for comparing several treatments with a control : %s \n", attr(out, "method"))

  return(out)

}






#' Hotelling's T2 Test
#' 
#' 
#' Hotelling's T2 test is the multivariate generlisation of the Student's t
#' test. A one-sample Hotelling's T2 test can be used to test if a set of
#' vectors of data (which should be a sample of a single statistical
#' population) has a mean equal to a hypothetical mean. A two-sample
#' Hotelling's T2 test may be used to test for significant differences between
#' the mean vectors (multivariate means) of two multivariate data sets are
#' different.
#' 
#' 
#' The classical test for testing the location of a multivariate population or
#' for testing the mean difference for two multivariate populations. When
#' \code{test = "f"} the F-distribution is used for the test statistic and it
#' is assumed that the data are normally distributed. If the chisquare
#' approximation is used, the normal assumption can be relaxed to existence of
#' second moments.  In the two sample case both populations are assumed to have
#' the same covariance matrix.
#' 
#' The formula interface is only applicable for the 2-sample tests.
#' 
#' @aliases HotellingsT2Test HotellingsT2Test.default HotellingsT2Test.formula
#' @param x a numeric data frame or matrix.
#' @param y an optional numeric data frame or matrix for the two sample test.
#' If \code{NULL} a one sample test is performed.
#' @param mu a vector indicating the hypothesized value of the mean (or
#' difference in means if a two sample test is performed). \code{NULL}
#' represents origin or no difference between the groups.
#' @param test if \code{"f"}, the decision is based on the F-distribution, if
#' \code{"chi"} a chi-squared approximation is used.
#' @param formula a formula of the form \code{x ~ g} where \code{x} is a
#' numeric matrix giving the data values and \code{g} a factor with two levels
#' giving the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see
#' \code{\link{model.frame}}) containing the variables in the formula
#' \code{formula}.  By default the variables are taken from
#' \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#' used.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. Defaults to \code{getOption("na.action")}.
#' @param \dots further arguments to be passed to or from methods.
#' @return A list with class 'htest' containing the following components:
#' \item{statistic }{the value of the T2-statistic. (That is the scaled value
#' of the statistic that has an F distribution or a chisquare distribution
#' depending on the value of \code{test}).}
#' \item{parameter}{the degrees of
#' freedom for the T2-statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{null.value}{the specified hypothesized value of the mean or mean
#' difference depending on whether it was a one-sample test or a two-sample
#' test.}
#' \item{alternative}{a character string with the value 'two.sided'.}
#' \item{method}{a character string indicating what type of test was
#' performed.}
#' \item{data.name}{a character string giving the name of the data
#' (and grouping vector).}
#' 
#' @author 
#' Klaus Nordhausen, <klaus.nordhausen@@uta.fi>
#' 
#' @references
#' Nordhausen K., Sirkia S., Oja H. and Tyler D. E. (2012)
#' \emph{ICSNP: Tools for Multivariate Nonparametrics}. R package version
#' 1.0-9.\cr
#' \url{https://cran.r-project.org/package=ICSNP}
#' 
#' Anderson, T.W. (2003), \emph{An introduction to multivariate analysis}, New
#' Jersey: Wiley.
#' 
#' @keywords htest multivariate
#'
#' @examples
#' 
#' math.teach <- data.frame(
#'   teacher = factor(rep(1:2, c(3, 6))),
#'   satis   = c(1, 3, 2, 4, 6, 6, 5, 5, 4),
#'   know    = c(3, 7, 2, 6, 8, 8, 10, 10, 6))
#' 
#' with(math.teach,
#'   HotellingsT2Test(cbind(satis, know) ~ teacher))
#' 
HotellingsT2Test <- function(x,...) {
  UseMethod("HotellingsT2Test")
}

HotellingsT2Test.default <- function(x, y=NULL, mu=NULL, test="f",...) {


  `HotellingsT.internal`  <-  function(x, y=NULL, mu, test) {
    n <- dim(x)[1]
    p <- dim(x)[2]

    if (is.null(y))     #one sample case
    {
      test.statistic <- n*as.numeric(t(colMeans(x)-mu)%*%solve(cov(x))%*%(colMeans(x)-mu))*switch(test,f=(n-p)/(p*(n-1)),chi=1)
      df.1 <- p
      df.2 <- switch(test,f=n-p,chi=NA)
      p.value <- 1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
      return(list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2))
    }

    # else two sample case
    n1 <- n
    n2 <- dim(y)[1]
    xmeans <- colMeans(x)
    ymeans <- colMeans(y)
    x.diff <- sweep(x,2,xmeans)
    y.diff <- sweep(y,2,ymeans)
    S.pooled <- 1/(n1+n2-2)*(t(x.diff)%*%x.diff+t(y.diff)%*%y.diff)
    test.statistic <- n1*n2/(n1+n2)*t(xmeans-ymeans-mu)%*%solve(S.pooled)%*%(xmeans-ymeans-mu)*switch(test,f=(n1+n2-p-1)/(p*(n1+n2-2)),chi=1)
    df.1 <- p
    df.2 <- switch(test,f=n1+n2-p-1,chi=NA)
    p.value <- 1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
    list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2)
  }


  if (is.null(y)) {
    DNAME <- deparse(substitute(x))
  } else {
    DNAME=paste(deparse(substitute(x)),"and",deparse(substitute(y)))
  }

  xok <- complete.cases(x)
  x <- x[xok,]
  if (!all(sapply(x, is.numeric))) stop("'x' must be numeric")
  x <- as.matrix(x)

  p <- dim(x)[2]

  if (!is.null(y)) {
    yok <- complete.cases(y)
    y <- y[yok,]

    if (!all(sapply(y, is.numeric))) stop("'y' must be numeric")
    if (p!=dim(y)[2]) stop("'x' and 'y' must have the same number of columns")
    y <- as.matrix(y)
  }

  if (is.null(mu)) mu <- rep(0,p)
  else if (length(mu)!=p) stop("length of 'mu' must equal the number of columns of 'x'")

  test <- match.arg(test,c("f","chi"))

  if (is.null(y) & test=="f") version <- "one.sample.f"
  if (is.null(y) & test=="chi") version <- "one.sample.chi"
  if (!is.null(y) & test=="f") version <- "two.sample.f"
  if (!is.null(y) & test=="chi") version <- "two.sample.chi"

  res1 <- switch(version,
                 "one.sample.f"={
                   result <- HotellingsT.internal(x,mu=mu,test=test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's one sample T2-test"
                   PARAMETER <- c(result$df.1,result$df.2)
                   names(PARAMETER) <- c("df1","df2")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)

                   RVAL}
                 ,
                 "one.sample.chi"={
                   result <- HotellingsT.internal(x,mu=mu,test=test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's one sample T2-test"
                   PARAMETER <- c(result$df.1)
                   names(PARAMETER) <- c("df")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)

                   RVAL}
                 ,
                 "two.sample.f"={
                   result <- HotellingsT.internal(x,y,mu,test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's two sample T2-test"
                   PARAMETER <- c(result$df.1,result$df.2)
                   names(PARAMETER) <- c("df1","df2")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)

                   RVAL}
                 ,
                 "two.sample.chi"={
                   result <- HotellingsT.internal(x,y,mu,test)
                   STATISTIC <- result$test.statistic
                   names(STATISTIC) <- "T.2"
                   PVAL <- result$p.value
                   METHOD <- "Hotelling's two sample T2-test"
                   PARAMETER <- c(result$df.1)
                   names(PARAMETER) <- c("df")
                   RVAL <- list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)

                   RVAL}
  )
  ALTERNATIVE="two.sided"
  NVAL <- paste("c(",paste(mu,collapse=","),")",sep="")
  if (is.null(y)) names(NVAL) <- "location" else names(NVAL) <- "location difference"
  res <- c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
  class(res) <- "htest"
  return(res)
}


HotellingsT2Test.formula <- function(formula, data, subset, na.action, ...) {

  if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
                                                                  "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if (nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  # DATA <- setNames(split(mf[[response]], g), c("x", "y"))
  DATA <- setNames(split(as.data.frame(mf[[response]]), g), c("x", "y"))

  y <- DoCall("HotellingsT2Test", c(DATA, list(...)))
  y$data.name <- DNAME
  y
}


###





#' Hosmer-Lemeshow Goodness of Fit Tests
#' 
#' The function computes Hosmer-Lemeshow goodness of fit tests for C and H
#' statistic as well as the le Cessie-van Houwelingen-Copas-Hosmer unweighted
#' sum of squares test for global goodness of fit.
#' 
#' Hosmer-Lemeshow goodness of fit tests are computed; see Lemeshow and Hosmer
#' (1982).
#' 
#' If \code{X} is specified, the le Cessie-van Houwelingen-Copas-Hosmer
#' unweighted sum of squares test for global goodness of fit is additionally
#' determined; see Hosmer et al. (1997). % A more general version of this test
#' is implemented in function % \code{\link[Design]{residuals.lrm}} in package
#' \pkg{Design}.
#' 
#' @param fit numeric vector with fitted probabilities.
#' @param obs numeric vector with observed values.
#' @param ngr number of groups for C and H statistic.
#' @param X covariate(s) for le Cessie-van Houwelingen-Copas-Hosmer global
#' goodness of fit test.
#' @param verbose logical, print intermediate results.
#' @return A list of tests.
#' 
#' @author 
#' Matthias Kohl <Matthias.Kohl@@stamats.de>
#' 
#' @seealso 
#' \code{\link{glm}}
#' 
#' @references 
#' Lemeshow, S. Hosmer, D.W., (1982): A review of goodness of fit
#' statistics for use in the development of logistic regression models.
#' \emph{American Journal of Epidemiology, \bold{115}(1), 92-106.}
#' 
#' Hosmer, D.W., Hosmer, T., le Cessie, S., Lemeshow, S. (1997). A comparison
#' of goodness-of-fit tests for the logistic regression model.
#' \emph{Statistics in Medicine}, \bold{16}, 965-980.
#' 
#' @keywords univar
#'
#' @examples
#' 
#' set.seed(111)
#' 
#' x1 <- factor(sample(1:3, 50, replace = TRUE))
#' x2 <- rnorm(50)
#' obs <- sample(c(0,1), 50, replace = TRUE)
#' 
#' fit <- glm(obs ~ x1+x2, family = binomial)
#' 
#' HosmerLemeshowTest(fit = fitted(fit), obs = obs, X = cbind(x1, x2))
#' 
HosmerLemeshowTest <- function(fit, obs, ngr = 10, X, verbose = FALSE) {

  # woher kommt das?? -> klaeren!
  # - > MKmisc

  ngr1 <- ngr
  # Hosmer-Lemeshow C statistic
  brks <- unique(quantile(fit, probs = seq(0, 1, by = 1/ngr)))
  cutfit <- cut(fit, breaks = brks, include.lowest = TRUE)
  if (length(brks) < ngr+1) {
    warning("Found only ", length(brks)-1, " different groups for Hosmer-Lemesho C statistic.")
    ngr <- length(brks)-1
  }
  if (verbose) {
    cat("Groups for Hosmer-Lemeshow C statistic:\n")
    print(table(cutfit))
  }
  Obs <- xtabs(cbind("0s" = 1 - obs, "1s" = obs) ~ cutfit)
  Exp <- xtabs(cbind("Os" = 1 - fit, "1s" = fit) ~ cutfit)
  chisq <- sum((Obs - Exp)^2/Exp, na.rm = TRUE)
  names(chisq) <- "X-squared"
  param <- ngr-2
  names(param) <- "df"
  P <- 1 - pchisq(chisq, param)

  # Hosmer-Lemeshow H statistic
  cutfit1 <- cut(fit, breaks = ngr1, include.lowest = TRUE)
  if (verbose) {
    cat("Groups for Hosmer-Lemeshow H statistic:\n")
    print(table(cutfit1))
  }
  Obs1 <- xtabs(cbind(1 - obs, obs) ~ cutfit1)
  Exp1 <- xtabs(cbind(1 - fit, fit) ~ cutfit1)
  chisq1 <- sum((Obs1 - Exp1)^2/Exp1, na.rm = TRUE)
  names(chisq1) <- "X-squared"
  param1 <- ngr1-2
  names(param1) <- "df"
  P1 <- 1 - pchisq(chisq1, param1)
  dname <- paste(deparse(substitute(fit)), "and", deparse(substitute(obs)))
  C <- structure(list(statistic = chisq, parameter = param,
                      p.value = P, method = "Hosmer-Lemeshow C statistic", data.name = dname,
                      observed = Obs, expected = Exp), class = "htest")
  H <- structure(list(statistic = chisq1, parameter = param1,
                      p.value = P1, method = "Hosmer-Lemeshow H statistic", data.name = dname,
                      observed = Obs1, expected = Exp1), class = "htest")


  if (!missing(X)) {
    # le Cessie-van Houwelingen-Copas-Hosmer unweighted sum of squares test for global goodness of fit
    #        X <- cbind(1, X)
    y <- obs == 1
    p <- fit
    sse <- sum((y - p)^2)
    wt <- p * (1 - p)
    d <- 1 - 2 * p
    z <- lm.wfit(X, d, wt, method = "qr")
    res <- z$residuals * sqrt(z$weights)
    sd <- sqrt(sum(res^2))
    ev <- sum(wt)
    z <- (sse - ev)/sd
    names(z) <- "z"
    P2 <- 2 * pnorm(abs(z), lower.tail = FALSE)
    stats <- c(sse, ev, sd, z, P)
    names(stats) <- c("Sum of squared errors", "Expected value|H0",
                      "SD", "Z", "P")
    gof <- structure(list(
      statistic = z, p.value = P2,
      method = "le Cessie-van Houwelingen-Copas-Hosmer global goodness of fit test",
      data.name = dname,
      observed = sse, expected = ev), class = "htest")
    
    return(list(C = C, H = H, gof = gof))
  }

  list(C = C, H = H)
}

