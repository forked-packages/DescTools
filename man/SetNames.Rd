\name{SetNames}
\alias{SetNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set the Names in an Object
%%  ~~function to do ... ~~
}
\description{This is a convenience function that sets the names of an object and returns it including the new names. It is most useful at the end of a function definition where one is creating the object to be returned and would prefer not to store it under a name just that the names can be assigned. In addition to the function \code{\link{setNames}} in base R the user can decide, whether rownames, colnames or simply the names are to be set. Names are recyled. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SetNames(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object for which a names attribute will be meaningful
%%     ~~Describe \code{object} here~~
}
  \item{...}{the names to be assigned to the object. This should be a character vector of names named \code{dimnames}, \code{rownames}, \code{colnames} or \code{names}. Setting \code{rownames=NULL} would remove existing rownames. All kind of names can be changed at the same time. Default would be \code{names}. Abbreviations are supported.
%%     ~~Describe \code{nm} here~~
}
}
\value{An object of the same sort as object with the new names assigned.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>\cr
}

\seealso{\code{\link{setNames}}, \code{\link{Rename}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
SetNames(1:5, names=letters[1:5])

# the default, if no argument names are provided, is "names"
SetNames(1:5, letters[1:5])

tab <- table(d.pizza$driver, d.pizza$wine_delivered)
# rownames and columnnames can be set at the same time
SetNames(BinomCI(tab[,1], rowSums(tab)), 
         rownames=rownames(tab), colnames=c("perc", "lci", "uci"))
         
# can also be used to set the names to an empty string
SetNames(diag(6), rownames="", colnames="")

# setting dimnames works as well
tab <- SetNames(
  as.table(rbind(c(84,43), c(10,92))), 
    dimnames= list(
       dipstick=c("positive","negative"),
       culture=c("positive","negative")))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ list }



