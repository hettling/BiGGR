\name{createLIMFromBiGG}
\encoding{latin1}
\Rdversion{1.1}
\alias{createLIMFromBiGG}

\title{Create a LIM model object from a BiGG database file}

\description{Creates a LIM model object from a file containing reactions
  extreacted from BiGG to be run for simulations of metabolic fluxes}

\usage{createLIMFromBiGG(reactions.filename, ...)}

\arguments{

  \item{reactions.filename}{file which contains the reactions extracted
    from the BiGG database}
  \item{...}{arguments passed to createLIMfromSBML}
  
}
\references{Soetaert K, van Oevelen D (2009). LIM: Linear Inverse Model examples and solution methods.
R package version 1.3}
\author{Anand K. Gavai <anand.gavai@bioinformatics.nl>, Hannes hettling <j.hettling@vu.nl>}
\note{none}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{createLIMFromSBML}}
}

\examples{

##maximize flux for reaction R_PYK
maximize <- "R_PYK"  

##setting equality constraint R_HEX = 1
equation_var <- "R_HEX1"  
equation_value <- 1  
eq <- list(equation_var, equation_value)

##range  of  possible  fluxes for R_PYK 
constraint <- list("R_PYK", 0, 1000)  
externals  <-  c("glc_c",  "pyr_c",  "h_c","nad_c",
                 "nadh_c",  "pi_c",  "fad_m",  "fadh2_m",
                 "o2_c",  "adp_c",  "atp_c",  "nadp_c",
                 "co2_c",  "o2_c",  "gdp_c",  "gtp_c")

##build LIM model from reactions file in package examples
path <- system.file("extdata", "Reactions.txt", package="BiGGR")
limfile.name <- tempfile()
createLIMFromBiGG(path, maximize, equations=eq, constraints=constraint,
                  externals=externals, file.name=limfile.name)

}
\keyword{Linear Inverse Models}
\keyword{Linear optimization model file}

