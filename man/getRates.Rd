\name{getRates}
\encoding{latin1}
\Rdversion{1.1}
\alias{getRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get Optimized Rates}
\description{
getRates takes the model file as the argument and based on the description of the model file generates flux values for "minimum" or "maximum" reaction rates
}
\usage{
getRates(modelFile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelFile}{A model file as generated from the createBiggModel function}
}

\value{The value returned is one dimentional numeric vector of flux rates for each reaction }

\author{Anand K. Gavai <anand.gavai@bioinformatics.nl>}

\examples{
	data("Glycolysis")
	rates <-getRates(Glycolysis)
	rates
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{flux rates}
\keyword{Linear Optimization}% __ONLY ONE__ keyword per line
