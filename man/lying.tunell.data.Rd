\name{lying.tunell.data}
\docType{data}
\alias{lying.tunell.data}
\title{Dataset of in vivo cerebral metabolite uptake and release rates in healthy humans (old subjects)}
\description{These data were taken from a publication of Lying-Tunell et al. (1980) reporting cerebral metabolic uptakes and release rates
             in older subjects (n=5).  The data were published as micromole/kg/min, but converted to mmole/min for this dataset (see details).

}
\usage{data(lying.tunell.data)}
\format{An object of class data.frame}

\details{Data were taken from table 2 (page 271) of the publication.
	 From the given median and range values, mean and standard deviation was estimated using a method
	 by   Hozo et al. (2005). Units were converted from micromole/kg/min to mmole/min assuming a brain mass of 1.4kg.}
	 

\source{http://www.ncbi.nlm.nih.gov/pubmed/7468149}
\references{
Lying-Tunell U, Lindblad BS, Malmlund HO, Persson B: Cerebral blood flow and metabolic rate of
oxygen, glucose, lactate, pyruvate, ketone bodies and amino acids. Acta Neurol Scand 1980, 62:265–75.

Hozo SP, Djulbegovic B, Hozo I: Estimating the mean and variance from the median, range, 
and the size of a sample. BMC Med Res Methodol 2005, 5:13.
}

\examples{
\dontrun{
##The dataset was generated as follows:

##Uptake rates given in micromole/kg/min from Lying-Tunell (1980), n=5 old patients
##converted to mmol/min and assuming a brain mass of 1.4 kg
brain.mass <- 1.4 ## in kg
oxygen.median <- 1679 * brain.mass / 1000
oxygen.range <- c(1184, 1872) * brain.mass / 1000
glucose.median <- 203 * brain.mass / 1000
glucose.range <- c(187, 321) * brain.mass / 1000
lactate.median <- -9.2 * brain.mass / 1000
lactate.range <- c(-68, 7.9) * brain.mass / 1000
pyruvate.median <- -2.4 * brain.mass / 1000
pyruvate.range <- c(-10, -brain.mass) * brain.mass / 1000
glutamine.median <- -11 * brain.mass / 1000
glutamine.range <- c(-61, 22) * brain.mass / 1000

##This implements eq 4 from Hozo et al. to estimate
##sample mean from median and range
##m: median, a: minimum, b: maximum, n: number of samples
estimate.sample.mean <- function(m, a, b, n)
  (a + 2*m + b)/4 + (a-2*m + b)/(4*n)

##This implements eq 16 from Hozo et al. to estimate
##sample standard deviation from median and range
##m: median, a: minimum, b: maximum, n: number of samples
estimate.sample.sd <- function(m, a, b, n)
  sqrt((((a - 2*m + b)^2)/4 + (b-a)^2)/12)

##Calculate mean and standard deviation from median and range values using the method of Hoxo et al. 
oxygen.mean <- estimate.sample.mean(oxygen.median, oxygen.range[1], oxygen.range[2], 5)
oxygen.sd <- estimate.sample.sd(oxygen.median, oxygen.range[1], oxygen.range[2], 5)

glucose.mean <- estimate.sample.mean(glucose.median, glucose.range[1], glucose.range[2], 5)
glucose.sd <- estimate.sample.sd(glucose.median, glucose.range[1], glucose.range[2], 5)

lactate.mean <- estimate.sample.mean(lactate.median, lactate.range[1], lactate.range[2], 5)
lactate.sd <- estimate.sample.sd(lactate.median, lactate.range[1], lactate.range[2], 5)

pyruvate.mean <- estimate.sample.mean(pyruvate.median, pyruvate.range[1], pyruvate.range[2], 5)
pyruvate.sd <- estimate.sample.sd(pyruvate.median, pyruvate.range[1], pyruvate.range[2], 5)

glutamine.mean <- estimate.sample.mean(glutamine.median, glutamine.range[1], glutamine.range[2], 5)
glutamine.sd <- estimate.sample.sd(glutamine.median, glutamine.range[1], glutamine.range[2], 5)


lying.tunell.data <- data.frame(median=c(oxygen.median, glucose.median, lactate.median, pyruvate.median, glutamine.median),
                                mean=c(oxygen.mean, glucose.mean, lactate.mean, pyruvate.mean, glutamine.mean),
                                sd=c(oxygen.sd, glucose.sd, lactate.sd, pyruvate.sd, glutamine.sd),
                                low=c(oxygen.range[1], glucose.range[1], lactate.range[1], pyruvate.range[1], glutamine.range[1]),
                                high=c(oxygen.range[2], glucose.range[2], lactate.range[2], pyruvate.range[2], glutamine.range[2]),
                                row.names=c("o2", "glucose", "lactate", "pyruvate", "glutamine"))

save(lying.tunell.data, file="lying.tunell.data.rda")


}

##load data
data(lying.tunell.data)
##get median value for glucose uptake
lying.tunell.data["glucose", "median"]
}

\keyword{cerebral metabolic uptake rates dataset}
