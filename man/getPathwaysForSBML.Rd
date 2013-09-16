\name{getPathwaysForSBML}
\alias{getPathwaysForSBML}

\title{
  Extract all pathways from a database that are relevant for a given SBML model 
}
\description{
  Extracts all pathway information from an rsbml document containing a
  metabolic reconstruction database (e.g. Recon2) and returns the subset
  of these pathways that is associated with the reactions in the given model.
  The pathway information is parsed from the "<notes>" tag of each reaction.
}
\usage{
getPathwaysForSBML(model, database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{an rsbml object of class \code{\linkS4class{SBMLDocument}}}
  \item{model}{an rsbml object of class \code{\linkS4class{Model}}}
}

\value{
  A vector of type \code{character} that contains all the pathways
  relevant for the given model according to the specified database. Note
  that duplicate pathways do not appear twice in the return value.
}
\references{
  Thiele, I. et al. Nat Biotech, 2013
}
\author{
Anand Gavai, Hannes Hettling
}
\note{
  If the reactions in the database document provided in the argument
  \code{database} do not contain any "<notes>" with tags with pathway information
  indicated by the string "SUBSYSTEM", no pathway information can be extracted.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{buildSBMLFromPathways}}
  \code{\link{buildSBMLFromGenes}}
  \code{\link{extractPathways}}
}
\examples{

##Build a model from query genes
data("Recon2")
database <- Recon2
query <- c("218.1", "223.1") ##query gene identifiers
m <- buildSBMLFromGenes(query, database)

##extract all pathways for that model
getPathwaysForSBML(m, database)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathway }

