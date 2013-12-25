\name{extractGeneAssociations}
\alias{extractGeneAssociations}

\title{
  Extract informations on genes from a given database
}
\description{
  Extracts all information on genes associated to reactions from an rsbml document containing a
  metabolic reconstruction database (e.g. Recon2).
  The associated information is parsed from the "<notes>" tag of each
  reaction's SBML representation.
}
\usage{
extractGeneAssociations(database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database}{an object of class \code{\linkS4class{SBMLDocument}}}
}

\value{
  a \code{list} with length being the number of reactions
    in the database passed as argument each entry containing a
    \code{character} containing the assiciated gene identifiers
    and the reaction IDs as names. For reactions without gene
    annotation, the list will contain \code{NA}.
}
\references{
  Thiele, I. et al. Nat Biotech, 2013
}
\author{
Anand Gavai, Hannes Hettling
}
\note{
  If the reactions in the database document provided in the argument
  \code{database} do not contain any "<notes>" with tags with gene information
  indicated by the string "GENE*ASSOCIATION" (the star stands for any
  character), 
  no gene association information can be
  extracted and thus the returned SBML mdel is empty..
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{buildSBMLFromGenes}}
}
\examples{

data("Recon2")
database <- Recon2
gene.info <- extractGeneAssociations(database)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene }

