\name{buildSBMLFromReactionIDs}
\alias{buildSBMLFromReactionIDs}

\title{
Build an SBML model for specific reactions in a given database
}
\description{
Creates an SBML model containing all species, reactions and
compartments that are associated with a number of reaction identifiers
in the database document (e.g. Recon2) passed as an argument.
}
\usage{
buildSBMLFromReactionIDs(reaction.ids, database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{reaction.ids}{a \code{character} or a \code{vector} or \code{list} of
\code{character} containing the names of the query reaction IDs}
\item{database}{an object of class \code{\linkS4class{SBMLDocument}}}
}
\value{
a rsbml \code{Model} object containing all reactions, species and
compartments that are present in the database for the query
reaction(s) or NULL if none of the query reaction IDs is found in the database.
}
\references{
Thiele, I. et al. Nat Biotech, 2013
}
\author{
Anand Gavai, Hannes Hettling
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{buildSBMLFromGenes}}
}
\examples{

##get list of reactions with Recon 2 identifiers from examples
path <- system.file("extdata", "Glycolysis_TCA_recon2_reactionIDs.txt", package="BiGGR")
reaction.ids <- scan(path, what=" ")

data("Recon2")
model <- buildSBMLFromReactionIDs(reaction.ids, Recon2)

}

\keyword{ reactions }

