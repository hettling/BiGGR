\name{buildSBMLFromPathways}
\alias{buildSBMLFromPathways}

\title{
Build an SBML model for specific pathway(s) in a given database
}
\description{
Creates an SBML model containing all species, reactions and
compartments that are part of (a) given pathway(s) in the database
document (e.g. Recon2) passed as an argument.
}
\usage{
buildSBMLFromPathways(query, database, match.exact=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{query}{a \code{character} or a \code{vector} or \code{list} of
\code{character} containing the names of the query pathways}
\item{database}{an object of class \code{\linkS4class{SBMLDocument}}}
\item{match.exact}{\code{logical} whether only the exact pathway name
should be matched or whether a pathway should match if one
keyword is in the pathway description in the database.}    
}

\value{
a rsbml \code{Model} object containing all reactions, species and
compartments that are present in the database for the query
pathway(s) or NULL if none of the pathways in the database match the query.
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
indicated by the string "SUBSYSTEM", no pathway information can be
extracted and thus the SBML model returned will be empty.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{extractPathways}}
}
\examples{

data("Recon2")
database <- Recon2

##Get Model for specific pathway
m1 <- buildSBMLFromPathways("Arginine and Proline Metabolism", database)

##Get Model for specific pathway "Metabolism": does not exist!
m2 <- buildSBMLFromPathways("Metabolism", database)

##Get model of all pathways which contain keyword "metabolism"
m3 <- buildSBMLFromPathways("Metabolism", database, match.exact=FALSE)

##Multi-query:
query <- c("Transport, endoplasmic reticular", "Arginine and Proline Metabolism")
m4 = buildSBMLFromPathways(query, database)
m5 = buildSBMLFromPathways(query[1], database)
length(m4@species)
length(m5@species)

##different database
data(H.pylori_ilT341)
database <- H.pylori_ilT341
m7 <- buildSBMLFromPathways("Metabolism", database, match.exact=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathways }

