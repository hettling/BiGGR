\name{extractPathways}
\alias{extractPathways}

\title{
Extract all pathways from given database
}
\description{
Extracts all pathway information from an rsbml document containing a
metabolic reconstruction database (e.g. Recon2).
The pathway information is parsed from the "<notes>" tag of each reaction.
}
\usage{
extractPathways(database)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{database}{an object of class \code{\linkS4class{SBMLDocument}}}
}

\value{
a \code{list} with length being the number of reactions
in the database passed as argument each entry containing a \code{character} with the pathway
information and the reaction IDs as names. For reactions without
pathway annotation, the list will contain \code{NA}.
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
\code{\link{getPathwaysForSBML}}
}
\examples{

data(Recon2)
pathways.recon2 <- extractPathways(Recon2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pathway }

