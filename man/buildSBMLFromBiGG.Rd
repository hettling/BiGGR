\name{buildSBMLFromBiGG}
\alias{buildSBMLFromBiGG}

\title{
  Build an SBML model from a given reactions file obtained from the BiGG database
}
\description{
  Creates an SBML model containing all species, reactions and
  compartments that occur in a reactions file obtained from the BiGG
  database.
}
\usage{
buildSBMLFromBiGG(reactions.filename, model.id=character(0), model.name=character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reactions.filename}{name of the file containing the reactions
    extracted from BiGG}
  \item{model.id}{id for the SBML model created by the
    function. Defaults to \code{reactions.filename}}    
  \item{model.name}{name for the SBML model created by the
    function. Defaults to \code{reactions.filename}}    
  
}

\value{
  a rsbml \code{Model} object containing all reactions, species and
  compartments that are associated with the reactions in the given input file.
}
\references{
  Schellenberger, J., Park, J. O., Conrad, T. C., and Palsson, B. , BiGG: a Biochemical Genetic and Genomic knowledgebase of large scale metabolic reconstructions, BMC Bioinformatics, 11:213, (2010).
    http://bigg.ucsd.edu/
}
\author{
Anand Gavai, Hannes Hettling
}
\note{
  Note that it can be the case that species can be present in multiple compartments.
  In order to avoid any ambiguities in the model returned by the
  function, each species identifier is composed of the
  species identifier given in the reactions file and the compartment
  identifier, joined by "_". Example: adp in compartment c (cytosol) has
  the id \code{adp_c}.
  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{createLIMFromSBML}}
}
\examples{

##build model from file Reactions.txt from the package examples  
path <- system.file("extdata", "Reactions.txt", package="BiGGR")
model <- buildSBMLFromBiGG(path, model.id="myid") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BiGG }

