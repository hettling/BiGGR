\name{buildSBMLFromGenes}
\alias{buildSBMLFromGenes}

\title{
Build an SBML model for specific genes in a given database
}
\description{
Creates an SBML model containing all species, reactions and
compartments that are associated with (a) specific gene(s) in the
database document (e.g. Recon2) passed as an argument. 
}
\usage{
buildSBMLFromGenes(query, database, logical.fun="any")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{query}{a \code{character} or a \code{vector} or \code{list} of
\code{character} containing the query genes with identifiers as
specified in the database.
}
\item{database}{an object of class \code{\linkS4class{SBMLDocument}}}
\item{logical.fun}{function which specifies the logical
relation of the query genes within the reactions
(e.g. \code{all} or \code{any}, see details).}
}
\details{
The function \code{all} as argument \code{logical.fun}
would mean that all genes in the query have to be associated with
a certain reaction from the database in order to be included in the returned
model. The default \code{any} means that a reaction is included if any
of the query genes are associated with it. Custom functions are
possible if they take a \code{vector} of type \code{logical} as an argument
and return a \code{logical}. The argument of logical.fun is a
\code{vector} of type \code{logical} having the same length as the
query and for each gene the value is \code{TRUE} if it is associated
with a specific reaction. See 'examples' section for an example of a
custom function as logical.fun.
}

\value{
a rsbml \code{Model} object containing all reactions, species and
compartments that are present in the database and are associated
with the query
gene(s) or NULL if none of the genes in the database match the query.
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

%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{buildSBMLFromPathways}}
\code{\link{extractGeneAssociations}}
}
\examples{


##Query genes in Recon 2 database
data("Recon2")
database <- Recon2
m1 <- buildSBMLFromGenes("8884.1", database)
m2 <- buildSBMLFromGenes(c("8884.1", "6509.1"), database)

##different databases
data(H.pylori_ilT341)
database <- H.pylori_ilT341
m3 <- buildSBMLFromGenes("HP0069", database)

data(M.barkeri_iAF692)
database <- M.barkeri_iAF692
m4 <- buildSBMLFromGenes(c("MBd0456", "MBd4814", "MBd4098"), database)


data(S.aureus_iSB619)
database <- S.aureus_iSB619
m5 <- buildSBMLFromGenes(c("SA0594", "SA1599", "SA0950", "SA0259"), database)

database <- Recon2
query <- c("218.1", "223.1")
m6 <- buildSBMLFromGenes(query, database)
m7 <- buildSBMLFromGenes(query, database, logical.fun="all")
##m6 has more reactions than m7
## because m7 has only reactions which match both genes in the query
length(m6@reactions) > length(m7@reactions)

##Custom logical function: Get model with all reactions
## which are not associated with the query gene
m8 <- buildSBMLFromGenes(query, database, logical.fun=function(x)!any(x))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gene }

