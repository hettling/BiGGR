\name{gprMapping}
\encoding{latin1}
\Rdversion{1.1}
\alias{gprMapping}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GPR mapping}
\description{
Continuous gene expression levels are mapped from genes to reactions using the gene-protein-reaction (GPR) association rules as found in ReconX databases. The expression level of reactions catalyzed by enzyme complexes (\code{and} operator) can be set to the minimum,maximum,mean and median functions. Similarly expression level of the associated genes, and the expression level of reactions catalyzed by isoenzymes (\code{or} operator) can also be set to either minimum,maximum,mean and median functions for the associated genes.
Operator Precedence: "AND" followed by "OR"
}
\usage{
gprMapping(gene_express,react_gene_map,OR=c("mean","median","min","max"),AND=c("min","max","mean","median"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gene_express}{The path to a gene expression file with three columns gene_symbol,entrez_id and foldchanges}
\item{react_gene_map}{Database file created from ReconX database using functions such as \code{rmvSpliceVariant}}
\item{OR}{Takes values from statistical functions such as \code{mean},\code{median},\code{min},\code{max}}
\item{AND}{Takes values from statistical functions such as \code{mean},\code{median},\code{min},\code{max}}
}

\value{Returns a dataframe with Reaction_id, GPR formulae and Calculated values }

\author{Anand K. Gavai <anand.gavai@bioinformatics.nl>, Hannes Hettling}

\examples{
# Read gene expression data
file <- system.file("extdata", "Gene_Symbol_Entrez_Foldchanges.csv", package="BiGGR")
gene_express<-read.csv(file,header=TRUE)
data(Recon2)
gene.info <- extractGeneAssociations(Recon2)

gene.info<-do.call(rbind.data.frame,gene.info)
colnames(gene.info)<-c("GPR")
gene.info$react_id<-row.names(gene.info)
gene.info<-gene.info[,c(2,1)]
rownames(gene.info)<-NULL
react_gene_map<-rmvSpliceVariant(gene.info)

gpr.map<-gprMapping(gene_express,react_gene_map,OR="mean",AND="min")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GPR Mapping}
\keyword{\code{AND} operator}% __ONLY ONE__ keyword per line
\keyword{\code{OR} operator}