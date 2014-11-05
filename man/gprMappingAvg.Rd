\name{gprMappingAvg}
\encoding{latin1}
\Rdversion{1.1}
\alias{gprMappingAvg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GPR mapping ignoring AND & OR operators}
\description{
Continuous gene expression levels are mapped from genes to reactions using the gene-protein-reaction (GPR) association rules as found in ReconX databases. These rules are comprised of \code{AND} and \code{OR} operators. This function ignores these rules and take average of all genes
}
\usage{
gprMappingAvg(gene_express,react_gene_map)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gene_express}{The path to a gene expression file with three columns gene_symbol,entrez_id and foldchanges}
\item{react_gene_map}{Database file created from ReconX database using functions such as \code{rmvSpliceVariant}}
}

\value{Returns a dataframe with Reaction_id, GPR formulae and average values }

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

gpr.map.avg<-gprMappingAvg(gene_express,react_gene_map)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GPR Mapping}
