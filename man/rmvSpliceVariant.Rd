\name{rmvSpliceVariant}
\alias{rmvSpliceVariant}

\title{
Remove splicing variants from the database.
}
\description{
Removes alternative splicing information from the database.
}
\usage{
rmvSpliceVariant(gene.info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{gene.info}{A reaction gene maping from the ReconX database created from functions \code{extractGeneAssociations}}
}

\value{
A n x 2 dimentional dataframe of Reaction-Gene(Entrez number) mapping from ReconX database
}
\references{
Thiele, I. et al. Nat Biotech, 2013
}
\author{
Anand Gavai <anand.gavai@bioinformatics.nl>, Hannes Hettling
}

\examples{
data(Recon2)
gene.info <- extractGeneAssociations(Recon2)

gene.info<-do.call(rbind.data.frame,gene.info)
colnames(gene.info)<-c("GPR")
gene.info$react_id<-row.names(gene.info)
gene.info<-gene.info[,c(2,1)]
rownames(gene.info)<-NULL

react_gene_map<-rmvSpliceVariant(gene.info)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{GPR Mapping}

