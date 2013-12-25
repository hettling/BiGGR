\name{sbml2hyperdraw}
\encoding{latin1}
\Rdversion{1.1}
\alias{sbml2hyperdraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns a graph representation of an SBML model}
\description{Convert an SBML model to a \code{RagraphBPH} using
  \code{hypergraph}. Metabolites are displayed as nodes and reactions are
  displayed as directed edges connecting the nodes. If a vector of rates
  is given, edge widths are weighted according to the rates. For
  negative rates, edges are drawn in red and the arrow between the
  metabolites is reversed to represent the correct direction of the flux.}


\usage{sbml2hyperdraw(sbml.model, rates ,relevant.species,
relevant.reactions,layoutType, lwd.max, lwd.min, plt.margins)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sbml.model}{an rsbml \code{Model} object}
  \item{rates}{a named \code{vector} with the rates of the reactions in
    the model. The names of the rates must agree with the reaction
    identifiers in the sbml.model}
  \item{relevant.species}{a \code{vector} of type \code{character}
    defining a subset of species in the sbml.model to be
    plotted. Defaults to all species identifiers in the sbml.model.}
  \item{relevant.reactions}{a \code{vector} of type \code{character}
    defining a subset of reactions in the sbml.model to be
    plotted. Defaults to all reactions identifiers in the sbml.model.}
  \item{layoutType}{is a character string representing the layout engine to
    be used for visualization. Current supported layouts are
    "dot", "twopi","neato","fdp","sfdp" and "circo". Defaults to
    "dot". See \code{?GraphvizLayouts} for further documentation.}
  \item{lwd.max}{a \code{numeric} given the maximum edge
    width. Defaults to 3.}
  \item{lwd.min}{a \code{numeric} given the minimum edge
    width. Defaults to 0.5.}
  \item{plt.margins}{A numerical vector of the form c(bottom, left,
    top, right) giving additional white space around the graph (in case
    long node or edge labels fall outside the plotting region). Defaults
    to c(150,150,150,150).}
}

\value{
  Object of class \code{RagraphBPH} with the hypergraph representation of the SBML object.
}

\author{Hannes Hettling <j.hettling@vu.nl>, Anand K. Gavai <anand.gavai@bioinformatics.nl>}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  RagraphBPH
  hyperdraw
  
}
\examples{

##Generate an example model
path <- system.file("extdata", "Glycolysis_TCA_recon2_reactionIDs.txt", package="BiGGR")
reaction.ids <- scan(path, what=" ")

data("Recon2")
model <- buildSBMLFromReactionIDs(reaction.ids, Recon2)

##Plot ATP and ADP in cytosol and mitochondrion in model without rates
rel.sp <- c("M_adp_c", "M_atp_c", "M_adp_m", "M_atp_m")
hd <- sbml2hyperdraw(model, relevant.species=rel.sp)
plot(hd)

##Plot model with random rates
rates <- rnorm(length(model@reactions))
names(rates) <- sapply(model@reactions, id)
hd <- sbml2hyperdraw(model, rates=rates, relevant.species=rel.sp, lwd.max=4)
plot(hd)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hypergraph}
\keyword{hyperdraw}
