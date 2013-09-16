\name{BiGGR-package}
\alias{BiGGR-package}
\alias{BiGGR}
\docType{package}
\title{
Creates an interface to the BiGG database, provides a framework for simulation and produces flux graphs
}
\description{
This package provides an interface to simulate metabolic reconstruction 
	     from the BiGG database(http://bigg.ucsd.edu/) and other metabolic reconstruction 
	     databases. The package aids in performing flux balance analysis (FBA).  
	     Metabolic networks and estimated fluxes can be visualized
	     using hypergraphs.
}
\details{
\tabular{ll}{
Package: \tab BiGGR\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2013-07-10\cr
Depends: \tab R (>= 2.14.0), rsbml, hyperdraw, LIM\cr
Imports: \tab hypergraph\cr
License: \tab GPL (>=2)\cr
URL: \tab http://www.bioconductor.org/\cr
Copyright: \tab see inst/COPYRIGHTS for the license of the BiGG database\cr
biocViews: \tab NetworkAnalysis, Visualization, Metabolomics\cr
LazyLoad: \tab yes\cr
Packaged: \tab 2013-08-05 12:44:08 UTC; hettling\cr
Built: \tab R 3.0.0; ; 2013-08-05 12:44:22 UTC; unix\cr
}

Index:
\preformatted{
E.coli_iAF1260                                              Ecoli dataset with ORFs and thermodynamic information
E.coli_iJR904                                               Ecoli genome-scale model
Glycolysis                                                  Metabolic reconstruction of Glycolysis pathway
H.pylori_ilT341                                             H.pylori in silico genome-scale characterization of single and double deletion mutants
H.sapiens_Recon_1                                           Reconstruction of human metabolism from the BiGG database
M.barkeri_iAF692                                            Metabolic reconstruction of M.barkeri from the BiGG database
Recon2                                                      Human metabolic reconstruction Recon2
S.aureus_iSB619                                             Metabolic reconstruction of S.aureus from the BiGG database
S.cerevisiae_iND750                                         Metabolic reconstruction of S.cerevisiae from the BiGG database
buildSBMLFromBiGG                                           Build an SBML model from a given reactions file obtained from the BiGG database
buildSBMLFromGenes                                          Build an SBML model for specific genes in a given database
buildSBMLFromPathways                                       Build an SBML model for specific pathway(s) in a given database
buildSBMLFromReactionIDs                                    Build an SBML model for specific reactions in a given database
createLIMFromBiGG                                           Create a LIM model object from a BiGG database file
createLIMFromSBML                                           Create a LIM model object from an SBML file
extractGeneAssociations                                     Extract informations on genes from a given database
extractPathways                                             Extract all pathways from given database
getPathwaysForSBML                                          Extract all pathways from a database that are relevant for a given SBML model
getRates                                                    Get Optimized Rates
sbml2hyperdraw                                              Returns a graph representation of an SBML model
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{BiGGR} \tab BiGGR (source, pdf)\cr
}
}
\author{
Anand K. Gavai, Hannes Hettling

Maintainer: Anand K. Gavai <anand.gavai@bioinformatics.nl>, Hannes Hettling <j.hettling@vu.nl>

}

\keyword{biggr}
\seealso{
  \code{rsbml} 
}
\examples{


library("BiGGR")


library(help="BiGGR")

##load reaction identifiers from package examples
file.name <- system.file("extdata", 
                         "Glycolysis_TCA_recon1_reactionIDs.txt", 
                         package="BiGGR")
reaction.ids <- scan(file.name, what=" ")

##load database
data("H.sapiens_Recon_1")

##build SBML model
sbml.model <- buildSBMLFromReactionIDs(reaction.ids, H.sapiens_Recon_1)

##following term is to be maximized
maximize <- "R_ATPS4m - R_NDPK1m - R_HEX1 - R_PFK - R_PGK + R_PYK"

##specify the external metabolites of the system
externals <- c("M_glc_DASH_D_e", "M_lac_DASH_L_e",
               "M_ala_DASH_L_e", "M_gln_DASH_L_c", "M_h2o_e",
               "M_co2_e", "M_o2_e", "M_h_e", "M_pi_c",
               "M_o2s_m", "M_nh4_m", "M_adp_c",
               "M_atp_c", "M_nadp_c", "M_nadph_c", "M_h_c")
##specify the values of following fluxes: 
##R_GLCt1r=0.4, R_O2t=2.4, R_L_LACt2r=R_GLNtm=0

equation.vars <- c("R_GLCt1r", "R_O2t", "R_L_LACt2r", "R_GLNtm")
equation.values <- c(0.4, 2.4, 0.0, 0.0)
eqns <- list(equation.vars, equation.values)

##create LIM file 
limfile.name <- tempfile()
createLIMFromSBML(sbml.model, maximize, equations=eqns,
                  externals=externals, file.name=limfile.name)

rates <- getRates(limfile.name)

relevant.species <- c("M_glc_DASH_D_c", "M_g6p_c", "M_f6p_c", 
                      "M_fdp_c", "M_dhap_c", "M_g3p_c", 
                      "M_13dpg_c", "M_3pg_c", "M_2pg_c", 
                      "M_pep_c", "M_pyr_c")
##generate graphical representation
hd <- sbml2hyperdraw(sbml.model, rates=rates, 
                     relevant.species=relevant.species, 
                     layoutType="dot", plt.margins=c(20, 0, 20, 0))

##plot hypergraph
plot(hd)


}

