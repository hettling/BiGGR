\name{createLIMFromSBML}
\encoding{latin1}
\Rdversion{1.1}
\alias{createLIMFromSBML}

\title{Create a LIM model object from an SBML file}

\description{creates a model file to be run for simulations of metabolic fluxes}

\usage{createLIMFromSBML(model, maximize, equations=NULL,
            constraints=NULL, externals=NULL, file.name="model.lim")}

\arguments{

  \item{model}{an SBML object of reactions/metabolites participating in a metabolic pathway.}
  
  \item{maximize}{
    a character vector consisting the tag of the reaction(s) to be maximized or minimized
  }
  \item{equations}{
    a \code{list} specifying equality constraints on the system.
    The list must have two entries, the first one being a
    vector of class \code{character} containing the left hand side(s) of
    the equation(s), the second one being a vector of type
    \code{character} or \code{mumeric} with the right hand side(s) of
    the equation(s). See also 'examples'.
  }
  \item{constraints}{
    a \code{list} specifying constrained on the solution space of the
    flux vector. The list must have three entries, the first one being
    a vector of class \code{character} with the reaction id(s) to be
    constrained, the second and third one a \code{numeric} vector with
    the lower and upper flux bounds, respectively,  for the reactions to be contrained.
    is a character vector specifying the minimum and maximum values(boundary) under which the solution for the 	\code{maximize} reaction should fall
  }
  \item{externals}{
    a character vector of metabolites as provided by the user for speficific pathways for which FBA (flux balance analysis needs to be performed)
  }
  \item{file.name}{
    a character string specifying the name of the LIM file created by
    the function.
  }
}

\value{A model file with with extension \code{".lim"} is created}

\references{Soetaert K, van Oevelen D (2009). LIM: Linear Inverse Model examples and solution methods.
R package version 1.3}
\author{Anand K. Gavai <anand.gavai@bioinformatics.nl>, Hannes Hettling <j.hettling@vu.nl>}
\note{none}



%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##Create a LIM model file from a reactions file in the examples
path <- system.file("extdata", "Glycolysis_TCA_recon1_reactionIDs.txt", package="BiGGR")
reaction.ids <- scan(path, what=" ")
data("H.sapiens_Recon_1")
sbml.model <- buildSBMLFromReactionIDs(reaction.ids, H.sapiens_Recon_1)
maximize <- c("R_ATPS4m - R_NDPK1m - R_HEX1 - R_PFK  - R_PGK + R_PYK")
externals <- c("M_glc_DASH_D_e", "M_lac_DASH_L_e",
               "M_ala_DASH_L_e", "M_gln_DASH_L_c", "M_h2o_e", 
               "M_co2_e", "M_o2_e", "M_h_e", "M_pi_c",
               "M_o2s_m", "M_nh4_m", "M_adp_c",
               "M_atp_c", "M_nadp_c", "M_nadph_c", "M_h_c")
equation.vars <- c("R_GLCt1r", "R_O2t", "R_L_LACt2r", "R_GLNtm")
equation.values <- c(0.4, 2.4, 0.0, 0.0)
eqns <- list(equation.vars, equation.values)
constraints <- list(c("R_GLCt1r", "R_CYOOm3"), c(-1000, -1000), c(1000, 1000)) 
limfile.name <- tempfile()
createLIMFromSBML(sbml.model, maximize, equations=eqns,
                  constraints=constraints, externals=externals, file.name=limfile.name)

}
\keyword{Linear Inverse Models}
\keyword{Linear optimization model file}
\keyword{SBML object}

