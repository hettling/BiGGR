\name{NEWS}
\title{NEWS file for the BiGGR package}
\encoding{UTF-8}

\section{Changes in version 1.1.5}{
  \itemize{
  Three new functions
    \item rmvSpliceVariant: remove splicing variants from entrez gene ids.
	\item gprMapping: maps gene expression values based on GPR rules as found in ReconX database. GPR rules are interpreted based on summary statistics.
	\item gprMappingAvg: maps gene expression values with GPR rules based on averages
	\item Recon2: Updated new version of Recon2 from www.humanmetabolism.org
  }
}

\section{Changes in version 1.1.4}{
  \itemize{
    \item Updated vignette. Equality constraints on GABA shunt and
    pentose phosphate pathway fluxes are now used in FBA and ensemble
    simulations.
  }
}

\section{Changes in version 1.1.3}{
  \itemize{
    \item Updated vignette and bibtex file. FBA and ensemble sampling
    are now done with integration of brain metabolite uptake measurement
    data.
    \item Added new function 'sampleFluxEnsemble' which can be used to
    generate a posterior distribution of flux vectors with respect to
    measurement error in certain fluxes. An example of
    'sampleFluxEnsemble' has been added to the package vignette. 
    \item Added LIM model file 'Glycolysis_TCA.LIM' to the examples
    directory. This file contains the model which is build in the
    examples.    
    \item Fixed warnings in functions buildSBMLFromReactionIDs
    and createLIMFromSBML to seperate missing identifiers with a comma.
  }
}

\section{Changes in version 0.99.3}{
  \itemize{
    \item Initial release of BiGGR in Bioconductor. Package was formerly
    available on CRAN (http://cran.r-project.org/src/contrib/Archive/BiGGR/).    
  }
}
