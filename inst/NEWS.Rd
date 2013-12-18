\name{NEWS}
\title{NEWS file for the BiGGR package}
\encoding{UTF-8}

\section{Changes in version 1.1.3}{
  \itemize{
    \item Updated vignette and bibtex file. FBA and ensemble sampling
    are now done with integration of brain metabolite uptake measurement
    data.
    \item Added new function 'sampleFluxEnsemble' which can be used to
    generate a posterior distribution of flux vectors with respect to
    measuerement error in certain fluxes. An example of
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
    available on CRAN.    
  }
}
