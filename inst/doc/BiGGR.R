### R code from vignette source 'BiGGR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: e (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("BiGGR")


###################################################
### code chunk number 2: f
###################################################
library("BiGGR")


###################################################
### code chunk number 3: g (eval = FALSE)
###################################################
## library(help="BiGGR")


###################################################
### code chunk number 4: h
###################################################
##load reaction identifiers from package examples
file.name <- system.file("extdata", 
                         "brainmodel_reactions.txt", 
                         package="BiGGR")
reaction.ids <- scan(file.name, what=" ")

##load database
data("H.sapiens_Recon_1")

##build SBML model
sbml.model <- buildSBMLFromReactionIDs(reaction.ids, H.sapiens_Recon_1)


###################################################
### code chunk number 5: i
###################################################
##following term is to be maximized
maximize <- "R_ATPS4m - R_NDPK1m -R_HEX1 - R_PFK - R_PGK + R_PYK"

##specify the external metabolites of the system
externals <- c("M_glc_DASH_D_e", "M_lac_DASH_L_e", "M_ala_DASH_L_e",
               "M_gln_DASH_L_e", "M_h2o_e", "M_co2_e",
               "M_o2_e", "M_h_e", "M_o2s_m",
               "M_adp_c", "M_atp_c", "M_pi_c",
               "M_h_c", "M_nadp_c", "M_nadph_c",
               "M_na1_c", "M_na1_e", "M_gln_DASH_L_c",
               "M_nh4_c", "M_pyr_e")


###################################################
### code chunk number 6: eq_ineq
###################################################

##load lying-tunell data
data(lying.tunell.data)

##set equality constraints
equation.vars <- c("R_GLCt1r", "R_L_LACt2r", "R_GLNtN1", 
                   "R_PYRt2r", "R_GLUDC", "R_G6PDH2r")

equation.values <- c(as.character(
    lying.tunell.data[c("glucose", "lactate", "glutamine", "pyruvate"), 
                      "median"]), 
                     "R_GLCt1r * 0.32", "R_GLCt1r * 0.069" )

eqns <- list(equation.vars, equation.values)

##write LIM file to system's temporary directory
limfile.path <- tempfile()
createLIMFromSBML(sbml.model, maximize, equations=eqns, 
                  externals=externals, file.name=limfile.path)


###################################################
### code chunk number 7: j
###################################################
rates <- getRates(limfile.path)


###################################################
### code chunk number 8: random_seed
###################################################
set.seed(111)


###################################################
### code chunk number 9: sample_flux_ensemble
###################################################
##specify the fluxes with uncertainty given as SD in a data frame
uncertain.vars <- data.frame(var=equation.vars[1:4],
                             value=equation.values[1:4],
                             sd=c(0.058,0.032,0.034,0.004))

uncertain.vars <- data.frame(var=c(equation.vars[c(1,2,3,4)]),
                             value=as.numeric(c(equation.values[c(1,2,3,4)])),
                             sd=lying.tunell.data[c("glucose", 
                                 "lactate", "glutamine", "pyruvate"), "sd"])


limfile.path.ens <- tempfile()

##Create new LIM model 
equations <- list(c("R_G6PDH2r", "R_GLUDC", "R_G3PD2m") , 
                  c("R_GLCt1r * 0.069", "R_GLCt1r * 0.32", "0"))
createLIMFromSBML(sbml.model, maximize, externals=externals, 
                  file.name=limfile.path.ens, equations=equations) 

##sample feasible flux distributions with MCMC
ensemble <- sampleFluxEnsemble(limfile.path.ens, uncertain.vars, 
                               x0=rates, iter=1e5, burninlength=1e4, 
                               outputlength=1e4, type="mirror", jmp=0.1)


###################################################
### code chunk number 10: plot_code (eval = FALSE)
###################################################
## par(mfrow=c(2,2))
## metab <- c(as.vector(uncertain.vars[1:2,1]), "R_SUCD1m")
## for (m in metab){
##   title <- paste(m, "\n(", sbml.model@reactions[[m]]@name, ")", sep="")
##   myhist <- hist(ensemble[,m], breaks=9, plot=FALSE)
##   plot(myhist, ylim=c(0, max(myhist$counts) + max(myhist$counts / 10)), 
##        xlab="flux (mmol/min)",main=title, col="cornflowerblue", cex.lab=1.3, 
##        xlim=c(min(myhist$breaks) - sd(myhist$breaks), 
##            max(myhist$breaks)+sd(myhist$breaks)))
##   text(mean(myhist$mids), max(myhist$counts) + max(myhist$counts / 10),
##        label=bquote(mu== ~.(round(mean(ensemble[,m]),3)) ~ 
##            "," ~ sigma== ~.(round(sd(ensemble[,m]),3))), cex=1.2)  
## }
## 
## ## get ensemble of net ATP production
## atp.prod.ens <- eval(parse(text=maximize), envir=data.frame(ensemble))
## 
## ##plot ensemble
## title <- paste("Net ATP production")
## myhist <- hist(atp.prod.ens, breaks=9, plot=FALSE)
## plot(myhist, ylim=c(0, max(myhist$counts) + 
##                  max(myhist$counts / 10)), xlab="flux (mmol/min)",
##      main=title, col="cornflowerblue", cex.lab=1.3, 
##      xlim=c(min(myhist$breaks) - sd(myhist$breaks), 
##          max(myhist$breaks)+sd(myhist$breaks)))
## text(mean(myhist$mids), max(myhist$counts) + max(myhist$counts / 10),
##      label=bquote(mu== ~.(round(mean(atp.prod.ens),3)) ~ 
##          "," ~ sigma== ~.(round(sd(atp.prod.ens),3))), cex=1.2)  


###################################################
### code chunk number 11: hists_fluxes
###################################################
## plot selected fluxes
par(mfrow=c(3,2))
metab <- c(as.vector(uncertain.vars[1:3,1]), "R_SUCD1m","R_ICDHxm")
for (m in metab){
  title <- paste(m, "\n(", sbml.model@reactions[[m]]@name, ")", sep="")
  myhist <- hist(ensemble[,m], breaks=9, plot=FALSE)
  plot(myhist, ylim=c(0, max(myhist$counts) + max(myhist$counts / 10)), 
       xlab="flux (mmol/min)",main=title, col="cornflowerblue", cex.lab=1.3, 
       xlim=c(min(myhist$breaks) - sd(myhist$breaks), 
           max(myhist$breaks)+sd(myhist$breaks)))
  text(mean(myhist$mids), max(myhist$counts) + max(myhist$counts / 10),
       label=bquote(mu== ~.(round(mean(ensemble[,m]),3)) ~ 
           "," ~ sigma== ~.(round(sd(ensemble[,m]),3))), cex=1.2)  
}

## get ensemble of net ATP production
atp.prod.ens <- eval(parse(text=maximize), envir=data.frame(ensemble))

##plot ensemble
title <- paste("Net ATP production")
myhist <- hist(atp.prod.ens, breaks=9, plot=FALSE)
plot(myhist, ylim=c(0, max(myhist$counts) + 
                 max(myhist$counts / 10)), xlab="flux (mmol/min)",
     main=title, col="cornflowerblue", cex.lab=1.3, 
     xlim=c(min(myhist$breaks) - sd(myhist$breaks), 
         max(myhist$breaks)+sd(myhist$breaks)))
text(mean(myhist$mids), max(myhist$counts) + max(myhist$counts / 10),
     label=bquote(mu== ~.(round(mean(atp.prod.ens),3)) ~ 
         "," ~ sigma== ~.(round(sd(atp.prod.ens),3))), cex=1.2)  



###################################################
### code chunk number 12: k
###################################################
relevant.species <- c("M_glc_DASH_D_c", "M_g6p_c", "M_f6p_c",
                      "M_fdp_c", "M_dhap_c", "M_g3p_c",
                      "M_13dpg_c", "M_3pg_c", "M_2pg_c",
                      "M_pep_c", "M_pyr_c",
                      "M_6pgl_c", "M_6pgc_c", "M_ru5p_DASH_D_c", 
                      "M_xu5p_DASH_D_c", "M_r5p_c", "M_g3p_c", "M_s7p_c")

relevant.reactions <- c("R_HEX1", "R_PGI", "R_PFK", "R_FBA", "R_TPI", 
                        "R_GAPD", "R_PGK", "R_PGM", "R_ENO", "R_PYK",
                        "R_G6PDH2r", "R_PGL", "R_GND", "R_RPE", "R_RPI", "R_TKT1")

hd <- sbml2hyperdraw(sbml.model, rates=rates, 
                     relevant.species=relevant.species, 
                     relevant.reactions=relevant.reactions,
                     layoutType="dot", plt.margins=c(20, 0, 20, 80))


###################################################
### code chunk number 13: l (eval = FALSE)
###################################################
## plot(hd)


###################################################
### code chunk number 14: fig_glyc
###################################################
plot(hd)


###################################################
### code chunk number 15: m (eval = FALSE)
###################################################
## relevant.species <- c("M_cit_m", "M_icit_m" , "M_akg_m",
##                       "M_succoa_m", "M_succ_m", "M_fum_m",
##                       "M_mal_DASH_L_m", "M_oaa_m")
## relevant.reactions <- c("R_CSm", "R_ACONTm", "R_ICDHxm",
##                         "R_AKGDm", "R_SUCOAS1m", "R_SUCD1m",
##                         "R_FUMm", "R_MDHm", "R_ICDHyrm", "R_ME1m",
##                         "R_ME2m", "R_ASPTAm","R_AKGMALtm", "R_GLUDym",
##                         "R_ABTArm", "R_SSALxm","R_CITtam")
## hd <- sbml2hyperdraw(sbml.model, rates=rates,
##                      relevant.reactions=relevant.reactions,
##                      relevant.species=relevant.species,
##                      layoutType="circo", plt.margins=c(150, 235, 150, 230))
## dev.new() ##Open a new plotting device
## plot(hd)


###################################################
### code chunk number 16: fig_tca
###################################################
relevant.species <- c("M_cit_m", "M_icit_m" , "M_akg_m",
                      "M_succoa_m", "M_succ_m", "M_fum_m",
                      "M_mal_DASH_L_m", "M_oaa_m")
relevant.reactions <- c("R_CSm", "R_ACONTm", "R_ICDHxm",
                        "R_AKGDm", "R_SUCOAS1m", "R_SUCD1m",
                        "R_FUMm", "R_MDHm", "R_ICDHyrm", "R_ME1m",
                        "R_ME2m", "R_ASPTAm","R_AKGMALtm", "R_GLUDym",
                        "R_ABTArm", "R_SSALxm","R_CITtam")
hd <- sbml2hyperdraw(sbml.model, rates=rates,
                     relevant.reactions=relevant.reactions,
                     relevant.species=relevant.species,
                     layoutType="circo", plt.margins=c(150, 235, 150, 230))
plot(hd)


