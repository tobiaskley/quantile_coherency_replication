###################################################################
## Code to replicate the results in the Barunik and Kley (2018)
###################################################################
##
## script is to install all required packages.
##
## Authors: Jozef Barunik, and Tobias Kley (2018)
###################################################################

install.packages(c("alphahull",
                   "astsa",
                   "copula",
                   "latex2exp",
                   "mAr",
                   "pbivnorm",
                   "quantreg",
                   "quantspec",
                   "rugarch",
                   "vars",
                   "zoo"),
               repos = "http://cran.us.r-project.org",
               dependencies = TRUE)

# package on which quantilogram depends
install.packages("np",
               repos = "http://cran.us.r-project.org",
               dependencies = TRUE)
           
install.packages("quantilogram_0.1.tar.gz",
                 dependencies = TRUE,
                 repos = NULL,
                 type = "source")
