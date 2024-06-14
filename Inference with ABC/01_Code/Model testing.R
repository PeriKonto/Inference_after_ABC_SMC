
# author: Periklis Kontoroupis
# date : 15/10/2018

# this is the script to generate data based on known transmission parameters which later on will be
# for model testing by the ABC-SMC to recover the parameters used

# Prep session
rm(list = ls())

library(data.table)
library(foreach)
library(doParallel)
library(XML)
library(ggplot2)
library(osDesign)


base.dir    <- "F:/to archive to X drive/2. ABC testing/Task 2.1/Inference with ABC"
code.dir    <- file.path(base.dir, "01_Code")
output.dir  <- file.path(base.dir, "02_Output")
source.dir  <- file.path(base.dir, "03_Source_code")
wormsim.dir <- file.path(base.dir, "04_wormsim-v2.58Ap25")

input.template <- file.path(source.dir, "template2MDA_all_6_10_18.xml")
input.schema   <- file.path(wormsim.dir, "wormsim.xsd")


# Load functions
setwd(source.dir)
source("xml_substitute_functions.r")  # Functions to adjust an xml file
source("create_xml_functions.r")      # Automated script to create xml file from a template
source("basic_functions.r")           # Functions to run WORMSIM and process output
source("par_grid_functions.r")        # Function to create a grid of parameter values
source("par_def.R")                   # Default values for subset of parameters


# Default parameters (if different from those specified in "par_def.R)

par.def <- within(par.def, {
  surveys <- list(start = list(year = 2000,
                               month = 0),
                  stop = list(year = 2010,
                              month = 1),  # make sure that at least two surveys take place for code to work
                  interval = list(year = 1,
                                  month = 0))})
  # mda <- list(nrounds = 5,  # Cannot exceed the number of treatments in template
  #             years = 2011:2013,
  #             months = c(0, 6) ) })
#------------------------------------------------------------------------------
setwd(wormsim.dir)
generate.inputfile(template = input.template, schema = input.schema,
                   xml.name = "input_default", pars = par.def)
#------------------------------------------------------------------------------
mbr         <- 600
exposure.p1 <- 0.4
foi         <- 0

#-----------------------------
df10 <- foi * 1#Foi$foi.DA.hyrl.80[1]
df11 <- foi * 0#Foi$foi.DA.hyrl.80[2]
df12 <- foi * 0#Foi$foi.DA.hyrl.80[3]
df13 <- foi * 0#Foi$foi.DA.hyrl.80[4]
#====================================
simnumber <- 1
# sample from a uniform distribution
seed <- 45345345;
set.seed(seed);


#====================================#====================================#====================================#====================================
par.alt <- list(mbr=mbr, exposure.p1=exposure.p1, foi=foi,
                df10=df10, df11=df11, df12=df12,  df13=df13)
#====================================#====================================#====================================#====================================

param <- function(parameter.alt) {

  n.alt.par.values <- sapply(parameter.alt, length)
  n.alt.par.comb <- n.alt.par.values[[1]]

  # Initialize list of combinations of parameter values and determine indices
  # needed to pull values from parameter.alt.
  index.alt.par <- data.frame(lapply(parameter.alt, function(x) {1:length(x)}))
  par.alt.combi.list <- as.list(rep(list(as.list(rep(NA, dim(index.alt.par)[2]))), n.alt.par.comb))

  # Fill list with unique combinations of parameter values
  for (i in 1:dim(index.alt.par)[1]) {
    for (j in 1:dim(index.alt.par)[2]) {

      par.alt.combi.list[[i]][[j]] <- parameter.alt[[j]][[index.alt.par[i,j]]]

    }

    names(par.alt.combi.list[[i]]) <- names(parameter.alt)

  }

  # Return list of lists of alternative parameter values
  return(par.alt.combi.list)

}

par.alt.combi.list <- param(par.alt)

#====================================#====================================#====================================#====================================

# Load functions
setwd(source.dir)
#---------------------------------------------------------------#---------------------------------------------------------------#---------------------------------------------------------------
start.time <- Sys.time()
# Run simulations and save in a list object ("output") with each element of the
# representing a parameter set
cluster <- makeCluster(16)

registerDoParallel(cluster)

output <- foreach(i = 1:length(par.alt.combi.list),
                  .inorder = TRUE,
                  .errorhandling = "remove",
                  .packages = c("XML", "data.table")) %dopar% {

                    # Set paths and load functions (re-execute for parallel sessions)
                    setwd(source.dir)
                    source("xml_substitute_functions.r")
                    source("create_xml_functions.r")
                    source("basic_functions.r")
                    source("par_grid_functions.r")


                    seed.start  <- 1
                    seed.end    <- 1000

                    # Create folder to work in
                    setwd(wormsim.dir)
                    input.file.name <- paste("par_set_", i, sep = "")
                    dir.create(file.path(wormsim.dir, input.file.name),
                               showWarnings = FALSE)
                    file.copy(dir(), file.path(wormsim.dir, input.file.name),
                              overwrite = TRUE)
                    setwd(file.path(wormsim.dir, input.file.name))

                    # Run model
                    generate.inputfile(template = "input_default.xml",
                                       schema = input.schema,
                                       xml.name = input.file.name,
                                       pars = par.alt.combi.list[[i]])

                    run.proc.sim(
                      input.file = input.file.name,
                      seed.start = seed.start,
                      seed.end = seed.end,
                      delete.txt = FALSE
                    )

                    # Read individual simulation files and remove failed simulations
                    summ <- read.output.ind(input.file.name,
                                            seed.start = seed.start,
                                            seed.end = seed.end,
                                            type = "")

                    age <- read.output.ind(input.file.name,
                                           seed.start = seed.start,
                                           seed.end = seed.end,
                                           type = "X")

                    intens <- read.output.ind(input.file.name,
                                             seed.start = seed.start,
                                             seed.end = seed.end,
                                             type = "Y")

                    fail <- summ[1,4,] == 0
                    summ <- data.table(apply(summ[,, !fail], 1:2, mean))
                    age <- data.table(apply(age[,, !fail], 1:2, mean))
                    intens <- data.table(apply(intens[,, !fail], 1:2, mean))

                    # Clean up
                    setwd(wormsim.dir)
                    unlink(input.file.name, recursive = TRUE)

                    # Return result
                    with(par.alt.combi.list[[i]],
                         list(par_set = i,
                              summ = summ,
                              age = age,
                              intens = intens))
                  }

stopCluster(cluster)

# Clean up WORMSIM folder
unlink(file.path(wormsim.dir, "input_default.xml"))

# Save output
setwd(output.dir)
save.image(file = file.path(output.dir, "Pruda_Wb_based_data_21018.RData"))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
#---------------------------------------------------------------#---------------------------------------------------------------#---------------------------------------------------------------
############################################################################################

par.alt.combi.table2 <- rbindlist(par.alt.combi.list)
max_value=dim(par.alt.combi.table2)
value=seq(from=1,  to = max_value[1], 1)
par.alt.combi.table2[ , id := value]

data <- output[[1]]$summ
y <- data[year =="2010"]



# collect.distances <- function(y) {
#   x <- y$par_set
#   y <- y$summ[year >"2009" & year < "2021"]
#   data = y[, list(year, mfPr, aNmf)]
#   zz <- dim(data)
#   data[, id:=  rep.int(x,  zz[1] ) ]
#   data.table(data)[, .(year, mfPr, id, aNmf)]
# }
# model.pred =  rbindlist(lapply(output, collect.distances))
# worm.data4list_x2MDA <- merge(x = model.pred, y = par.alt.combi.table2, by = "id", incomparables = NA)
# names(worm.data4list_x2MDA)[names(worm.data4list_x2MDA) == "exposure.p1"] = "k"

dataX <- output[[1]]$age


  y1 <- dataX[year=="2010" & age > 5]
  ydata = y1[, list(year, age,  M,   F)]
  ydata [, N := round((M+F))  ]
  ydata[, age := factor(x = age,
                       levels  =c(10, 20, 30, 40, 50, 99),
                       labels = c("05-9", "10-19","20-29","30-39","40-49","50+"),
                       ordered = TRUE)]

  ydata <- ydata[, c("age", "N")]
  write.csv(ydata, file = "ydata.csv")

data <- output[[1]]$intens
names(data)[3:6] <- paste("M", names(data)[3:6], sep = "")
names(data)[7:10] <- paste("F", names(data)[7:10], sep = "")
data[, None := get("M-0") + get("F-0")]
data[, Light := get("M-2") + get("F-2")]
data[, Moderate := get("M-16") + get("F-16")]
data[, Heavy := get("M16+") + get("F16+")]
data[, N := None + Light + Moderate + Heavy]
data1 <- melt(data = data,
          measure.vars = c("None", "Light", "Moderate", "Heavy"),
          variable = "intens",
          value.name = "cases")


  data1  <- data.table(subset(data1, year==2010))

  profile1 <- aggregate(cases ~ intens, data1, mean)
  profile1$village <- "TestVillage"
  profile1$year <- 2011
  profile1$N    <-  round(sum(profile1$cases))
  profile1$prev <- profile1$cases / sum(profile1$cases)

  syn.obs <-  profile1[, c("village", "year", "intens", "cases", "N", "prev") ]
  syn.obs$cases <-  round(syn.obs$cases)
  write.csv(syn.obs, file = "syn.obs.csv")
  # Village Year   intens cases    N       prev
  # 1:   Pruda 2011    Heavy    11 1027 0.01071081
  # 2:   Pruda 2011 Moderate    40 1027 0.03894839
  # 3:   Pruda 2011    Light    38 1027 0.03700097
  # 4:   Pruda 2011     None   938 1027 0.91333982

