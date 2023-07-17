# reference:  https://doi.org/10.5061/dryad.3bk3j9khj
# original publication: Pioneering polyploids: the impact of whole-genome duplication on biome shifting in New Zealand Coprosma (Rubiaceae) and Veronica (Plantaginaceae)
# by L. G. Liddell, W. G. Lee, E. E. Dale, H. M. Meudt and N. J. Matzke
# Biology Letters 2021 Vol. 17 Issue 9 Pages 20210297
# DOI: doi:10.1098/rsbl.2021.0297

#Remove all data

rm(list=ls())

# Load the package (after installation, see above).

library(GenSA)  
library(FD)      
library(snow)     
library(parallel)
library(optimx)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

# Set the work directory

setwd("...")

# read the tree

trfn = "#.newick"
moref(trfn)
tr = read.tree(trfn)
tr

# read the distribution data

geogfn = "#.txt"
moref(geogfn)

# Look at your geographic range data:

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size, set the maximum range size (e.g., = 2 here)

max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 2


#######################################################
# Traits-only model -- 1 rate
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = 2
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

geog_values = getranges_from_LagrangePHYLIP("#.txt")
geog_values

trait_fn = "#.txt"
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)

# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0

# No multipliers on geog (set m1 and m2 to 1)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object$on_NaN_error = -1000000

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "sim_traitsOnly_1rate_v1_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resTrait_1rate = res
} else {
  # Loads to "res"
  load(resfn)
  resTrait_1rate = res
} # END if (runslow)

#######################################################
# Traits-only model -- 2 rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = 2
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

# Set up DEC model, but set all rates to 0 (data are 1 invariant area)
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["a","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = 0.0

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999


tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()

trait_fn = "#.txt"
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = 0.001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0

# No multipliers on geog (set m1 and m2 to 1)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = 1.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "desc"] = "trait-based dispersal rate multipliers m1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "desc"] = "trait-based dispersal rate multipliers m2"

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object$on_NaN_error = -1000000

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.

resfn = "sim_traitsOnly_2rates_v1_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resTrait_2rates = res
} else {
  # Loads to "res"
  load(resfn)
  resTrait_2rates = res
} # END if (runslow)



#######################################################
# Run DEC (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

#tr = read.tree(BioGeoBEARS_run_object$trfn)
BioGeoBEARS_run_object$timesfn = "#.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "#.txt"

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


resfn = "DEC_inf_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}


#######################################################
# Run DEC+J (on geography only)
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

#tr = read.tree(BioGeoBEARS_run_object$trfn)
BioGeoBEARS_run_object$timesfn = "#.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "#.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)


dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
xstart = resDEC$outputs@params_table["x","est"]

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

print("Printing warnings: 'warnings()':")
print(warnings())


resfn = "DECj_inf_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

#######################################################
# Run DEC + t12 + t21 + m2, starting from DEC-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()
BioGeoBEARS_run_object$timesfn = "#.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "#.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)


trait_fn = "#.txt"
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)


# Look at the params table
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table


# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resDEC$outputs@params_table["d","est"]
estart = max(c(resDEC$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001
xstart = resDEC$outputs@params_table["x","est"]


BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 100

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DEC+t12+t21+m2_inf_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC_t12_t21_m2 = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC_t12_t21_m2 = res
}


#######################################################
# Run DECj + t12 + t21 + m2, starting from DECj-geog and 2rates
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000

tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()
BioGeoBEARS_run_object$timesfn = "#.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "#.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)


trait_fn = "#.txt"
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resTrait_2rates$outputs@params_table["t12","est"]
t21_start = resTrait_2rates$outputs@params_table["t21","est"]
m2_start = 1
dstart = resDECj$outputs@params_table["d","est"]
estart = max(c(resDECj$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = resDECj$outputs@params_table["j","est"]
xstart = resDECj$outputs@params_table["x","est"]


# Set up DEC+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Crash fix
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 100

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DECJ+t12+t21+m2_inf_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj_t12_t21_m2 = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj_t12_t21_m2 = res
}

#######################################################
# Run DECj + t12 + t21 + m2, starting from DEC + t12 + t21 + m2
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$print_optim = TRUE
BioGeoBEARS_run_object$calc_ancprobs=TRUE        # get ancestral states from optim run
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$num_cores_to_use=8
BioGeoBEARS_run_object$use_optimx="GenSA"
BioGeoBEARS_run_object$speedup=TRUE
BioGeoBEARS_run_object$geogfn = "#.txt"
BioGeoBEARS_run_object$trfn = "#.newick"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE
BioGeoBEARS_run_object$on_NaN_error = -1000000



tr = read.tree(BioGeoBEARS_run_object$trfn)
#plot(tr); axisPhylo()
BioGeoBEARS_run_object$timesfn = "#.txt"
BioGeoBEARS_run_object$areas_allowed_fn = "#.txt"
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$force_sparse = FALSE  
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)


trait_fn = "#.txt"
trait_values = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn)
trait_values

# Add the traits data and model
BioGeoBEARS_run_object = add_trait_to_BioGeoBEARS_run_object(BioGeoBEARS_run_object, trait_fn=trait_fn)



# Starting values from ML results of simpler run
t12_start = resDEC_t12_t21_m2$outputs@params_table["t12","est"]
t21_start = resDEC_t12_t21_m2$outputs@params_table["t21","est"]
m2_start = resDEC_t12_t21_m2$outputs@params_table["m2","est"]
dstart = resDEC_t12_t21_m2$outputs@params_table["d","est"]
estart = max(c(resDEC_t12_t21_m2$outputs@params_table["e","est"], 1.1*BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"]))
jstart = 0.0001
xstart = resDEC_t12_t21_m2$outputs@params_table["x","est"]

# Set up DEC+J model
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Crash fix
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 1e-13
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 1e-13

#######################################################
# Manual modifications of trait-based model
#######################################################
# Edit t12 and t21 rates
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "init"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "est"] = t12_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t12", "max"] = 0.1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["t21", "max"] = 0


# Set 0/1 multipliers on dispersal rate
# For flightlessness (m2), max multiplier is 1, and
# fix to a small value, or estimate
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "init"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "est"] = 1
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "min"] = 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m1", "max"] = 1

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "init"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "est"] = m2_start
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2", "max"] = 100

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "DECJ+t12+t21+m2_rep2_inf_PL.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj_t12_t21_m2_rep2 = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj_t12_t21_m2_rep2 = res
}



param_names = c("lnL", "d", "e", "j", "t12", "t21", "m1", "m2")

Trait_1rate_results = c(
  resTrait_1rate$total_loglikelihood,
  resTrait_1rate$output@params_table["d", "est"], 
  resTrait_1rate$output@params_table["e", "est"], 
  resTrait_1rate$output@params_table["j", "est"], 
  resTrait_1rate$output@params_table["t12", "est"], 
  resTrait_1rate$output@params_table["t21", "est"], 
  resTrait_1rate$output@params_table["m1", "est"], 
  resTrait_1rate$output@params_table["m2", "est"]
)
names(Trait_1rate_results) = paste("Trait_1rate_", param_names, sep="")

Trait_2rates_results = c(
  resTrait_2rates$total_loglikelihood,
  resTrait_2rates$output@params_table["d", "est"], 
  resTrait_2rates$output@params_table["e", "est"], 
  resTrait_2rates$output@params_table["j", "est"], 
  resTrait_2rates$output@params_table["t12", "est"], 
  resTrait_2rates$output@params_table["t21", "est"], 
  resTrait_2rates$output@params_table["m1", "est"], 
  resTrait_2rates$output@params_table["m2", "est"]
)
names(Trait_2rates_results) = paste("Trait_2rates_", param_names, sep="")

DEC_results = c(
  resDEC$total_loglikelihood,
  resDEC$output@params_table["d", "est"], 
  resDEC$output@params_table["e", "est"], 
  resDEC$output@params_table["j", "est"], 
  resDEC$output@params_table["t12", "est"], 
  resDEC$output@params_table["t21", "est"], 
  resDEC$output@params_table["m1", "est"], 
  resDEC$output@params_table["m2", "est"]
)
names(DEC_results) = paste("DEC_", param_names, sep="")


DECj_results = c(
  resDECj$total_loglikelihood,
  resDECj$output@params_table["d", "est"], 
  resDECj$output@params_table["e", "est"], 
  resDECj$output@params_table["j", "est"], 
  resDECj$output@params_table["t12", "est"], 
  resDECj$output@params_table["t21", "est"], 
  resDECj$output@params_table["m1", "est"], 
  resDECj$output@params_table["m2", "est"]
)
names(DECj_results) = paste("DECj_", param_names, sep="")

DEC_t12_t21_m2_results = c(
  resDEC_t12_t21_m2$total_loglikelihood,
  resDEC_t12_t21_m2$output@params_table["d", "est"], 
  resDEC_t12_t21_m2$output@params_table["e", "est"], 
  resDEC_t12_t21_m2$output@params_table["j", "est"], 
  resDEC_t12_t21_m2$output@params_table["t12", "est"], 
  resDEC_t12_t21_m2$output@params_table["t21", "est"], 
  resDEC_t12_t21_m2$output@params_table["m1", "est"], 
  resDEC_t12_t21_m2$output@params_table["m2", "est"]
)
names(DEC_t12_t21_m2_results) = paste("DEC_t12_t21_m2_", param_names, sep="")

DECj_t12_t21_m2_results = c(
  resDECj_t12_t21_m2$total_loglikelihood,
  resDECj_t12_t21_m2$output@params_table["d", "est"], 
  resDECj_t12_t21_m2$output@params_table["e", "est"], 
  resDECj_t12_t21_m2$output@params_table["j", "est"], 
  resDECj_t12_t21_m2$output@params_table["t12", "est"], 
  resDECj_t12_t21_m2$output@params_table["t21", "est"], 
  resDECj_t12_t21_m2$output@params_table["m1", "est"], 
  resDECj_t12_t21_m2$output@params_table["m2", "est"]
)
names(DECj_t12_t21_m2_results) = paste("DECj_t12_t21_m2_", param_names, sep="")


DECj_t12_t21_m2_rep2_results = c(
  resDECj_t12_t21_m2_rep2$total_loglikelihood,
  resDECj_t12_t21_m2_rep2$output@params_table["d", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["e", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["j", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["t12", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["t21", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["m1", "est"], 
  resDECj_t12_t21_m2_rep2$output@params_table["m2", "est"]
)
names(DECj_t12_t21_m2_rep2_results) = paste("DECj_t12_t21_m2_rep2_", param_names, sep="")


tmp_results = c(Trait_1rate_results, Trait_2rates_results, DEC_results, DECj_results, DEC_t12_t21_m2_results, DECj_t12_t21_m2_results, DECj_t12_t21_m2_rep2_results)
tmp_results_mat = as.matrix(tmp_results, nrow=7, byrow=TRUE)
tmp_results_mat = as.data.frame(tmp_results_mat, stringsAsFactors=FALSE)

outfn = slashslash("dec_params_inferred_PL.txt")
write.table(x=tmp_results, file=outfn, append=FALSE, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
moref(outfn)


#######################################################
# Likelihood profile of m2
#######################################################
# Vary m2 and get the log-likelihood, for approximate
# uncertainty
m2_vals = c(0.01, 0.1, 0.5, 0.75, seq(1,50,0.1))

lnLs = NULL
for (i in 1:length(m2_vals))
{
  m2_val = m2_vals[i]
  
  BioGeoBEARS_run_object = resDEC_t12_t21_m2$inputs
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","max"] = 100.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","init"] = m2_val
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","est"] = m2_val
  lnL = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_loglike")
  txt = paste0(i, ": When m2=", m2_val, ", lnL=", lnL)
  cat("\n")
  cat("\n")
  cat(txt)
  cat("\n")
  cat("\n")
  lnLs = c(lnLs, lnL)
}

ivals = 1:length(m2_vals)
lnL_table = cbind(ivals, m2_vals, lnLs)
lnL_table = adf2(lnL_table)
write.table(lnL_table, file="LINEAGES_m2_DEC_vs_lnL_table_v1_PL.txt", sep="\t", quote=FALSE)

# Confidence interval
likelihood = exp(resDEC_t12_t21_m2$total_loglikelihood)
lower95_likelihood = (1/20) * likelihood
lower95_lnL = log(lower95_likelihood)

pdffn = "LINEAGES_m2_DEC_vs_lnL_v2.pdf"
pdf(file=pdffn, width=9, height=12)

par(mfrow=c(2,1))

plot(m2_vals, lnLs, main="Profile log-likelihood of m2", xlab="m2", ylab="lnL")
lines(m2_vals, lnLs)
abline(h=lower95_lnL, lty="dashed")

plot(m2_vals, exp(lnLs), main="Profile likelihood of m2", xlab="m2", ylab="likelihood")
lines(m2_vals, exp(lnLs))

abline(h=lower95_likelihood, lty="dashed")
dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

#make table with m2 values with less than <5% removed
DEC_95_table = lnL_table[lnL_table$lnLs > lower95_lnL, ]
DEC_95_table = subset(DEC_95_table, select = -c(ivals))
DEC_95_table$likelihoods = exp(DEC_95_table$lnLs)
DEC_95_table$model = "DEC"
DEC_95_table$slikelihoods = DEC_95_table$likelihoods-lower95_likelihood
DEC_95_table$AICc = 0.2735
write.table(DEC_95_table, file = "DEC_95_table_PL.txt", sep="\t", quote=FALSE)


# 95% CI
m2_ML_est = resDEC_t12_t21_m2$outputs@params_table["m2","est"]

TF = lnLs > lower95_lnL
lower95_m2 = min(m2_vals[TF])
upper95_m2 = max(m2_vals[TF])

m2_ML_est
lower95_m2
upper95_m2

#######################################################
# Likelihood profile of m2+j
#######################################################
# Vary m2 and get the log-likelihood, for approximate
# uncertainty
m2_vals = c(0.01, 0.1, 0.5, 0.75, seq(1,50,0.1))

lnLs = NULL
for (i in 1:length(m2_vals))
{
  m2_val = m2_vals[i]
  
  BioGeoBEARS_run_object = resDECj_t12_t21_m2_rep2$inputs
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","max"] = 100.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","init"] = m2_val
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["m2","est"] = m2_val
  lnL = bears_optim_run(BioGeoBEARS_run_object, skip_optim=TRUE, skip_optim_option="return_loglike")
  txt = paste0(i, ": When m2=", m2_val, ", lnL=", lnL)
  cat("\n")
  cat("\n")
  cat(txt)
  cat("\n")
  cat("\n")
  lnLs = c(lnLs, lnL)
}

ivals = 1:length(m2_vals)
lnL_table = cbind(ivals, m2_vals, lnLs)
lnL_table = adf2(lnL_table)
write.table(lnL_table, file="LINEAGES_m2_DECj_vs_lnL_table_v1_PL.txt", sep="\t", quote=FALSE)

# Confidence interval
likelihood = exp(resDECj_t12_t21_m2_rep2$total_loglikelihood)
lower95_likelihood = (1/20) * likelihood
lower95_lnL = log(lower95_likelihood)

pdffn = "LINEAGES_m2_DECj_vs_lnL_v2_PL.pdf"
pdf(file=pdffn, width=9, height=12)

par(mfrow=c(2,1))

plot(m2_vals, lnLs, main="Profile log-likelihood of m2", xlab="m2", ylab="lnL")
lines(m2_vals, lnLs)
abline(h=lower95_lnL, lty="dashed")

plot(m2_vals, exp(lnLs), main="Profile likelihood of m2", xlab="m2", ylab="likelihood")
lines(m2_vals, exp(lnLs))

abline(h=lower95_likelihood, lty="dashed")
dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)

#make table with m2 values with less than <5% removed
DECj_95_table = lnL_table[lnL_table$lnLs > lower95_lnL, ]
DECj_95_table = subset(DECj_95_table, select = -c(ivals))
DECj_95_table$likelihoods = exp(DECj_95_table$lnLs)
DECj_95_table$model = "DEC+J"
DECj_95_table$slikelihoods = DECj_95_table$likelihoods-lower95_likelihood
DECj_95_table$AICc = 52.1097
write.table(DECj_95_table, file = "DECj_95_table_PL.txt", sep="\t", quote=FALSE)


# 95% CI
m2_ML_est = resDECj_t12_t21_m2_rep2$outputs@params_table["m2","est"]

TF = lnLs > lower95_lnL
lower95_m2 = min(m2_vals[TF])
upper95_m2 = max(m2_vals[TF])

m2_ML_est
lower95_m2
upper95_m2

