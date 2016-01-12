rm(list = ls())
source('~/Documents/R Codes/Offsets/offsets_functions_parceled.R')
source('~/Documents/R Codes/Offsets/Params_regional.R')
# 
# graphics.off()

load_object <- readRDS('~/R_Data/load_object.rds')
global_params <- initialise_global_params()
region_params <- initialise_region_params(global_params$region_num, load_object$development_vector)

# development_vector <- split_vector(global_params$time_steps, 60, sd = 1, min_width = -1)
# parcels <- initialise_shape_parcels(global_params)
# index_object <- initialise_index_object(parcels, global_params)
# initial_ecology <- initialise_ecology(global_params, parcels)
# decline_rates_initial <- build_decline_rates(parcels, region_params)


 # load_object$development_vector = development_vector
# load_object$parcels = parcels
# load_object$index_object = index_object
# load_object$initial_ecology = initial_ecology
# load_object$decline_rates_initial = decline_rates_initial


outputs <- calc_trajectories_by_parcel(global_params, region_params, current_ecology = load_object$initial_ecology, decline_rates = load_object$decline_rates_initial, 
                                       load_object$parcels, load_object$index_object, perform_offsets = TRUE, record_parcel_sets = TRUE)

#counterfactuals <- initialise_counterfactuals_vectorised(global_params, region_params, parcels$land_parcels, initial_ecology, decline_rates = initial_decline_rates)

assess_object <- assess_parcel_sets(outputs$parcel_sets, load_object$parcels$land_parcels, outputs$trajectories, load_object$global_params$time_steps, load_object$decline_rates_initial)

assess_no_avoided_clearing <- readRDS('~/R_Data/assess_object.rds')
#saveRDS(assess_object, '~/R_Data/assess_object.rds')

#plot_parcel_sets(assess_object, assess_type = 'set', parcel_set_num = 10)

plot_net_trajectories(assess_object)