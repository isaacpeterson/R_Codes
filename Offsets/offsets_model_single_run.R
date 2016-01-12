rm(list = ls())
source('~/Documents/R Codes/Offsets/offsets_functions_parceled.R')
source('~/Documents/R Codes/Offsets/Params_regional.R')

graphics.off()
global_params <- initialise_global_params()
development_vector <- split_vector(global_params$time_steps, 60, sd = 1, min_width = -1)
region_params <- initialise_region_params(global_params$region_num, development_vector)
parcels <- initialise_shape_parcels(global_params)
index_object <- initialise_index_object(parcels, global_params)
initial_ecology <- initialise_ecology(global_params, parcels)
decline_rates_initial <- build_decline_rates(parcels, region_params)
outputs <- calc_trajectories_by_parcel(global_params, region_params, current_ecology = initial_ecology, decline_rates = decline_rates_initial, 
                                         parcels, index_object, perform_offsets = TRUE, record_parcel_sets = TRUE)

#counterfactuals <- initialise_counterfactuals_vectorised(global_params, region_params, parcels$land_parcels, initial_ecology, decline_rates = initial_decline_rates)

assess_object <- assess_parcel_sets(outputs$parcel_sets, parcels$land_parcels, outputs$trajectories, global_params$time_steps, decline_rates_initial)

plot_parcel_sets(assess_object, assess_type = 'set', parcel_set_num = 10)

#plot_net_trajectories(assess_object)