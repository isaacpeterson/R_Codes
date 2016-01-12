rm(list = ls())
source('~/Documents/R Codes/Offsets/offsets_functions_parceled.R')
source('~/Documents/R Codes/Offsets/Params_regional.R')

graphics.off()
global_params <- initialise_global_params()
region_params <- initialise_region_params(global_params$region_num)
parcels <- initialise_shape_parcels(global_params)
index_object <- initialise_index_object(parcels, global_params)
initial_ecology <- initialise_ecology(global_params, parcels)
decline_rates <- build_decline_rates(parcels, region_params)
counterfactuals_object <- initialise_counterfactuals_vectorised(global_params, region_params, parcels, initial_ecology, decline_rates$decline_rates_array)
outputs <- initialise_outputs(counterfactuals_object, parcels, global_params)

outputs <- calc_trajectories_by_parcel(outputs, global_params, region_params, current_ecology = initial_ecology, decline_rates$decline_rates_array, 
                                       parcels, index_object, counterfactuals_object, perform_offsets = FALSE)
                                    

traj_sums <- find_sums(1:length(parcels$land_parcels), outputs$trajectories, parcels, global_params$time_steps)
counter_sums <- find_sums(1:length(parcels$land_parcels), counterfactuals_object$counterfactuals, parcels, global_params$time_steps)
#offset_sums <- find_sums(outputs$offset_list, outputs$trajectories, parcels, global_params$time_steps)
dev_sums <- find_sums(outputs$development_list, outputs$trajectories, parcels, global_params$time_steps)

#true_offset_sums <- true_sums(outputs$offset_list, offset_sums, counter_sums)

#a_object <- assess_parcel_sets(outputs$parcel_sets, parcels$land_parcels, outputs$trajectories, counterfactuals_object$counterfactuals, global_params$time_steps)
# plot_outs(counter_sums$region_sums, traj_sums$region_sums, c(1, 2))
# 
# plot_outs(a_object$dev_true_sum[, 1:6], a_object$dev_traj_rel_initial_sum[, 1:6], a_object$dev_initial_rel_counter_sum[, 1:6], c(2, 3))
# plot_outs(a_object$offset_true_sum[, 1:6], a_object$offset_traj_rel_initial_sum[, 1:6], a_object$offset_initial_rel_counter_sum[, 1:6], c(2, 3))


