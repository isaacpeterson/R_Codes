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
counterfactuals_object <- initialise_counterfactuals_vectorised(global_params, region_params, parcels$land_parcels, initial_ecology, decline_rates$array)

outputs <- calc_trajectories_by_parcel(outputs, global_params, region_params, current_ecology = initial_ecology, decline_rates, 
                                       parcels, index_object, counterfactuals_object, perform_offsets = FALSE)

# .libPaths('~/Library/R/3.2/library')
# library(foreach)
# library(doParallel)
# cl<-makeCluster(4)
# registerDoParallel(cl)
# 
# strt<-Sys.time()
# run_num = 100
# traj_runs <- foreach(run_ind=1:run_num) %dopar%{
#   outputs <- calc_trajectories_by_parcel(outputs, global_params, region_params, current_ecology = initial_ecology, decline_rates$array, 
#                                          parcels, index_object, counterfactuals_object, perform_offsets = FALSE)
# }
# 
# stopCluster(cl)
# 
# p_reals = parcel_realisations(traj_runs, parcel_ind = 1, parcels$land_parcels, global_params$time_steps)
# plot(rowSums(p_reals$sums))
# overlay_plots(p_reals$sums)

#mean_traj <- Reduce('+', traj_runs)
#mean_traj = mean_traj/run_num

# print(Sys.time()-strt)

traj_sums <- find_sums(1:length(parcels$land_parcels), outputs$trajectories, parcels, global_params$time_steps)
counter_sums <- find_sums(1:length(parcels$land_parcels), counterfactuals_object$counterfactuals, parcels, global_params$time_steps)

# offset_sums <- find_sums(outputs$offset_list, outputs$trajectories, parcels, global_params$time_steps)
# dev_sums <- find_sums(outputs$development_list, outputs$trajectories, parcels, global_params$time_steps)
# 
# true_offset_sums <- true_sums(outputs$offset_list, offset_sums, counter_sums)
# true_dev_sums <- true_sums(outputs$development_list, dev_sums, counter_sums)
# a_object <- assess_parcel_sets(outputs$parcel_sets, parcels$land_parcels, outputs$trajectories, counterfactuals_object$counterfactuals, global_params$time_steps)
# 
#plot_outs(traj_sums$parcel_sums[, 1:6], traj_sums$parcel_sums[, 1:6], c(2, 3))
