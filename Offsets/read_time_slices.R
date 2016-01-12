source('~/Documents/R Codes/Offsets/offsets_functions.R')
params = readRDS('~/R Data/params.rds')
parcels = readRDS('~/R Data/parcels.rds')
counterfactuals = readRDS('~/R Data/counterfactuals.rds')

mean_trajectories = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))

trial_num = 1000

for (t in 1:params$time_steps){
  mean_time_slice = array(0, c(params$ecology_size, params$ecology_size))
  for (s in 1:trial_num){  
    filename = paste('~/R Data/test_', s, '_time_slice_', t, '.rds', sep = '')    
    time_slice <- readRDS(filename)
    mean_time_slice <- mean_time_slice + time_slice
  }
  mean_trajectories[, , t] = mean_time_slice/trial_num
  print(t)
}

counterfactual_sums <- find_sums(params$time_steps, parcels, counterfactuals$elements)
mean_trajectory_sums <- find_sums(params$time_steps, parcels, mean_trajectories)

plot_params = c(params$region_num_x, params$region_num_y, parcels$region_num)
plot_outs(mean_trajectory_sums$region_sums, counterfactual_sums$region_sums, plot_params)