split_vector <- function(N, M, sd) {
  min_width = 3;
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  while (any(vec <= min_width)) {
    negs <- vec <= min_width
    pos  <- vec > min_width
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}


initialise_shape_parcels <- function(params){
  parcels = list()
  parcel_num_x = params$parcel_num_x
  parcel_num_y = params$parcel_num_y
  parcel_vx = split_vector(parcel_num_x, params$ecology_size, 5)
  parcel_vy = split_vector(parcel_num_y, params$ecology_size, 5)
  
  pixel_indexes = 1:(params$ecology_size*params$ecology_size)
  dim(pixel_indexes) = c(params$ecology_size, params$ecology_size)
  land_parcels = mcell(pixel_indexes, parcel_vx, parcel_vy) #split the ecology array into a series of subarrays with dimensions sz_x by sz_y
  land_parcel_num = length(land_parcels$elements) #total number of parcels
  parcel_indexes = 1:land_parcel_num #index all parcels
  dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
  region_vx = split_vector(params$region_num_x, parcel_num_x, 1) 
  region_vy = split_vector(params$region_num_y, parcel_num_y, 1)
  
  regions = mcell(parcel_indexes, region_vx, region_vy)
  
  region_num = length(regions$elements)
  parcels$land_parcel_num = land_parcel_num
  parcels$land_parcels = land_parcels$elements
  parcels$land_parcel_dims = land_parcels$dims
  parcels$regions = regions$elements
  parcels$region_dims = regions$dims
  parcels$region_num = region_num
  parcels$parcel_vx = parcel_vx
  parcels$parcel_vy = parcel_vy
  
  return(parcels)
}


initialise_index_object <- function(parcels){
  index_object = list()
  index_object$ind_available = parcels$regions
  index_object$break_flag = FALSE
  return(index_object)
}

initialise_ecology <- function(params, parcels){
  land_parcels = parcels$land_parcels
  land_parcel_num = parcels$land_parcel_num
  ecology_initial = matrix(1,params$ecology_size,params$ecology_size)
  
  for (s in 1:land_parcel_num){
    initial_parcel_value = params$min_initial_eco_val + (params$max_initial_eco_val - params$min_initial_eco_val - params$initial_eco_noise)*runif(1) 
    inds = parcels$land_parcels[[s]]
    dim(inds) = c(1, length(inds))
    ecology_initial[inds] = ecology_initial[inds]*initial_parcel_value
  }
  ecology_initial = ecology_initial + params$initial_eco_noise*matrix(runif(params$ecology_size*params$ecology_size),params$ecology_size,params$ecology_size)
  return(ecology_initial)
}

initialise_offset_object <- function(){
  object = list()
  object$current_parcel = list()
  object$current_parcel_vals = list()
  return(object)
}


mcell <- function(x, vx, vy){
  
  rowsizes = vy;
  colsizes = vx;
  rows = length(rowsizes);
  cols = length(colsizes);
  
  a = 1
  B = vector('list', rows*cols)
  colStart = 0
  for (i in 1:cols){
    rowStart = 0
    for (j in 1:rows){
      B[[a]] = x[rowStart+(1:rowsizes[j]), colStart+(1:colsizes[i])]
      rowStart = rowStart + rowsizes[j]
      a = a + 1
    }
    colStart = colStart + colsizes[i]
  }
  
  parcel = list()
  parcel$dims = c(length(vy), length(vx))
  parcel$elements = B
  return(parcel)
  
}  

ind2sub <- function(rows, ind){
  rw = ((ind-1) %% rows) + 1 
  cl = floor((ind-1) / rows) + 1
  loc = c(rw, cl)
  return(loc)
}

ecology_change_function <- function(pos_y, params, rate, t){
  mn = params$min_eco_val
  mx = params$max_eco_val
  t_sh = -1/rate * log( ((pos_y - mn)/(mx - pos_y)))
  ecology_change_function = mn + (mx - mn)/(1 + exp(-rate*( t - t_sh)))
}

ecology_change_function_vectorised <- function(parcel_vals, params, curve_rate, time_remaining){
  time_array = matrix(rep(time_remaining, length(parcel_vals)), ncol = length(parcel_vals))
  mn = params$min_eco_val
  mx = params$max_eco_val
  t_sh = -1/curve_rate * log( ((parcel_vals - mn)/(mx - parcel_vals)))
  t_sh_array = matrix(rep(t_sh, length(time_remaining)), ncol = length(time_remaining))
  t_sh_array <- t(t_sh_array)
  ecology_change_function = mn + (mx - mn)/(1 + exp(-curve_rate*( time_array - t_sh_array)))
}




build_decline_rates <- function(region_num, regions, params){
  decline_rates = vector('list', region_num)
  
  for (region_ind in 1:region_num){ 
    current_region = regions[[region_ind]]
    current_parcel_num = length(current_region)    
    decline_params = c(length(current_region), params$mean_decline_rates[region_ind], params$decline_rate_std) #params$decline_rate_std[region_ind])
    current_decline_rates = matrix(rnorm(decline_params[1], mean = decline_params[2], sd = decline_params[3]), ncol = ncol(current_region))
    decline_rates[[region_ind]] = current_decline_rates
  }
  return(decline_rates)
}


build_counterfactuals <- function(region_num, regions, params, decline_rates, land_parcels, ecology_initial, time_vec){
  elements = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))
  for (region_ind in 1:region_num){ 
    current_parcel_pool = regions[[region_ind]]
    current_parcel_num = length(current_parcel_pool)    
    current_decline_rates = decline_rates[[region_ind]]
    
    for (parcel_ind in 1:current_parcel_num){
      current_parcel_ind = current_parcel_pool[parcel_ind]
      current_parcel = land_parcels[[current_parcel_ind]]
      current_dec_rate = current_decline_rates[parcel_ind]      
      predicted_traj = ecology_change_function_vectorised(ecology_initial[current_parcel], params, current_dec_rate, time_vec)
      for (pix_ind in 1:length(current_parcel)){
        loc = ind2sub(params$ecology_size, current_parcel[pix_ind])
        elements[loc[1], loc[2], 1:params$time_steps] = predicted_traj[, pix_ind]        
      }
    }
  }
  return(elements)
}

build_counterfactuals_blur <- function(region_num, regions, params, decline_rates, land_parcels, ecology_initial, time_vec){
  elements = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))
  for (yr in 1:params$time_steps){
    
    for (region_ind in 1:region_num){ 
      current_parcel_pool = regions[[region_ind]]
      current_parcel_num = length(current_parcel_pool)    
      current_decline_rates = decline_rates[[region_ind]]
      
      for (parcel_ind in 1:current_parcel_num){
        current_parcel_ind = current_parcel_pool[parcel_ind]
        current_parcel = land_parcels[[current_parcel_ind]]
        current_dec_rate = current_decline_rates[parcel_ind]      
        predicted_traj = ecology_change_function_vectorised(ecology_initial[current_parcel], params, current_dec_rate, time_vec)
        for (pix_ind in 1:length(current_parcel)){
          loc = ind2sub(params$ecology_size, current_parcel[pix_ind])
          elements[loc[1], loc[2], yr] = predicted_traj[, pix_ind]        
        }
      }
    }
    
  }
  return(elements)
}


initialise_counterfactuals_vectorised <- function(params, parcels){
  counterfactuals = list() 
  ecology_initial <- initialise_ecology(params, parcels)
  
  region_num = parcels$region_num
  regions = parcels$regions
  land_parcels = parcels$land_parcels
  time_vec = 1:params$time_steps 
  
  decline_rates = build_decline_rates(parcels$region_num, parcels$regions, params)
  elements = build_counterfactuals(region_num, regions, params, decline_rates, land_parcels, ecology_initial, time_vec)
  final_ecology = elements[, , params$time_steps]
  final_parcel_vals = find_current_parcel_vals(parcels, final_ecology, 1:length(land_parcels))
  
  counterfactuals$final_parcel_vals = final_parcel_vals
  counterfactuals$final_ecology = final_ecology
  counterfactuals$elements = elements
  counterfactuals$decline_rates = decline_rates
  return(counterfactuals)
}




initialise_outputs <- function(params, counterfactuals, parcels){
  parcel_num_x = length(parcels$parcel_vx)
  parcel_num_y = length(parcels$parcel_vy)
  outputs = list()
  outputs$trajectories = counterfactuals$elements
  outputs$adjusted_counterfactuals = counterfactuals$elements
  outputs$parcel_sets = list()
  outputs$parcel_sets$dev_count = 1
  return(outputs)
  
}


update_development_object <- function(params, ind_available, parcels, yr, current_ecology, final_ecology){
  development_object = list()
  land_parcels = parcels$land_parcels
  parcel_index = ind_available[sample(1:length(ind_available), 1)]
  current_parcel = land_parcels[[parcel_index]]
  current_parcel_val = sum(current_ecology[current_parcel])
  if (params$dev_val_type == 'current'){
    parcel_val_to_use = current_parcel_val
  } else if (params$dev_val_type == 'predicted'){
    parcel_val_to_use= sum(final_ecology[current_parcel])
  }
  
  development_object$parcel_index = parcel_index
  development_object$current_parcel_val = current_parcel_val
  development_object$parcel_val_to_use = parcel_val_to_use
  development_object$current_parcel = current_parcel
  return(development_object)
}


update_index <- function(index_object, current_parcel_indexes, region_ind, params){
  
  ind_available = index_object$ind_available[[region_ind]]
  ind_available <- setdiff(ind_available, current_parcel_indexes) #remove development parcel from available list   
  if (length(ind_available) < params$region_dev_num[[region_ind]] ){   
    break_flag = TRUE
  } else {break_flag = FALSE}
  
  index_object$ind_available[[region_ind]] = ind_available
  index_object$break_flag = break_flag
  
  return(index_object)
  
}

write_development_trajectory <- function(trajectories, counterfactuals, development_object, current_ecology, params, yr){
  
  current_parcel = development_object$current_parcel
  dim(current_parcel) = c(1, length(current_parcel))
  for (pix_ind in current_parcel){
    loc = ind2sub(params$ecology_size, pix_ind)
    trajectories[loc[1], loc[2], (yr + 1):params$time_steps] = 0
  }   
  return(trajectories)
}


find_offset_pool <- function(index_object, region_ind, params, parcels, current_ecology, time_remaining, final_ecology){
  
  offset_pool = list()
  current_offset_pool = index_object$ind_available[[region_ind]]
  current_parcel_vals <- find_current_parcel_vals(parcels, current_ecology, current_offset_pool)
  
  if (params$offset_selection_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
    parcel_vals_to_use <- predict_final_parcel_vals(parcels, current_offset_pool, current_ecology, params, params$offset_rate, time_remaining)
    if (params$predict_offset_mode == 'gains'){
      parcel_vals_to_use = parcel_vals_to_use - counterfactuals$final_parcel_vals[current_offset_pool]
    }
    
  } else if (params$offset_selection_type == 'current'){
    parcel_vals_to_use = current_parcel_vals
  }
  
  offset_pool$current_parcel_pool = current_offset_pool
  offset_pool$current_parcel_vals = current_parcel_vals
  offset_pool$parcel_vals_to_use = parcel_vals_to_use
  
  return(offset_pool)
}


predict_final_parcel_vals <- function(parcels, parcel_pool, current_ecology, params, curve_rate, time_remaining){
  
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(parcel_pool))  
  predicted_parcel_trajs = predict_parcel_trajectories(parcels, parcel_pool, params, current_ecology, curve_rate, time_remaining)
  predicted_parcel_vals = predicted_parcel_trajs[, length(time_remaining)]                                  
  return(predicted_parcel_vals)  
  
}

predict_parcel_trajectories <- function(parcels, parcel_pool, params, current_ecology, curve_rate, time_remaining){
  
  land_parcels = parcels$land_parcels
  predicted_parcel_trajs = array(0, c(length(parcel_pool), length(time_remaining)))
  
  for (parcel_ind in 1:length(parcel_pool)){
    current_parcel_ind = parcel_pool[parcel_ind]
    current_parcel = land_parcels[[current_parcel_ind]]
    predicted_pixel_trajs = ecology_change_function_vectorised(current_ecology[current_parcel], params, curve_rate, time_remaining)
    predicted_parcel_trajs[parcel_ind, ] = rowSums(predicted_pixel_trajs)
  }
  
  return(predicted_parcel_trajs)
  
}

update_offset_object <- function(offset_pool, params, parcels, development_val, yr, current_ecology){      
  
  offset_object = list()
  current_offset_pool = offset_pool$current_parcel_pool
  parcel_vals_to_use = offset_pool$parcel_vals_to_use  
  land_parcels = parcels$land_parcels
  val_to_match = params$offset_multiplier*development_val
  
  parcel_vals_used = vector()
  parcel_indexes = vector()
  
  if (sum(parcel_vals_to_use)<val_to_match){
    offset_object$break_flag = TRUE
    return(offset_object)
  }  else {offset_object$break_flag = FALSE}
  
  #  err = abs(val_to_match - parcel_vals_to_use)
  #  best_ind = which(err == min(err))
  #  parcel_indexes = current_offset_pool[best_ind]
  
  while (val_to_match > 0){
    err = val_to_match - parcel_vals_to_use
    if ( all(err > 0) ){
      max_ind = which(parcel_vals_to_use == max(parcel_vals_to_use))
      parcel_indexes = c(parcel_indexes, current_offset_pool[max_ind])
      parcel_vals_used = c(parcel_vals_used, parcel_vals_to_use[max_ind])
      val_to_match = val_to_match - parcel_vals_to_use[max_ind]     
      current_offset_pool = current_offset_pool[-max_ind]
      parcel_vals_to_use = parcel_vals_to_use[-max_ind]
    } else{
      max_ind = which(err == max(err[which(err <= 0)]))
      parcel_indexes = c(current_offset_pool[max_ind], parcel_indexes)
      parcel_vals_used = c(parcel_vals_used, parcel_vals_to_use[max_ind])
      break
    }
  }
  
  parcel_num = length(parcel_indexes)
  offset_object$current_parcels = list()
  offset_object$current_parcel_vals = array(0, parcel_num)
  
  for (parcel_ind in 1:parcel_num){
    current_parcel = land_parcels[[parcel_indexes[parcel_ind]]]
    offset_object$current_parcels[[parcel_ind]] = current_parcel
    offset_object$current_parcel_vals[parcel_ind] = sum(current_ecology[current_parcel])
  }  
  
  offset_object$current_parcel_indexes = parcel_indexes
  offset_object$parcel_vals_used = parcel_vals_used
  return(offset_object)  
}

write_null_offset_object <- function(){
  offset_object = list()
  offset_object$current_parcel_indexes = list()
  offset_object$current_parcel_vals = list()
  offset_object$parcel_vals_used = list()
  return(offset_object)
}



write_offset_trajectory_vectorised <- function(trajectories, offset_object, current_ecology, params, time_remaining, yr){
  current_parcels = offset_object$current_parcels
  offset_rate = params$offset_rate
  
  parcel_num = length(offset_object$current_parcels)
  for (parcel_ind in 1:parcel_num){
    current_parcel = current_parcels[[parcel_ind]]
    predicted_traj = ecology_change_function_vectorised(current_ecology[current_parcel], params, offset_rate, time_remaining)
    for (pix_ind in 1:length(current_parcel)){
      loc = ind2sub(params$ecology_size, current_parcel[pix_ind])
      trajectories[loc[1], loc[2], yr:params$time_steps] = predicted_traj[, pix_ind]        
    }
  }
  return(trajectories)
}


find_current_parcel_vals <- function(parcels, current_ecology, parcel_indexes){
  
  land_parcels = parcels$land_parcels
  parcel_num = length(parcel_indexes)
  current_parcel_vals = rep(0, parcel_num)
  
  for (s in 1:parcel_num){
    parcel_ind = parcel_indexes[s]
    current_parcel = land_parcels[[parcel_ind]]
    current_parcel_vals[s] = sum(current_ecology[current_parcel])
  }
  return(current_parcel_vals)
}





predict_parcels_vectorised <- function(region_ind, parcels, ind_available, current_ecology, params, time_remaining, final_ecology){
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(ind_available))
  
  for (s in 1:length(ind_available)){
    current_ind = ind_available[s]
    current_parcel = land_parcels[[current_ind]]
    predicted_trajs = ecology_change_function_vectorised(current_ecology[current_parcel], params, offset_rate, time_remaining)
    predicted_parcel_vals[s] = sum(predicted_trajs[nrow(predicted_trajs), ])
  }
  
  return(predicted_parcel_vals)
  
}




find_parcel_trajectory <- function(parcels, parcel_ind, params, trajectories){
  land_parcels = parcels$land_parcels
  current_parcel = land_parcels[[parcel_ind]]
  parcel_traj = array(0, params$time_steps)
  
  for (yr in 1:params$time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    parcel_traj[yr] = sum(current_slice[current_parcel])
  }
  return(parcel_traj)
}


find_parcel_trajectories <- function(parcels, parcel_indexes, params, trajectories){
  land_parcels = parcels$land_parcels
  parcel_trajs = array(0, c(length(parcel_indexes), params$time_steps))
  
  for (parcel_ind in 1:length(parcel_indexes)){
    current_ind = parcel_indexes[parcel_ind]
    parcel_trajs[parcel_ind, ] = find_parcel_trajectory(parcels, current_ind, params, trajectories)
  }
  
  if (length(parcel_indexes) == 1){
    dim(parcel_trajs) = c(length(parcel_trajs), 1)
  }
  
  return(parcel_trajs)
}

write_parcel_sets <- function(parcel_sets, dev_count, yr, region_ind, development_object, offset_object){
  parcel_set = list()
  parcel_set$yr = yr
  parcel_set$region = region_ind
  parcel_set$developed_parcel = development_object$parcel_index
  parcel_set$dev_val = development_object$current_parcel_val
  parcel_set$offset_parcels = offset_object$current_parcel_indexes
  parcel_set$offset_vals = offset_object$current_parcel_vals
  
  parcel_sets[[dev_count]] = parcel_set
  #print(c(parcel_sets[[dev_count]], yr, dev_count))  
  return(parcel_sets)
}



plot_outs <- function(...){
  
  dots = list(...)
  plot_params = dots[[length(dots)]]
  graphics.off()
  plot_num = length(dots) - 1
  A = 1:plot_params[3]
  dim(A) = c(plot_params[1], plot_params[2])
  layout(A)
  ymaxs = array(0, plot_num)
  ymins = array(0, plot_num)
  for (s in 1:(length(dots) - 1)){
    ymaxs[s] = max(dots[[s]])
    ymins[s] = min(dots[[s]])
  }
  ymax = max(ymaxs)
  ymin = min(ymins)
  for (s in 1:plot_params[3]){
    plot(dots[[1]][, s], type = 'l', col = 'red', ylim = c(ymin, ymax))
    if (plot_num > 1){
      for (t in 2:plot_num){
        lines(dots[[t]][, s],  ylim = c(ymin, ymax))
      }
    } 
  }
}





calc_diffs <- function(parcel_list, trajectories, params, parcels, counterfactuals){
  
  outputs = list()
  true_change = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  traj_rel_to_initial = true_change
  initial_rel_to_counter = true_change
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  
  for (yr in 1:params$time_steps){
    
    for (parcel_ind in (parcel_list[[yr]])){
      current_parcel = land_parcels[[parcel_ind]]
      for (pix_ind in current_parcel){
        loc = ind2sub(params$ecology_size, pix_ind)
        initial_val = trajectories[loc[1], loc[2], yr]
        current_trajectory = trajectories[loc[1], loc[2], yr:params$time_steps]
        true_change[loc[1], loc[2], yr:params$time_steps] = current_trajectory - counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]       
        traj_rel_to_initial[loc[1], loc[2], yr:params$time_steps] = current_trajectory - initial_val
        initial_rel_to_counter[loc[1], loc[2], yr:params$time_steps] = initial_val - counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]     
      }
    }
  }
  
  outputs$true_change = true_change
  outputs$traj_rel_to_initial = traj_rel_to_initial
  outputs$initial_rel_to_counter = initial_rel_to_counter
  return(outputs)
  
}









sum_parcels <- function(time_steps, parcels, trajectories){
  
  land_parcel_num = length(parcels$land_parcels)
  land_parcels = parcels$land_parcels
  parcel_sums = array(0, c(time_steps, land_parcel_num))
  
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    
    for (parcel_ind in 1:land_parcel_num){   
      current_parcel = land_parcels[[parcel_ind]]
      parcel_sums[yr, parcel_ind] = sum(current_slice[current_parcel])
    }
  }
  
  return(parcel_sums)
  
}


sum_regions <- function(time_steps, parcel_sums, parcels){
  
  region_num = parcels$region_num
  regions = parcels$regions
  region_sums = array(0, c(time_steps, region_num))
  for (yr in 1:time_steps){ 
    current_slice = parcel_sums[yr, ]
    for (region_ind in 1:region_num){ #number of regions  
      current_region = regions[[region_ind]]
      region_sums[yr, region_ind] = sum(current_slice[current_region])    
    }
  }
  return(region_sums)
}



find_sums <- function(time_steps, parcels, trajectories){
  object = list()
  object$parcel_sums = sum_parcels(time_steps, parcels, trajectories)
  object$region_sums = sum_regions(time_steps, object$parcel_sums, parcels)
  object$net_sums = rowSums(object$region_sums)
  return(object)
}





view_individual_offsets <- function(parcel_set, parcels, params, trajectories, counterfactuals){
  assess_object = list()
  offset_yr = parcel_set$yr
  developed_parcel = parcel_set$developed_parcel
  offset_parcels = parcel_set$offset_parcels
  dev_trajectory = find_parcel_trajectory(parcels, developed_parcel, params, trajectories)
  dev_counterfactual = find_parcel_trajectory(parcels, developed_parcel, params, counterfactuals$elements)
  
  offset_trajectories = find_parcel_trajectories(parcels, offset_parcels, params, trajectories)
  
  offset_counterfactuals = find_parcel_trajectories(parcels, offset_parcels, params, counterfactuals$elements)
  
  if (length(offset_parcels) > 1){
    offset_sum = colSums(offset_trajectories)
    offset_counter_sum = colSums(offset_counterfactuals)
  } else {
    offset_sum = offset_trajectories
    offset_counter_sum = offset_counterfactuals
  }
  
  assess_object$dev_trajectory = dev_trajectory
  assess_object$offset_trajectories = offset_trajectories
  
  graphics.off()
  mx = max(c(dev_trajectory, offset_sum))
  plot(assess_object$dev_trajectory, type = 'l', ylim = c(0, mx))
  lines((offset_sum - offset_counter_sum), col = 'red', ylim = c(0, mx))
  lines((offset_sum), ylim = c(0, mx))
  lines(dev_counterfactual, ylim = c(0, mx))
  lines(offset_counter_sum, ylim = c(0, mx))
  chk = (offset_sum[params$time_steps] - offset_counter_sum[params$time_steps])/(params$offset_multiplier*dev_counterfactual[params$time_steps])
  print(chk)
  return(assess_object)
}



perform_offset_routine <- function(){
  
  offset_pool <- find_offset_pool(index_object, region_ind, params, parcels, current_ecology, time_remaining, counterfactuals$final_ecology)
  offset_object <- update_offset_object(offset_pool, params, parcels, development_object$parcel_val_to_use, yr, current_ecology)                  
  if (offset_object$break_flag == TRUE){
    break
  }
  
  outputs$trajectories = write_offset_trajectory_vectorised(outputs$trajectories, offset_object, current_ecology, params, time_remaining, yr)         
  index_object = update_index(index_object, offset_object$current_parcel_indexes, region_ind, params)           
  if (index_object$break_flag == TRUE){
    break
  }
  outputs$parcel_sets <- write_parcel_sets(outputs$parcel_sets, dev_count, yr, region_ind, development_object, offset_object)
  
}

find_current_dev_nums <- function(params){
  current_dev_nums = params$region_dev_num
  return(current_dev_nums)
}

calc_trajectories <- function(params, outputs, parcels, index_object, counterfactuals, perform_offsets){
  dev_count = 1;
  for (yr in 1:params$time_steps){      #main time loop    
    time_remaining = 0:(params$time_steps - yr) 
    current_ecology = outputs$trajectories[ , , yr]
    
    if (yr%%params$develop_every == 0 && yr < params$time_steps){
      current_dev_nums <- find_current_dev_nums(params)
      
      for (region_ind in 1:parcels$region_num){
        
        current_develop_num = current_dev_nums[region_ind]
        outputs$adjusted_counterfactuals <- probability_counterfactuals(outputs$adjusted_counterfactuals, index_object$ind_available[[region_ind]], current_develop_num, parcels$land_parcels, parcels$regions[[region_ind]], yr)
        
        if (current_develop_num > 0){
          for (dev_index in 1:current_develop_num){    #loop over the number of possible developments           
            if (length(index_object$ind_available[[region_ind]]) < 1){
              break
            }
            
            development_object <- update_development_object(params, index_object$ind_available[[region_ind]], parcels, yr, current_ecology, counterfactuals$final_ecology)
            outputs$trajectories <- write_development_trajectory(outputs$trajectories, counterfactuals, development_object, current_ecology, params, yr)     
            index_object <- update_index(index_object, development_object$parcel_index, region_ind, params)      
            if (index_object$break_flag == TRUE){
              break
            } 
            
            if (perform_offsets == TRUE){
              offset_pool <- find_offset_pool(index_object, region_ind, params, parcels, current_ecology, time_remaining, counterfactuals$final_ecology)
              offset_object <- update_offset_object(offset_pool, params, parcels, development_object$parcel_val_to_use, yr, current_ecology)                  
              if (offset_object$break_flag == TRUE){
                break
              }
              if (params$offset_policy == 'managed'){
                outputs$trajectories = write_offset_trajectory_vectorised(outputs$trajectories, offset_object, current_ecology, params, time_remaining, yr)
              }
              
              index_object = update_index(index_object, offset_object$current_parcel_indexes, region_ind, params)           
              if (index_object$break_flag == TRUE){
                break
              }               
            } else {offset_object <- write_null_offset_object()}
            
            outputs$parcel_sets <- write_parcel_sets(outputs$parcel_sets, dev_count, yr, region_ind, development_object, offset_object)
            dev_count = dev_count + 1;
            
          }
          
        } 
      }
    }
    #print(yr)
  }
  
  return(outputs)
}


calc_trajectories_blur <- function(params, outputs, parcels, index_object, counterfactuals, perform_offsets){
  dev_count = 1;
  for (yr in 1:params$time_steps){      #main time loop    
    time_remaining = 0:(params$time_steps - yr) 
    current_ecology = outputs$trajectories[ , , yr]
    
    if (yr%%params$develop_every == 0 && yr < params$time_steps){
      current_dev_nums <- find_current_dev_nums(params)
      
      for (region_ind in 1:parcels$region_num){
        
        current_develop_num = current_dev_nums[region_ind]
        outputs$adjusted_counterfactuals <- probability_counterfactuals(outputs$adjusted_counterfactuals, index_object$ind_available[[region_ind]], current_develop_num, parcels$land_parcels, parcels$regions[[region_ind]], yr)
        
        if (current_develop_num > 0){
          for (dev_index in 1:current_develop_num){    #loop over the number of possible developments           
            if (length(index_object$ind_available[[region_ind]]) < 1){
              break
            }
            
            development_object <- update_development_object(params, index_object$ind_available[[region_ind]], parcels, yr, current_ecology, counterfactuals$final_ecology)
            outputs$trajectories <- write_development_trajectory(outputs$trajectories, counterfactuals, development_object, current_ecology, params, yr)     
            index_object <- update_index(index_object, development_object$parcel_index, region_ind, params)      
            if (index_object$break_flag == TRUE){
              break
            } 
            
            if (perform_offsets == TRUE){
              offset_pool <- find_offset_pool(index_object, region_ind, params, parcels, current_ecology, time_remaining, counterfactuals$final_ecology)
              offset_object <- update_offset_object(offset_pool, params, parcels, development_object$parcel_val_to_use, yr, current_ecology)                  
              if (offset_object$break_flag == TRUE){
                break
              }
              if (params$offset_policy == 'managed'){
                outputs$trajectories = write_offset_trajectory_vectorised(outputs$trajectories, offset_object, current_ecology, params, time_remaining, yr)
              }
              
              index_object = update_index(index_object, offset_object$current_parcel_indexes, region_ind, params)           
              if (index_object$break_flag == TRUE){
                break
              }               
            } else {offset_object <- write_null_offset_object()}
            
            outputs$parcel_sets <- write_parcel_sets(outputs$parcel_sets, dev_count, yr, region_ind, development_object, offset_object)
            dev_count = dev_count + 1;
            
          }
          
        } 
      }
    }
    #print(yr)
  }
  
  return(outputs)
}

# update_counterfactuals <- function(adjusted_counterfactuals, current_parcel_set, current_develop_num, land_parcels, region, yr){
#   current_parcel_num = length(current_parcel_set)
#   current_prob = current_develop_num/current_parcel_num
#   current_slice = adjusted_counterfactuals[, , yr]
#   for (parcel_ind in current_parcel_set){
#     current_parcel = land_parcels[[parcel_ind]]
#     current_slice[current_parcel] = (1 - current_prob)*current_slice[current_parcel]
#   }
#   for (parcel_ind in setdiff(region, current_parcel_set)){
#     current_parcel = land_parcels[[parcel_ind]]
#     current_slice[current_parcel] = 0
#   }
#   adjusted_counterfactuals[, , yr] = current_slice
#   return(adjusted_counterfactuals)
# }

probability_counterfactuals <- function(adjusted_counterfactuals, current_parcel_set, current_develop_num, land_parcels, region, yr){
  current_parcel_num = length(current_parcel_set)
  current_prob = 1 - current_develop_num/current_parcel_num
  current_slice = adjusted_counterfactuals[, , yr]
  for (parcel_ind in length(region)){
    current_parcel = land_parcels[[parcel_ind]]
    current_slice[current_parcel] = current_prob*current_slice[current_parcel]
  }
  adjusted_counterfactuals[, , yr] = current_slice
  return(adjusted_counterfactuals)
}

cartesian_mesh <- function(N, M){
  mesh = list()
  xx = seq(-floor(M/2), (ceiling(M/2) - 1))
  yy = seq(-floor(N/2), (ceiling(N/2) - 1))
  x = matrix(rep(xx,each=N),nrow=N);
  y = matrix(rep(yy,M),nrow=N)
  mesh$x = x
  mesh$y = y
  return(mesh)
#  m=length(x); n=length(y);
#  X=matrix(rep(x,each=n),nrow=n);
#  Y=matrix(rep(y,m),nrow=n)
}
  
gauss <- function(x, y, sig_x, sig_y){
  g = exp(-x^2/(sig_x^2)) * exp(-y^2/(sig_y^2))
  return(g)
}


fftshift <- function(fft_array){
  numDims = 2
  idx <- vector('list', numDims)
  for (k in 1:numDims){
    m = dim(fft_array)[k];
    p = ceiling(m/2);
    idx[[k]] = c((p+1):m, (1:p));
  }
  fft_array = fft_array[idx[[1]], ]
  fft_array = fft_array[, idx[[2]]]
  return(fft_array)
}


Blur_2D <- function(A, sig_x, sig_y){
  dims = dim(A)
  M = dims[1]
  N = dims[2]
  mesh = cartesian_mesh(M, N)
  knl = gauss(mesh$x, mesh$y, sig_x, sig_y)
  knl = knl / sum(knl)
  convolve = fftshift( fft( fftshift(knl), inverse = TRUE )) * fftshift( fft( fftshift(A), inverse = TRUE ) )
  B = (fftshift( fft( fftshift( convolve ) ) ))/(M*N)
  B = Re(B)
  return(A)
}



collate_parcel_sets <- function(parcel_sets, params){
  parcels_object = list()
  offset_parcels = vector('list', params$time_steps)
  offset_vals = vector('list', params$time_steps)
  offset_parcels_total = vector()
  offset_vals_total = vector()
  developed_parcels = vector('list', params$time_steps)
  for (dev_ind in 1:length(parcel_sets)){
    yr = parcel_sets[[dev_ind]]$yr
    offset_parcels[[yr]] = c(offset_parcels[[yr]], parcel_sets[[dev_ind]]$offset_parcels)
    offset_vals[[yr]] = c(offset_vals[[yr]], parcel_sets[[dev_ind]]$offset_vals)
    offset_parcels_total = c(offset_parcels_total, parcel_sets[[dev_ind]]$offset_parcels)
    offset_vals_total = c(offset_vals_total, parcel_sets[[dev_ind]]$offset_vals)
    developed_parcels[[yr]] = c(developed_parcels[[yr]], parcel_sets[[dev_ind]]$developed_parcel)
  }
  parcels_object$offset_parcels = offset_parcels
  parcels_object$offset_vals = offset_vals
  parcels_object$developed_parcels = developed_parcels
  parcels_object$offset_vals_total = offset_vals_total
  parcels_object$offset_parcels_total = offset_parcels_total
  return(parcels_object)
}