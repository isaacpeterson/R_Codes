initialise_parameters <- function(){
  params = list()
  params$offset_multiplier = 1
  params$parcel_num_x = 20 #number of parcels in x
  params$parcel_num_y = 30 #number of parcels in y
  params$region_num_x = 2 
  params$region_num_y = 3
  params$region_num = params$region_num_x*params$region_num_y
  params$ecology_size = 250 #array size to be broken up into sub arrays
  params$time_steps = 100 #number of years etc
  params$max_developments = 1 #maximum number of developments per year 
  params$update_every = 10
  params$min_eco_val = 0
  params$max_eco_val = 100
  params$min_initial_eco_val = 40
  params$max_initial_eco_val = 80
  params$max_decline_rate = 0.02 #max regional decline rate
  params$offset_rate = 0.05 #rate of improvement in land parcel after parcel is offset
  params$eco_noise = 10
  params$parcel_selection_type = 'regional' #regional or national
  params$offset_selection_type = 'predicted' 
  params$dev_val_type = 'current'
  params$predict_offset_mode = 'gains'
  params$offset_value_type = 'match' #match - match absolute value of offset parcel to the development value (predicted or current)
  params$region_num = params$region_num_x*params$region_num_y 
  params$development_rates = sample(1:params$max_developments, params$region_num, replace = T)
  #params$development_rates = array(1, params$region_num)
  params$decline_rates = -params$max_decline_rate*runif((params$region_num - 1))
  params$decline_rates = c(params$decline_rates, params$max_decline_rate)
  #params$decline_rate_std = 0.1*params$max_decline_rate*runif((params$region_num))
  params$decline_rate_std = 0.5*params$max_decline_rate
  return(params)
}

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
    initial_parcel_value = params$min_initial_eco_val + (params$max_initial_eco_val - params$min_initial_eco_val - params$eco_noise)*runif(1) 
    inds = parcels$land_parcels[[s]]
    dim(inds) = c(1, length(inds))
    ecology_initial[inds] = ecology_initial[inds]*initial_parcel_value
  }
  ecology_initial = ecology_initial + params$eco_noise*matrix(runif(params$ecology_size*params$ecology_size),params$ecology_size,params$ecology_size)
  return(ecology_initial)
}

initialise_offset_object <- function(params){
  object = list()
  object$parcel_list = vector('list', params$time_steps)
  object$parcel_vals_list = vector('list', params$time_steps)
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

offset_rule <- function(pos_y, params, rate, t){
  mn = params$min_eco_val
  mx = params$max_eco_val
  t_sh = -1/rate * log( ((pos_y - mn)/(mx - pos_y)))
  offset_rule = mn + (mx - mn)/(1 + exp(-rate*( t - t_sh)))
}

offset_rule_vectorised <- function(parcel_vals, params, rate, time_remaining){
  time_array = matrix(rep(time_remaining, length(parcel_vals)), ncol = length(parcel_vals))
  mn = params$min_eco_val
  mx = params$max_eco_val
  t_sh = -1/rate * log( ((parcel_vals - mn)/(mx - parcel_vals)))
  t_sh_array = matrix(rep(t_sh, length(time_remaining)), ncol = length(time_remaining))
  t_sh_array <- t(t_sh_array)
  offset_rule = mn + (mx - mn)/(1 + exp(-rate*( time_array - t_sh_array)))
}

initialise_counterfactuals <- function(params, parcels){
  counterfactuals = list() 
  ecology_initial <- initialise_ecology(params, parcels)
  
  region_num = parcels$region_num
  regions = parcels$regions
  land_parcels = parcels$land_parcels
  current_ecology = ecology_initial
  time_vec = 1:params$time_steps 
  elements = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))
  decline_rates = vector('list', region_num)
  
  for (region_ind in 1:region_num){
    
    current_region = regions[[region_ind]]
    current_parcel_num = length(current_region)    
    decline_params = c(length(current_region), params$decline_rates[region_ind], params$decline_rate_std) #params$decline_rate_std[region_ind])
    current_decline_rates = matrix(rnorm(decline_params[1], mean = decline_params[2], sd = decline_params[3]), ncol = ncol(current_region))
    decline_rates[[region_ind]] = current_decline_rates
    parcel_inds = current_region
    dim(parcel_inds) = c(1, current_parcel_num)
    
    for (parcel_count in 1:current_parcel_num){
      parcel_ind = parcel_inds[parcel_count]
      current_parcel = land_parcels[[parcel_ind]]
      current_dec_rate = current_decline_rates[parcel_count]      
      for (pix_ind in current_parcel){      
        loc <- ind2sub(params$ecology_size, pix_ind)
        elements[loc[1], loc[2], ] = offset_rule(current_ecology[pix_ind], params, current_dec_rate, time_vec)
      }
    }
  }
  counterfactuals$final_ecology = elements[, , params$time_steps]
  counterfactuals$elements = elements
  counterfactuals$decline_rates = decline_rates
  return(counterfactuals)
}

initialise_parcel_sets <- function(params){
  parcel_sets = vector('list', params$time_steps)
  for (yr in 1:params$time_steps){
    parcel_sets[[yr]] = list()
  }
  return(parcel_sets)
}


initialise_outputs <- function(params, counterfactuals, parcels){
  parcel_num_x = length(parcels$parcel_vx)
  parcel_num_y = length(parcels$parcel_vy)
  outputs = list()
  outputs$trajectories = counterfactuals$elements
  outputs$offset_object <- initialise_offset_object(params)
  outputs$development_object <- initialise_offset_object(params)
  outputs$parcel_sets = list() #initialise_parcel_sets(params)
  return(outputs)
  
}



select_parcel_to_develop <- function(params, development_object, ind_available, parcels, yr, current_ecology, final_ecology){

  land_parcels = parcels$land_parcels
  parcel_index = ind_available[sample(1:length(ind_available), 1)]
   #determine ecological value of selected parcel   
  if (params$dev_val_type == 'current'){
    parcel_val = find_current_parcel_value(parcels, parcel_index, current_ecology)
  } else if (params$dev_val_type == 'predicted'){
    parcel_val = sum(final_ecology[land_parcels[[parcel_index]]])
  }
  current_parcel = land_parcels[[parcel_index]]
  
  development_object$parcel_index = parcel_index
  development_object$current_parcel_val = parcel_val
  development_object$current_parcel = current_parcel

  development_object$parcel_list[[yr]] = c(development_object$parcel_list[[yr]], parcel_index)
  development_object$parcel_vals_list[[yr]] = c(development_object$parcel_vals_list[[yr]], parcel_val)
  
  return(development_object)
}


update_index <- function(index_object, current_parcel_indexes, region_ind, params){
  
  ind_available = index_object$ind_available[[region_ind]]
  ind_available <- setdiff(ind_available, current_parcel_indexes) #remove development parcel from available list   
  if (length(ind_available) < params$development_rates[[region_ind]] ){   
    break_flag = TRUE
  } else {break_flag = FALSE}
  
  index_object$ind_available[[region_ind]] = ind_available
  index_object$break_flag = break_flag
  
  return(index_object)
  
}


select_offset_parcel <- function(offset_pool, offset_object, params, parcels, development_val, yr){      
 
  offset_ind_to_use = offset_pool$offset_ind_to_use
  parcel_vals_to_use = offset_pool$current_parcel_vals
  #current_parcel_vals = offset_pool$current_parcel_vals
  #parcel_vals_to_use = current_parcel_vals[offset_ind_to_use] 
  parcel_indexes = vector()
  parcel_vals = vector()
  
  land_parcels = parcels$land_parcels
  
  val_to_match = params$offset_multiplier*development_val
 
  parcel_indexes = vector()
  
  if (sum(parcel_vals_to_use)<val_to_match){
    offset_object$break_flag = TRUE
    return(offset_object)
  }  else {offset_object$break_flag = FALSE}
  
  err = abs(val_to_match - parcel_vals_to_use)
  best_ind = which(err == min(err))
  parcel_indexes = offset_ind_to_use[best_ind]
#   while (val_to_match > 0){
#     err = val_to_match - parcel_vals_to_use
#     if ( all(err > 0) ){
#       max_ind = which(parcel_vals_to_use == max(parcel_vals_to_use))
#       parcel_indexes = c(parcel_indexes, offset_ind_to_use[max_ind])
#       parcel_vals = c(parcel_vals, parcel_vals_to_use[max_ind])
#       val_to_match = val_to_match - parcel_vals_to_use[max_ind]     
#       offset_ind_to_use = offset_ind_to_use[-max_ind]
#       #parcel_vals_to_use = current_parcel_vals[offset_ind_to_use]
#       parcel_vals_to_use = parcel_vals_to_use[-max_ind]
#     } else{
#       max_ind = which(err == max(err[which(err <= 0)]))
#       parcel_indexes = c(offset_ind_to_use[max_ind], parcel_indexes)
#       parcel_vals = c(parcel_vals, parcel_vals_to_use[max_ind])
#       break
#     }
#   }
  
  offset_object$current_parcels = list()
  for (parcel_ind in 1:length(parcel_indexes)){
    parcel_index = parcel_vals
    offset_object$current_parcels[[parcel_ind]] = land_parcels[[parcel_indexes[parcel_ind]]]
  }
  
  offset_object$current_parcel_indexes = parcel_indexes
  offset_object$parcel_list[[yr]] = c(offset_object$parcel_list[[yr]], parcel_indexes)
  offset_object$parcel_vals_list[[yr]] = c(offset_object$parcel_vals_list[[yr]], parcel_vals)
  
  return(offset_object)  
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

write_offset_trajectory <- function(trajectories, offset_object, current_ecology, params, time_remaining, yr){
  current_parcels = offset_object$current_parcels
  offset_rate = params$offset_rate

  parcel_num = length(offset_object$current_parcels)
  for (parcel_ind in 1:parcel_num){
    current_parcel = current_parcels[[parcel_ind]]
    dim(current_parcel) = c(1, length(current_parcel))
    for (pix_ind in current_parcel){
      loc = ind2sub(params$ecology_size, pix_ind)
      current_val = current_ecology[pix_ind]
      current_trajectory = offset_rule(current_val, params, offset_rate, time_remaining)
      trajectories[loc[1], loc[2], yr:params$time_steps] = current_trajectory        
    }
  }
  return(trajectories)
}

write_offset_trajectory_vectorised <- function(trajectories, offset_object, current_ecology, params, time_remaining, yr){
  current_parcels = offset_object$current_parcels
  offset_rate = params$offset_rate
  
  parcel_num = length(offset_object$current_parcels)
  for (parcel_ind in 1:parcel_num){
    current_parcel = current_parcels[[parcel_ind]]
    predicted_traj = offset_rule_vectorised(current_ecology[current_parcel], params, offset_rate, time_remaining)
    for (pix_ind in 1:length(current_parcel)){
      loc = ind2sub(params$ecology_size, current_parcel[pix_ind])
      trajectories[loc[1], loc[2], yr:params$time_steps] = predicted_traj[, pix_ind]        
    }
  }
  return(trajectories)
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
    plot(dots[[1]][s, ], type = 'l', col = 'red', ylim = c(ymin, ymax))
    if (plot_num > 1){
      for (t in 2:plot_num){
        lines(dots[[t]][s, ],  ylim = c(ymin, ymax))
      }
    } 
  }
}


determine_devs <- function(development_object, params, parcels){
  
  outputs = list()
  parcel_list = development_object$parcel_list
  development_losses = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  static_losses = development_losses
  dev_initial_rel_to_counter = development_losses
  current_parcels = development_object$current_parcels
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  for (yr in 1:params$time_steps){
    
    for (parcel_ind in (parcel_list[[yr]])){
      current_parcel = land_parcels[[parcel_ind]]
      for (pix_ind in current_parcel){
        loc = ind2sub(params$ecology_size, pix_ind)
        current_val = counterfactuals$elements[loc[1], loc[2], yr]
        development_losses[loc[1], loc[2], yr:params$time_steps] = counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]
        static_losses[loc[1], loc[2], yr:params$time_steps] = current_val
        dev_initial_rel_to_counter[loc[1], loc[2], yr:params$time_steps] = current_val - counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]
      }
    }
  }
  outputs$development_losses = development_losses
  outputs$static_losses = static_losses
  outputs$dev_initial_rel_to_counter = dev_initial_rel_to_counter
  return(outputs)
}


determine_offs <- function(object, trajectories, params, parcels){
  
  outputs = list()
  parcel_list = object$parcel_list
  offset_trajectories = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  true_offset_gains = offset_trajectories
  offset_gain_rel_to_current = offset_trajectories
  avoided_degredation = offset_trajectories
  current_parcels = object$current_parcels
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  
  for (yr in 1:params$time_steps){
    
    for (parcel_ind in (parcel_list[[yr]])){
      current_parcel = land_parcels[[parcel_ind]]
      for (pix_ind in current_parcel){
        loc = ind2sub(params$ecology_size, pix_ind)
        current_val = trajectories[loc[1], loc[2], yr]
        current_trajectory = trajectories[loc[1], loc[2], yr:params$time_steps]
        offset_trajectories[loc[1], loc[2], yr:params$time_steps] = current_trajectory
        true_offset_gains[loc[1], loc[2], yr:params$time_steps] = current_trajectory - counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]
        offset_gain_rel_to_current[loc[1], loc[2], yr:params$time_steps] = current_trajectory - current_val
        avoided_degredation[loc[1], loc[2], yr:params$time_steps] = current_val - counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]

      }
    }
  }
  
  outputs$offset_trajectories = offset_trajectories
  outputs$true_offset_gains = true_offset_gains
  outputs$offset_gain_rel_to_current = offset_gain_rel_to_current
  outputs$avoided_degredation = avoided_degredation
  return(outputs)
  
}


calc_diffs <- function(object, trajectories, params, parcels, counterfactuals){
  
  outputs = list()
  parcel_list = object$parcel_list
  true_change = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  traj_rel_to_initial = true_change
  initial_rel_to_counter = true_change
  current_parcels = object$current_parcels
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



find_parcel_vals <- function(parcels, current_ecology, parcel_indexes){
  
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

find_current_parcel_value <- function(parcels, parcel_ind, current_ecology){
  inds = parcels$land_parcels[[parcel_ind]]
  parcel_val = sum(current_ecology[inds])
  return(parcel_val)
}

find_predicted_parcel_value <- function(parcels, parcel_ind, trajectories){
  inds = parcels$land_parcels[[parcel_ind]]
  final_ecology = trajectories[, , params$time_steps]
  parcel_val = sum(final_ecology[inds])
  return(parcel_val)
}


sum_parcels <- function(time_steps, parcels, trajectories){
  
  land_parcel_num = length(parcels$land_parcels)
  land_parcels = parcels$land_parcels
  parcel_sums = array(0, c(land_parcel_num, time_steps))
  
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    
    for (parcel_ind in 1:land_parcel_num){   
      current_parcel = land_parcels[[parcel_ind]]
      parcel_sums[parcel_ind, yr] = sum(current_slice[current_parcel])
    }
  }

return(parcel_sums)

}


sum_regions <- function(time_steps, parcel_sums, parcels){
  
  region_num = parcels$region_num
  regions = parcels$regions
  region_sums = array(0, c(region_num, time_steps))
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = parcel_sums[, yr]
    for (region_ind in 1:region_num){ #number of regions  
      current_region = regions[[region_ind]]
      region_sums[region_ind, yr] = sum(current_slice[current_region])    
    }
  }
  return(region_sums)
}



find_sums <- function(time_steps, parcels, trajectories){
  object = list()
  object$parcel_sums = sum_parcels(time_steps, parcels, trajectories)
  object$region_sums = sum_regions(time_steps, object$parcel_sums, parcels)
  object$net_sums = colSums(object$region_sums)
  return(object)
}
  

find_offset_pool <- function(index_object, region_ind, params, parcels, current_ecology, time_remaining, final_ecology){
  
  offset_pool = list()
  offset_ind_to_use = index_object$ind_available[[region_ind]]
  
  if (params$offset_selection_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
    current_parcel_vals <- predict_parcels_vectorised(region_ind, parcels, offset_ind_to_use, current_ecology, params, time_remaining, final_ecology)
  } else if (params$offset_selection_type == 'current'){
    current_parcel_vals <- find_parcel_vals(parcels, current_ecology, offset_ind_to_use)
  }
  
  offset_pool$offset_ind_to_use = offset_ind_to_use
  offset_pool$current_parcel_vals = current_parcel_vals
  return(offset_pool)
}



predict_parcels <- function(region_ind, parcels, ind_available, current_ecology, params, time_remaining, final_ecology){
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(ind_available))
  
  for (s in 1:length(ind_available)){
    current_ind = ind_available[s]
    current_parcel = land_parcels[[current_ind]]
    predicted_offset_parcel = array(0, length(current_parcel))
    parcel_pix_ind = 1
    for (pix_ind in current_parcel){
      current_val = current_ecology[pix_ind]
      predicted_pixel_trajectory = offset_rule(current_val, params, offset_rate, time_remaining)
      if (params$predict_offset_mode == 'gains'){ 
        counterfactual_val = final_ecology[pix_ind]
        predicted_offset_parcel[parcel_pix_ind] = predicted_pixel_trajectory[length(predicted_pixel_trajectory)] - counterfactual_val
      } else if (params$predict_offset_mode == 'saved'){
        predicted_offset_parcel[parcel_pix_ind] = predicted_pixel_trajectory[length(predicted_pixel_trajectory)]
      }
      parcel_pix_ind <- parcel_pix_ind + 1
    }     
    predicted_parcel_vals[s] = sum(predicted_offset_parcel)
  }
  return(predicted_parcel_vals)
  
}

predict_parcels_vectorised <- function(region_ind, parcels, ind_available, current_ecology, params, time_remaining, final_ecology){
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(ind_available))
  
  for (s in 1:length(ind_available)){
    current_ind = ind_available[s]
    current_parcel = land_parcels[[current_ind]]
    predicted_trajs = offset_rule_vectorised(current_ecology[current_parcel], params, offset_rate, time_remaining)
    predicted_parcel_vals[s] = sum(predicted_trajs[nrow(predicted_trajs), ])
  }
  
  return(predicted_parcel_vals)
  
}

parcel_trajectory <- function(parcels, parcel_ind, params, trajectories){
  land_parcels = parcels$land_parcels
  current_parcel = land_parcels[[parcel_ind]]
  parcel_traj = array(0, params$time_steps)
  
  for (yr in 1:params$time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    parcel_traj[yr] = sum(current_slice[current_parcel])
  }
  return(parcel_traj)
}


write_parcel_sets <- function(dev_counter, parcel_sets, yr, region_ind, dev_index, dev_parcel, offset_parcels){  
  parcel_set = list()
  parcel_set$yr = yr
  parcel_set$region = region_ind
  parcel_set$dev_parcel = dev_parcel
  parcel_set$offset_parcels = offset_parcels
  #parcel_sets = c(parcel_sets, parcel_set)
  parcel_sets[[dev_counter]] = parcel_set
  #print(c(parcel_sets[[dev_counter]], yr, dev_counter))
  #parcel_sets$dev_counter[[1]] = parcel_sets$dev_counter[[1]] + 1
  
  return(parcel_sets)
}

view_individual_offsets <- function(parcel_set, parcels, params, trajectories, counterfactuals){
  assess_object = list()
  
  offset_yr = parcel_set$yr
  dev_parcel = parcel_set$dev_parcel
  offset_parcels = parcel_set$offset_parcels
  dev_trajectory = parcel_trajectory(parcels, dev_parcel, params, trajectories)
  counterfactual_trajectory = parcel_trajectory(parcels, dev_parcel, params, counterfactuals$elements)
  offset_trajectories = array(0, c(length(offset_parcels), params$time_steps))
  for (parcel_ind in 1:length(offset_parcels)){
    offset_trajectories[parcel_ind, ] = parcel_trajectory(parcels, offset_parcels[parcel_ind], params, trajectories)
  }
  if (length(offset_parcels) > 1){
    offset_sum = colSums(offset_trajectories)
  } else (offset_sum = offset_trajectories[1, ])
  
  assess_object$dev_trajectory = dev_trajectory
  assess_object$offset_trajectories = offset_trajectories
  
  graphics.off()
  mx = max(c(dev_trajectory, offset_sum))
  plot(assess_object$dev_trajectory, type = 'l', ylim = c(0, mx))
  lines((offset_sum - offset_sum[offset_yr]), col = 'red', ylim = c(0, mx))
  lines((offset_sum), col = 'red', ylim = c(0, mx))
  lines(counterfactual_trajectory, ylim = c(0, mx))
  return(assess_object)
}

calc_trajectories <- function(params, outputs, parcels, index_object, counterfactuals){
  dev_counter = 1;
  for (yr in 1:params$time_steps){      #main time loop    
    time_remaining = 0:(params$time_steps - yr) 
    current_ecology = outputs$trajectories[ , , yr]
    
    #if (yr%%params$update_every == 0 && yr < params$time_steps){
    if (yr == 4){
      develop_nums = params$development_rates
      #for (region_ind in 1:parcels$region_num){
      for (region_ind in 1:1){
        develop_num = develop_nums[region_ind]
        
        if (develop_num > 0){
          for (dev_index in 1:develop_num){    #loop over the number of possible developments           
            if (length(index_object$ind_available[[region_ind]]) < 1){
              break
            }

            outputs$development_object <- select_parcel_to_develop(params, outputs$development_object, index_object$ind_available[[region_ind]], parcels, yr, current_ecology, counterfactuals$final_ecology)
            outputs$trajectories <- write_development_trajectory(outputs$trajectories, counterfactuals, outputs$development_object, current_ecology, params, yr)     
            index_object <- update_index(index_object, outputs$development_object$parcel_index, region_ind, params)      
            if (index_object$break_flag == TRUE){
              break
            } 
            
            offset_pool <- find_offset_pool(index_object, region_ind, params, parcels, current_ecology, time_remaining, counterfactuals$final_ecology)
            outputs$offset_object <- select_offset_parcel(offset_pool, outputs$offset_object, params, parcels, outputs$development_object$current_parcel_val, yr)                  
            if (outputs$offset_object$break_flag == TRUE){
              break
            }     
            
            outputs$trajectories = write_offset_trajectory_vectorised(outputs$trajectories, outputs$offset_object, current_ecology, params, time_remaining, yr)         
            index_object = update_index(index_object, outputs$offset_object$current_parcel_indexes, region_ind, params)           
            if (index_object$break_flag == TRUE){
              break
            }
            outputs$parcel_sets <- write_parcel_sets(dev_counter, outputs$parcel_sets, yr, region_ind, dev_index, outputs$development_object$parcel_index, outputs$offset_object$current_parcel_indexes)
            dev_counter = dev_counter + 1;
            }
          #image(outputs$trajectories[ , , yr],zlim = c(0, 100))
        } 
      }
    }
   print(yr)
  }
  return(outputs)
}





graphics.off()
params <- initialise_parameters()
parcels <- initialise_shape_parcels(params)
index_object <- initialise_index_object(parcels)
counterfactuals <- initialise_counterfactuals(params, parcels)
outputs <- initialise_outputs(params, counterfactuals, parcels)

outputs <- calc_trajectories(params, outputs, parcels, index_object, counterfactuals)

offset_outs = calc_diffs(outputs$offset_object, outputs$trajectories, params, parcels, counterfactuals)
development_outs = calc_diffs(outputs$development_object, outputs$trajectories, params, parcels, counterfactuals)

net_outs = list();
net_outs$traj_rel_to_initial = development_outs$traj_rel_to_initial + offset_outs$traj_rel_to_initial
net_outs$initial_rel_to_counter = development_outs$initial_rel_to_counter + offset_outs$initial_rel_to_counter
net_outs$true_change = development_outs$true_change + offset_outs$true_change

trajectory_sums <- find_sums(params$time_steps, parcels, outputs$trajectories)
offset_absolute_sums <- find_sums(params$time_steps, parcels, offset_outs$offset_trajectories)
counterfactual_sums <- find_sums(params$time_steps, parcels, counterfactuals$elements)

net_impacts <- find_sums(params$time_steps, parcels, (outputs$trajectories - counterfactuals$elements))

true_offset_gains <- find_sums(params$time_steps, parcels, offset_outs$true_change)
offset_traj_rel_to_initial <- find_sums(params$time_steps, parcels, offset_outs$traj_rel_to_initial)
avoided_degradation <- find_sums(params$time_steps, parcels, offset_outs$initial_rel_to_counter)

development_true_losses <- find_sums(params$time_steps, parcels, development_outs$true_change)
development_traj_rel_to_initial <- find_sums(params$time_steps, parcels, development_outs$traj_rel_to_initial)
dev_initial_rel_to_counter <- find_sums(params$time_steps, parcels, development_outs$initial_rel_to_counter)

net_true_change <- find_sums(params$time_steps, parcels, net_outs$true_change)
net_traj_rel_to_initial <- find_sums(params$time_steps, parcels, net_outs$traj_rel_to_initial)
net_initial_rel_to_current <- find_sums(params$time_steps, parcels, net_outs$initial_rel_to_counter)
 
plot_params = c(params$region_num_x, params$region_num_y, parcels$region_num)

#plot_outs(true_offset_gains$region_sums[1:6, ], avoided_degradation$region_sums[1:6, ], offset_traj_rel_to_initial$region_sums[1:6, ], plot_params)

#plot_outs(net_true_change$region_sums[1:6, ], net_traj_rel_to_initial$region_sums[1:6, ], net_initial_rel_to_current$region_sums[1:6, ], plot_params)

#plot_outs(net_true_change$parcel_sums[1:6, ], net_traj_rel_to_initial$parcel_sums[1:6, ], net_initial_rel_to_current$parcel_sums[1:6, ], plot_params)

a <- view_individual_offsets(outputs$parcel_sets[[1]], parcels, params, outputs$trajectories, counterfactuals)

#plot_outs(true_offset_gains$parcel_sums[1:6, ], avoided_degradation$parcel_sums[1:6, ], offset_traj_rel_to_initial$parcel_sums[1:6, ], plot_params)

# 
# parcel_sums = sum_parcels(params$time_steps, parcels, outputs$trajectories)
# parcel_avoided_degradation = sum_parcels(params$time_steps, parcels, offset_outs$avoided_degredation)
# parcel_offset_gain_rel_to_current = sum_parcels(params$time_steps, parcels, offset_outs$offset_gain_rel_to_current)
# parcel_true_offset_gains = sum_parcels(params$time_steps, parcels, offset_outs$true_offset_gains)








# parcel_sum <- function(time_steps, parcels, parcel_ind, trajectories){
#   
#   land_parcels = parcels$land_parcels
#   parcel_trajectory = array(0, time_steps)
#   
#   current_parcel = land_parcels[[parcel_ind]]
#   
#   for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
#     current_trajectory = trajectories[ , , yr]
#     parcel_trajectory[yr] = sum(current_trajectory[current_parcel])
#   }
#   return(parcel_sums)
# }










#par(mfrow = c(2, 3))
#par(cex = 0.6)
#par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
#par(tcl = -0.25)
#par(mgp = c(2, 0.6, 0))
#for (i in 1:6) {
#plot(1, axes = FALSE, type = "n")
#mtext(letters[i], side = 3, line = -1, adj = 0.1, cex = 0.6,
#col = "grey40")
#if (i %in% c(4, 5, 6))
#axis(1, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
#if (i %in% c(1, 4))
#    axis(2, col = "grey40", col.axis = "grey20", at = seq(0.6,1.2, 0.2))
#box(col = "grey60")
#}
#mtext("x axis", side = 1, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")
#mtext("y axis", side = 2, outer = TRUE, cex = 0.7, line = 2.2,col = "grey20")



# 
# initialise_square_parcels <- function(params){
#   parcels = list()
#   pixel_indexes = 1:(params$ecology_size*params$ecology_size)
#   dim(pixel_indexes) = c(params$ecology_size, params$ecology_size)
#   land_parcels = matcell(pixel_indexes, c(params$parcel_num_y, params$parcel_num_x)) #split the ecology array into a series of subarrays with dimensions sz_x by sz_y
#   parcel_num_y = land_parcels$dims[1]
#   parcel_num_x = land_parcels$dims[2]
#   land_parcel_num = length(land_parcels$elements) #total number of parcels
#   parcel_indexes = 1:land_parcel_num #index all parcels
#   dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
#   region_size_x = ceiling(parcel_num_x/params$region_num_x)
#   region_size_y = ceiling(parcel_num_y/params$region_num_y)
#   regions = matcell(parcel_indexes, c(region_size_y, region_size_x))
#   region_num = length(regions$elements)
#   parcels$land_parcel_num = land_parcel_num
#   parcels$land_parcels = land_parcels$elements
#   parcels$land_parcel_dims = land_parcels$dims
#   parcels$regions = regions$elements
#   parcels$region_dims = regions$dims
#   parcels$region_num = region_num
#   return(parcels)
# }
# 
# 
# 
# 
# find_region <- function(parcels, parcel_ind){
#   region_num = parcels$region_num
#   regions = parcels$regions
#   for (s in 1:region_num){
#     if (any(regions[[s]] == parcel_ind)){
#       region_index = s
#       break
#     }
#     
#   }
#   return(region_index)
# }
# 
# write_trajectory <- function(outputs, current_parcel, current_ecology, params, time_remaining, yr){
#   
#   min_eco_val = params$min_eco_val
#   max_eco_val = params$max_eco_val
#   offset_rate = params$offset_rate
#   dim(current_parcel) = c(1, length(current_parcel))
#   
#   for (pix_ind in current_parcel){
#     loc = ind2sub(params$ecology_size, pix_ind)  
#     current_trajectory = offset_rule(current_ecology[pix_ind], params, offset_rate, time_remaining)
#     outputs$avoided_losses[loc[1], loc[2], yr:params$time_steps] = current_ecology[pix_ind]
#     counterfactual = counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]            
#   }
#   return(outputs)
#   
# }
# 
# 
# write_current_value <- function(A, current_parcel, current_ecology, params, yr){
#   current_array = matrix(0, params$ecology_size, params$ecology_size)
#   current_array[current_parcel] = current_ecology[current_parcel]
#   A[, , yr] = A[, , yr] + current_array
#   return(A)
# }
# 
# plot_outputs <- function(sums1, sums2, params){
#   
#   graphics.off()
#   A = 1:parcels$region_num
#   dim(A) = c(params$region_num_x, params$region_num_y)
#   layout(A)
#   ymax = max(cbind(sums1, sums2))
#   ymin = min(cbind(sums1, sums2))
#   
#   for (s in 1:params$region_num){
#     plot(sums1[ , s], type = 'l', ylim = c(ymin, ymax))
#     lines(sums2[ , s], col = 'red', ylim = c(ymin, ymax))
#   }
#   
#   #mtext(plot_type, side = 3, line = -2, outer = TRUE)
# }
# 
# matcell <- function(x, tilesize){
#   
#   sz = c(nrow(x), ncol(x))
#   A = vector('list', 2)
#   for (s in 1:2){
#     dm = sz[s]
#     T = min(dm, tilesize[s])
#     nn = floor( dm / T ) 
#     resid = dm %% tilesize[s]
#     if (resid == 0 ) {
#       A[[s]]=c(rep(1, nn)*T)
#     } else A[[s]]=c(rep(1, nn)*T,resid)
#   }
#   
#   rowsizes = A[[1]]
#   colsizes = A[[2]]
#   rows = length(rowsizes)
#   cols = length(colsizes)
#   B = vector('list', rows*cols)
#   
#   rowStart = 0
#   
#   a = 1  
#   for (i in 1:rows){
#     colStart = 0
#     for (j in 1:cols){
#       B[[a]] = x[colStart+(1:colsizes[j]), rowStart+(1:rowsizes[i])]
#       colStart = colStart + colsizes[j]
#       a = a + 1
#     }
#     rowStart = rowStart + rowsizes[i]
#   }
#   
#   parcel = list()
#   parcel$dims = c(length(A[[1]]), length(A[[2]]))
#   parcel$elements = B
#   return(parcel)
# }
# 


# sum_regions <- function(time_steps, parcels, trajectories){
#   
#   region_num = parcels$region_num
#   land_parcels = parcels$land_parcels
#   regions = parcels$regions
#   region_sums = array(0, c(time_steps, region_num))
#   
#   for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
#     current_slice = trajectories[ , , yr]
#     
#     for (region_ind in 1:region_num){ #number of regions  
#       
#       current_region = regions[[region_ind]]
#       current_parcel_num = length(current_region)
#       current_region_sum = 0
#       
#       for (parcel_ind in current_region){ 
#         current_parcel = land_parcels[[parcel_ind]]
#         current_region_sum = current_region_sum + sum(current_slice[current_parcel])
#       }
#       region_sums[yr, region_ind] = current_region_sum
#     }
#   }
#   return(region_sums)
# }

# select_offset_pool <- function(region_ind, offset_object, params, parcels, development_object, index_object, current_parcel_vals, predicted_parcel_vals, predict_flag){
#   
#   ind_available = index_object$ind_available
#   
#   region_num = parcels$region_num
#   regions = parcels$regions
#   parcel_selection_type = params$parcel_selection_type 
#   offset_selection_type = params$offset_selection_type
#   offset_value_type = params$offset_value_type
#   development_parcel_ind = development_object$parcel_index
#   development_parcel_val = development_object$parcel_val
#   
#   #if (parcel_selection_type == 'regional'){        
#   #  for (s in 1:region_num){
#   #    if (any(regions[[s]] == development_parcel_ind)) {
#   #      offset_ind_to_use = intersect(regions[[s]], ind_available)                                        
#   #      break
#   #    }
#   #  }
#   #} else if (parcel_selection_type == 'national'){ 
#   #  offset_ind_to_use = ind_available
#   #}
#   
#   if (length(offset_ind_to_use) < 1 ){   
#     offset_object$break_flag = TRUE
#     return(offset_object)
#   } else {offset_object$break_flag = FALSE}
#   
#   return(offset_ind_to_use)
# }


# determine_offset_index <- function(params, parcels, development_object, ind_available){
#   region_num = parcels$region_num
#   regions = parcels$regions
#   if (params$parcel_selection_type == 'regional') {
#     for (s in 1:region_num){
#       if ( any(regions[[s]] == development_object$parcel_ind )) {
#         offset_ind_to_use = intersect(regions[[s]], ind_available)
#         break
#       }
#     }
#   } else if (params$parcel_selection_type == 'national'){
#     offset_ind_to_use = ind_available
#   }
#   return(offset_ind_to_use)  
# }
