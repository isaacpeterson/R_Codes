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


initialise_shape_parcels <- function(global_params){
  parcels = list()
  parcel_num_x = global_params$parcel_num_x
  parcel_num_y = global_params$parcel_num_y
  parcel_vx = split_vector(parcel_num_x, global_params$ecology_size, 5)
  parcel_vy = split_vector(parcel_num_y, global_params$ecology_size, 5)
  
  pixel_indexes = 1:(global_params$ecology_size*global_params$ecology_size)
  dim(pixel_indexes) = c(global_params$ecology_size, global_params$ecology_size)
  land_parcels = mcell(pixel_indexes, parcel_vx, parcel_vy) #lit the ecology array into a series of subarrays with dimensions sz_x by sz_y
  land_parcel_num = length(land_parcels$elements) #total number of parcels
  parcel_indexes = 1:land_parcel_num #index all parcels
  dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
  region_vx = split_vector(global_params$region_num_x, parcel_num_x, 1) 
  region_vy = split_vector(global_params$region_num_y, parcel_num_y, 1)
  
  regions = mcell(parcel_indexes, region_vx, region_vy)
  
  region_num = length(regions$elements)
  parcels$parcel_indexes = parcel_indexes
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


initialise_index_object <- function(parcels, global_params){
  index_object = list()
  index_object$ind_available = parcels$regions
  index_object$developments = vector()
  index_object$offsets = vector()
  index_object$parcel_sets = list()
  index_object$dev_count = 1
  index_object$break_flag = FALSE
  return(index_object)
}

initialise_ecology <- function(global_params, parcels){
  land_parcels = parcels$land_parcels
  land_parcel_num = parcels$land_parcel_num
  initial_ecology = matrix(1,global_params$ecology_size,global_params$ecology_size)
  
  for (s in 1:land_parcel_num){
    initial_parcel_value = global_params$min_initial_eco_val + (global_params$max_initial_eco_val - global_params$min_initial_eco_val - global_params$initial_eco_noise)*runif(1) 
    inds = parcels$land_parcels[[s]]
    dim(inds) = c(1, length(inds))
    initial_ecology[inds] = initial_ecology[inds]*initial_parcel_value
  }
  initial_ecology = initial_ecology + global_params$initial_eco_noise*matrix(runif(global_params$ecology_size*global_params$ecology_size),global_params$ecology_size,global_params$ecology_size)
  return(initial_ecology)
}

initialise_offset_object <- function(){
  object = list()
  object$current_parcels = list()
  object$current_parcel_sums = list()
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


eco_change <- function(parcel_vals, min_eco_val, max_eco_val, decline_rate, time_step){
  
  t_sh = -1/decline_rate * log( ((parcel_vals - min_eco_val)/(max_eco_val - parcel_vals)))
  eco_shift = min_eco_val + (max_eco_val - min_eco_val)/(1 + exp(-decline_rate*(time_step - t_sh)))
  return(eco_shift)
}


build_counterfactuals_sliced <- function(region_num, regions, ecology_size, time_steps, decline_rates, land_parcels, initial_ecology){
  current_ecology = initial_ecology
  counterfactuals = array(0, c(ecology_size, ecology_size, time_steps))
  counterfactuals[, , 1] = initial_ecology
  parcel_num = length(land_parcels)
  
  for (yr in 1:time_steps){
    for (region_ind in 1:region_num){ 
      current_parcel_pool = regions[[region_ind]]  
      current_decline_rates = decline_rates[[region_ind]]
      for (parcel_ind in seq_along(current_parcel_pool)){
        current_parcel_ind = current_parcel_pool[parcel_ind]
        current_parcel = land_parcels[[current_parcel_ind]]
        current_dec_rate = current_decline_rates[parcel_ind]      
        updated_parcel <- sapply(current_ecology[current_parcel], eco_change, min_eco_val = 0, max_eco_val = 100, decline_rate = current_dec_rate, time_step = 1)
        current_ecology[current_parcel] = updated_parcel 
      }
    }
    #current_ecology = Blur_2D(current_ecology, 0.5, 0.5)  
    counterfactuals[, , yr] = current_ecology
    print(yr)
  }
  
  return(counterfactuals)
}


build_counterfactuals_by_parcel <- function(ecology_size, time_steps, decline_rates, land_parcels, initial_ecology){
  current_ecology = initial_ecology
  counterfactuals = array(0, c(ecology_size, ecology_size, time_steps))
  counterfactuals[, , 1] = initial_ecology
  parcel_num = length(land_parcels)
  
  for (yr in 1:time_steps){
    for (parcel_ind in 1:parcel_num){
      current_parcel = land_parcels[[parcel_ind]]
      current_dec_rate = decline_rates[parcel_ind]      
      updated_parcel <- sapply(current_ecology[current_parcel], eco_change, min_eco_val = 0, max_eco_val = 100, decline_rate = current_dec_rate, time_step = 1)
      current_ecology[current_parcel] = updated_parcel 
    }
    #current_ecology = Blur_2D(current_ecology, 0.5, 0.5)  
    counterfactuals[, , yr] = current_ecology
    print(yr)
  }
  
  return(counterfactuals)
}

initialise_counterfactuals_vectorised <- function(global_params, region_params, parcels, initial_ecology, decline_rates){
  counterfactuals_object = list() 
  land_parcels = parcels$land_parcels
  time_steps = global_params$time_steps 
  ecology_size = global_params$ecology_size
  counterfactuals = build_counterfactuals_by_parcel(ecology_size, time_steps, decline_rates$decline_rates_array, land_parcels, initial_ecology)
  final_ecology = counterfactuals[, , global_params$time_steps]
  final_parcel_vals = find_current_parcel_sums(parcels, final_ecology, 1:length(land_parcels))
  
  counterfactuals_object$final_parcel_vals = final_parcel_vals
  counterfactuals_object$final_ecology = final_ecology
  counterfactuals_object$counterfactuals = counterfactuals
  counterfactuals_object$initial_decline_rates = decline_rates
  return(counterfactuals_object)
}

build_decline_rates <- function(parcels, region_params){
  regions = parcels$regions
  region_num = length(regions)
  dec_rates = vector('list', region_num)
  decline_rates = list()
  decline_rates_array = array(0, dim(parcels$parcel_indexes))
  
  for (region_ind in 1:region_num){ 
    current_region = regions[[region_ind]]
    current_parcel_num = length(current_region)    
    decline_params = c(length(current_region), region_params[[region_ind]]$mean_decline_rate, region_params[[region_ind]]$decline_rate_std) #params$decline_rate_std[region_ind])
    current_decline_rates = matrix(rnorm(decline_params[1], mean = decline_params[2], sd = decline_params[3]), ncol = ncol(current_region))
    dec_rates[[region_ind]] = current_decline_rates
    decline_rates_array[current_region] = current_decline_rates
  }
  decline_rates$dec_rates = dec_rates
  decline_rates$decline_rates_array = decline_rates_array
  return(decline_rates)
}


initialise_outputs <- function(counterfactuals, parcels, global_params){
  parcel_num_x = length(parcels$parcel_vx)
  parcel_num_y = length(parcels$parcel_vy)
  outputs = list()
  outputs$trajectories = array(0, c(global_params$ecology_size, global_params$ecology_size, global_params$time_steps))
  outputs$adjusted_counterfactuals = counterfactuals_object$counterfactuals
  return(outputs)
  
}


update_development_object <- function(dev_score_type, ind_available, parcels, yr, current_ecology, final_ecology){
  development_object = list()
  land_parcels = parcels$land_parcels
  current_parcel_index = ind_available[sample(1:length(ind_available), 1)]
  current_parcel = land_parcels[[current_parcel_index]]
  current_parcel_val = sum(current_ecology[current_parcel])
  
  if (dev_score_type == 'current'){
    parcel_vals_to_use = current_parcel_val
  } else if (dev_score_type == 'predicted'){
    parcel_vals_to_use = sum(final_ecology[current_parcel])
  }
  
  current_parcel_ecology = current_ecology[current_parcel]
  dim(current_parcel_ecology) = dim(current_parcel)
  
  development_object$current_parcel_indexes = current_parcel_index
  development_object$current_parcel_sums = current_parcel_val
  development_object$parcel_vals_to_use = parcel_vals_to_use
  development_object$current_parcels = current_parcel
  development_object$current_parcel_ecology = current_parcel_ecology
  return(development_object)
}







find_offset_pool <- function(index_object, region_ind, min_eco_val, max_eco_val, time_steps, region_params, parcels, current_ecology, time_remaining, counterfactuals){
  
  offset_pool = list()
  current_offset_pool = index_object$ind_available[[region_ind]]
  current_parcel_sums <- find_current_parcel_sums(parcels, current_ecology, current_offset_pool)
  
  if (region_params[[region_ind]]$offset_selection_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
    parcel_vals_to_use <- predict_final_parcel_vals(parcels, current_offset_pool, current_ecology, min_eco_val, max_eco_val, region_params[[region_ind]]$offset_rate, time_remaining)
    if (region_params[[region_ind]]$predict_offset_mode == 'gains'){
      predicted_slice = counterfactuals[, , time_steps]
      parcel_vals_to_use = parcel_vals_to_use - counterfactuals_object$final_parcel_vals[current_offset_pool]
    }
    
  } else if (region_params[[region_ind]]$offset_selection_type == 'current'){
    parcel_vals_to_use = current_parcel_sums
  }
  
  offset_pool$current_parcel_pool = current_offset_pool
  offset_pool$current_parcel_sums = current_parcel_sums
  offset_pool$parcel_vals_to_use = parcel_vals_to_use
  
  return(offset_pool)
}                 

update_offset_object <- function(offset_pool, offset_multiplier, parcels, development_val, yr, current_ecology){      
  
  offset_object = list()
  current_offset_pool = offset_pool$current_parcel_pool
  parcel_vals_to_use = offset_pool$parcel_vals_to_use  
  land_parcels = parcels$land_parcels
  val_to_match = offset_multiplier*development_val
  
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
      parcel_indexes = c(parcel_indexes, current_offset_pool[max_ind])
      parcel_vals_used = c(parcel_vals_used, parcel_vals_to_use[max_ind])
      break
    }
  }
  
  parcel_num = length(parcel_indexes)
  offset_object$current_parcels = list()
  offset_object$current_parcel_sums = array(0, parcel_num)
  offset_object$current_parcel_ecology = vector('list', parcel_num)
  
  for (parcel_ind in 1:parcel_num){
    current_parcel = land_parcels[[parcel_indexes[parcel_ind]]]
    offset_object$current_parcels[[parcel_ind]] = current_parcel
    current_parcel_ecology = current_ecology[current_parcel]
    dim(current_parcel_ecology) = dim(current_parcel)
    offset_object$current_parcel_ecology[[parcel_ind]] = current_parcel_ecology
    offset_object$current_parcel_sums[parcel_ind] = sum(current_ecology[current_parcel])
  }  
  offset_object$current_parcel_indexes = parcel_indexes
  offset_object$parcel_vals_used = parcel_vals_used
  return(offset_object)  
}


write_null_offset_object <- function(){
  offset_object = list()
  offset_object$current_parcel_indexes = list()
  offset_object$current_parcel_sums = list()
  offset_object$parcel_vals_used = list()
  return(offset_object)
}


write_development <- function(development_object, current_ecology){
  
  current_parcel = development_object$current_parcels
  current_ecology[current_parcel] = 0
  
  return(current_ecology)
  
}

write_offset <- function(offset_object, current_ecology, min_eco_val, max_eco_val, ecology_size, offset_rate, yr){
  current_parcels = offset_object$current_parcels
  parcel_num = length(offset_object$current_parcels)
  for (parcel_ind in 1:parcel_num){
    current_parcel = current_parcels[[parcel_ind]]
    updated_parcel <- sapply(current_ecology[current_parcel], eco_change, min_eco_val = 0, max_eco_val = 100, decline_rate = offset_rate, time_step = 1)
    current_ecology[current_parcel] = updated_parcel 
  }
  return(trajectories)
}


find_current_parcel_sums <- function(parcels, current_ecology, parcel_indexes){
  
  land_parcels = parcels$land_parcels
  parcel_num = length(parcel_indexes)
  current_parcel_sums = rep(0, parcel_num)
  
  for (s in 1:parcel_num){
    parcel_ind = parcel_indexes[s]
    current_parcel = land_parcels[[parcel_ind]]
    current_parcel_sums[s] = sum(current_ecology[current_parcel])
  }
  return(current_parcel_sums)
}


predict_parcels_vectorised <- function(region_ind, parcels, ind_available, current_ecology, global_params, time_remaining, final_ecology){
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(ind_available))
  
  for (s in 1:length(ind_available)){
    current_ind = ind_available[s]
    current_parcel = land_parcels[[current_ind]]
    predicted_trajs = ecology_change_function_vectorised(current_ecology[current_parcel], global_params$min_eco_val, global_params$max_eco_val, offset_rate, time_remaining)
    predicted_parcel_vals[s] = sum(predicted_trajs[nrow(predicted_trajs), ])
  }
  
  return(predicted_parcel_vals)
  
}




find_parcel_trajectory <- function(parcels, parcel_ind, time_steps, trajectories){
  land_parcels = parcels$land_parcels
  current_parcel = land_parcels[[parcel_ind]]
  parcel_traj = array(0, time_steps)
  
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    parcel_traj[yr] = sum(current_slice[current_parcel])
  }
  return(parcel_traj)
}


find_parcel_trajectories <- function(parcels, parcel_indexes, time_steps, trajectories){
  land_parcels = parcels$land_parcels
  parcel_trajs = array(0, c(length(parcel_indexes), time_steps))
  
  for (parcel_ind in 1:length(parcel_indexes)){
    current_ind = parcel_indexes[parcel_ind]
    parcel_trajs[parcel_ind, ] = find_parcel_trajectory(parcels, current_ind, time_steps, trajectories)
  }
  
  if (length(parcel_indexes) == 1){
    dim(parcel_trajs) = c(length(parcel_trajs), 1)
  }
  
  return(parcel_trajs)
}


write_developments_offsets <- function(index_object, development_object, offset_object, region_ind, yr){
  index_object$developments = c(index_object$developments, development_object$current_parcel_indexes)
  index_object$offsets = c(index_object$offsets, offset_object$current_parcel_indexes)
  index_object$parcel_sets[[index_object$dev_count]] <- write_parcel_sets(development_object, offset_object, yr, region_ind)
  index_object$dev_count <- index_object$dev_count + 1
  return(index_object)
}

write_parcel_sets <- function(development_object, offset_object, yr, region_ind){
  parcel_set = list()
  parcel_set$yr = yr
  parcel_set$development_object = development_object
  parcel_set$offset_object = offset_object
  return(parcel_set)
}


# write_parcel_sets <- function(development_object, offset_object, yr, region_ind){
#   parcel_set = list()
#   parcel_set$yr = yr
#   parcel_set$region = region_ind
#   parcel_set$developed_parcel = development_object$current_parcel_indexes
#   parcel_set$dev_val = development_object$current_parcel_sums
#   parcel_set$offset_parcels = offset_object$current_parcel_indexes
#   parcel_set$offset_vals = offset_object$current_parcel_sums
#   return(parcel_set)
# }


update_ind_available <- function(index_object, current_parcel_indexes, region_ind, max_region_dev_num){
  
  ind_available = index_object$ind_available[[region_ind]]
  ind_available <- setdiff(ind_available, current_parcel_indexes) #remove development parcel from available list   
  if (length(ind_available) < max_region_dev_num ){   
    break_flag = TRUE
  } else {break_flag = FALSE}
  
  index_object$ind_available[[region_ind]] = ind_available
  index_object$break_flag = break_flag
  
  return(index_object)
  
}


update_decline_rates <- function(decline_rates, parcel_indexes, offset_rate){
  for (parcel_ind in 1:length(parcel_indexes)){
    current_parcel_ind = parcel_indexes[parcel_ind]
    decline_rates[current_parcel_ind] = offset_rate
  }
  return(decline_rates)
}

find_current_dev_nums <- function(region_params){
  region_num = length(region_params)
  current_dev_nums = array(0, region_num)
  for (region_ind in 1:region_num){
    current_dev_nums[region_ind] = region_params[[region_ind]]$max_region_dev_num
  }
  return(current_dev_nums)
}

calc_trajectories_by_parcel <- function(outputs, global_params, region_params, current_ecology, decline_rates, parcels, index_object, counterfactuals_object, perform_offsets){
  for (yr in 1:global_params$time_steps){      #main time loop    
    time_remaining = 0:(global_params$time_steps - yr)     
    if (yr%%global_params$develop_every == 0 && yr < global_params$time_steps){
      current_dev_nums <- find_current_dev_nums(region_params) 
      for (region_ind in 1:parcels$region_num){       
        current_develop_num = current_dev_nums[region_ind]
        if (current_develop_num > 0){
          for (dev_index in 1:current_develop_num){    #loop over the number of possible developments           
            if (length(index_object$ind_available[[region_ind]]) < 1){
              break
            }           
            development_object <- update_development_object(region_params[[region_ind]]$development_value_type, index_object$ind_available[[region_ind]], parcels, yr, current_ecology, counterfactuals_object$final_ecology)   
            index_object <- update_ind_available(index_object, development_object$current_parcel_indexes, region_ind, region_params[[region_ind]]$max_region_dev_num)      
            decline_rates <- update_decline_rates(decline_rates, development_object$current_parcel_indexes, offset_rate = 0)
            if (index_object$break_flag == TRUE){
              break
            }                      
            if (perform_offsets == TRUE){
              offset_pool <- find_offset_pool(index_object, region_ind, global_params$min_eco_val, global_params$max_eco_val, global_params$time_steps, region_params, parcels, current_ecology, time_remaining, counterfactuals_object$counterfactuals)
              offset_object <- update_offset_object(offset_pool, region_params[[region_ind]]$offset_multiplier, parcels, development_object$parcel_vals_to_use, yr, current_ecology)                  
              if (offset_object$break_flag == TRUE){
                break
              }             
              index_object = update_ind_available(index_object, offset_object$current_parcel_indexes, region_ind, region_params[[region_ind]]$max_region_dev_num)              
              decline_rates <- update_decline_rates(decline_rates, offset_object$current_parcel_indexes, offset_rate = region_params[[region_ind]]$offset_rate)
              if (index_object$break_flag == TRUE){
                break
              }               
            } else {offset_object <- write_null_offset_object()}
            
            index_object <- write_developments_offsets(index_object, development_object, offset_object, region_ind, yr)
          }
        }
        if (global_params$display_object == TRUE){
          image(current_ecology, zlim = c(global_params$min_eco_val, global_params$max_eco_val))
        }
      }
    }
    #current_ecology = Blur_2D(current_ecology, 1.0, 1.0)  
    current_ecology <- write_ecology_by_parcel(current_ecology, parcels$land_parcels, decline_rates, global_params$min_eco_val, global_params$max_eco_val, global_params$ecology_size, yr)
    outputs$trajectories[, , yr] = current_ecology
    print(yr)
  }
  
  outputs$offset_list = index_object$offsets
  outputs$development_list = index_object$developments
  outputs$parcel_sets = index_object$parcel_sets
  
  return(outputs)
}


write_ecology_by_parcel <- function(current_ecology, land_parcels, decline_rates, min_eco_val, max_eco_val, ecology_size, yr){
  parcel_num = length(land_parcels)
  for (parcel_ind in 1:parcel_num){  
    current_parcel = land_parcels[[parcel_ind]]
    decline_rate = decline_rates[parcel_ind]
    if (decline_rate == 0){
      updated_parcel = 0
    } else {      
      updated_parcel <- sapply(current_ecology[current_parcel], eco_change, min_eco_val, max_eco_val, decline_rate, time_step = 1)
    }
    current_ecology[current_parcel] = updated_parcel 
  }
  return(current_ecology) 
}

write_ecology <- function(current_ecology, write_type, land_parcels, current_parcel_indexes, decline_rates, min_eco_val, max_eco_val, ecology_size, offset_rate, yr){
  
  for (region_ind in regions){
    parcel_num = length(current_parcel_indexes)
    for (parcel_ind in 1:parcel_num){
      current_parcel_ind = current_parcel_indexes[[parcel_ind]]
      current_parcel = land_parcels[[current_parcel_ind]]
      if (write_type == 'development'){
        updated_parcel = 0
      } else {
        if (write_type == 'offset'){
          decline_rate = offset_rate
        } else if(write_type == 'remaining'){
          decline_rate = decline_rates[current_parcel_ind]       
        }
        updated_parcel <- sapply(current_ecology[current_parcel], eco_change, min_eco_val, max_eco_val, decline_rate, time_step = 1)
      }
      current_ecology[current_parcel] = updated_parcel 
    }
  }
  return(current_ecology)
  
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





calc_diffs <- function(parcel_list, trajectories, ecology_size, time_steps, parcels, counterfactuals_){
  graphics.off()
  outputs = list()
  true_change = array(0, c(ecology_size, ecology_size, time_steps)) 
  traj_rel_to_initial = true_change
  initial_rel_to_counter = true_change
  
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  
  for (yr in 1:time_steps){
    
    for (parcel_ind in (parcel_list[[yr]])){
      current_parcel = land_parcels[[parcel_ind]]
      for (pix_ind in current_parcel){
        loc = ind2sub(global_ecology_size, pix_ind)
        initial_val = trajectories[loc[1], loc[2], yr]
        current_trajectory = trajectories[loc[1], loc[2], yr:time_steps]
        true_change[loc[1], loc[2], yr:time_steps] = current_trajectory - counterfactuals_object$counterfactuals[loc[1], loc[2], yr:time_steps]       
        traj_rel_to_initial[loc[1], loc[2], yr:time_steps] = current_trajectory - initial_val
        initial_rel_to_counter[loc[1], loc[2], yr:time_steps] = initial_val - counterfactuals_object$counterfactuals[loc[1], loc[2], yr:time_steps]     
      }
    }
  }
  
  outputs$true_change = true_change
  outputs$traj_rel_to_initial = traj_rel_to_initial
  outputs$initial_rel_to_counter = initial_rel_to_counter
  return(outputs)
  
}

# 
# calc_diffs_by_parcel <- function(parcel_list, trajectories, ecology_size, time_steps, land_parcels, counterfactuals){
#   outs = list()
#   outs$traj_rel_initial = list()
#   outs$initial_rel_counter = list()
#   for (yr in 1:time_steps){
#     current_slice = counterfactuals[, , yr]
#     for (parcel_ind in (parcel_list[[yr]])){
#       current_parcel = land_parcels[[parcel_ind]]
#       traj_rel_initial = array(0, c(dim(current_parcel), time_steps))
#       initial_rel_counter = array(0, c(dim(current_parcel), time_steps))
#       current_parcel_sums = current_slice[current_parcel]
#       loc_1 = ind2sub(ecology_size, current_parcel[1])
#       loc_2 = ind2sub(ecology_size, current_parcel[length(current_parcel)])
#       initial_rel_counter[, , yr:time_steps] = current_parcel_sums - counterfactuals[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], yr:time_steps]
#       traj_rel_initial[, , yr:time_steps] = trajectories[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], yr:time_steps] - current_parcel_sums
#     }
#   }
#   
#   outs$traj_rel_initial = traj_rel_initial
#   outs$initial_rel_counter = initial_rel_counter
#   outs$true_change = traj_rel_initial + initial_rel_counter
#   return(outputs)
#   
# }
# 

assess_development_set <- function(parcel_set, trajectories, ecology_size, time_steps, land_parcels, counterfactuals){
  development_parcel_inds = parcel_set$development_parcels
  
  offset_parcel_inds = parcel_set$offset_parcels
  offset_parcel_sums = parcel_sums(, offset_parcels)
  counter_sums = counter_parcel_sums(, offset_parcels)
  offset_sums = colSums()
  
}
  


# collate_parcel_sets <- function(parcel_sets, parcel_type, trajectories, ecology_size, time_steps, land_parcels, counterfactuals){
#   parcel_set_num = length(parcel_sets)
#   dev_parcels = vector('list', parcel_set_num)
#   
#   for (dev_index in 1:parcel_set_num){
#     parcel_set = parcel_sets[[dev_index]] 
#     yr = parcel_set$yr
#     if (parcel_type == development){
#       parcel_to_assess = parcel_set$developed_parcel
#     } else if (parcel_type == offset){
#       parcel_to_assess = parcel_set$offset_parcels
#     }
#     assessed_parcel_set = assess_parcel_set(parcel_to_assess, trajectories, ecology_size, time_steps, land_parcels, counterfactuals, yr)
#     assessed_parcels[[dev_index]] = assessed_parcel_set
#     traj_rel_initial[, dev_index] = assessed_parcel_set$traj_rel_initial_sum
#     initial_rel_counter[, dev_index] = assessed_parcel_set$initial_rel_counter_sum
#   }
#   outs$traj_rel_initial = traj_rel_initial
#   outs$initial_rel_counter = initial_rel_counter
#   return(outs)    
# }
# 
# 
# assess_parcel_set <- function(parcel_inds, trajectories, ecology_size, time_steps, land_parcels, counterfactuals, yr){
#   outs = list()
#   traj_rel_initial_sum = vector()
#   initial_rel_counter_sum = vector()
#   for (parcel_ind in 1:length(parcel_inds)){
#     current_parcel_ind = parcel_inds[parcel_ind]
#     current_parcel = land_parcels[[current_parcel_ind]]
#     traj_rel_initial = array(0, c(dim(current_parcel), time_steps))
#     initial_rel_counter = array(0, c(dim(current_parcel), time_steps))
#     current_parcel_sums = current_slice[current_parcel]
#     loc_1 = ind2sub(ecology_size, current_parcel[1])
#     loc_2 = ind2sub(ecology_size, current_parcel[length(current_parcel)])
#     initial_rel_counter[, , yr:time_steps] = current_parcel_sums - counterfactuals[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], yr:time_steps]
#     traj_rel_initial[, , yr:time_steps] = trajectories[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], yr:time_steps] - current_parcel_sums
#     outs$initial_rel_counter[[parcel_ind]] = initial_rel_counter
#     outs$traj_rel_initial[[parcel_ind]] = traj_rel_initial
#     outs$traj_rel_initial_sum = traj_rel_initial_sum + apply(traj_rel_initial, 3, sum)
#     outs$initial_rel_counter_sum = initial_rel_counter_sum + apply(initial_rel_counter, 3, sum)
#   }
#   return(outs)
# }

# collate_parcel_sets <- function(parcel_sets, time_steps){
#   parcels_object = list()
#   offset_parcels = vector('list', time_steps)
#   offset_vals = vector('list', time_steps)
#   offset_parcels_total = vector()
#   offset_vals_total = vector()
#   developed_parcels = vector('list', time_steps)
#   for (dev_ind in 1:length(parcel_sets)){
#     yr = parcel_sets[[dev_ind]]$yr
#     offset_parcels[[yr]] = c(offset_parcels[[yr]], parcel_sets[[dev_ind]]$offset_parcels)
#     offset_vals[[yr]] = c(offset_vals[[yr]], parcel_sets[[dev_ind]]$offset_vals)
#     offset_parcels_total = c(offset_parcels_total, parcel_sets[[dev_ind]]$offset_parcels)
#     offset_vals_total = c(offset_vals_total, parcel_sets[[dev_ind]]$offset_vals)
#     developed_parcels[[yr]] = c(developed_parcels[[yr]], parcel_sets[[dev_ind]]$developed_parcel)
#   }
#   parcels_object$offset_parcels = offset_parcels
#   parcels_object$offset_vals = offset_vals
#   parcels_object$developed_parcels = developed_parcels
#   parcels_object$offset_vals_total = offset_vals_total
#   parcels_object$offset_parcels_total = offset_parcels_total
#   return(parcels_object)
# }




parcel_3D <- function(current_parcel, trajectories){
  eco_size = dim(trajectories)[1]
  loc_1 = ind2sub(eco_size, current_parcel[1])
  loc_2 = ind2sub(eco_size, current_parcel[length(current_parcel)])
  current_parcel_trajectory = trajectories[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], ]
  return(current_parcel_trajectory)
}


sum_parcels <- function(trajectories, parcels_to_use, land_parcels, time_steps){
  
  parcel_num = length(parcels_to_use)
  parcel_sums = array(0, c(time_steps, parcel_num))
  
  for (parcel_ind in 1:parcel_num){   
    current_parcel_ind = parcels_to_use[parcel_ind]
    current_parcel = land_parcels[[current_parcel_ind]]
    current_parcel_trajectory = parcel_3D(current_parcel, trajectories)
    parcel_sums[, parcel_ind] = apply(current_parcel_trajectory, 3, sum)
  }  
  return(parcel_sums)  
}

sum_regions <- function(parcel_sums, parcel_index_list, regions, time_steps){
  region_num = length(regions)
  region_sums = array(0, c(time_steps, region_num))
  parcel_num = length(parcel_index_list)
  
  for (parcel_ind in 1:parcel_num){
    current_parcel_ind = parcel_index_list[parcel_ind]
    region_ind = find_region(current_parcel_ind, regions)
    region_sums[, region_ind] = region_sums[, region_ind] + parcel_sums[, parcel_ind]
  }
  
  return(region_sums)
}

find_region <- function(parcel, regions){
  
  for (region_ind in 1:length(regions)){
    if (is.element(parcel, regions[[region_ind]])){
      return(region_ind)
    }
  }
}


true_sums <- function(parcel_sums, counter_sums){
  true_sums = parcel_sums - counter_sums
}

# 
# sum_regions <- function(time_steps, parcel_sums, parcels){
#   
#   region_num = parcels$region_num
#   regions = parcels$regions
#   region_sums = array(0, c(time_steps, region_num))
#   for (yr in 1:time_steps){ 
#     current_slice = parcel_sums[yr, ]
#     for (region_ind in 1:region_num){ #number of regions  
#       current_region = regions[[region_ind]]
#       region_sums[yr, region_ind] = sum(current_slice[current_region])    
#     }
#   }
#   return(region_sums)
# }
# 
# 
# 
# find_sums <- function(time_steps, parcels, trajectories){
#   object = list()
#   object$parcel_sums = sum_parcels(time_steps, parcels, trajectories)
#   object$region_sums = sum_regions(time_steps, object$parcel_sums, parcels)
#   object$net_sums = rowSums(object$region_sums)
#   return(object)
# }
# 
# 
# 
# 
# 
# view_individual_offsets <- function(parcel_set, parcels, time_steps, trajectories, counterfactuals){
#   assess_object = list()
#   offset_yr = parcel_set$yr
#   developed_parcel = parcel_set$developed_parcel
#   offset_parcels = parcel_set$offset_parcels
#   dev_trajectory = find_parcel_trajectory(parcels, developed_parcel, time_steps, trajectories)
#   dev_counterfactual = find_parcel_trajectory(parcels, developed_parcel, time_steps, counterfactuals_object$counterfactuals)
#   
#   offset_trajectories = find_parcel_trajectories(parcels, offset_parcels, time_steps, trajectories)
#   
#   offset_counterfactuals = find_parcel_trajectories(parcels, offset_parcels, time_steps, counterfactuals_object$counterfactuals)
#   
#   if (length(offset_parcels) > 1){
#     offset_sum = colSums(offset_trajectories)
#     offset_counter_sum = colSums(offset_counterfactuals)
#   } else {
#     offset_sum = offset_trajectories
#     offset_counter_sum = offset_counterfactuals
#   }
#   
#   assess_object$dev_trajectory = dev_trajectory
#   assess_object$offset_trajectories = offset_trajectories
#   
#   graphics.off()
#   mx = max(c(dev_trajectory, offset_sum))
#   plot(assess_object$dev_trajectory, type = 'l', ylim = c(0, mx))
#   lines((offset_sum - offset_counter_sum), col = 'red', ylim = c(0, mx))
#   lines((offset_sum), ylim = c(0, mx))
#   lines(dev_counterfactual, ylim = c(0, mx))
#   lines(offset_counter_sum, ylim = c(0, mx))
#   chk = (offset_sum[params$time_steps] - offset_counter_sum[params$time_steps])/(params$offset_multiplier*dev_counterfactual[params$time_steps])
#   print(chk)
#   return(assess_object)
# }
# 
# 
# 
# 




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
  A = (fftshift( fft( fftshift( convolve ) ) ))/(M*N)
  A = Re(A)
  return(A)
}



# collate_parcel_sets <- function(parcel_sets, time_steps){
#   parcels_object = list()
#   offset_parcels = vector('list', time_steps)
#   offset_vals = vector('list', time_steps)
#   offset_parcels_total = vector()
#   offset_vals_total = vector()
#   developed_parcels = vector('list', time_steps)
#   for (dev_ind in 1:length(parcel_sets)){
#     yr = parcel_sets[[dev_ind]]$yr
#     offset_parcels[[yr]] = c(offset_parcels[[yr]], parcel_sets[[dev_ind]]$offset_parcels)
#     offset_vals[[yr]] = c(offset_vals[[yr]], parcel_sets[[dev_ind]]$offset_vals)
#     offset_parcels_total = c(offset_parcels_total, parcel_sets[[dev_ind]]$offset_parcels)
#     offset_vals_total = c(offset_vals_total, parcel_sets[[dev_ind]]$offset_vals)
#     developed_parcels[[yr]] = c(developed_parcels[[yr]], parcel_sets[[dev_ind]]$developed_parcel)
#   }
#   parcels_object$offset_parcels = offset_parcels
#   parcels_object$offset_vals = offset_vals
#   parcels_object$developed_parcels = developed_parcels
#   parcels_object$offset_vals_total = offset_vals_total
#   parcels_object$offset_parcels_total = offset_parcels_total
#   return(parcels_object)
# }