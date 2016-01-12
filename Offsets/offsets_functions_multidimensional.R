split_vector <- function(N, M, sd, min_width) {

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
  parcel_vx = split_vector(parcel_num_x, global_params$ecology_size, sd = 5, min_width = 3)
  parcel_vy = split_vector(parcel_num_y, global_params$ecology_size, sd = 5, min_width = 3)
  
  pixel_indexes = 1:(global_params$ecology_size*global_params$ecology_size)
  dim(pixel_indexes) = c(global_params$ecology_size, global_params$ecology_size)
  land_parcels = mcell(pixel_indexes, parcel_vx, parcel_vy) #lit the ecology array into a series of subarrays with dimensions sz_x by sz_y
  land_parcel_num = length(land_parcels$elements) #total number of parcels
  parcel_indexes = 1:land_parcel_num #index all parcels
  dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
  region_vx = split_vector(global_params$region_num_x, parcel_num_x, 1, min_width = 3) 
  region_vy = split_vector(global_params$region_num_y, parcel_num_y, 1, min_width = 3)
  
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

initialise_ecologies <- function(global_params, parcels){
  initial_ecologies <- array(0, c(global_params$ecology_size,global_params$ecology_size, global_params$eco_dims))
  
  for (eco_num in 1:global_params$eco_dims){
    initial_ecologies[, , eco_num] = initialise_ecology(global_params, parcels)
  }
  return(initial_ecologies)
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

build_counterfactuals_by_parcel <- function(global_params, decline_rates, land_parcels, initial_ecology){
  current_ecology = initial_ecology
  counterfactuals = array(0, c(global_params$ecology_size, global_params$ecology_size, global_params$time_steps))
  counterfactuals[, , 1] = initial_ecology
  parcel_num = length(land_parcels)
  
  for (yr in 1:global_params$time_steps){
    for (parcel_ind in 1:parcel_num){
      current_parcel = land_parcels[[parcel_ind]]
      current_dec_rate = decline_rates[parcel_ind]      
      updated_parcel = sapply(current_ecology[current_parcel], eco_change, min_eco_val = 0, max_eco_val = 100, decline_rate = current_dec_rate, time_step = 1)
      current_ecology[current_parcel] = updated_parcel 
    }
    if (global_params$blur == TRUE){
      current_ecology = Blur_2D(current_ecology, 0.5, 0.5)
    } 
    counterfactuals[, , yr] = current_ecology
  }
  
  return(counterfactuals)
}

initialise_counterfactuals_multi <- function(global_params, region_params, land_parcels, initial_ecologies, decline_rates){
  counterfactuals = vector('list', global_params$eco_dims)
  time_steps = global_params$time_steps 
  ecology_size = global_params$ecology_size
  for (eco_num in 1:global_params$eco_dims){
    counterfactuals[[eco_num]] = build_counterfactuals_by_parcel(global_params, decline_rates, land_parcels, initial_ecologies[, , eco_num])
  }
  return(counterfactuals)
}

build_decline_rates <- function(parcels, region_params){
  regions = parcels$regions
  region_num = length(regions)
  dec_rates = vector('list', region_num)
  decline_rates = list()
  decline_rates = array(0, dim(parcels$parcel_indexes))
  
  for (region_ind in 1:region_num){ 
    current_region = regions[[region_ind]]
    current_parcel_num = length(current_region)    
    decline_params = c(length(current_region), region_params[[region_ind]]$mean_decline_rate, region_params[[region_ind]]$decline_rate_std) #params$decline_rate_std[region_ind])
    current_decline_rates = matrix(rnorm(decline_params[1], mean = decline_params[2], sd = decline_params[3]), ncol = ncol(current_region))
    dec_rates[[region_ind]] = current_decline_rates
    decline_rates[current_region] = current_decline_rates
  }
  decline_rates_region_list = dec_rates
  
  return(decline_rates)
}

build_decline_rates_multi <- function(parcels, region_params){
  
}

initialise_outputs <- function(counterfactuals, global_params){
  outputs = list()
  outputs$trajectories = array(0, c(global_params$ecology_size, global_params$ecology_size, global_params$time_steps))

  return(outputs)
  
}


select_development_index <- function(ind_available, parcel_num){
  parcel_inds = ind_available[sample(1:length(ind_available), parcel_num)]
  return(parcel_inds)
}



update_development_object <- function(region_params, global_params, region_ind, ind_available, current_ecology, decline_rates, land_parcels, yr){
  development_object = list()
  
  dev_value_type = region_params[[region_ind]]$dev_value_type
  offset_rate = region_params[[region_ind]]$offset_rate
  
  min_eco_val = global_params$min_eco_val
  max_eco_val = global_params$max_eco_val
  time_steps = global_params$time_steps
  
  parcel_inds = select_development_index(ind_available, 1)
  
  current_parcel = land_parcels[[parcel_inds]]
  current_parcel_val = sum(current_ecology[current_parcel])
  
  if (dev_value_type == 'current'){
    parcel_vals_used = current_parcel_val
  } else if (dev_value_type == 'predicted'){
    parcel_vals_used = predict_parcel_vals(predict_type = 'development', current_ecology, parcel_inds, land_parcels, decline_rates, offset_rate, 
                                                 min_eco_val, max_eco_val, time_step = (time_steps - yr))
  }
  development_object = record_parcel_info(development_object, parcel_inds, land_parcels, current_ecology, parcel_vals_used)
  return(development_object)
  
}


predict_parcel_vals <- function(predict_type, current_ecology, parcel_inds, land_parcels, decline_rates, offset_rate, min_eco_val, max_eco_val, time_step){
  
  parcel_num = length(parcel_inds)
  predicted_parcel_vals = array(0, parcel_num)
  current_decline_rates = decline_rates
  
  for (parcel_ind in 1:parcel_num){
    current_parcel_ind = parcel_inds[parcel_ind]
    current_parcel = land_parcels[[current_parcel_ind]]
    current_parcel_ecology = current_ecology[current_parcel]
    if (predict_type == 'development'){
      decline_rate = current_decline_rates[current_parcel_ind]      
    } else (decline_rate = offset_rate)
    predicted_parcel = sapply(current_parcel_ecology, eco_change, min_eco_val, max_eco_val, decline_rate, time_step)
    predicted_parcel_vals[parcel_ind] = sum(predicted_parcel)
  }
  
  return(predicted_parcel_vals)
  
}


find_offset_pool <- function(index_object, region_ind, decline_rates, min_eco_val, max_eco_val, time_steps, region_params, land_parcels, current_ecology, time_remaining, counterfactuals){
  
  offset_pool = list()
  current_offset_pool = index_object$ind_available[[region_ind]]
  current_parcel_sums = find_current_parcel_sums(land_parcels, current_ecology, current_offset_pool)
  
  if (region_params[[region_ind]]$offset_selection_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
      parcel_vals_pool = predict_parcel_vals(predict_type = 'offset', current_ecology, parcel_inds = current_offset_pool, land_parcels, decline_rates, region_params[[region_ind]]$offset_rate,
                                                 min_eco_val, max_eco_val, time_step = length(time_remaining))
      
      
      
    if (region_params[[region_ind]]$predict_offset_mode == 'gains'){
      predicted_slice = counterfactuals[, , time_steps]
      parcel_vals_pool = parcel_vals_pool - counterfactuals_object$final_parcel_vals[current_offset_pool]
    }
    
  } else if (region_params[[region_ind]]$offset_selection_type == 'current'){
    parcel_vals_pool = current_parcel_sums
  }
  
  offset_pool$current_parcel_pool = current_offset_pool
  offset_pool$current_parcel_sums = current_parcel_sums
  offset_pool$parcel_vals_pool = parcel_vals_pool
  
  return(offset_pool)
}                 

select_offset_index <-function(offset_object, val_to_match, parcel_vals_pool, current_offset_pool, offset_multiplier, land_parcels, development_val, current_ecology){
  outs = list()
  parcel_vals_used = vector()
  parcel_indexes = vector()

  while (val_to_match > 0){
    err = val_to_match - parcel_vals_pool
    if ( all(err > 0) ){
      max_ind = which(parcel_vals_pool == max(parcel_vals_pool))
      parcel_indexes = c(parcel_indexes, current_offset_pool[max_ind])
      parcel_vals_used = c(parcel_vals_used, parcel_vals_pool[max_ind])
      val_to_match = val_to_match - parcel_vals_pool[max_ind]     
      current_offset_pool = current_offset_pool[-max_ind]
      parcel_vals_pool = parcel_vals_pool[-max_ind]
    } else{
      max_ind = which(err == max(err[which(err <= 0)]))
      parcel_indexes = c(parcel_indexes, current_offset_pool[max_ind])
      parcel_vals_used = c(parcel_vals_used, parcel_vals_pool[max_ind])
      break
    }
  }
  outs$parcel_indexes = parcel_indexes
  outs$parcel_vals_used = parcel_vals_used
  return(outs)
}
  
record_parcel_info <- function(out_object, parcel_indexes, land_parcels, current_ecology,  parcel_vals_used){

  parcel_num = length(parcel_indexes)
  out_object$current_parcels = list()
  out_object$current_parcel_sums = array(0, parcel_num)
  out_object$current_parcel_ecologies = vector('list', parcel_num)
  
  for (parcel_ind in 1:parcel_num){
    current_parcel = land_parcels[[parcel_indexes[parcel_ind]]]
    out_object$current_parcels[[parcel_ind]] = current_parcel
    current_parcel_ecologies = current_ecology[current_parcel]
    dim(current_parcel_ecologies) = dim(current_parcel)
    out_object$current_parcel_ecologies[[parcel_ind]] = current_parcel_ecologies
    out_object$current_parcel_sums[parcel_ind] = sum(current_ecology[current_parcel])
  }  
  out_object$parcel_indexes = parcel_indexes
  out_object$parcel_vals_used = parcel_vals_used
  return(out_object)
  
}


update_offset_object <- function(offset_pool, offset_multiplier, land_parcels, development_val, yr, current_ecology){      
  
  offset_object = list()
  parcel_vals_pool = offset_pool$parcel_vals_pool
  val_to_match = offset_multiplier*development_val[1]
  if (sum(parcel_vals_pool)<val_to_match){
    offset_object$break_flag = TRUE
    return(offset_object)
  }  else {offset_object$break_flag = FALSE}

  selected_parcel_object = select_offset_index(offset_object, val_to_match, parcel_vals_pool, offset_pool$current_parcel_pool, offset_multiplier, land_parcels, development_val, current_ecology)  
  offset_object = record_parcel_info(offset_object, selected_parcel_object$parcel_indexes, land_parcels, current_ecology,  selected_parcel_object$parcel_vals_used)
    
  return(offset_object)
}


write_null_offset_object <- function(){
  offset_object = list()
  offset_object$parcel_indexes = list()
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
    updated_parcel = sapply(current_ecology[current_parcel], eco_change, min_eco_val = 0, max_eco_val = 100, decline_rate = offset_rate, time_step = 1)
    current_ecology[current_parcel] = updated_parcel 
  }
  return(trajectories)
}


find_current_parcel_sums <- function(land_parcels, current_ecology, parcel_indexes){
  
  parcel_num = length(parcel_indexes)
  current_parcel_sums = rep(0, parcel_num)
  
  for (s in 1:parcel_num){
    parcel_ind = parcel_indexes[s]
    current_parcel = land_parcels[[parcel_ind]]
    current_parcel_sums[s] = sum(current_ecology[current_parcel])
  }
  return(current_parcel_sums)
}






find_parcel_trajectory <- function(land_parcels, parcel_ind, time_steps, trajectories){

  current_parcel = land_parcels[[parcel_ind]]
  parcel_traj = array(0, time_steps)
  
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_slice = trajectories[ , , yr]
    parcel_traj[yr] = sum(current_slice[current_parcel])
  }
  return(parcel_traj)
}


find_parcel_trajectories <- function(land_parcels, parcel_indexes, time_steps, trajectories){

  parcel_trajs = array(0, c(length(parcel_indexes), time_steps))
  
  for (parcel_ind in 1:length(parcel_indexes)){
    current_ind = parcel_indexes[parcel_ind]
    parcel_trajs[parcel_ind, ] = find_parcel_trajectory(land_parcels, current_ind, time_steps, trajectories)
  }
  
  if (length(parcel_indexes) == 1){
    dim(parcel_trajs) = c(length(parcel_trajs), 1)
  }
  
  return(parcel_trajs)
}


write_developments_offsets <- function(index_object, development_object, offset_object, region_ind, yr){
  index_object$developments = c(index_object$developments, development_object$parcel_indexes)
  index_object$offsets = c(index_object$offsets, offset_object$parcel_indexes)
  index_object$parcel_sets[[index_object$dev_count]] = write_parcel_sets(development_object, offset_object, yr, region_ind)
  index_object$dev_count = index_object$dev_count + 1
  return(index_object)
}

write_parcel_sets <- function(development_object, offset_object, yr, region_ind){
  parcel_set = list()
  parcel_set$yr = yr
  parcel_set$development_object = development_object
  parcel_set$offset_object = offset_object
  return(parcel_set)
}



update_ind_available <- function(index_object, parcel_indexes, region_ind, current_develop_num){
  
  ind_available = index_object$ind_available[[region_ind]]
  ind_available = setdiff(ind_available, parcel_indexes) #remove development parcel from available list   
  if (length(ind_available) < current_develop_num ){   
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


find_current_dev_nums <- function(region_params, yr){
  region_num = global_params$region_num
  dev_nums = array(0, region_num)
  for (region_ind in 1:region_num){
    current_dev_num = region_params[[region_ind]]$dev_nums[yr]
    dev_nums[region_ind] = current_dev_num
  }
  return(dev_nums)
}

initialise_trajectories <- function(eco_dims, ecology_size, time_steps){
  trajectories = vector('list', eco_dims)
  for (dm in 1:eco_dims){
    trajectories[[dm]] = array(0, c(ecology_size, ecology_size, time_steps))
  }
  return(trajectories)
}
  
calc_trajectories_multi <- function(outputs, trajectories, global_params, region_params, current_ecology, decline_rates, parcels, index_object, counterfactuals_object, perform_offsets, record_parcel_sets){
  
  for (yr in 1:global_params$time_steps){      #main time loop    
    time_remaining = 0:(global_params$time_steps - yr)     
    if (yr%%global_params$develop_every == 0 && yr < global_params$time_steps){
      current_dev_nums <- find_current_dev_nums(region_params, yr) 
      for (region_ind in 1:parcels$region_num){       
        current_develop_num = current_dev_nums[region_ind]
        if (current_develop_num > 0){
          
          for (dev_index in 1:current_develop_num){    #loop over the number of possible developments           
            if (length(index_object$ind_available[[region_ind]]) < current_develop_num){
              break
            }                                             
        
            development_object <- update_development_object(region_params, global_params, region_ind, index_object$ind_available[[region_ind]], current_ecology, decline_rates, parcels$land_parcels, yr)  
            index_object <- update_ind_available(index_object, development_object$parcel_indexes, region_ind, current_develop_num)      
            decline_rates <- update_decline_rates(decline_rates, development_object$parcel_indexes, offset_rate = 0)
            if (index_object$break_flag == TRUE){
              break
            }                      
            if (perform_offsets == TRUE){
              offset_pool <- find_offset_pool(index_object, region_ind, decline_rates, global_params$min_eco_val, global_params$max_eco_val, global_params$time_steps, region_params, parcels$land_parcels, current_ecology, time_remaining, counterfactuals_object$counterfactuals)
              offset_object <- update_offset_object(offset_pool, region_params[[region_ind]]$offset_multiplier, parcels$land_parcels, development_object$parcel_vals_used, yr, current_ecology)                  
              if (offset_object$break_flag == TRUE){
                break
              }             
              index_object = update_ind_available(index_object, offset_object$parcel_indexes, region_ind, current_develop_num)              
              decline_rates <- update_decline_rates(decline_rates, offset_object$parcel_indexes, offset_rate = region_params[[region_ind]]$offset_rate)
              if (index_object$break_flag == TRUE){
                break
              }               
            } else {
              offset_object <- write_null_offset_object() 
            }
            
            index_object <- write_developments_offsets(index_object, development_object, offset_object, region_ind, yr)
          }
        }
        
        if (global_params$display_object == TRUE){
          image(current_ecology, zlim = c(global_params$min_eco_val, global_params$max_eco_val))
        }
      }
    }
    
    if (global_params$blur == TRUE){
      current_ecology = Blur_2D(current_ecology, 0.5, 0.5)
    }  
    current_ecology <- update_ecology_by_parcel(current_ecology, parcels$land_parcels, decline_rates, global_params$min_eco_val, global_params$max_eco_val, time_step = 1)
    #outputs$trajectories[, , yr] = current_ecology
    trajectories[, , yr] = current_ecology
  }
  
  if (record_parcel_sets == TRUE){
    outputs = list()
    outputs$offset_list = index_object$offsets
    outputs$development_list = index_object$developments
    outputs$parcel_sets = index_object$parcel_sets
    outputs$trajectories = trajectories
    return(outputs)
  } else{
    return(trajectories)
  }
}


update_ecology_by_parcel <- function(current_ecology, land_parcels, decline_rates, min_eco_val, max_eco_val, time_step){
  parcel_num = length(land_parcels)
  
  for (parcel_ind in 1:parcel_num){  
    current_parcel = land_parcels[[parcel_ind]]
    decline_rate = decline_rates[parcel_ind]
    if (decline_rate == 0){
      updated_parcel = 0
    } else {      
      updated_parcel = sapply(current_ecology[current_parcel], eco_change, min_eco_val, max_eco_val, decline_rate, time_step)
    }
    current_ecology[current_parcel] = updated_parcel 
  }
  
  return(current_ecology) 
  
}





plot_outs <- function(...){
  
  dots = list(...)
  plot_params = dots[[length(dots)]]
  graphics.off()
  plot_num = length(dots) - 1
  sub_plot_num = plot_params[1]*plot_params[2]
  A = 1:sub_plot_num
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
  for (s in 1:sub_plot_num){
    plot(dots[[1]][, s], type = 'l', col = 'red', ylim = c(ymin, ymax))
    if (plot_num > 1){
      for (t in 2:plot_num){
        lines(dots[[t]][, s],  ylim = c(ymin, ymax))
      }
    } 
  }
}






initialise_parcel_set_list <- function(parcel_set_num){
  object = list()
  object$traj_rel_initial = vector('list', parcel_set_num)
  object$initial_rel_counter = vector('list', parcel_set_num)
  object$initial_sum = vector('list', parcel_set_num)
  object$net_parcel = vector('list', parcel_set_num)
  object$net_counter = vector('list', parcel_set_num)
  
  object$parcel_list = vector()
  
  object$traj_rel_initial_set = vector()
  object$initial_rel_counter_set = vector()
  object$net_parcel_set = vector()
  object$net_counter_set = vector()
  object$initial_sum_set = vector()
  return(object)
  
}




assess_parcel_sets <- function(parcel_sets, land_parcels, trajectories, counterfactuals, time_steps){
  
  parcel_sets_object = list()
  parcel_set_num = length(parcel_sets)
  
  parcel_sets_object$developments = initialise_parcel_set_list(parcel_set_num)
  parcel_sets_object$offsets = initialise_parcel_set_list(parcel_set_num)
  
  for (parcel_set_ind in 1:parcel_set_num){
    parcel_set = parcel_sets[[parcel_set_ind]]
    yr = parcel_set$yr
    parcel_sets_object$developments = assess_parcel_set(parcel_sets_object$developments, parcel_set$development_object, parcel_set_ind, land_parcels, trajectories, counterfactuals, yr, time_steps)
    
    if (length(parcel_set$offset_object$parcel_indexes) > 0){
      parcel_sets_object$offsets = assess_parcel_set(parcel_sets_object$offsets, parcel_set$offset_object, parcel_set_ind, land_parcels, trajectories, counterfactuals, yr, time_steps)
    } else offset_set = list()
  }
  
  return(parcel_sets_object)
  
}  




assess_parcel_set <- function(parcel_sets_object, current_parcel_set_object, parcel_set_index, land_parcels, trajectories, counterfactuals, yr, time_steps){
  
  parcel_inds = current_parcel_set_object$parcel_indexes
  initial_parcel_ecologies = current_parcel_set_object$current_parcel_ecologies
  
  parcel_num = length(parcel_inds)
  
  traj_rel_initial_list = vector('list', parcel_num)
  initial_rel_counter_list = vector('list', parcel_num)
  
  initial_rel_counter_current = array(0, c(time_steps, parcel_num))
  traj_rel_initial_current = array(0, c(time_steps, parcel_num))
  initial_sum_current = array(0, parcel_num)
  
  net_counter_current = array(0, c(time_steps, parcel_num))
  net_parcel_current = array(0, c(time_steps, parcel_num))
  
  
  for (parcel_ind in 1:parcel_num){
    
    current_parcel_ind = parcel_inds[parcel_ind]
    current_parcel = land_parcels[[current_parcel_ind]]
    current_parcel_trajectory = parcel_through_time(current_parcel, trajectories)
    current_counter_trajectory = parcel_through_time(current_parcel, counterfactuals)
    initial_parcel_ecology = initial_parcel_ecologies[[parcel_ind]]
    
    traj_rel_initial_array = array(0, c(dim(current_parcel), time_steps))
    initial_rel_counter_array = array(0, c(dim(current_parcel), time_steps))
    
    traj_rel_initial_array[, , yr:time_steps] = current_parcel_trajectory[, , yr:time_steps] - as.vector(initial_parcel_ecology)
    initial_rel_counter_array[, , yr:time_steps] = as.vector(initial_parcel_ecology) - current_counter_trajectory[, , yr:time_steps]
    
    traj_rel_initial_current[, parcel_ind] = apply(traj_rel_initial_array, 3, sum)
    initial_rel_counter_current[, parcel_ind] = apply(initial_rel_counter_array, 3, sum)
    initial_sum_current[parcel_ind] = sum(initial_parcel_ecology)
    
    net_parcel_current[, parcel_ind] = apply(current_parcel_trajectory, 3, sum)
    net_counter_current[, parcel_ind] = apply(current_counter_trajectory, 3, sum)
    
    #traj_rel_initial_list[[parcel_ind]] = traj_rel_initial
    #initial_rel_counter_list[[parcel_ind]] = initial_rel_counter
    
  }
  
  parcel_sets_object$traj_rel_initial[[parcel_set_index]] = traj_rel_initial_current
  parcel_sets_object$initial_rel_counter[[parcel_set_index]] = initial_rel_counter_current
  parcel_sets_object$net_parcel[[parcel_set_index]] = net_parcel_current
  parcel_sets_object$net_counter[[parcel_set_index]] = net_counter_current
  parcel_sets_object$initial_sum[[parcel_set_index]] = initial_sum_current
  
  parcel_sets_object$traj_rel_initial_set = cbind(parcel_sets_object$traj_rel_initial_set, rowSums(traj_rel_initial_current))
  parcel_sets_object$initial_rel_counter_set = cbind(parcel_sets_object$initial_rel_counter_set, rowSums(initial_rel_counter_current))
  parcel_sets_object$net_parcel_set = cbind(parcel_sets_object$net_parcel_set, rowSums(net_parcel_current))
  parcel_sets_object$net_counter_set = cbind(parcel_sets_object$net_counter_set, rowSums(net_counter_current))
  parcel_sets_object$initial_sum_set = c(parcel_sets_object$initial_sum_set, sum(initial_sum_current))
  
  parcel_sets_object$parcel_list = c(parcel_sets_object$parcel_list, parcel_inds)
  return(parcel_sets_object)
  
}





plot_parcel_sums <- function(plot_type, assess_type, assess_object, parcel_set_num){
  
  if (plot_type == 'offsets'){
    plot_object = assess_object$offsets
  } else if(plot_type == 'developments') {
    plot_object = assess_object$developments
  }
  
  
  if (assess_type == 'set'){
    plot_1a = plot_object$net_parcel_set[, parcel_set_num]
    plot_1b = plot_object$net_counter_set[, parcel_set_num]
    plot_1c = plot_object$initial_sum_set[parcel_set_num]*array(1, length(plot_1a))
    plot_2a = plot_object$traj_rel_initial_set[, parcel_set_num]
    plot_2b = plot_object$initial_rel_counter_set[, parcel_set_num]

  } else if(assess_type == 'individual'){
    plot_1a = plot_object$net_parcel[, parcel_set_num]
    plot_1b = plot_object$net_counter[, parcel_set_num]
    plot_1c = plot_object$initial_sum[parcel_set_num]*array(1, length(plot_1a))
    plot_2a = plot_object$traj_rel_initial[, parcel_set_num]
    plot_2b = plot_object$initial_rel_counter[, parcel_set_num]
  }
 
  mx = max(cbind(plot_1a, plot_1b))
  mn = min(cbind(plot_1a, plot_1b))
  plot(plot_1a, ylim = c(mn, mx), type = 'l', xlab="time", ylab="trajectory", col = 'red')
  lines(plot_1b, ylim = c(mn, mx), col = 'blue')
  lines(plot_1c, ylim = c(mn, mx), lty = 2)
  
  
  
  plot_2c = plot_2a + plot_2b
  mx = max(cbind(plot_2a, plot_2b, plot_2c))
  mn = min(cbind(plot_2a, plot_2b, plot_2c))
  plot(plot_2a, ylim = c(mn, mx), type = 'l', xlab="time", ylab="split")
  lines(plot_2b, ylim = c(mn, mx), col = 'blue')
  lines(plot_2c, ylim = c(mn, mx), col = 'red')
  
}




plot_parcel_sets <- function(assess_object, assess_type, parcel_set_num){
 
  graphics.off()
  sub_plots = 1:6
  dim(sub_plots) = c(2, 3)
  layout(sub_plots)
  plot_parcel_sums(plot_type = 'developments', assess_type, assess_object, parcel_set_num)
  plot_parcel_sums(plot_type = 'offsets', assess_type, assess_object, parcel_set_num)
  
  plot_a = assess_object$developments$traj_rel_initial_set[, parcel_set_num] + assess_object$developments$initial_rel_counter_set[, parcel_set_num]
  plot_b = assess_object$offsets$traj_rel_initial_set[, parcel_set_num] + assess_object$offsets$initial_rel_counter_set[, parcel_set_num]
  two_plot(plot_a, plot_b, cols = c('red', 'black'))
  lines((plot_a + plot_b), col = 'blue')
  grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted")
  
}



plot_net_parcel_sets <- function(assess_object){
  graphics.off()
  sub_plots = 1:3
  dim(sub_plots) = c(1, 3)
  layout(sub_plots)
  
  plot_a = rowSums(assess_object$developments$traj_rel_initial_set) 
  plot_b = rowSums(assess_object$developments$initial_rel_counter_set)
  plot_c = rowSums(assess_object$offsets$traj_rel_initial_set) 
  plot_d = rowSums(assess_object$offsets$initial_rel_counter_set)
  
  two_plot(plot_a, plot_b, cols = c('red', 'black'))
  lines((plot_a + plot_b), col = 'blue')
  grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted")
  
  two_plot(plot_c, plot_d, cols = c('red', 'black'))
  lines((plot_c + plot_d), col = 'blue')
  grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted")

  two_plot(plot_a + plot_b, plot_c + plot_d, cols = c('red', 'black'))
  lines((plot_a + plot_b) + (plot_c + plot_d), col = 'blue')
  grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted")
  
}



two_plot <- function(plot_a, plot_b, cols){
  mx = max(c(plot_a, plot_b))
  mn = min(c(plot_a, plot_b))
  plot(plot_a, type = 'l', col = cols[1], ylim = c(mn, mx))
  lines(plot_b, col = cols[2], ylim = c(mn, mx))
}










parcel_through_time <- function(current_parcel, trajectories){
  eco_size = dim(trajectories)[1]
  loc_1 = ind2sub(eco_size, current_parcel[1])
  loc_2 = ind2sub(eco_size, current_parcel[length(current_parcel)])
  current_parcel_trajectory = trajectories[loc_1[1]:loc_2[1], loc_1[2]:loc_2[2], ]
  return(current_parcel_trajectory)
}


parcel_realisations <- function(traj_runs, parcel_ind, land_parcels, time_steps){
  realisation_num = length(traj_runs)
  current_parcel = land_parcels[[parcel_ind]]
  
  parcel_realisations = list()
  parcel_realisations$list = vector('list', realisation_num)
  parcel_realisations$sums = array(0, c(time_steps, realisation_num))
  
  for (realisation_ind in 1:realisation_num){
    parcel_3D = parcel_through_time(current_parcel, traj_runs[[realisation_ind]])
    parcel_realisations$list[[realisation_ind]] = parcel_3D
    parcel_realisations$sums[, realisation_ind] = apply(parcel_3D, 3, sum)
  }
  return(parcel_realisations)
}

sum_parcels <- function(trajectories, parcels_to_use, land_parcels, time_steps){
  
  parcel_num = length(parcels_to_use)
  parcel_sums = array(0, c(time_steps, parcel_num))
  
  for (parcel_ind in 1:parcel_num){   
    current_parcel_ind = parcels_to_use[parcel_ind]
    current_parcel = land_parcels[[current_parcel_ind]]
    current_parcel_trajectory = parcel_through_time(current_parcel, trajectories)
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


true_sums <- function(parcel_inds, object_sums, counter_sums){
  true_sums_object = list()
  true_sums_object$parcel_sums = object_sums$parcel_sums - counter_sums$parcel_sums[, parcel_inds]
  true_sums_object$net_sums = object_sums$net_sums - counter_sums$net_sums
  return(true_sums_object)
}


find_sums <- function(parcel_inds, trajectories, parcels, time_steps){
  object = list()
  object$parcel_sums = sum_parcels(trajectories, parcel_inds, parcels$land_parcels, time_steps)
  object$region_sums = sum_regions(object$parcel_sums, parcel_inds, parcels$regions, time_steps)
  object$net_sums = rowSums(object$region_sums)
  return(object)
}


sum_dev_vec <- function(dev_vec){
  sum_array = array(0, length(dev_vec))
  for (ind in 1:100){
    sum_array[ind] = sum(dev_vec[1:ind])
  }
  return(sum_array)
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


overlay_plots <- function(plot_array){
  graphics.off()
  mx = max(plot_array)
  mn = min(plot_array)
  plot(plot_array[, 1], type = 'l', ylim = c(mn, mx))
  for (s in 2:dim(plot_array)[2]){
    lines(plot_array[, s],  ylim = c(mn, mx))
  }
}

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
