initialise_parameters <- function(){
  params = list()
  params$offset_multiplier = 4
  params$parcel_num_x = 20 #number of parcels in x
  params$parcel_num_y = 30 #number of parcels in y
  params$region_num_x = 2 
  params$region_num_y = 3
  params$region_num = params$region_num_x*params$region_num_y
  params$ecology_size = 250 #array size to be broken up into sub arrays
  params$time_steps = 100 #number of years etc
  params$max_developments = 3 #maximum number of developments per year 
  params$min_eco_val = 30
  params$max_eco_val = 100
  params$min_initial_eco_val = 30
  params$max_initial_eco_val = 80
  params$max_decline_rate = 0.02 #max regional decline rate
  params$offset_rate = 0.05 #rate of improvement in land parcel after parcel is offset
  params$eco_noise = 20
  params$parcel_selection_type = 'regional' #regional or national
  params$offset_selection_type = 'current' 
  params$offset_value_type = 'match' #match - match absolute value of offset parcel to the developed value (predicted or current)
  params$region_num = params$region_num_x*params$region_num_y 
  params$development_rates = sample(1:params$max_developments, params$region_num, replace = T)
  #params$development_rates = array(1, params$region_num)
  params$decline_rates = -params$max_decline_rate*runif((params$region_num - 1))
  params$decline_rates = c(params$decline_rates, params$max_decline_rate)
  #params$decline_rate_std = 0.1*params$max_decline_rate*runif((params$region_num))
  params$decline_rate_std = 0.5*params$max_decline_rate
  return(params)
}

rand_vec <- function(N, M, sd) {
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
  parcel_vx = rand_vec(parcel_num_x, params$ecology_size, 5)
  parcel_vy = rand_vec(parcel_num_y, params$ecology_size, 5)
  
  pixel_indexes = 1:(params$ecology_size*params$ecology_size)
  dim(pixel_indexes) = c(params$ecology_size, params$ecology_size)
  land_parcels = mcell(pixel_indexes, parcel_vx, parcel_vy) #split the ecology array into a series of subarrays with dimensions sz_x by sz_y
  land_parcel_num = length(land_parcels$elements) #total number of parcels
  parcel_indexes = 1:land_parcel_num #index all parcels
  dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
  region_vx = rand_vec(params$region_num_x, parcel_num_x, 1) 
  region_vy = rand_vec(params$region_num_y, parcel_num_y, 1)
  
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

initialise_counterfactuals <- function(params, parcels, ecology_initial){
  counterfactuals = list()
  
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
  counterfactuals$elements = elements
  counterfactuals$decline_rates = decline_rates
  return(counterfactuals)
}

initialise_outputs <- function(params, counterfactuals, parcels){
  parcel_num_x = length(parcels$parcel_vx)
  parcel_num_y = length(parcels$parcel_vy)
  outputs = list()
  outputs$trajectories = counterfactuals$elements
  return(outputs)
  
}






select_parcel_to_develop <- function(developed_object, ind_available, current_parcel_vals, parcels){

  land_parcels = parcels$land_parcels
  parcel_index = ind_available[sample(1:length(ind_available), 1)]
  parcel_val = current_parcel_vals[parcel_index] #determine ecological value of selected parcel   
  current_parcel = land_parcels[[parcel_index]]
  
  developed_object$parcel_index = parcel_index
  developed_object$parcel_val = parcel_val
  developed_object$current_parcel = current_parcel
  developed_object$parcel_list[[yr]] = c(developed_object$parcel_list[[yr]], parcel_index)
  developed_object$parcel_vals_list[[yr]] = c(developed_object$parcel_vals_list[[yr]], parcel_val)
  return(developed_object)
}


update_index <- function(index_object, current_parcel_indexes, region_ind, params){
  
  ind_available = index_object$ind_available[[region_ind]]
  ind_available <- setdiff(ind_available, current_parcel_indexes) #remove developed parcel from available list   
  if (length(ind_available) < params$development_rates[[region_ind]] ){   
    break_flag = TRUE
  } else {break_flag = FALSE}
  
  index_object$ind_available[[region_ind]] = ind_available
  index_object$break_flag = break_flag
  
  return(index_object)
  
}


select_offset_pool <- function(region_ind, offset_object, params, parcels, developed_object, index_object, current_parcel_vals, predicted_parcel_vals, predict_flag){
  
  ind_available = index_object$ind_available
  
  region_num = parcels$region_num
  regions = parcels$regions
  parcel_selection_type = params$parcel_selection_type 
  offset_selection_type = params$offset_selection_type
  offset_value_type = params$offset_value_type
  developed_parcel_ind = developed_object$parcel_index
  developed_parcel_val = developed_object$parcel_val
  
  #if (parcel_selection_type == 'regional'){        
  #  for (s in 1:region_num){
  #    if (any(regions[[s]] == developed_parcel_ind)) {
  #      offset_ind_to_use = intersect(regions[[s]], ind_available)                                        
  #      break
  #    }
  #  }
  #} else if (parcel_selection_type == 'national'){ 
  #  offset_ind_to_use = ind_available
  #}
  
  if (length(offset_ind_to_use) < 1 ){   
    offset_object$break_flag = TRUE
    return(offset_object)
  } else {offset_object$break_flag = FALSE}
    
  return(offset_ind_to_use)
}



select_offset_parcel <- function(offset_ind_to_use, offset_object, params, parcels, developed_object, index_object, current_parcel_vals, predicted_parcel_vals, predict_flag){      
  
  if (predict_flag == 'true'){
    current_parcel_vals = predicted_parcel_vals
  }
  
  parcel_vals_to_use = current_parcel_vals[offset_ind_to_use]
  
  land_parcels = parcels$land_parcels
  region_num = parcels$region_num
  regions = parcels$regions
  parcel_selection_type = params$parcel_selection_type 
  offset_selection_type = params$offset_selection_type
  offset_value_type = params$offset_value_type
  developed_parcel_ind = developed_object$parcel_index
  developed_parcel_val = developed_object$parcel_val
  
  parcel_indexes = vector()
  parcel_vals = vector()
  val_to_match = params$offset_multiplier*developed_parcel_val
 
  parcel_indexes = vector()
  
  if (sum(parcel_vals_to_use)<val_to_match){
    offset_object$break_flag = TRUE
    return(offset_object)
  }  else {offset_object$break_flag = FALSE}
  
  while (val_to_match > 0){
    err = val_to_match - parcel_vals_to_use
    if ( all(err > 0) ){
      max_ind = which(parcel_vals_to_use == max(parcel_vals_to_use))
      parcel_indexes = c(parcel_indexes, offset_ind_to_use[max_ind])
      parcel_vals = c(parcel_vals, parcel_vals_to_use[max_ind])
      val_to_match = val_to_match - parcel_vals_to_use[max_ind]     
      offset_ind_to_use = offset_ind_to_use[-max_ind]
      parcel_vals_to_use = current_parcel_vals[offset_ind_to_use]      
    } else{
      max_ind = which(err == max(err[which(err < 0)]))
      parcel_indexes = c(offset_ind_to_use[max_ind], parcel_indexes)
      parcel_vals = c(parcel_vals, parcel_vals_to_use[max_ind])
      break
    }
  }
  
  offset_object$current_parcels = list()
  for (parcel_ind in 1:length(parcel_indexes)){
    offset_object$current_parcels[[parcel_ind]] = land_parcels[[parcel_indexes[parcel_ind]]]
  }
  
  offset_object$current_parcel_indexes = parcel_indexes
  offset_object$parcel_list[[yr]] = c(offset_object$parcel_list[[yr]], parcel_indexes)
  offset_object$parcel_vals_list[[yr]] = c(offset_object$parcel_vals_list[[yr]], parcel_vals)
  
  return(offset_object)  
}


determine_offset_index <- function(params, parcels, developed_object, ind_available){
  region_num = parcels$region_num
  regions = parcels$regions
  if (params$parcel_selection_type == 'regional') {
    for (s in 1:region_num){
      if ( any(regions[[s]] == developed_object$parcel_ind )) {
        offset_ind_to_use = intersect(regions[[s]], ind_available)
        break
      }
    }
  } else if (params$parcel_selection_type == 'national'){
    offset_ind_to_use = ind_available
  }
  return(offset_ind_to_use)  
}

 
write_developed_trajectory <- function(outputs, counterfactuals, developed_object, current_ecology, params, yr){
  
  current_parcel = developed_object$current_parcel
  dim(current_parcel) = c(1, length(current_parcel))
  for (pix_ind in current_parcel){
    loc = ind2sub(params$ecology_size, pix_ind)
    outputs$trajectories[loc[1], loc[2], yr:params$time_steps] = 0
  }
    
  return(outputs)
}

write_offset_trajectory <- function(outputs, counterfactuals, offset_object, current_ecology, params, time_remaining, yr){
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
      outputs$trajectories[loc[1], loc[2], yr:params$time_steps] = current_trajectory        
    }
  }
  return(outputs)
}




sum_regions <- function(time_steps, parcels, trajectories){
  
  region_num = parcels$region_num
  land_parcels = parcels$land_parcels
  regions = parcels$regions
  region_sums = array(0, c(time_steps, region_num))
  
  for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
    current_trajectory = trajectories[ , , yr]

    for (region_ind in 1:region_num){ #number of regions  

      current_region = regions[[region_ind]]
      current_parcel_num = length(current_region)
      current_region_sum = 0

      for (parcel_ind in current_region){ 
        current_parcel = land_parcels[[parcel_ind]]
        current_region_sum = current_region_sum + sum(current_trajectory[current_parcel])
      }
      region_sums[yr, region_ind] = current_region_sum
    }
  }
  return(region_sums)
}

find_parcel_vals <- function(parcels, current_ecology){
  land_parcels = parcels$land_parcels
  land_parcel_num = length(land_parcels)
  current_parcel_vals = rep(0, land_parcel_num)
  
  for (s in 1:land_parcel_num){
    inds = land_parcels[[s]]
    current_parcel_vals[s] = sum(current_ecology[inds])
  }
  return(current_parcel_vals)
}

predict_parcels <- function(region_ind, parcels, index_object, current_ecology, params, time_remaining){
  ind_available = index_object$ind_available[[region_ind]]
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  predicted_parcel_vals = array(0, length(ind_available))

  for (s in ind_available){

    current_parcel = land_parcels[[s]]
    predicted_offset_parcel = array(0, length(current_parcel))
    parcel_pix_ind = 1
    for (pix_ind in current_parcel){
      predicted_pixel_trajectory = offset_rule(current_ecology[pix_ind], params, offset_rate, time_remaining)
      predicted_offset_parcel[parcel_pix_ind] = predicted_pixel_trajectory[length(predicted_pixel_trajectory)]
      parcel_pix_ind <- parcel_pix_ind + 1
    }     
    predicted_parcel_vals[s] = sum(predicted_offset_parcel)
  }
  return(predicted_parcel_vals)

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
    plot(dots[[1]][ , s], type = 'l', col = 'red', ylim = c(ymin, ymax))
    if (plot_num > 1){
      for (t in 2:plot_num){
        lines(dots[[t]][ , s],  ylim = c(ymin, ymax))
      }
    } 
  }
}

determine_devs <- function(developed_object, params, parcels){
  
  outputs = list()
  parcel_list = developed_object$parcel_list
  developed_losses = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  static_losses = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  current_parcels = developed_object$current_parcels
  offset_rate = params$offset_rate
  land_parcels = parcels$land_parcels
  for (yr in 1:params$time_steps){
    
    for (parcel_ind in (parcel_list[[yr]])){
      current_parcel = land_parcels[[parcel_ind]]
      for (pix_ind in current_parcel){
        loc = ind2sub(params$ecology_size, pix_ind)
        developed_losses[loc[1], loc[2], yr:params$time_steps] = counterfactuals$elements[loc[1], loc[2], yr:params$time_steps]
        static_losses[loc[1], loc[2], yr:params$time_steps] = counterfactuals$elements[loc[1], loc[2], yr]
      }
    }
  }
  outputs$developed_losses = developed_losses
  outputs$static_losses = static_losses
  return(outputs)
  
}

determine_offs <- function(object, trajectories, params, parcels){
  
  outputs = list()
  parcel_list = object$parcel_list
  offset_trajectories = array(0, c(params$ecology_size, params$ecology_size, params$time_steps)) 
  true_offset_gains = offset_trajectories
  offset_gain_rel_to_current = offset_trajectories
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
      }
    }
  }
  
  outputs$offset_trajectories = offset_trajectories
  outputs$true_offset_gains = true_offset_gains
  outputs$offset_gain_rel_to_current = offset_gain_rel_to_current
  return(outputs)
  
}


params <- initialise_parameters()
parcels <- initialise_shape_parcels(params)
index_object <- initialise_index_object(parcels)
ecology_initial <- initialise_ecology(params, parcels)
offset_object <- initialise_offset_object(params)
developed_object <- initialise_offset_object(params)
counterfactuals <- initialise_counterfactuals(params, parcels, ecology_initial)
outputs <- initialise_outputs(params, counterfactuals, parcels)

graphics.off()

image(outputs$trajectories[ , , yr],zlim = c(0, 100))

t1=proc.time() 




for (yr in 1:params$time_steps){      #main time loop

  time_remaining = 0:(params$time_steps - yr)    #remaining time vector
  current_ecology = outputs$trajectories[ , , yr]
  
  if (yr%%4 == 0){
    
    current_parcel_vals <- find_parcel_vals(parcels, current_ecology)
   
    for (region_ind in 1:parcels$region_num){
    
    develop_num = params$development_rates[region_ind]
    
    if (develop_num > 0){
      
      for (dev_index in 1:develop_num){    #loop over the number of possible developments           
        if (length(index_object$ind_available[[region_ind]]) < 1){
          break
        }
        developed_object <- select_parcel_to_develop(developed_object, index_object$ind_available[[region_ind]], current_parcel_vals, parcels)      
        outputs <- write_developed_trajectory(outputs, counterfactuals, developed_object, current_ecology, params, yr)     
           
        index_object = update_index(index_object, developed_object$parcel_index, region_ind, params)      
        if (index_object$break_flag == TRUE){
          break
        }  
      
        if (params$offset_selection_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
          predict_flag = 'true'
          predicted_parcel_vals <- predict_parcels(region_ind, parcels, index_object, current_ecology, params, time_remaining)
        } else (predict_flag = 'false')
        
        offset_ind_to_use = index_object$ind_available[[region_ind]]
        
        offset_object <- select_offset_parcel(offset_ind_to_use, offset_object, params, parcels, developed_object, index_object, current_parcel_vals, predicted_parcel_vals, predict_flag)      
        
        if (offset_object$break_flag == TRUE){
          break
        } 
        
        outputs = write_offset_trajectory(outputs, counterfactuals, offset_object, current_ecology, params, time_remaining, yr)         
        index_object = update_index(index_object, offset_object$current_parcel_indexes, region_ind, params)
      
        if (index_object$break_flag == TRUE){
          break
        }             
      }
      #image(outputs$trajectories[ , , yr],zlim = c(0, 100))
    } 

 #   if (length(index_object$ind_available[[region_ind]]) < 1){ #if no parcels are available break out of inner loop
#      print('all parcels used')        
#      break     
#    }
  }
  }
  
print(yr)

}
offset_outs = determine_offs(offset_object, outputs$trajectories, params, parcels)
development_outs <- determine_devs(developed_object, params, parcels)

(proc.time()-t1)[3]
#image(outputs$trajectories[ , , yr],zlim = c(0, 100))

region_sums <- sum_regions(params$time_steps, parcels, outputs$trajectories)
offset_absolute_sums <- sum_regions(params$time_steps, parcels, offset_outs$offset_trajectories)
counterfactual_sums <- sum_regions(params$time_steps, parcels, counterfactuals$elements)

net_impacts <- sum_regions(params$time_steps, parcels, (outputs$trajectories - counterfactuals$elements))
static_losses <- sum_regions(params$time_steps, parcels, development_outs$static_losses)
true_offset_gains <- sum_regions(params$time_steps, parcels, offset_outs$true_offset_gains)
true_losses <- sum_regions(params$time_steps, parcels, development_outs$developed_losses)
offset_gains_rel_to_current <- sum_regions(params$time_steps, parcels, offset_outs$offset_gain_rel_to_current)


plot_params = c(params$region_num_x, params$region_num_y, parcels$region_num)
plot_outs(net_impacts, true_offset_gains, true_losses, plot_params)

plot_outs(net_impacts, true_offset_gains, -true_losses, offset_absolute_sums, plot_params)

























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
# sum_parcels <- function(time_steps, parcels, trajectories){
#   
#   region_num = parcels$region_num
#   land_parcels = parcels$land_parcels
#   regions = parcels$regions
#   region_sums = array(0, c(time_steps, region_num))
#   
#   for (yr in 1:time_steps){ #determine net regional offsets, net regional development_losses
#     current_trajectory = trajectories[ , , yr]
#     
#     for (region_ind in 1:region_num){ #number of regions  
#       
#       current_region = regions[[region_ind]]
#       current_parcel_num = length(current_region)
#       current_region_sum = 0
#       
#       for (parcel_ind in current_region){ 
#         current_parcel = land_parcels[[parcel_ind]]
#         current_region_sum = current_region_sum + sum(current_trajectory[current_parcel])
#       }
#       region_sums[yr, region_ind] = current_region_sum
#     }
#   }
#   return(region_sums)
# }
# 
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


