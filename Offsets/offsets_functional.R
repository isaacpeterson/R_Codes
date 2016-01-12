initialise_parameters <- function(){
  params = list()
  params$ecology_size = 25 #array size to be broken up into sub arrays
  params$time_steps = 5 #number of years etc

  params$max_developments = 1 #maximum number of developments per year
  params$parcel_sz_x = 10 #width of basic land parcel (eg for N = 500 and parcel_sz_x = 10 the ecology is composed of 50 parcels in the x direction)
  params$parcel_sz_y = 10 #width of basic land parcel
  params$min_eco_val = 0
  params$max_eco_val = 100

  params$min_initial_eco_val = 40
  params$max_initial_eco_val = 60

  params$max_decline_rate = 0.002 #max regional decline rate

  params$offset_rate = 0.002 #rate of improvement in land parcel after parcel is offset
  params$region_dimensions = c(2,2) #number of regions in [y x] 
  params$eco_noise = 10

  params$decline_type = 'individual'
  params$parcel_selection_type = 'national' #regional or national
  params$offset_value_type = 'current' #current - use current developed values to select parcel to offset) 
  #predicted - predict the eventual offset value after all timesteps used and select the best package
  params$offset_type = 'match' #match - just match the offset value to the developed value (predicted or current)
  params$decline_rates = c(-0.02,   -0.01,   -0.015,    0.005)
  #params$decline_rates = -p$max_decline_rate*ones(region_num, 1)
  return(params)
}

initialise_ecology <- function(params){
  ecology_initial = matrix(1,params$ecology_size,params$ecology_size)
  pixel_indexes = 1:length(ecology_initial)
  dim(pixel_indexes) = c(nrow(ecology_initial), ncol(ecology_initial))
  land_parcels = matcell(pixel_indexes, c(params$parcel_sz_y, params$parcel_sz_x)) #split the ecology array into a series of subarrays with dimensions sz_x by sz_y
  parcel_num_y = land_parcels$dims[1]
  parcel_num_x = land_parcels$dims[2]
  parcel_num = length(land_parcels$elements) #total number of parcels
  parcel_indexes = 1:parcel_num #index all parcels
  dim(parcel_indexes) = c(parcel_num_y, parcel_num_x) #arrange indicies into array with dimensions of land parcels
  region_size_x = ceiling(parcel_num_x/params$region_dimensions[2])
  region_size_y = ceiling(parcel_num_y/params$region_dimensions[1])

  regions = matcell(parcel_indexes, c(region_size_y, region_size_x))
  region_num = length(regions$elements)
  ind_available = 1:parcel_num #parcel indexes available for development or offsetting

  for (s in 1:parcel_num){
    initial_parcel_value = params$min_initial_eco_val + (params$max_initial_eco_val - params$min_initial_eco_val - params$eco_noise)*runif(1) 
    inds = land_parcels$elements[[s]]
    ecology_initial[inds] = ecology_initial[inds]*initial_parcel_value
  }

  ecology_initial = ecology_initial + params$eco_noise*matrix(runif(params$ecology_size*params$ecology_size),params$ecology_size,params$ecology_size)
  eco_list = list()
  eco_list$ecology_initial = ecology_initial
  eco_list$land_parcels = land_parcels
  eco_list$regions = regions
  eco_list$region_num = region_num 
  eco_list$ind_available = ind_available
  return(eco_list)
}

matcell <- function(x, tilesize){
  
  sz = c(nrow(x), ncol(x))
  A = vector('list', 2)
  for (s in 1:2){
    dm = sz[s]
    T = min(dm, tilesize[s])
    nn = floor( dm / T ) 
    resid = dm %% tilesize[s]
    if (resid == 0 ) {
      A[[s]]=c(rep(1, nn)*T)
    } else A[[s]]=c(rep(1, nn)*T,resid)
  }
  
  rowsizes = A[[1]]
  colsizes = A[[2]]
  rows = length(rowsizes)
  cols = length(colsizes)
  B = vector('list', rows*cols)
  
  rowStart = 0
  
  a = 1
  
  for (i in 1:rows){
    colStart = 0
    for (j in 1:cols){
      B[[a]] = x[colStart+(1:colsizes[j]), rowStart+(1:rowsizes[i])]
      colStart = colStart + colsizes[j]
      a = a + 1
    }
    rowStart = rowStart + rowsizes[i]
  }
  parcel = list()
  parcel$dims = c(length(A[[1]]), length(A[[2]]))
  parcel$elements = B
  return(parcel)
}

determine_current_pixel <- function(pix_ind, ecology){
  
  loc <- ind2sub(ecology, pix_ind)
  val = ecology[pix_ind]   
  
  current_pixel = list()
  current_pixel$loc = loc
  current_pixel$val = val
  return(current_pixel)
}

ind2sub <- function(A, ind){
  m = nrow(A)
  rw = ((ind-1) %% m) + 1 
  cl = floor((ind-1) / m) + 1
  loc = c(rw, cl)
  return(loc)
}

offset_rule <- function(pos_y, params, rate, t){
  mn = params$min_eco_val
  mx = params$max_eco_val
  t_sh = -1/rate * log( ((pos_y - mn)/(mx - pos_y)))
  offset_rule = mn + (mx - mn)/(1 + exp(-rate*( t - t_sh)))
}

initialise_counterfactuals <- function(params, eco_list){
  region_num = eco_list$region_num
  regions = eco_list$regions
  land_parcels = eco_list$land_parcels
  
  time_vec = 1:params$time_steps 
  counterfactuals = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))
  current_ecology = eco_list$ecology_initial
  
  for (region_ind in 1:region_num){

    current_region = regions$elements[[region_ind]]
    current_parcel_num = length(current_region)    

    decline_rate = params$decline_rates[region_ind]
    parcel_inds = current_region
    dim(parcel_inds) = c(1, current_parcel_num)
    for (parcel_ind in parcel_inds){
      current_parcel = land_parcels$elements[[parcel_ind]]
      current_pixel_num = length(current_parcel)
      pixel_inds = current_parcel
      dim(pixel_inds) = c(1, current_pixel_num)
      for (pix_ind in pixel_inds){
        current_pixel <- determine_current_pixel(pix_ind, current_ecology)
        counterfactuals[current_pixel$loc[1], current_pixel$loc[2], ] = offset_rule(current_pixel$val, params, decline_rate, time_vec)
      }
    }
  }
  return(counterfactuals)
}


initialise_outputs <- function(params, eco_list, counterfactuals){
  outputs = list()
  outputs$developed_indexes = vector('list', params$time_steps)
  outputs$offset_indexes = vector('list', params$time_steps)
  offset_trajectories = array(0, c(params$ecology_size, params$ecology_size, params$time_steps))
  outputs$offset_trajectories = offset_trajectories 
  outputs$developed_losses = offset_trajectories
  outputs$avoided_losses = offset_trajectories
  outputs$true_gains = offset_trajectories
  outputs$true_losses = offset_trajectories
  outputs$trajectories = counterfactuals
  outputs$counterfactuals = counterfactuals
  return(outputs)
}

find_parcel_vals <- function(land_parcels, current_ecology){

  parcel_num = length(land_parcels$elements)
  current_parcel_vals = rep(0, parcel_num)
  
  for (s in 1:parcel_num){
    inds = land_parcels$elements[[s]]
    current_parcel_vals[s] = sum(current_ecology[inds])
  }
  return(current_parcel_vals)
}

select_parcel_to_develop <- function(ind_available, current_parcel_vals, land_parcels){
  parcel_ind = sample(1:length(ind_available), 1) #randomly select an index from those available
  parcel_val = current_parcel_vals[parcel_ind] #determine ecological value of selected parcel
  developed_ind = ind_available[parcel_ind]
  ind_available <- ind_available[-parcel_ind] #remove developed parcel from available list      
  current_parcel = land_parcels$elements[[developed_ind]]
  
  developed_list = list()
  developed_list$developed_ind = developed_ind
  developed_list$parcel_val = parcel_val
  developed_list$ind_available = ind_available
  developed_list$current_parcel = current_parcel
  return(developed_list)
}

find_region <- function(eco_list, parcel_ind){
  region_num = eco_list$region_num
  regions = eco_list$regions
  for (s in 1:region_num){
    if (any(regions$elements[[s]] == parcel_ind)){
      region_index = s
      break
    }
    
  }
  return(region_index)
}

write_developed_trajectory <- function(outputs, developed_list, current_ecology, yr, params){
  current_parcel = developed_list$current_parcel
  dim(current_parcel) = c(1, length(current_parcel))
  for (pix_ind in current_parcel){
    current_pixel = determine_current_pixel(pix_ind, current_ecology)         
    current_trajectory = 0
    outputs$trajectories[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_trajectory        
    outputs$developed_losses[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_pixel$val #rep(current_pixel$val, length(yr:params$time_steps))
  }
  return(outputs)
}

select_parcel_to_offset <- function(params, eco_list, developed_list, current_parcel_vals){
  offset_list = list()
  region_num = eco_list$region_num
  regions = eco_list$regions$elements
  parcel_selection_type = params$parcel_selection_type 
  offset_value_type = params$offset_value_type
  offset_type = params$offset_type
  developed_parcel_ind = developed_list$parcel_ind
  developed_parcel_val = developed_list$parcel_val
  ind_available = developed_list$ind_available
  if (parcel_selection_type == 'regional'){        
    for (s in 1:region_num){
      if (any(regions[[s]] == developed_parcel_ind)) {
        offset_ind_available = intersect(regions[[s]], ind_available)                                        
        break
      }
    }
  } else if (parcel_selection_type == 'national'){ 
    offset_ind_available = ind_available
  }
  
  if (length(offset_ind_available) == 0){   
    return(offset_list)
  }
  
  if (offset_value_type == 'current') { 
    parcel_vals_to_use = current_parcel_vals[offset_ind_available]
  } else if (offset_value_type == 'predicted'){
    parcel_vals_to_use = predicted_offset_vals[offset_ind_available]
  }
  
  if (offset_type == 'match') {
    err = abs(parcel_vals_to_use - developed_parcel_val)
    err_ind = which(err == min(err))
    parcel_offset_ind = offset_ind_available[err_ind]
    if (length(parcel_offset_ind) > 1){   #if more than one parcel is of equal value chose at random
      parcel_offset_ind = sample(parcel_offset_ind)
    } 
  }
  
  ind_available = setdiff(ind_available, parcel_offset_ind) #remove offset parcel(s) from available list
  
  offset_list$parcel_offset_ind = parcel_offset_ind
  offset_list$ind_available = ind_available
  return(offset_list)
}

write_offset_trajectory <- function(outputs, current_parcel, current_ecology, params, time_remaining, yr){

  min_eco_val = params$min_eco_val
  max_eco_val = params$max_eco_val
  offset_rate = params$offset_rate
  dim(current_parcel) = c(1, length(current_parcel))

  for (pix_ind in current_parcel){
    current_pixel = determine_current_pixel(pix_ind, current_ecology)   
    current_trajectory = offset_rule(current_pixel$val, params, offset_rate, time_remaining)
    outputs$trajectories[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_trajectory        
    outputs$avoided_losses[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_pixel$val
    outputs$offset_trajectories[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_trajectory
    counterfactual = counterfactuals[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps]
    outputs$true_gains[current_pixel$loc[1], current_pixel$loc[2], yr:params$time_steps] = current_trajectory - counterfactual             
  }
  return(outputs)
  
}

determine_offset_index <- function(params, eco_list, developed_list){
  ind_available = developed_list$ind_available
  region_num = eco_list$region_num
  regions = eco_list$regions$elements
  if (params$parcel_selection_type == 'regional') {
    for (s in 1:region_num){
      if ( any(regions[[s]] == developed_list$parcel_ind )) {
        offset_ind_available = intersect(regions[[s]], ind_available)
        break
      }
    }
  } else if (params$parcel_selection_type == 'national'){
    offset_ind_available = ind_available
  }
  return(offset_ind_available)  
}

  







params <- initialise_parameters()
eco_list <- initialise_ecology(params)
counterfactuals <- initialise_counterfactuals(params, eco_list)
outputs <- initialise_outputs(params, eco_list, counterfactuals)

ind_available = eco_list$ind_available

for (yr in 1:params$time_steps){      #main time loop

  if (length(ind_available) == 0){ #if no parcels are available break out of inner loop
    #disp(['all parcels used yr = ', num2str(yr)])        
    break     
  }

  time_remaining = 0:(params$time_steps - yr)    #remaining time vector
  current_ecology = outputs$trajectories[ , , yr]
  current_parcel_vals <- find_parcel_vals(eco_list$land_parcels, current_ecology)

  if (params$offset_value_type == 'predicted'){   #determine predicted offset parcel values for each remaining parcel
    predicted_parcel_vals = predict_parcels(current_parcel_vals, land_parcels, current_ecology, params, time_remaining)
  }

  current_max_developments = min(params$max_developments, length(ind_available))

#  if (yr%%12 == 0){ 
#    if (current_max_developments > 1){ 
#      develop_num = sample(current_max_developments, 1)
#    } else {develop_num = current_max_developments
#    }
#  } else develop_num = 0
 develop_num = 1  

  if (develop_num > 0){
    for (dev_index in 1:develop_num){    #loop over the number of possible developments           

      developed_list <- select_parcel_to_develop(ind_available, current_parcel_vals, eco_list$land_parcels)   
      ind_available = developed_list$ind_available
      outputs <- write_developed_trajectory(outputs, developed_list, current_ecology, yr, params)
    
      print(developed_list$developed_ind)
      print(ind_available)
      #offset_list <- select_parcel_to_offset(params, eco_list, developed_list, current_parcel_vals)
      #ind_available = offset_list$ind_available
      #current_parcel = land_parcels$elements[[offset_list$parcel_offset_ind]]
      #outputs = write_offset_trajectory(outputs, current_parcel, current_ecology, params, time_remaining, yr)
      
      #current_offset_inds = [current_offset_inds parcel_offset_ind]
      #current_developed_inds = [current_developed_inds developed_parcel_ind]
      
      
      #if (length(ind_available) == 0){ #if no parcels are available break out of inner loop        
      #break 
      #}
    }

  } 
  
  #developed_indexes[yr] = current_developed_inds
  #offset_indexes[yr] = current_offset_inds

#print(yr)
image(outputs$trajectories[, , yr], zlim = c(1,100))
}


