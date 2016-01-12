initialise_parameters <- function(){
  params = list()
  params$time_steps = 100 #number of timesteps simulation will be run
  params$ecology_size = 50 #ecology array size (square dimension) to be broken up into regions and land parcels  
  params$region_num_x = 1 #numnber of regions in x
  params$region_num_y = 1 #numnber of regions in y
  params$region_num = params$region_num_x*params$region_num_y
  params$parcel_num_x = 10 #numnber of parcels in x
  params$parcel_num_y = 12 #numnber of parcels in x
  params$min_eco_val = 0  #minimum allowable ecological value of smallest ecological element (pixel)
  params$max_eco_val = 100 #maximum "   "     "           "
  params$min_initial_eco_val = 40 #minimum allowable initial ecological value of smallest ecological element (pixel)
  params$max_initial_eco_val = 80 #maximum "   "     "           "
  params$initial_eco_noise = 10 #how much initial variation in pixels per land parcel
  
  params$offset_multiplier = 1 #multiply developed parcel value by this value to determine value to offset
  params$max_developments = 1 #maximum number of developments per year 
  params$develop_every = 1  #how often the policy is implemented
  
  params$max_decline_rate = 0.02 #max regional decline rate parameter used in logisitic equation
  params$offset_rate = 1e-10 #rate of improvement in land parcel after parcel is offset according to logisitic equation
  
  params$perform_offsets = TRUE
  params$offset_policy = 'save'
  params$parcel_selection_type = 'regional' #'regional' - select offset parcel from same region as developed parcel
  #'national' - select offset parcel from any region
  params$offset_selection_type = 'current' #'current' - use current value of offset parcel as offset metric
  #'predicted' - use predicted value (i.e. the value the parcel would attain if the parcel was offset) as offset metric
  
  params$predict_offset_mode = 'gains'       #'saved' - use final predicted value of saved parcel as offset metric
  #'gains' - use final predicted value of saved parcel relative to counterfactual as offset metric
  params$dev_val_type = 'current'         #'current' - use current value of development parcel to determine offset 
  #'predicted' - use predicted value (i.e. the value the parcel would attain if the parcel was NOT developed) to determine offset
  
  params$offset_action_type = 'match' #match - try to match offset parcel to the development value
  #sum - find the highest offset parcel then if there is any left over find the next highest etc. 
  
  params$region_dev_num = sample(1:params$max_developments, params$region_num, replace = T)   #number of developments that occur per region per development cycle
  if (params$region_num > 1){
    params$mean_decline_rates = -params$max_decline_rate*runif((params$region_num - 1))     #set up regions composed of land parcels with decline rates distributed according to the 
    params$mean_decline_rates = c(params$mean_decline_rates, params$max_decline_rate)
  } else {params$mean_decline_rates = -params$max_decline_rate}
  
  #specified mean (params$mean_decline_rates) and std dev (params$decline_rate_std)
  #params$decline_rate_std = 0.1*params$max_decline_rate*runif((params$region_num))
  params$decline_rate_std = 0.5*params$max_decline_rate
  return(params)
}