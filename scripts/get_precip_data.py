#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# set up

# The data contains daily ensemble forecasts for 2017 and 2018
# and 20 years of reforecasts (209 forecasts per year). 
# The 2017-18 data and the reforecasts need to be extracted separately.
# The forecasts and observations need to be extracted separately. 
# For station data, data for each country needs to be extracted separately.

# For more details:
# https://github.com/EUPP-benchmark/climetlab-eumetnet-postprocessing-benchmark
# https://eupp-benchmark.github.io/EUPPBench-doc/files/EUPPBench_datasets.html


import xarray as xr
import climetlab as cml

# specify which folder to save the data in
path_save = 'Data\\EUMetNet\\'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# load station reforecasts

# specify path: '-processed' should be added if looking at processed/accumulated variables, e.g. precip
path = 'EUPPBench-training-data-stations-reforecasts-surface-processed'

ds_aus = cml.load_dataset(path, 'austria')
ds_bel = cml.load_dataset(path, 'belgium')
ds_fra = cml.load_dataset(path, 'france')
ds_ger = cml.load_dataset(path, 'germany')
ds_net = cml.load_dataset(path, 'netherlands')


#%% save station reforecasts as netcdf

# specify variable of interest: 'tp6' is total precipitation
# other codes can be found in the online documentation
vari = 'tp6' 

# convert to xarray
fcs_aus = ds_aus.to_xarray()[vari]
fcs_bel = ds_bel.to_xarray()[vari]
fcs_fra = ds_fra.to_xarray()[vari]
fcs_ger = ds_ger.to_xarray()[vari]
fcs_net = ds_net.to_xarray()[vari]

fcs = xr.merge([fcs_aus, fcs_bel, fcs_fra, fcs_ger, fcs_net])
fcs.to_netcdf(path = path_save + vari + '_station_refo_fc.ncdf4') # change path to desired folder


#%% save station reforecast observations as netcdf

obs_aus = ds_aus.get_observations_as_xarray()[vari]
obs_bel = ds_bel.get_observations_as_xarray()[vari]
obs_fra = ds_fra.get_observations_as_xarray()[vari]
obs_ger = ds_ger.get_observations_as_xarray()[vari]
obs_net = ds_net.get_observations_as_xarray()[vari]

obs = xr.merge([obs_aus, obs_bel, obs_fra, obs_ger, obs_net])
obs.to_netcdf(path = path_save + vari + '_station_refo_obs.ncdf4')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# load station forecasts (2017-2018)

path = 'EUPPBench-training-data-stations-forecasts-surface-processed'

ds_aus = cml.load_dataset(path, 'ensemble', 'austria')
ds_bel = cml.load_dataset(path, 'ensemble', 'belgium')
ds_fra = cml.load_dataset(path, 'ensemble', 'france')
ds_ger = cml.load_dataset(path, 'ensemble', 'germany')
ds_net = cml.load_dataset(path, 'ensemble', 'netherlands')

#%% save station forecasts (2017-2018) as netcdf

# convert to xarray
fcs_aus = ds_aus.to_xarray()[vari]
fcs_bel = ds_bel.to_xarray()[vari]
fcs_fra = ds_fra.to_xarray()[vari]
fcs_ger = ds_ger.to_xarray()[vari]
fcs_net = ds_net.to_xarray()[vari]

fcs = xr.merge([fcs_aus, fcs_bel, fcs_fra, fcs_ger, fcs_net])
fcs.to_netcdf(path = path_save + vari + '_station_1718_fc.ncdf4')

#%% save station observations (2017-2018) as netcdf

obs_aus = ds_aus.get_observations_as_xarray()[vari]
obs_bel = ds_bel.get_observations_as_xarray()[vari]
obs_fra = ds_fra.get_observations_as_xarray()[vari]
obs_ger = ds_ger.get_observations_as_xarray()[vari]
obs_net = ds_net.get_observations_as_xarray()[vari]

obs = xr.merge([obs_aus, obs_bel, obs_fra, obs_ger, obs_net])
obs.to_netcdf(path = path_save + vari + '_station_1718_obs.ncdf4')

