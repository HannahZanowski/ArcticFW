import xarray as xr
import numpy as np
import os
import get_models as gm

##############Arctic_external_fw_fluxes_CNRM-CM6-1.py#################
#Compute river and p-e fw fluxes into the Arctic Ocean
#for the ensemble members of CNRM-CM6-1
institution='CNRM-CERFACS'
model='CNRM-CM6-1'
experiment='historical' #any of the ssps as well
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,20),'ssp126':np.arange(1,7),
                'ssp585':np.arange(1,6)}
nmax=int(len(ens_num_dict[experiment]))
ens_nums=ens_num_dict[experiment]

#Set constants and get model output
#Constants
s_ref=34.8 #Reference salinity [PSU]
s_pr=0 #Precip salinity [PSU], zero because it's pure freshwater
s_riv=0 #River salinity [PSU], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Cell areas
direca_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f2/fx/'
xda=xr.open_dataset(direca_static+'areacella/gr/v20180814/areacella/areacella_fx_'+model+'_piControl_r1i1p1f2_gr.nc')
aarea=xda['areacella'].sel({'lat':slice(50,None)}) #Atmo grid cell area [m^2] north of 50N
direco_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f2/Ofx/'
xdo=xr.open_dataset(direco_static+'areacello/gn/v20180814/areacello/areacello_Ofx_'+model+'_piControl_r1i1p1f2_gn.nc')
#This [:-1,1:-1] gets rid of masked edge points
oarea=xdo['areacello'][:-1,1:-1] #Ocean tracer grid cell area [m^2]
oarea.shape
#Arctic ocean fraction on atmo grid and
#Arctic mask for river flux into ocean
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+model+'_Arctic_Ocean_fraction.nc')
arctic_ocean_fraction=xdf['ocean_fraction'].sel({'lat':slice(50,None)})
xdm=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xdm['arctic_mask_with_davis']
oarea=oarea.where(amask==1,drop=True).load()

#NOTE:version numbers change with each ensemble member,
#so to account for this the code below sorts
#the contents of os.listdir() at the version level
#and chooses the last element, which should be the latest version date.

for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f2'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)
	simulation_a=simulation[:-2]+'r_' #atmos is not on the native grid
	#Atmo directories for CNRM-CM6-1 for CMIP6
	direc_a='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Amon/'
	evs_path=direc_a+'evspsbl/gr/'+sorted(os.listdir(direc_a+'evspsbl/gr/'))[-1]+'/evspsbl/evspsbl_Amon_'+simulation_a+time_period+'.nc'
	pr_path=direc_a+'pr/gr/'+sorted(os.listdir(direc_a+'pr/gr/'))[-1]+'/pr/pr_Amon_'+simulation_a+time_period+'.nc'
	#Ocean directories for CNRM-CM6-1 for CMIP6
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	riv_path=direc_o+'friver/gn/'+sorted(os.listdir(direc_o+'friver/gn/'))[-1]+'/friver/friver_Omon_'+simulation+'*.nc'

	#This model leaves extra edge points for the cyclic boundaries and the northern boundary
	#so get rid of these by slicing like this: [:,:,:-1,1:-1]
	#Ocean variables
	print('Getting ocean variables')
	#Get river flux
	xdo=xr.open_mfdataset(riv_path,combine='by_coords',data_vars='minimal')
	riv=xdo['friver'][:,:-1,1:-1] #river flux [kg/m2/s]

	#Atmo variables
	print('Getting atmo variables')
	#Get sea ice thickness, snow thickness, and sea ice velocities
	xda=xr.open_dataset(evs_path,chunks={})
	evs=xda['evspsbl'].sel({'lat':slice(50,None)}) #evapotranspiration flux into atmo [kg/m2/s], N of 50N
	xda=xr.open_dataset(pr_path,chunks={})
	pr=xda['pr'].sel({'lat':slice(50,None)}) #precip flux (solid and liquid) into atmo [kg/m2/s], N of 50N

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop. Will have to do this here for
	#ice variables anyway because they are in one file
	if experiment=='historical':
		evs=evs.sel({'time':slice('1950-01-01',None)})
		pr=pr.sel({'time':slice('1950-01-01',None)})
		riv=riv.sel({'time':slice('1950-01-01',None)})

	print('Model output read in and ready for computation')

	####-------------Compute the river fluxes-------------####
	riv_arctic=riv.where(amask==1,drop=True).load()
	riv_fw=rho_fw_inv*riv_arctic*oarea
	riv_fw_flux=km3yr*riv_fw.sum(dim=('y','x'))

	print('River Flux Done!')

	####-------------Compute the Precip, Evap and P-E fluxes-------------####
	pr_fw=rho_fw_inv*pr*aarea*arctic_ocean_fraction
	pr_fw_flux=km3yr*pr_fw.sum(dim=('lat','lon')).compute()
	evs_fw=rho_fw_inv*evs*aarea*arctic_ocean_fraction
	evs_fw_flux=km3yr*evs_fw.sum(dim=('lat','lon')).compute()
	pme_fw_flux=pr_fw_flux-evs_fw_flux

	print('Atmo Fluxes Done!')

	####-------------Save everything-------------####

	#Save everything as a netcdf
	#Make the time series an xarray Dataset
	sv_dims=['time'] #['ensemble','time']
	dvs={'river_fw_flux':(sv_dims,riv_fw_flux),'pme_fw_flux':(sv_dims,pme_fw_flux)}
	ds=xr.Dataset(data_vars=dvs, coords={'time':pr.coords['time']}) #you can use any variable here that has a time coord
	#Change units attributes to be km3/yr, check the encoding params and the attributes
	for a in [v for v in ds.variables if 'flux' in v]:
		ds[a].attrs['units']='km3 yr-1'
	ds.attrs['Description']='Fluxes are the total input to the Arctic Ocean as defined by the boundaries of the Barents Sea Opening and Bering, Davis, and Fram straits'
	#Save it as a netcdf
	svdirec='/glade/u/home/zanowski/ArcticFW/'
	#Opening in 'a' mode overwrites exsiting variables and 'w' overwrites the whole file
	ds.to_netcdf(path=svdirec+'Arctic_fw_external_ts_'+model+'_'+experiment+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')           

