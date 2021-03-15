import xarray as xr
import numpy as np
import get_models as gm

#####Arctic_external_fw_fluxes_CESM2.py#######
#Code to compute the river and pme fluxes into the 
#Arctic Ocean in CESM2 (CAM6)
institution='NCAR'
model='CESM2'
experiment='ssp585' #historical, ssp126, or ssp585
#Current number of ensemble members:
#Historical: CAM-11
#SSP126: CAM-3
#SSP585: CAM-3
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,12),'ssp126':[4,10,11],'ssp585':[4,10,11]}
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
direca_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/fx/'
xda=xr.open_dataset(direca_static+'areacella/gn/latest/areacella_fx_'+model+'_piControl_r1i1p1f1_gn.nc')
aarea=xda['areacella'].sel({'lat':slice(50,None)}) #Atmo grid cell area [m^2] north of 50N
direco_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/Ofx/'
xdo=xr.open_dataset(direco_static+'areacello/gn/latest/areacello_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
#Arctic ocean fraction on atmo grid and
#Arctic mask for river flux into ocean
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+model+'_Arctic_Ocean_fraction.nc')
arctic_ocean_fraction=xdf['ocean_fraction'].sel({'lat':slice(50,None)})
xdm=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xdm['arctic_mask_with_davis']
oarea=oarea.where(amask==1,drop=True).load()

#Read in the atmo/ocean output for each ensemble member
for ens_num in ens_nums:
	print('Starting Ensemble Member '+str('%i' %ens_num)) #+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Atmo directory for CESM2 for CMIP6
	direc_a='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Amon/'
	xda=xr.open_mfdataset(direc_a+'evspsbl/gn/latest/evspsbl_Amon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	evs=xda['evspsbl'].sel({'lat':slice(50,None)}) #evapotranspiration flux INTO atmo [kg/m2s], N of 50N
	if ens_num==9: #dumb shit because of the extra 1850-2014 file in the 'latest' folder, that doesn't actually exist
	    file_list=[direc_a+'pr/gn/latest/pr_Amon_'+simulation+'185001-189912.nc',
		       direc_a+'pr/gn/latest/pr_Amon_'+simulation+'190001-194912.nc',
		      direc_a+'pr/gn/latest/pr_Amon_'+simulation+'195001-199912.nc',
		       direc_a+'pr/gn/latest/pr_Amon_'+simulation+'200001-201412.nc']
	    xda=xr.open_mfdataset(file_list,combine='by_coords',data_vars='minimal')
	    pr=xda['pr'].sel({'lat':slice(50,None)}) #atmo precipitation flux (solid and liquid) [kg/m2s], N of 50N
	else:
	    xda=xr.open_mfdataset(direc_a+'pr/gn/latest/pr_Amon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	    pr=xda['pr'].sel({'lat':slice(50,None)}) #atmo precipitation flux (solid and liquid) [kg/m2s], N of 50N
	#Ocean directory
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	xdo=xr.open_mfdataset(direc_o+'vsfriver/gn/latest/vsfriver_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	riv=xdo['vsfriver'] #virtual salt flux from rivers [kg/m2s]

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop.
	if experiment=='historical':
		riv=riv.sel({'time':slice('1950-01-01',None)})
		evs=evs.sel({'time':slice('1950-01-01',None)})
		pr=pr.sel({'time':slice('1950-01-01',None)})

	print('Model output read in and ready for computation')

	####-------------Compute the river fluxes-------------####
	riv_arctic=riv.where(amask==1,drop=True).load()
	riv_fw=rho_fw_inv*riv_arctic*oarea
	riv_fw_flux=km3yr*riv_fw.sum(dim=('nlat','nlon'))

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
