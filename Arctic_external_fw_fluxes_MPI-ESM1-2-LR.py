import xarray as xr
import numpy as np
import os
import get_models as gm


##############Arctic_external_fw_fluxes_MPI-ESM1-2-LR.py#################
#Compute p-e fw fluxes into Arctic Ocean
#for the ensemble members of MPI-ESM1-2-LR
#NOTE: No river fluxes (virtual or real water) available
institution='MPI-M'
model='MPI-ESM1-2-LR'
experiment='historical' #any of the ssps as well
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
#Ensemble numbers are also discontinuous (like this on ESGF)
ens_num_dict={'historical':np.arange(1,11),'ssp126':np.arange(1,11),
                'ssp585':np.arange(1,11)}
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
direca_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/fx/'
xda=xr.open_dataset(direca_static+'areacella/gn/v20190710/areacella/areacella_fx_'+model+'_piControl_r1i1p1f1_gn.nc')
aarea=xda['areacella'].sel({'lat':slice(50,None)}) #Atmo grid cell area [m^2] north of 50N
direco_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/Ofx/'
xdo=xr.open_dataset(direco_static+'areacello/gn/v20190710/areacello/areacello_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
#This [:-1,1:-1] gets rid of masked edge points
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
oarea.shape
#Arctic ocean fraction on atmo grid and
#Arctic mask for river flux into ocean
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+model+'_Arctic_Ocean_fraction.nc')
arctic_ocean_fraction=xdf['ocean_fraction'].sel({'lat':slice(50,None)})
xdm=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xdm['arctic_mask_with_davis']
oarea=oarea.where(amask==1,drop=True).load()

####IMPORTANT:
#MPI's ocean grid is upside down s.t. index zero corresponds to the
#farthest point in the NH and index -1 corresponds to the
#farthest point in the SH. Everything plots upside down, etc.
#The mask and indices for the ocean code are all
#in the frame of this original indexing convention--I have
#NOT flipped the lats so that the NH points are last instead
#of first. The atmo is the right way up.

#NOTE: MPI-ESM1-2-LR has no river fluxes available, so
#this paart of the code is commented out and ignored for now
for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)) #+'/'+str('%i' %nmax))
	print('Getting ocean and atmo variables')
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean variables
	#Ocean directory for MPI for CMIP6
	#direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	#Ocean variables
	#Get river fluxes
	#xdo=xr.open_mfdataset(direc_o+'friver/gn/'+sorted(os.listdir(direc_o+'friver/gn/'))[-1]+'/friver/friver_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	#riv=xdo['friver'] #river flux into ocean [kg/m2/s]

	#Atmo variables
	#Atmo directory for MPI for CMIP6
	direc_a='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Amon/'
	#Get precip and evap fluxes
	xda=xr.open_mfdataset(direc_a+'evspsbl/gn/'+sorted(os.listdir(direc_a+'evspsbl/gn/'))[-1]+'/evspsbl/evspsbl_Amon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	evs=xda['evspsbl'].sel({'lat':slice(50,None)}) #evapotranspiration flux into atmo [kg/m2/s], N of 50N
	xda=xr.open_mfdataset(direc_a+'pr/gn/'+sorted(os.listdir(direc_a+'pr/gn/'))[-1]+'/pr/pr_Amon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	pr=xda['pr'].sel({'lat':slice(50,None)}) #precip flux (solid+liquid) into atmo [kg/m2/s], N of 50N

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop. 
	if experiment=='historical':
		evs=evs.sel({'time':slice('1950-01-01',None)})
		pr=pr.sel({'time':slice('1950-01-01',None)})
		#riv=riv.sel({'time':slice('1950-01-01',None)})

	print('Model output read in and ready for computation')

	####-------------Compute the river fluxes-------------####
	#riv_arctic=riv.where(amask==1,drop=True).load()
	#riv_fw=rho_fw_inv*riv_arctic*oarea
	#riv_fw_flux=km3yr*riv_fw.sum(dim=('j','i'))

	#print('River Flux Done!')

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
	dvs={'pme_fw_flux':(sv_dims,pme_fw_flux)} #'river_fw_flux':(sv_dims,riv_fw_flux)
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

