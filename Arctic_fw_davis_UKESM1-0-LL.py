import xarray as xr
import numpy as np
import os
import get_models as gm

##############Arctic_fw_davis_UKESM1-0-LL.py#################
#Compute liquid and solid fw fluxes through Davis Strait
#for the ensemble members of UKESM1-0-LL
institution='MOHC'
model='UKESM1-0-LL'
experiment='ssp585' #any of the ssps as well
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
#Ensemble numbers are also discontinuous (like this on ESGF)
ens_num_dict={'historical':[1,2,3,4,8,9,10,11,12,16,17,18,19],'ssp126':[1,2,3,4,8],
		'ssp585':[1,2,3,4,8]}
nmax=int(len(ens_num_dict[experiment]))
ens_nums=ens_num_dict[experiment]

#Constants. UKESM uses CICE so I am assuming some constants are as in CESM2
s_ref=34.8 #Reference salinity [PSU]
s_si=4 #Sea ice salinity [PSU]. #NO CLUE IF THIS IS RIGHT. NOT IN ESDOCS
s_sn=0 #Snow salinity [PSU], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3] 
rho_si=917.0 #Sea ice density [kg/m3] #NO CLUE IF THIS IS RIGHT. NOT IN ESDOCS
rho_sn=330.0 #Snow denisty [kg/m3] (from the UKESM website)
rho_si_over_rho_fw=rho_si/rho_fw #Scale factor for converting ice volume to liquid freshwater volume [unitless]
rho_sn_over_rho_fw=rho_sn/rho_fw #Scale factor for converting snow volume to liquid freshwater volume [unitless]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Loop over all ensemble members
#NOTE:version numbers change with each ensemble member,
#so to account for this the code below sorts
#the contents of os.listdir() at the version level
#and chooses the last element, which should be the latest version date.

#Read in grid info outside the loop
#Ocean and sea ice grid variables
#Grid distance file
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+model+'_grid_distances.nc')
#get x and y distances: tracer x and y, x,y distances at v points and x,y distances at u points
DYV=xdf['dyv'] #y-spacing, v points
DXU=xdf['dxu'] #x-spacing, u points
DXT=xdf['dxt'] #x-spacing, tracer pts
DYT=xdf['dyt'] #y-spacing, tracer pts
DXV=xdf['dxv'] #x-spacing, v points
DYU=xdf['dyu'] #y-spacing, u points

#Arctic mask,etc for storage calculations
#Arctic mask
xda=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xda['arctic_mask_with_davis']
#Read in areacello
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f2/Ofx/'
xdo=xr.open_dataset(direc_static+'areacello/gn/v20190705/areacello/areacello_Ofx_'+model+'_piControl_r1i1p1f2_gn.nc')
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
#Read in ocean depth
xdo=xr.open_dataset(direc_static+'deptho/gn/v20190705/deptho/deptho_Ofx_'+model+'_piControl_r1i1p1f2_gn.nc')
odepth=xdo['deptho'] #Ocean depth (tracer) [m]
odepth=odepth.where(amask==1,drop=True)
oarea=oarea.where(amask==1,drop=True)
#Compute the volume of the Arctic
arctic_vol=m3tokm3*(oarea*odepth).sum(dim=('j','i')).compute()

for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)) #+'/'+str('%i' %nmax))
	print('Getting ocean and sea ice variables')
	variant='r'+str('%i' %ens_num)+'i1p1f2'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean variables
	#Ocean directory for UKESM for CMIP6
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	#Ocean variables
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_o+'so/gn/'+sorted(os.listdir(direc_o+'so/gn/'))[-1]+'/so/so_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	salt=xdo['so'] #ocean salinity in ppt
	xdo=xr.open_mfdataset(direc_o+'uo/gn/'+sorted(os.listdir(direc_o+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	u=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_o+'vo/gn/'+sorted(os.listdir(direc_o+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	v=xdo['vo'] #ocean meridional velocity in [m/s]
	xdo=xr.open_mfdataset(direc_o+'thkcello/gn/'+sorted(os.listdir(direc_o+'thkcello/gn/'))[-1]+'/thkcello/thkcello_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	dz=xdo['thkcello'] #ocean cell thicknesses in [m]

	#Ice variables
	#Ice directory for UKESM for CMIP6
	direc_i='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/SImon/'
	#Get sea ice thickness and concentration, snow thickness, and sea ice velocities
	xdi=xr.open_mfdataset(direc_i+'siconc/gn/'+sorted(os.listdir(direc_i+'siconc/gn/'))[-1]+'/siconc/siconc_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sic=xdi['siconc']*0.01 #sea ice concentration [%], converted to fraction [1]
	xdi=xr.open_mfdataset(direc_i+'sithick/gn/'+sorted(os.listdir(direc_i+'sithick/gn/'))[-1]+'/sithick/sithick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sit=xdi['sithick'] #sea ice thickness [m]
	xdi=xr.open_mfdataset(direc_i+'sisnthick/gn/'+sorted(os.listdir(direc_i+'sisnthick/gn/'))[-1]+'/sisnthick/sisnthick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	snt=xdi['sisnthick'] #sea ice snow thickness [m]
	xdi=xr.open_mfdataset(direc_i+'siu/gn/'+sorted(os.listdir(direc_i+'siu/gn/'))[-1]+'/siu/siu_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	siu=xdi['siu'] #sea ice zonal velocity [m/s]
	xdi=xr.open_mfdataset(direc_i+'siv/gn/'+sorted(os.listdir(direc_i+'siv/gn/'))[-1]+'/siv/siv_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	siv=xdi['siv'] #sea ice meridional velocity [m/s]

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs.
	if experiment=='historical':
		salt=salt.sel({'time':slice('1950-01-01',None)})
		u=u.sel({'time':slice('1950-01-01',None)})
		v=v.sel({'time':slice('1950-01-01',None)})
		dz=dz.sel({'time':slice('1950-01-01',None)})
		sic=sic.sel({'time':slice('1950-01-01',None)})
		sit=sit.sel({'time':slice('1950-01-01',None)})
		snt=snt.sel({'time':slice('1950-01-01',None)})
		siu=siu.sel({'time':slice('1950-01-01',None)})
		siv=siv.sel({'time':slice('1950-01-01',None)})

	print('Model output read in and ready for computation')
	

	####-------------Compute the freshwater storage-------------####
	salt_arctic=salt.where(amask==1,drop=True).load()
	sit_arctic=rho_si_over_rho_fw*(1-(s_si/s_ref))*sit.where(amask==1,drop=True).load()
	sic_arctic=sic.where(amask==1,drop=True).load()
	snt_arctic=rho_sn_over_rho_fw*snt.where(amask==1,drop=True).load()
	dz_arctic=dz.where(amask==1,drop=True).load()

	#Compute fw storage (total that can be negative, 
	#and anywhere s is < s_ref; ensures positive volumes only, also s>=s_ref)
	liq_fw=((s_ref-salt_arctic)/s_ref)*dz_arctic*oarea
	liq_fw_storage=m3tokm3*liq_fw.sum(dim=('lev','j','i'))
	#keep track of the water with salinities above/below sref too
	liq_fw_storage_below_sref=m3tokm3*liq_fw.where(salt_arctic<=s_ref,other=0).sum(dim=('lev','j','i'))
	liq_fw_storage_above_sref=m3tokm3*liq_fw.where(salt_arctic>s_ref,other=0).sum(dim=('lev','j','i'))

	#Solid fw storage--no need for snow concentration correction
	solid_storage=(sit_arctic+snt_arctic)*sic_arctic*oarea #s_sn is zero so the salinity scale factor for snow is 1
	solid_fw_storage=m3tokm3*solid_storage.sum(dim=('i','j'))

	print('Storage done!')

	####-------------Compute the solid flux through Davis Strait-------------####

	#-------Davis Strait-------#
	# x=[233,239], y=[289] Indices on full world map
	# v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. 
	#Flow out of Davis is southward so the sign convention is correct.
	j=289
	istart=233
	iend=239
	sivdav=siv[:,j-1:j+1,istart:iend+1]
	sivdav=sivdav.where(xr.ufuncs.isfinite(sivdav),other=0).compute()
	sicdav=sic[:,j,istart:iend+1]
	sicdav=sicdav.where(xr.ufuncs.isfinite(sicdav),other=0).compute()
	sitdav=sit[:,j,istart:iend+1]
	sitdav=sitdav.where(xr.ufuncs.isfinite(sitdav),other=0).compute()
	sntdav=snt[:,j,istart:iend+1]
	sntdav=sntdav.where(xr.ufuncs.isfinite(sntdav),other=0).compute()
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitdav),coords=sitdav.coords,dims=sitdav.dims)
	vmid=(sivdav[:,0]*DYVdav[0,:]+sivdav[:,1]*DYVdav[1,:])/DYVdav.sum(dim='j')
	#Compute the fluxes
	dav_vel=vmid*DXTdav
	dav_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitdav
	#snow is pure fw (s_sn=0) so the salinity scale factor is 1
	dav_snow=rho_sn_over_rho_fw*sntdav
	dav_solid=(dav_ice+dav_snow)*sicdav*dav_vel
	dav_solid_flux=km3yr*dav_solid.sum(dim=('i'))

	print('Solid Flux Done!')


	####-------------Compute the liquid flux through Davis Strait-------------####

	#-------Davis Strait-------#
	#Liquid flux through Davis Strait: x=[233,239], y=[289] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. 
	#Flow out of Davis is southward so the sign convention is correct.
	j=289
	istart=233
	iend=239
	kmax=47 #maximum numdav of depth levels across the transect
	vdav=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again 
	#although it takes basically no time
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(saltdav),coords=saltdav.coords,dims=saltdav.dims)
	vmid=(vdav[:,:,0]*DYVdav[0,:]+vdav[:,:,1]*DYVdav[1,:])/DYVdav.sum(dim='j')
	#Compute the fluxes
	#Volume flux
	dav_vol=vmid*dzdav*DXTdav
	dav_vol_flux=km3yr*dav_vol.sum(dim=('lev','i')) #only one lat pt so no sum over lat
	#FW
	dav_liq=((s_ref-saltdav)/s_ref)*dav_vol
	dav_liq_flux=km3yr*dav_liq.sum(dim=('lev','i')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	dav_liq_flux_below_sref=km3yr*dav_liq.where(saltdav<=s_ref,other=0).sum(dim=('lev','i'))
	dav_liq_flux_above_sref=km3yr*dav_liq.where(saltdav>s_ref,other=0).sum(dim=('lev','i'))

	print('Liquid flux done!')

	####-------------Save everything-------------####

	#Save everything as a netcdf
	#Make the time series an xarray Dataset
	sv_dims=['time'] #['ensemble','time']
	dvs={'davis_solid_flux':(sv_dims,dav_solid_flux),'davis_liq_flux':(sv_dims,dav_liq_flux),
	     'davis_liq_flux_above_sref':(sv_dims,dav_liq_flux_above_sref),'davis_liq_flux_below_sref':(sv_dims,dav_liq_flux_below_sref),
	     'davis_vol_flux':(sv_dims,dav_vol_flux),
	     'solid_fw_storage_with_davis':(sv_dims,solid_fw_storage),'liq_fw_storage_with_davis':(sv_dims,liq_fw_storage),
	     'liq_fw_storage_above_sref_with_davis':(sv_dims,liq_fw_storage_above_sref),
	     'liq_fw_storage_below_sref_with_davis':(sv_dims,liq_fw_storage_below_sref),
	     'arctic_vol_with_davis':(arctic_vol)}
	ds=xr.Dataset(data_vars=dvs, coords={'time':u.coords['time']}) #you can use any variable here that has a time coord
	#Change units attributes to be km3/yr or km3, check the encoding params and the attributes
	for a in [v for v in ds.variables if 'flux' in v]:
		ds[a].attrs['units']='km3 yr-1'
	for a in [v for v in ds.variables if 'storage' in v]:
		ds[a].attrs['units']='km3'
	ds['arctic_vol_with_davis'].attrs['units']='km3'
	ds['arctic_vol_with_davis'].attrs['long_name']='Total Arctic Ocean volume within the boundaries of Bering, Davis, Fram and the Baarents Sea Opening'
	#Save it as a netcdf
	svdirec='/glade/u/home/zanowski/ArcticFW/'
	#Opening in 'a' mode overwrites exsiting variables and 'w' overwrites the whole file
	ds.to_netcdf(path=svdirec+'Arctic_fw_davis_ts_'+model+'_'+experiment+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')
