import xarray as xr
import numpy as np
import get_models as gm


#####Arctic_fw_davis_CESM2.py#######
#Code to compute the solid and liquid fluxes through Davis Strait for
#CESM2 (CAM6 version).
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
s_si=4 #Sea ice salinity [PSU]. CESM2 has variable sea ice salinity but the model's fw budget is not impacted by it
s_sn=0 #Snow salinity [PSU], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3]
rho_si=917.0 #Sea ice density [kg/m3]
rho_sn=330.0 #Snow density [kg/m3]
rho_si_over_rho_fw=rho_si/rho_fw #Scale factor for converting ice volume to liquid freshwater volume [unitless]
rho_sn_over_rho_fw=rho_sn/rho_fw #Scale factor for converting snow volume to liquid freshwater volume [unitless]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Read in ocean and sea ice grid info outside of the for loop--get this from local direcs
simulation, time_period=gm.choose_model_and_sim(model,'historical',ens_mem='001')
direc_o='/glade/collections/cdg/timeseries-cmip6/'+simulation+'/ocn/proc/tseries/month_1/'
xdo=xr.open_mfdataset(direc_o+simulation+'.pop.h.SALT.*.nc',combine='by_coords',data_vars='minimal')
#Ocean
DYU=xdo['DYU']*0.01 #U cell y-spacing [cm], converted to [m]
DXU=xdo['DXU']*0.01 #U cell x-spacing, converted to [m]
DXT=xdo['DXT']*0.01 #T cell x-spacing, converted to [m]
DYT=xdo['DYT']*0.01 #T cell y-spacing, converted to [m]
HUW=xdo['HUW']*0.01 #U cell widths on west sides [cm], converted to [m]
HUS=xdo['HUS']*0.01 #U cell widths on south sides [cm], converted to [m]
HTE=xdo['HTE']*0.01 #T cell widths on east sides [cm], converted to [m]
HTN=xdo['HTN']*0.01 #T cell widths on north sides [cm], converted to [m]
zt=xdo['z_t']*0.01 #Depth of center of each ocean layer [cm], converted to [m]
dz=xdo['dz']*0.01 #Ocean layer thickness [cm], converted to [m]
#Rename zt and dz vertical dimension from 'z_t' to 'lev'
#to be consistent with CMIP6 output
zt=zt.rename({'z_t':'lev'})
dz=dz.rename({'z_t':'lev'})
#Sea Ice
#Set to ocean grid; no weird masking in the ocean fields
#Makes a negligible difference
iDYU=DYU.rename({'nlat':'nj','nlon':'ni'})
iDXU=DXU.rename({'nlat':'nj','nlon':'ni'})
iDXT=DXT.rename({'nlat':'nj','nlon':'ni'})
iDYT=DYT.rename({'nlat':'nj','nlon':'ni'})
iHUW=HUW.rename({'nlat':'nj','nlon':'ni'})
iHUS=HUS.rename({'nlat':'nj','nlon':'ni'})
iHTE=HTE.rename({'nlat':'nj','nlon':'ni'})
iHTN=HTN.rename({'nlat':'nj','nlon':'ni'})
#Arctic mask, area for storage
#Read in areacello
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/Ofx/'
xdo=xr.open_dataset(direc_static+'areacello/gn/latest/areacello_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
xdo=xr.open_dataset(direc_static+'deptho/gn/latest/deptho_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
odepth=xdo['deptho'] #Ocean depth (tracer) [m]
direc_f2='/glade/work/zanowski/ArcticFW/'
xda=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xda['arctic_mask_with_davis']
iamask=amask.rename({'nlat':'nj','nlon':'ni'})
oarea=oarea.where(amask==1,drop=True).load()
odepth=odepth.where(amask==1,drop=True).load()
iarea=oarea.rename({'nlat':'nj','nlon':'ni'})
arctic_vol=m3tokm3*(oarea*odepth).sum(dim=('nlat','nlon'))


#Read in ice/ocean output for each ensemble member and compute the fluxes
for ens_num in ens_nums:
	print('Starting Ensemble Member '+str('%i' %ens_num)) #+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ice directory for CESM2 for CMIP6
	direc_i='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/SImon/'
	xdi=xr.open_mfdataset(direc_i+'siconc/gn/latest/siconc_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sic=xdi['siconc']*0.01 #sea ice concentration [%], converted to fraction [1]
	xdi=xr.open_mfdataset(direc_i+'sithick/gn/latest/sithick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sit=xdi['sithick'] #sea ice thickness in [m]
	xdi=xr.open_mfdataset(direc_i+'sisnthick/gn/latest/sisnthick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	snt=xdi['sisnthick'] #snow thickness [m]
	xdi=xr.open_mfdataset(direc_i+'siu/gn/latest/siu_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	siu=xdi['siu'] #sea ice zonal velocity in [m/s]
	#Dumb catch for siv in r2i1p1f1 and rr3i1p1f1 which have repeated filesets.GAHHHHHHHHHH. 
	#This is only for the historrical. Fortunately r2 and r3i1p1f1
	#Don't exist in the SSPs, but the code below should fail if they ever do
	if variant in ['r2i1p1f1','r3i1p1f1']:
		xdi=xr.open_mfdataset(direc_i+'siv/gn/latest/siv_SImon_'+simulation+time_period+'.nc',combine='by_coords',data_vars='minimal')
		siv=xdi['siv'] #sea ice meridional velocity in [m/s]
	else:
		xdi=xr.open_mfdataset(direc_i+'siv/gn/latest/siv_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
		siv=xdi['siv'] #sea ice meridional velocity in [m/s]		
	#Ocean directory
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	#Ocean variables
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_o+'so/gn/latest/so_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	salt=xdo['so'] #ocean salinity in PSU
	xdo=xr.open_mfdataset(direc_o+'uo/gn/latest/uo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	u=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_o+'vo/gn/latest/vo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	v=xdo['vo'] #ocean meridional velocity in [m/s]

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop.
	if experiment=='historical':
		salt=salt.sel({'time':slice('1950-01-01',None)})
		u=u.sel({'time':slice('1950-01-01',None)})
		v=v.sel({'time':slice('1950-01-01',None)})
		sic=sic.sel({'time':slice('1950-01-01',None)})
		sit=sit.sel({'time':slice('1950-01-01',None)})
		snt=snt.sel({'time':slice('1950-01-01',None)})
		siu=siu.sel({'time':slice('1950-01-01',None)})
		siv=siv.sel({'time':slice('1950-01-01',None)})

	print('Model output read in and ready for computation')

	#########1. Arctic fw storage [km3]
	salt_arctic=salt.where(amask==1,drop=True).load()
	sit_arctic=rho_si_over_rho_fw*(1-(s_si/s_ref))*sit.where(iamask==1,drop=True).load()
	sic_arctic=sic.where(iamask==1,drop=True).load()
	snt_arctic=rho_sn_over_rho_fw*snt.where(iamask==1,drop=True).load()

	#Compute fw storage (total and anywhere s is >, or <= s_ref; 
	#ensures positive volumes only)
	liq_fw=((s_ref-salt_arctic)/s_ref)*dz*oarea
	liq_fw_storage=m3tokm3*liq_fw.sum(dim=('lev','nlat','nlon')) 
	#keep track of the water with salinities above/below sref 
	liq_fw_storage_below_sref=m3tokm3*liq_fw.where(salt_arctic<=s_ref,other=0).sum(dim=('lev','nlat','nlon'))
	liq_fw_storage_above_sref=m3tokm3*liq_fw.where(salt_arctic>s_ref,other=0).sum(dim=('lev','nlat','nlon'))

	#Solid fw storage--no need for snow concentration correction
	solid_storage=(sit_arctic+snt_arctic)*sic_arctic*iarea 
	solid_fw_storage=m3tokm3*solid_storage.sum(dim=('ni','nj'))

	print('Storage Done!')

	#########2. Solid flux through Davis [km3/yr]

	#Solid flux through Davis Strait: x=[293:304], y=364 #Indices on the full world map
	#Compute the fw volume relative to 34.8. 
	#Flow out of Davis is southward so the sign convention is correct
	j=364
	istart=293
	iend=304
	irange=np.arange(istart,iend+1)
	sicdav=sic[:,j,istart:iend+1]
	sicdav=sicdav.where(xr.ufuncs.isfinite(sicdav),other=0).compute()
	sitdav=sit[:,j,istart:iend+1]
	sitdav=sitdav.where(xr.ufuncs.isfinite(sitdav),other=0).compute()
	sntdav=snt[:,j,istart:iend+1]
	sntdav=sntdav.where(xr.ufuncs.isfinite(sntdav),other=0).compute()
	#Average meridional velocities in y first
	siv2=siv[:,j-1:j+1,istart-1:iend+1]
	siv2=siv2.where(xr.ufuncs.isfinite(siv2),other=0).compute()
	DYUdav=iDYU[j-1:j+1,istart-1:iend+1].load()
	sivdav=(siv2[:,0]*DYUdav[0]+siv2[:,1]*DYUdav[1])/(DYUdav[0]+DYUdav[1])
	#Average in x
	DXUdav=iHUS[j,istart-1:iend+1].load()
	DXTdav=iDXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitdav),coords=sitdav.coords,dims=sitdav.dims)
	for i in range(0,len(irange)):
		vmid[:,i]=(sivdav[:,i]*DXUdav[i]+sivdav[:,i+1]*DXUdav[i+1])/sum(DXUdav[i:i+2])
	#Compute the fluxes
	dav_vel=vmid*DXTdav
	dav_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitdav
	dav_snow=rho_sn_over_rho_fw*sntdav
	dav_solid=(dav_ice+dav_snow)*sicdav*dav_vel
	dav_solid_flux=km3yr*dav_solid.sum(dim=('ni'))


	print('Solid flux through Davis done!')

	#########3. Liquid flux through Davis [km3/yr]

	#Liquid flux through Davis Strait: x=[293,304], y=364 #Indices on full world
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	j=364
	istart=293
	iend=304
	kmax=36 #maximum numdav of depth levels across the transect
	irange=np.arange(istart,iend+1)
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	vdav=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	DYUdav=DYU[j-1:j+1,istart-1:iend+1].load()
	DXUdav=HUS[j,istart-1:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[0:kmax].load()
	#Average meridional velocities in y first
	vdav=(vdav[:,:,0]*DYUdav[0]+vdav[:,:,1]*DYUdav[1])/(DYUdav[0]+DYUdav[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltdav),coords=saltdav.coords,dims=saltdav.dims)
	for i in range(0,len(irange)):
		vmid[:,:,i]=(vdav[:,:,i]*DXUdav[i]+vdav[:,:,i+1]*DXUdav[i+1])/sum(DXUdav[i:i+2])

	#Compute the fluxes
	#Volume flux
	dav_vol=vmid*dzdav*DXTdav
	dav_vol_flux=km3yr*dav_vol.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#FW flux
	dav_liq=((s_ref-saltdav)/s_ref)*dav_vol
	dav_liq_flux=km3yr*dav_liq.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	dav_liq_flux_below_sref=km3yr*dav_liq.where(saltdav<=s_ref,other=0).sum(dim=('lev','nlon'))
	dav_liq_flux_above_sref=km3yr*dav_liq.where(saltdav>s_ref,other=0).sum(dim=('lev','nlon'))

	print('Liquid flux through Davis done!')

	#########4.Save everything

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
	ds['arctic_vol_with_davis'].attrs['long_name']='Total Arctic Ocean volume within the boundaries of Bering, Davis, Fram and the Barents Sea Opening'
	#Save it as a netcdf
	svdirec='/glade/u/home/zanowski/ArcticFW/'
	#Opening in 'a' mode overwrites exsiting variables and 'w' overwrites the whole file
	ds.to_netcdf(path=svdirec+'Arctic_fw_davis_ts_'+model+'_'+experiment+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')

