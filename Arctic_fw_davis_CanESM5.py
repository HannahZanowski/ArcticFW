import xarray as xr
import numpy as np
import get_models as gm

##############Arctic_fw_davis_CanESM5.py#################
#Compute liquid and solid fw fluxes through Davis Strait
#for the ensemble members of CanESM5
institution='CCCma'
model='CanESM5'
experiment='ssp585' #any of the ssps as well
nmax=10 #Current number of ensemble members for each experiment
n_start=1
n_end=nmax #nmax,normally

#Constants
s_ref=34.8 #Reference salinity [PSU]
s_si=6 #Sea ice salinity [PSU]. From A. Shao
s_sn=0 #Snow salinity [PSU], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3]
rho_si=900.0 #Sea ice density [kg/m3]. From A. Shao
rho_sn=330.0 #Snow denisty [kg/m3]. From A. Shao
rho_si_over_rho_fw=rho_si/rho_fw #Scale factor for converting ice volume to liquid freshwater volume [unitless]
rho_sn_over_rho_fw=rho_sn/rho_fw #Scale factor for converting snow volume to liquid freshwater volume [unitless]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Ocean and sea ice grid variables
#Grid files from Andrew
direc_f2='/glade/work/zanowski/canESM5_info/'
xdf=xr.open_dataset(direc_f2+'nemo_grid.nc',decode_times=False) #xarray gets pissy otherwise
#get x and y distances: x distances ('1') at v points and y distances ('2') at u points
#Change the dimension names from 'x' and 'y' to 'i' and 'j' for these otherwise it confuses xarray and
#gives 2 dimensions ('i',and 'x' for example) when I sum
xdf2=xdf.rename({'x':'i','y':'j'})
DYV=xdf2['e2v'][:-1,1:-1] #y-spacing, v points, "
DXU=xdf2['e1u'][:-1,1:-1] #x-spacing, u points, "
DXV=xdf2['e1v'][:-1,1:-1] #x-spacing, v points, "
DYU=xdf2['e2u'][:-1,1:-1] #y-spacing, u points, "
DXT=xdf2['e1t'][:-1,1:-1] #x-spacing, tracer pts, "
DYT=xdf2['e2t'][:-1,1:-1] #y-spacing, tracer pts, "

#Cell areas, thicknesses
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p2f1/Ofx/'
xdo=xr.open_dataset(direc_static+'thkcello/gn/v20190429/thkcello/thkcello_Ofx_'+model+'_piControl_r1i1p2f1_gn.nc')
dz=xdo['thkcello'] #Ocean layer thickness [m] (this is 3D--z,y,x)
xdo=xr.open_dataset(direc_static+'areacello/gn/v20190429/areacello/areacello_Ofx_'+model+'_piControl_r1i1p2f1_gn.nc')
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
xdo=xr.open_dataset(direc_static+'deptho/gn/v20190429/deptho/deptho_Ofx_'+model+'_piControl_r1i1p2f1_gn.nc')
odepth=xdo['deptho'] #Ocean depth (tracer) [m]
#Arctic mask for fw storage
direc_f2='/glade/work/zanowski/ArcticFW/'
xda=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xda['arctic_mask_with_davis']
oarea=oarea.where(amask==1,drop=True).load()
odepth=odepth.where(amask==1,drop=True).load()
arctic_vol=m3tokm3*(oarea*odepth).sum(dim=('i','j'))


for ens_num in range(n_start,n_end+1):
	####-------------Read in the output-------------####
	#Loops over ensemble members. Technically I can skip the for loop for ice variables
	#by providing a list of the paths to the sea ice file to open_mfdataset instead
	#but easier to just do one member at a time

	print('Getting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p2f1'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ice directory for CanESM5 for CMIP6
	direc_i='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/SImon/'
	sit_path=direc_i+'sithick/gn/v20190429/sithick/sithick_SImon_'+simulation+time_period+'.nc'
	snt_path=direc_i+'sisnthick/gn/v20190429/sisnthick/sisnthick_SImon_'+simulation+time_period+'.nc'
	siu_path=direc_i+'siu/gn/v20190429/siu/siu_SImon_'+simulation+time_period+'.nc'
	siv_path=direc_i+'siv/gn/v20190429/siv/siv_SImon_'+simulation+time_period+'.nc'
	sic_path=direc_i+'siconc/gn/v20190429/siconc/siconc_SImon_'+simulation+time_period+'.nc'
	#A. Shao says that 'sisnthick' needs to be scaled by snow concentration as well
	snc_path=direc_i+'sisnconc/gn/v20190429/sisnconc/sisnconc_SImon_'+simulation+time_period+'.nc'

	#Ocean directory for CanESM5 for CMIP6
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	salt_path=direc_o+'so/gn/v20190429/so/so_Omon_'+simulation+'*.nc'
	u_path=direc_o+'uo/gn/v20190429/uo/uo_Omon_'+simulation+'*.nc'
	v_path=direc_o+'vo/gn/v20190429/vo/vo_Omon_'+simulation+'*.nc'

	print('Getting ocean variables')
	#Ocean variables
	#Get salinity and velocities
	xdo=xr.open_mfdataset(salt_path,combine='by_coords',data_vars='minimal')
	salt=xdo['so'] #ocean salinity in PSU
	xdo=xr.open_mfdataset(u_path,combine='by_coords',data_vars='minimal')
	u=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(v_path,combine='by_coords',data_vars='minimal')
	v=xdo['vo'] #ocean meridional velocity in [m/s]

	#Ice variables
	print('Getting sea ice variables')
	#Get sea ice thickness and concentration, snow thickness, and sea ice velocities
	xdi=xr.open_dataset(sit_path,chunks={})
	sit=xdi['sithick'] #sea ice thickness [m]
	xdi=xr.open_dataset(sic_path,chunks={})
	sic=xdi['siconc']*0.01 #sea ice concentration [%], converted to fraction [1]
	xdi=xr.open_dataset(snt_path,chunks={})
	snt=xdi['sisnthick'] #sea ice snow thickness [m]
	xdi=xr.open_dataset(siu_path,chunks={})
	siu=xdi['siu'] #sea ice zonal velocity [m/s]
	xdi=xr.open_dataset(siv_path,chunks={})
	siv=xdi['siv'] #sea ice meridional velocity [m/s]
	xdi=xr.open_dataset(snc_path,chunks={})
	snc=xdi['sisnconc']*0.01 #snow concentration [%], converted to fraction [1]


	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop. Will have to do this here for
	#ice variables anyway because they are in one file
	if experiment=='historical':
		salt=salt.sel({'time':slice('1950-01-01',None)})
		u=u.sel({'time':slice('1950-01-01',None)})
		v=v.sel({'time':slice('1950-01-01',None)})
		sit=sit.sel({'time':slice('1950-01-01',None)})
		sic=sic.sel({'time':slice('1950-01-01',None)})
		snt=snt.sel({'time':slice('1950-01-01',None)})
		snc=snc.sel({'time':slice('1950-01-01',None)})
		siu=siu.sel({'time':slice('1950-01-01',None)})
		siv=siv.sel({'time':slice('1950-01-01',None)})
		
	print('Model output read in and ready for computation')

	####-------------Compute the freshwater storage-------------####	
	salt_arctic=salt.where(amask==1,drop=True).load()
	sit_arctic=rho_si_over_rho_fw*(1-(s_si/s_ref))*sit.where(amask==1,drop=True).load()
	sic_arctic=sic.where(amask==1,drop=True).load()
	snt_arctic=rho_sn_over_rho_fw*snt.where(amask==1,drop=True).load()
	snc_arctic=snc.where(amask==1,drop=True).load()

	#Compute fw storage total and where s is >, <= s_ref
	#--ensures positive volumes only)
	liq_fw=((s_ref-salt_arctic)/s_ref)*dz*oarea
	liq_fw_storage=m3tokm3*liq_fw.sum(dim=('lev','j','i'))
	#keep track of the water with salinities above/below sref too
	liq_fw_storage_below_sref=m3tokm3*liq_fw.where(salt_arctic<=s_ref,other=0).sum(dim=('lev','j','i'))
	liq_fw_storage_above_sref=m3tokm3*liq_fw.where(salt_arctic>s_ref,other=0).sum(dim=('lev','j','i'))

	#Solid fw storage--no need for snow concentration correction
	solid_storage=(sit_arctic+snt_arctic*snc_arctic)*sic_arctic*oarea 
	solid_fw_storage=m3tokm3*solid_storage.sum(dim=('i','j'))

	print('Storage Done!')

	####-------------Compute the solid and liquid fluxes through Davis Strait-------------####

	#-------Solid-------#
	# x=[232,239], y=[249] Indices on full world map
	# v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=249
	istart=232
	iend=239
	sivdav=siv[:,j-1:j+1,istart:iend+1]
	sivdav=sivdav.where(xr.ufuncs.isfinite(sivdav),other=0).compute()
	sitdav=sit[:,j,istart:iend+1]
	sitdav=sitdav.where(xr.ufuncs.isfinite(sitdav),other=0).compute()
	sntdav=snt[:,j,istart:iend+1]
	sntdav=sntdav.where(xr.ufuncs.isfinite(sntdav),other=0).compute()
	sicdav=sic[:,j,istart:iend+1]
	sicdav=sicdav.where(xr.ufuncs.isfinite(sicdav),other=0).compute()
	sncdav=snc[:,j,istart:iend+1]
	sncdav=sncdav.where(xr.ufuncs.isfinite(sncdav),other=0).compute()
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sivdav[:,1]),coords=sivdav[:,1].coords,dims=sivdav[:,1].dims)
	vmid=(sivdav[:,0]*DYVdav[0,:]+sivdav[:,1]*DYVdav[1,:])/DYVdav.sum(dim='j')
	#Compute the fluxes
	dav_vel=vmid*DXTdav
	dav_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitdav
	#snow is pure fw (s_sn=0) so the salinity scale factor is 1
	dav_snow=rho_sn_over_rho_fw*sntdav*sncdav
	dav_solid=(dav_ice+dav_snow)*sicdav*dav_vel
	dav_solid_flux=km3yr*dav_solid.sum(dim='i')


	print('Solid Fluxes Done!')

	#-------Liquid-------#
	#Liquid flux through Davis Strait: x=[232,239], y=[249] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=249
	istart=232
	iend=239
	kmax=21 #maximum number of depth levels across the transect
	vdav=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again 
	#although it takes basically no time
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vdav[:,:,1]),coords=vdav[:,:,1].coords,dims=vdav[:,:,1].dims)
	vmid=(vdav[:,:,0]*DYVdav[0,:]+vdav[:,:,1]*DYVdav[1,:])/DYVdav.sum(dim='j')
	#Compute the fluxes
	#Volume flux
	dav_vol=vmid*dzdav*DXTdav
	dav_vol_flux=km3yr*dav_vol.sum(dim=('lev','i')) #only one lat pt so no sum over lat
	#FW flux
	dav_liq=((s_ref-saltdav)/s_ref)*dav_vol
	dav_liq_flux=km3yr*dav_liq.sum(dim=('lev','i')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	dav_liq_flux_below_sref=km3yr*dav_liq.where(saltdav<=s_ref,other=0).sum(dim=('lev','i'))
	dav_liq_flux_above_sref=km3yr*dav_liq.where(saltdav>s_ref,other=0).sum(dim=('lev','i'))	

	print('Liquid fluxes done!')

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
print('Ensemble members '+str('%02i' %n_start)+' to '+str('%02i' %n_end)+' done!')
