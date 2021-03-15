import xarray as xr
import numpy as np
import get_models as gm

#####Arctic_fw_gateways_CESM2.py#######
#Code to compute the solid and liquid fluxes through the Arctic gateways for
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
direc_i='/glade/collections/cdg/timeseries-cmip6/'+simulation+'/ice/proc/tseries/month_1/'
direc_o='/glade/collections/cdg/timeseries-cmip6/'+simulation+'/ocn/proc/tseries/month_1/'
xdi=xr.open_mfdataset(direc_i+simulation+'.cice.h.siv.*.nc',combine='by_coords',data_vars='minimal')
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
amask=xda['arctic_mask']
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
	#Dumb catch for siv in r2i1p1f1 and r3i1p1f1 which have repeated filesets.GAHHHHHHHHHH. 
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

	#########2. Solid fluxes through Arctic Gateways [km3/yr]
	#u,v are not on same grid as tracers, etc so you need to compute velocities at tracer points, and you need to average
	#twice because velocities are catty-corner to tracer points on a B grid. This undoubtedly introduces errors.
	#For a given tracer index, u[index] is at the northeast corner

	##--------Bering Strait
	#Solid flux through Bering Strait: x=[199:200], y=332 #Indices on the full world map
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct
	j=332
	istart=199
	iend=200
	sicber=sic[:,j,istart:iend+1]
	sicber=sicber.where(xr.ufuncs.isfinite(sicber),other=0).compute()
	sitber=sit[:,j,istart:iend+1]
	sitber=sitber.where(xr.ufuncs.isfinite(sitber),other=0).compute()
	sntber=snt[:,j,istart:iend+1]
	sntber=sntber.where(xr.ufuncs.isfinite(sntber),other=0).compute()
	#Average meridional velocities in y first
	siv2=siv[:,j-1:j+1,istart-1:iend+1]
	siv2=siv2.where(xr.ufuncs.isfinite(siv2),other=0).compute()
	DYUber=iDYU[j-1:j+1,istart-1:iend+1].load()
	sivber=(siv2[:,0]*DYUber[0]+siv2[:,1]*DYUber[1])/(DYUber[0]+DYUber[1])
	#Average in x
	DXUber=iHUS[j,istart-1:iend+1].load()
	DXTber=iDXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitber),coords=sitber.coords,dims=sitber.dims)
	vmid[:,0]=(sivber[:,0]*DXUber[0]+sivber[:,1]*DXUber[1])/sum(DXUber[0:2])
	vmid[:,1]=(sivber[:,1]*DXUber[1]+sivber[:,2]*DXUber[2])/sum(DXUber[1:])
	#Compute the fluxes
	ber_vel=vmid*DXTber
	ber_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitber
	ber_snow=rho_sn_over_rho_fw*sntber
	ber_solid=(ber_ice+ber_snow)*sicber*ber_vel
	ber_solid_flux=km3yr*ber_solid.sum(dim=('ni'))


	##--------Nares Strait
	#Solid flux through Nares: x=237, y=[376,377]. Indices on the full world map
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is east (positive) so to keep with the Arctic in=positive sign convention, multiply by -1
	jstart=376
	jend=377
	i=237
	sicnar=sic[:,jstart:jend+1,i]
	sicnar=sicnar.where(xr.ufuncs.isfinite(sicnar),other=0).compute()
	sitnar=sit[:,jstart:jend+1,i]
	sitnar=sitnar.where(xr.ufuncs.isfinite(sitnar),other=0).compute()
	sntnar=snt[:,jstart:jend+1,i]
	sntnar=sntnar.where(xr.ufuncs.isfinite(sntnar),other=0).compute()
	#Average zonal velocities in x first
	siu2=siu[:,jstart-1:jend+1,i-1:i+1]
	siu2=siu2.where(xr.ufuncs.isfinite(siu2),other=0).compute()
	DXUnar=iDXU[jstart-1:jend+1,i-1:i+1].load()
	siunar=(siu2[:,:,0]*DXUnar[:,0]+siu2[:,:,1]*DXUnar[:,1])/(DXUnar[:,0]+DXUnar[:,1])
	#Average in y
	DYUnar=iHUW[jstart-1:jend+1,i].load()
	DYTnar=iDYT[jstart:jend+1,i].load()
	umid=xr.DataArray(data=np.zeros_like(sitnar),coords=sitnar.coords,dims=sitnar.dims)
	umid[:,0]=-1.0*(siunar[:,0]*DYUnar[0]+siunar[:,1]*DYUnar[1])/sum(DYUnar[0:2])
	umid[:,1]=-1.0*(siunar[:,1]*DYUnar[1]+siunar[:,2]*DYUnar[2])/sum(DYUnar[1:])
	#Compute the fluxes
	nar_vel=umid*DYTnar
	nares_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitnar
	#Snow is pure fw so the salinity scale factor is 1
	nares_snow=rho_sn_over_rho_fw*sntnar
	nares_solid=(nares_ice+nares_snow)*sicnar*nar_vel
	nares_solid_flux=km3yr*nares_solid.sum(dim=('nj')) #sum over latitude


	##--------Barrow Strait
	#Solid flux through Barrow Strait: x=[233,237], y=360 #Indices on full world map
	#boundary out of Barrow
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are OUT OF the Arctic so multiply by -1 to get the right sign convention
	j=360
	istart=233
	iend=237
	irange=np.arange(istart,iend+1)
	sicbrw=sic[:,j,istart:iend+1]
	sicbrw=sicbrw.where(xr.ufuncs.isfinite(sicbrw),other=0).compute()
	sitbrw=sit[:,j,istart:iend+1]
	sitbrw=sitbrw.where(xr.ufuncs.isfinite(sitbrw),other=0).compute()
	sntbrw=snt[:,j,istart:iend+1]
	sntbrw=sntbrw.where(xr.ufuncs.isfinite(sntbrw),other=0).compute()
	#Average meridional velocities in y first
	siv2=siv[:,j-1:j+1,istart-1:iend+1]
	siv2=siv2.where(xr.ufuncs.isfinite(siv2),other=0).compute()
	DYUbrw=iDYU[j-1:j+1,istart-1:iend+1].load()
	sivbrw=(siv2[:,0]*DYUbrw[0]+siv2[:,1]*DYUbrw[1])/(DYUbrw[0]+DYUbrw[1])
	#Average in x
	DXUbrw=iHUS[j,istart-1:iend+1].load()
	DXTbrw=iDXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitbrw),coords=sitbrw.coords,dims=sitbrw.dims)
	for i in range(0,len(irange)):
		vmid[:,i]=-1.0*(sivbrw[:,i]*DXUbrw[i]+sivbrw[:,i+1]*DXUbrw[i+1])/sum(DXUbrw[i:i+2])
	#Compute the fluxes
	brw_vel=vmid*DXTbrw
	brw_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbrw
	#Snow is pure fw so the salinity scale factor is 1
	brw_snow=rho_sn_over_rho_fw*sntbrw
	brw_solid=(brw_ice+brw_snow)*sicbrw*brw_vel
	brw_solid_flux=km3yr*brw_solid.sum(dim=('ni')) #sum over longitude


	##--------Fram Strait
	#Solid flux through Fram Strait, x=[98,106], y=[370,378]. Indices on full world. This is 
	#is more or less a 45˚ line so it doesn't have a complex pattern for the grid cells it passes through
	#Fram is oriented such that eastward (positive) velocities are into the Arctic, whereas northward (positive)
	#velocities are out of the Arcitc, so multiply v by -1
	#Compute the fw volume relative to 34.8
	istart=98
	iend=106
	jstart=370
	jend=378
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jstart,jend+1)

	siufram=siu[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	siufram=siufram.where(xr.ufuncs.isfinite(siufram),other=0).compute()
	sivfram=siv[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	sivfram=sivfram.where(xr.ufuncs.isfinite(sivfram),other=0).compute()
	sicfram=sic[:,jstart:jend+1,istart:iend+1]
	sicfram=sicfram.where(xr.ufuncs.isfinite(sicfram),other=0).compute()
	sitfram=sit[:,jstart:jend+1,istart:iend+1]
	sitfram=sitfram.where(xr.ufuncs.isfinite(sitfram),other=0).compute()
	sntfram=snt[:,jstart:jend+1,istart:iend+1]
	sntfram=sntfram.where(xr.ufuncs.isfinite(sntfram),other=0).compute()
	DYUfram=iHUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNfram=iHTN[jstart:jend+1,istart:iend+1].load()
	DXUfram=iDXU[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEfram=iHTE[jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(sitfram),coords=sitfram.coords,dims=sitfram.dims)
	vmid=xr.DataArray(data=np.zeros_like(sitfram),coords=sitfram.coords,dims=sitfram.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,i_range_enum):
		#Avg in x first
		umidx=(siufram[:,j:j+2,i]*DXUfram[j:j+2,i]+siufram[:,j:j+2,i+1]*DXUfram[j:j+2,i+1])/DXUfram[j:j+2,i:i+2].sum(dim='ni')
		vmidx=(sivfram[:,j:j+2,i]*DXUfram[j:j+2,i]+sivfram[:,j:j+2,i+1]*DXUfram[j:j+2,i+1])/DXUfram[j:j+2,i:i+2].sum(dim='ni')
		#Then in y
		umid[:,j,i]=(umidx[:,0]*DYUfram[j,i]+umidx[:,1]*DYUfram[j+1,i])/DYUfram[j:j+2,i].sum(dim='nj')
		vmid[:,j,i]=-1.0*(vmidx[:,0]*DYUfram[j,i]+vmidx[:,1]*DYUfram[j+1,i])/DYUfram[j:j+2,i].sum(dim='nj')

	#Compute the fluxes
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Fram line makes a 45˚ angle with u and v
	#Can only do this when u and v are the same size and in the same place!
	fram_vel=(umid+vmid)*vel_scale*np.sqrt(HTEfram*HTEfram+HTNfram*HTNfram)
	fram_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitfram
	#snow is pure fw so salinity factor is 1
	fram_snow=rho_sn_over_rho_fw*sntfram
	fram_solid=(fram_ice+fram_snow)*sicfram*fram_vel
	fram_solid_flux=km3yr*fram_solid.sum(dim=('ni','nj')) #sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...

	##--------Barents Sea Opening (BSO)
	#Solid flux through the BSO, x=[79,90], y=[355,366]. Indices on full world. This is 
	#is more or less a 45˚ line so it doesn't have a complex pattern for the grid cells it passes through
	#The BSO is oriented such that eastward (positive) velocities are into the Arctic, whereas northward (positive)
	#velocities are out of the Arcitc, so multiply v by -1
	#Compute the fw volume relative to 34.8
	istart=79
	iend=90
	jstart=355
	jend=366
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jstart,jend+1)

	siubso=siu[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	siubso=siubso.where(xr.ufuncs.isfinite(siubso),other=0).compute()
	sivbso=siv[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	sivbso=sivbso.where(xr.ufuncs.isfinite(sivbso),other=0).compute()
	sicbso=sic[:,jstart:jend+1,istart:iend+1]
	sicbso=sicbso.where(xr.ufuncs.isfinite(sicbso),other=0).compute()
	sitbso=sit[:,jstart:jend+1,istart:iend+1]
	sitbso=sitbso.where(xr.ufuncs.isfinite(sitbso),other=0).compute()
	sntbso=snt[:,jstart:jend+1,istart:iend+1]
	sntbso=sntbso.where(xr.ufuncs.isfinite(sntbso),other=0).compute()
	DYUbso=iHUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNbso=iHTN[jstart:jend+1,istart:iend+1].load()
	DXUbso=iDXU[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEbso=iHTE[jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(sitbso),coords=sitbso.coords,dims=sitbso.dims)
	vmid=xr.DataArray(data=np.zeros_like(sitbso),coords=sitbso.coords,dims=sitbso.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,i_range_enum):
		#Avg in x first
		umidx=(siubso[:,j:j+2,i]*DXUbso[j:j+2,i]+siubso[:,j:j+2,i+1]*DXUbso[j:j+2,i+1])/DXUbso[j:j+2,i:i+2].sum(dim='ni')
		vmidx=(sivbso[:,j:j+2,i]*DXUbso[j:j+2,i]+sivbso[:,j:j+2,i+1]*DXUbso[j:j+2,i+1])/DXUbso[j:j+2,i:i+2].sum(dim='ni')
		#Then in y
		umid[:,j,i]=(umidx[:,0]*DYUbso[j,i]+umidx[:,1]*DYUbso[j+1,i])/DYUbso[j:j+2,i].sum(dim='nj')
		vmid[:,j,i]=-1.0*(vmidx[:,0]*DYUbso[j,i]+vmidx[:,1]*DYUbso[j+1,i])/DYUbso[j:j+2,i].sum(dim='nj')

	#Compute the fluxes
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes a 45˚ angle with u and v
	#Can only do this if u,v are same size and in the same place!
	bso_vel=(umid+vmid)*vel_scale*np.sqrt(HTEbso*HTEbso+HTNbso*HTNbso)
	bso_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbso
	#snow is pure fw so salinity factor is 1
	bso_snow=rho_sn_over_rho_fw*sntbso
	bso_solid=(bso_ice+bso_snow)*sicbso*bso_vel
	bso_solid_flux=km3yr*bso_solid.sum(dim=('ni','nj')) #sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
		

	print('Solid fluxes through gateways done!')

	#########3. Liquid fluxes through Arctic Gateways [km3/yr]

	##--------Bering Strait
	#Liquid flux through Bering Strait: x=[199,200], y=332 #Indices on full world
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	j=332
	istart=199
	iend=200
	kmax=4 #maximum number of depth levels across the transect
	saltber=salt[:,0:kmax,j,istart:iend+1].load()
	vber=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vber=vber.where(xr.ufuncs.isfinite(vber),other=0).compute()
	DYUber=DYU[j-1:j+1,istart-1:iend+1].load()
	DXUber=HUS[j,istart-1:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	dzber=dz[0:kmax].load()
	#Average meridional velocities in y first
	vber=(vber[:,:,0]*DYUber[0]+vber[:,:,1]*DYUber[1])/(DYUber[0]+DYUber[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltber),coords=saltber.coords,dims=saltber.dims)
	vmid[:,:,0]=(vber[:,:,0]*DXUber[0]+vber[:,:,1]*DXUber[1])/sum(DXUber[0:2])
	vmid[:,:,1]=(vber[:,:,1]*DXUber[1]+vber[:,:,2]*DXUber[2])/sum(DXUber[1:])

	#Compute the fluxes
	#Volume flux
	ber_vol=vmid*dzber*DXTber
	ber_vol_flux=km3yr*ber_vol.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#FW flux
	ber_liq=((s_ref-saltber)/s_ref)*ber_vol
	ber_liq_flux=km3yr*ber_liq.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	ber_liq_flux_below_sref=km3yr*ber_liq.where(saltber<=s_ref,other=0).sum(dim=('lev','nlon'))
	ber_liq_flux_above_sref=km3yr*ber_liq.where(saltber>s_ref,other=0).sum(dim=('lev','nlon'))


	##--------Nares Strait
	#Liquid flux through Nares: x=237, y=[376,377] #Indices on full world
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is east (positive) so to keep with the Arctic sign covention, multiply by -1
	#Nares is only 2 grid cells wide, and the walls are vertical
	jstart=376
	jend=377
	i=237
	kmax=11  #Maximum number of deptth levels across the transect
	saltnar=salt[:,0:kmax,jstart:jend+1,i].load()
	unar=u[:,0:kmax,jstart-1:jend+1,i-1:i+1]
	unar=unar.where(xr.ufuncs.isfinite(unar),other=0).compute()
	DXUnar=DXU[jstart-1:jend+1,i-1:i+1].load()
	DYUnar=HUW[jstart-1:jend+1,i].load()
	DYTnar=DYT[jstart:jend+1,i].load()
	dznar=dz[0:kmax].load()
	#Average zonal velocities in x first
	unar=(unar[:,:,:,0]*DXUnar[:,0]+unar[:,:,:,1]*DXUnar[:,1])/(DXUnar[:,0]+DXUnar[:,1])
	#Average in y
	umid=xr.DataArray(data=np.zeros_like(saltnar),coords=saltnar.coords,dims=saltnar.dims)
	umid[:,:,0]=-1.0*(unar[:,:,0]*DYUnar[0]+unar[:,:,1]*DYUnar[1])/sum(DYUnar[0:2])
	umid[:,:,1]=-1.0*(unar[:,:,1]*DYUnar[1]+unar[:,:,2]*DYUnar[2])/sum(DYUnar[1:])

	#Compute the fluxes
	#Volume flux
	nares_vol=umid*dznar*DYTnar
	nares_vol_flux=km3yr*nares_vol.sum(dim=('lev','nlat')) #only one lon pt so no sum over nlon
	#FW flux
	nares_liq=((s_ref-saltnar)/s_ref)*nares_vol
	nares_liq_flux=km3yr*nares_liq.sum(dim=('lev','nlat')) #only one lon pt so no sum over nlon
	#Fluxes above and below sref
	nares_liq_flux_below_sref=km3yr*nares_liq.where(saltnar<=s_ref,other=0).sum(dim=('lev','nlat'))
	nares_liq_flux_above_sref=km3yr*nares_liq.where(saltnar>s_ref,other=0).sum(dim=('lev','nlat'))


	##--------Barrow Strait
	#Liquid fw flux through Barrow Strait: x=[233,237], y=360 #Indices on full world map
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are OUT OF the Arctic so multiply by -1 to get the right sign convention
	j=360
	istart=233
	iend=237
	kmax=25 #maximum number of depth levels across the transect
	irange=np.arange(istart,iend+1)
	saltbrw=salt[:,0:kmax,j,istart:iend+1].load()
	vbrw=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vbrw=vbrw.where(xr.ufuncs.isfinite(vbrw),other=0).compute()
	DYUbrw=DYU[j-1:j+1,istart-1:iend+1].load()
	DXUbrw=HUS[j,istart-1:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	dzbrw=dz[0:kmax].load()
	#Average meridional velocities in y first
	vbrw=(vbrw[:,:,0]*DYUbrw[0]+vbrw[:,:,1]*DYUbrw[1])/(DYUbrw[0]+DYUbrw[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltbrw),coords=saltbrw.coords,dims=saltbrw.dims)
	for i in range(0,len(irange)):
		vmid[:,:,i]=-1.0*(vbrw[:,:,i]*DXUbrw[i]+vbrw[:,:,i+1]*DXUbrw[i+1])/sum(DXUbrw[i:i+2])

	#Compute the fluxes
	#Volume flux
	brw_vol=vmid*dzbrw*DXTbrw
	brw_vol_flux=km3yr*brw_vol.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#FW flux
	brw_liq=((s_ref-saltbrw)/s_ref)*brw_vol
	brw_liq_flux=km3yr*brw_liq.sum(dim=('lev','nlon')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	brw_liq_flux_below_sref=km3yr*brw_liq.where(saltbrw<=s_ref,other=0).sum(dim=('lev','nlon'))
	brw_liq_flux_above_sref=km3yr*brw_liq.where(saltbrw>s_ref,other=0).sum(dim=('lev','nlon'))

	
	##--------Fram Strait
	#Liquid flux through Fram Strait, x=[98,106], y=[370,378]. Indices on full world. This is
	#is a 45˚ line so it doesn't have a complex pattern for the grid cells it passes through
	#Fram is oriented such that eastward (positive) velocities are into the Arctic, whereas northward (positive)
	#velocities are out of the Arcitc, so multiply v by -1
	#Compute the fw volume relative to 34.8
	istart=98
	iend=106
	jstart=370
	jend=378
	kmax=47 #maximum number of depth levels across the transect
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jstart,jend+1)

	saltfram=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	ufram=u[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	ufram=ufram.where(xr.ufuncs.isfinite(ufram),other=0).compute()
	vfram=v[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	vfram=vfram.where(xr.ufuncs.isfinite(vfram),other=0).compute()
	DYUfram=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNfram=HTN[jstart:jend+1,istart:iend+1].load()
	DXUfram=DXU[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEfram=HTE[jstart:jend+1,istart:iend+1].load()
	dzfram=dz[0:kmax].load()
	umid=xr.DataArray(data=np.zeros_like(saltfram),coords=saltfram.coords,dims=saltfram.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltfram),coords=saltfram.coords,dims=saltfram.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,i_range_enum):
		#Avg in x first
		umidx=(ufram[:,:,j:j+2,i]*DXUfram[j:j+2,i]+ufram[:,:,j:j+2,i+1]*DXUfram[j:j+2,i+1])/DXUfram[j:j+2,i:i+2].sum(dim='nlon')
		vmidx=(vfram[:,:,j:j+2,i]*DXUfram[j:j+2,i]+vfram[:,:,j:j+2,i+1]*DXUfram[j:j+2,i+1])/DXUfram[j:j+2,i:i+2].sum(dim='nlon')
		#Then in y
		umid[:,:,j,i]=(umidx[:,:,0]*DYUfram[j,i]+umidx[:,:,1]*DYUfram[j+1,i])/DYUfram[j:j+2,i].sum(dim='nlat')
		vmid[:,:,j,i]=-1.0*(vmidx[:,:,0]*DYUfram[j,i]+vmidx[:,:,1]*DYUfram[j+1,i])/DYUfram[j:j+2,i].sum(dim='nlat')
	#Compute the fluxes
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Fram line makes a 45˚ angle with u and v
	#Can only do this if u,v are same size and in the same place!
	#Volume flux
	fram_vol=dzfram*(umid+vmid)*vel_scale*np.sqrt(HTEfram*HTEfram+HTNfram*HTNfram)
	fram_vol_flux=km3yr*fram_vol.sum(dim=('lev','nlat','nlon'))#sum over latitude and lon.
	#Fw flux
	fram_liq=((s_ref-saltfram)/s_ref)*fram_vol
	fram_liq_flux=km3yr*fram_liq.sum(dim=('lev','nlat','nlon'))#sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
	#Fluxes above and below sref
	fram_liq_flux_below_sref=km3yr*fram_liq.where(saltfram<=s_ref,other=0).sum(dim=('lev','nlat','nlon')) 
	fram_liq_flux_above_sref=km3yr*fram_liq.where(saltfram>s_ref,other=0).sum(dim=('lev','nlat','nlon')) 

	##--------Barents Sea Opening (BSO)
	#Liquid flux through the BSO, x=[79,90], y=[355,366]. Indices on full world. This is
	#is a 45˚ line so it doesn't have a complex pattern for the grid cells it passes through
	#The BSO is oriented such that eastward (positive) velocities are into the Arctic, whereas northward (positive)
	#velocities are out of the Arcitc, so multiply v by -1
	#Compute the fw volume relative to 34.8
	istart=79
	iend=90
	jstart=355
	jend=366
	kmax=47 #maximum number of depth levels across the transect
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jstart,jend+1)

	saltbso=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	ubso=u[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	ubso=ubso.where(xr.ufuncs.isfinite(ubso),other=0).compute()
	vbso=v[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	vbso=vbso.where(xr.ufuncs.isfinite(vbso),other=0).compute()
	DYUbso=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNbso=HTN[jstart:jend+1,istart:iend+1].load()
	DXUbso=DXU[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEbso=HTE[jstart:jend+1,istart:iend+1].load()
	dzbso=dz[0:kmax].load()
	umid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,i_range_enum):
		#Avg in x first
		umidx=(ubso[:,:,j:j+2,i]*DXUbso[j:j+2,i]+ubso[:,:,j:j+2,i+1]*DXUbso[j:j+2,i+1])/DXUbso[j:j+2,i:i+2].sum(dim='nlon')
		vmidx=(vbso[:,:,j:j+2,i]*DXUbso[j:j+2,i]+vbso[:,:,j:j+2,i+1]*DXUbso[j:j+2,i+1])/DXUbso[j:j+2,i:i+2].sum(dim='nlon')
		#Then in y
		umid[:,:,j,i]=(umidx[:,:,0]*DYUbso[j,i]+umidx[:,:,1]*DYUbso[j+1,i])/DYUbso[j:j+2,i].sum(dim='nlat')
		vmid[:,:,j,i]=-1.0*(vmidx[:,:,0]*DYUbso[j,i]+vmidx[:,:,1]*DYUbso[j+1,i])/DYUbso[j:j+2,i].sum(dim='nlat')
	#Compute the fluxes
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes a 45˚ angle with u and v
	#Volume flux
	bso_vol=dzbso*(umid+vmid)*vel_scale*np.sqrt(HTEbso*HTEbso+HTNbso*HTNbso)
	bso_vol_flux=km3yr*bso_vol.sum(dim=('lev','nlat','nlon'))#sum over latitude and lon
	#FW flux
	bso_liq=((s_ref-saltbso)/s_ref)*bso_vol
	bso_liq_flux=km3yr*bso_liq.sum(dim=('lev','nlat','nlon'))#sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
	#Fluxes above and below sref
	bso_liq_flux_below_sref=km3yr*bso_liq.where(saltbso<=s_ref,other=0).sum(dim=('lev','nlat','nlon')) 
	bso_liq_flux_above_sref=km3yr*bso_liq.where(saltbso>s_ref,other=0).sum(dim=('lev','nlat','nlon')) 
	print('Liquid fluxes through gateways done!')

	#########4.Compute flux totals [km3/yr]
	#Compute total solid flux
	arctic_solid_flux_total=ber_solid_flux+nares_solid_flux+brw_solid_flux+fram_solid_flux+bso_solid_flux
	#Compute total liquid flux
	arctic_liq_flux_total=ber_liq_flux+nares_liq_flux+brw_liq_flux+fram_liq_flux+bso_liq_flux
	#Compute total liquid flux above and below sref
	arctic_liq_flux_total_below_sref=(ber_liq_flux_below_sref+nares_liq_flux_below_sref+
		brw_liq_flux_below_sref+fram_liq_flux_below_sref+bso_liq_flux_below_sref)
	arctic_liq_flux_total_above_sref=(ber_liq_flux_above_sref+nares_liq_flux_above_sref+
		brw_liq_flux_above_sref+fram_liq_flux_above_sref+bso_liq_flux_above_sref)
	#Compute total volume flux
	arctic_vol_flux_total=ber_vol_flux+nares_vol_flux+brw_vol_flux+fram_vol_flux+bso_vol_flux
	#Compute total storage--ONLY INCLUDE POSITIVE VOLUMES FOR OCEAN STORAGE (salinity below s_ref)
	arctic_storage_total=solid_fw_storage+liq_fw_storage_below_sref


	#########5. Save the output
	#Save everything as a netcdf
	#Make the time series an xarray Dataset
	sv_dims=['time'] #['ensemble','time']
	dvs={'ber_solid_flux':(sv_dims,ber_solid_flux),'nares_solid_flux':(sv_dims,nares_solid_flux),'brw_solid_flux':(sv_dims,brw_solid_flux),
	     'fram_solid_flux':(sv_dims,fram_solid_flux),'bso_solid_flux':(sv_dims,bso_solid_flux),
	     'arctic_solid_flux_total':(sv_dims,arctic_solid_flux_total),'ber_liq_flux':(sv_dims,ber_liq_flux),
	     'ber_liq_flux_above_sref':(sv_dims,ber_liq_flux_above_sref),'ber_liq_flux_below_sref':(sv_dims,ber_liq_flux_below_sref),
	     'nares_liq_flux':(sv_dims,nares_liq_flux), 'nares_liq_flux_above_sref':(sv_dims,nares_liq_flux_above_sref),
	     'nares_liq_flux_below_sref':(sv_dims,nares_liq_flux_below_sref), 'brw_liq_flux':(sv_dims,brw_liq_flux),
	     'brw_liq_flux_above_sref':(sv_dims,brw_liq_flux_above_sref),'brw_liq_flux_below_sref':(sv_dims,brw_liq_flux_below_sref),
	     'fram_liq_flux':(sv_dims,fram_liq_flux), 'fram_liq_flux_above_sref':(sv_dims,fram_liq_flux_above_sref),
	     'fram_liq_flux_below_sref':(sv_dims,fram_liq_flux_below_sref), 'bso_liq_flux':(sv_dims,bso_liq_flux),
	     'bso_liq_flux_above_sref':(sv_dims,bso_liq_flux_above_sref),'bso_liq_flux_below_sref':(sv_dims,bso_liq_flux_below_sref),
	     'ber_vol_flux':(sv_dims,ber_vol_flux),'nares_vol_flux':(sv_dims,nares_vol_flux), 'brw_vol_flux':(sv_dims,brw_vol_flux),
	     'fram_vol_flux':(sv_dims,fram_vol_flux),'bso_vol_flux':(sv_dims,bso_vol_flux),'arctic_liq_flux_total':(sv_dims,arctic_liq_flux_total),
	     'arctic_liq_flux_total_below_sref':(sv_dims,arctic_liq_flux_total_below_sref),
	     'arctic_liq_flux_total_above_sref':(sv_dims,arctic_liq_flux_total_above_sref),'arctic_vol_flux_total':(sv_dims,arctic_vol_flux_total),
	     'solid_fw_storage':(sv_dims,solid_fw_storage),'liq_fw_storage':(sv_dims,liq_fw_storage),
	     'liq_fw_storage_above_sref':(sv_dims,liq_fw_storage_above_sref),'liq_fw_storage_below_sref':(sv_dims,liq_fw_storage_below_sref),
	     'arctic_vol':(arctic_vol)}
	     #'ber_salt_flux':(sv_dims,ber_salt_flux),
	     #'nares_salt_flux':(sv_dims,nares_salt_flux), 'brw_salt_flux':(sv_dims,brw_salt_flux),'fram_salt_flux':(sv_dims,fram_salt_flux),
	     #'bso_salt_flux':(sv_dims,bso_salt_flux), 'arctic_salt_flux_total':(sv_dims,arctic_salt_flux_total),'total_arctic_flux':(sv_dims,total_arctic_flux)}
	ds=xr.Dataset(data_vars=dvs, coords={'time':u.coords['time']}) #you can use any variable here that has a time coord
	#Change units attributes to be km3/yr or km3, check the encoding params and the attributes
	for a in [v for v in ds.variables if 'flux' in v]:
		ds[a].attrs['units']='km3 yr-1'
	for a in [v for v in ds.variables if 'storage' in v]:
		ds[a].attrs['units']='km3'
	ds['arctic_vol'].attrs['units']='km3'
	ds['arctic_vol'].attrs['long_name']='Total Arctic Ocean volume within the gateways boundaries'
	#Save it as a netcdf
	svdirec='/glade/u/home/zanowski/ArcticFW/'
	#Opening in 'a' mode overwrites exsiting variables and 'w' overwrites the whole file
	ds.to_netcdf(path=svdirec+'Arctic_fw_gateways_ts_'+model+'_'+experiment+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')
