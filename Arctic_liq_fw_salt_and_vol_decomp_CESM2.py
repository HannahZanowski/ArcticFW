import xarray as xr
import numpy as np
import get_models as gm


#####Arctic_liq_fw_salt_and_vol_decomp_CESM2.py#######
#Decomposes the anomalous salt and volume contributions to the ocean
#fw flux through the Arctic gateways following Holland et al. 2006
#s=s_mean + s', u=u_mean+u' where u_mean,s_mean are the 2000-2100
#mean salinity and velocity in each grid cell and u',s' are temporal
#departures from that mean. The relevant terms in the integrals for the
#fluxes are then s'*u_mean (salinity contribution) and s_mean*u' (velocity contribution)
institution='NCAR'
model='CESM2'
experiment='ssp585' #ssp126, or ssp585
#Current number of ensemble members:
#Historical: CAM-11
#SSP126: CAM-3
#SSP585: CAM-3
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,12),'ssp126':[4,10,11],
	'ssp585':[4,10,11]}
ens_nums=ens_num_dict[experiment]
nmax=len(ens_nums)

#Set constants and get model output
#Constants
s_ref=34.8 #Reference salinity [PSU]
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Read in ocean and sea ice grid info outside of the for loop--get this from local direcs
simulation, time_period=gm.choose_model_and_sim(model,'historical',ens_mem='001')
direc_o='/glade/collections/cdg/timeseries-cmip6/'+simulation+'/ocn/proc/tseries/month_1/'
xdo=xr.open_dataset(direc_o+simulation+'.pop.h.SALT.185001-201412.nc')
#Ocean
oarea=xdo['TAREA']*0.0001 #T cell area [cm^2], converted to [m^2]
odepth=xdo['HT']*0.01 #T cell ocean depth [cm], converted to [m]
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

#Read in ice/ocean output for each ensemble member and compute the fluxes
for ens_num in ens_nums:
	print('Starting Ensemble Member '+str('%i' %ens_num)) # +'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	#For getting the historical and ssp files
	simhist,miphist,tp_hist=gm.choose_model_and_sim_cmip6(model,'historical',variant)
	simssp, mipssp, tp_ssp=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean directory for CanESM5 for CMIP6
	direc_hist='/glade/collections/cmip/CMIP6/'+miphist+'/'+institution+'/'+model+'/historical/'+variant+'/Omon/'
	direc_ssp='/glade/collections/cmip/CMIP6/'+mipssp+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'

	#Ocean variables
	print('Getting ocean variables')
	#Historical
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_hist+'so/gn/latest/so_Omon_'+simhist+'*201412.nc',combine='by_coords',data_vars='minimal')
	salthist=xdo['so'].sel({'time':slice('2000-01-01',None)}) #ocean salinity in PSU
	xdo=xr.open_mfdataset(direc_hist+'uo/gn/latest/uo_Omon_'+simhist+'*201412.nc',combine='by_coords',data_vars='minimal')
	uhist=xdo['uo'].sel({'time':slice('2000-01-01',None)}) #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_hist+'vo/gn/latest/vo_Omon_'+simhist+'*201412.nc',combine='by_coords',data_vars='minimal')
	vhist=xdo['vo'].sel({'time':slice('2000-01-01',None)}) #ocean meridional velocity in [m/s]
	#ssp
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_ssp+'so/gn/latest/so_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	saltssp=xdo['so'] #ocean salinity in PSU
	xdo=xr.open_mfdataset(direc_ssp+'uo/gn/latest/uo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	ussp=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_ssp+'vo/gn/latest/vo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	vssp=xdo['vo'] #ocean meridional velocity in [m/s]

	#Concatenate
	salt=xr.concat([salthist,saltssp],dim='time')
	u=xr.concat([uhist,ussp],dim='time')
	v=xr.concat([vhist,vssp],dim='time')

	print('Model output read in and ready for computation')

	#########1. Liquid flux decomposition through Arctic Gateways [km3/yr]

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
	
	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltber.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltber-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	ber_area=dzber*DXTber
	ber_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*ber_area
	ber_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*ber_area
	ber_liq_flux_vol_cont=km3yr*ber_vol_cont.sum(dim=('lev','nlon'))
	ber_liq_flux_salt_cont=km3yr*ber_salt_cont.sum(dim=('lev','nlon'))
	#Split for salinities above and below s_ref
	ber_liq_flux_vol_cont_below_sref=km3yr*ber_vol_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','nlon'))
	ber_liq_flux_vol_cont_above_sref=km3yr*ber_vol_cont.where(saltber>s_ref,other=0).sum(dim=('lev','nlon'))
	ber_liq_flux_salt_cont_below_sref=km3yr*ber_salt_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','nlon'))
	ber_liq_flux_salt_cont_above_sref=km3yr*ber_salt_cont.where(saltber>s_ref,other=0).sum(dim=('lev','nlon'))


	##--------Nares Strait
	#Liquid flux through Nares: x=237, y=[376,377] #Indices on full world
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is east (positive) so to keep with the Arctic sign covention, multiply by -1
	#Nares is only 2 grid cells wide, and the walls are vertical
	jstart=376
	jend=377
	i=237
	kmax=11  #Maximum number of depth levels across the transect
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

	#Compute mean and anomalies for salinity and velocity
	umid_mn=umid.mean('time')
	salt_mn=saltnar.mean('time')
	umid_anom=umid-umid_mn
	salt_anom=saltnar-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	nar_area=dznar*DYTnar
	nar_vol_cont=((s_ref-salt_mn)/s_ref)*umid_anom*nar_area
	nar_salt_cont=(-1.0*salt_anom/s_ref)*umid_mn*nar_area
	nar_liq_flux_vol_cont=km3yr*nar_vol_cont.sum(dim=('lev','nlat'))
	nar_liq_flux_salt_cont=km3yr*nar_salt_cont.sum(dim=('lev','nlat'))
	#Split for salinities above and below s_ref
	nar_liq_flux_vol_cont_below_sref=km3yr*nar_vol_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','nlat'))
	nar_liq_flux_vol_cont_above_sref=km3yr*nar_vol_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','nlat'))
	nar_liq_flux_salt_cont_below_sref=km3yr*nar_salt_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','nlat'))
	nar_liq_flux_salt_cont_above_sref=km3yr*nar_salt_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','nlat'))


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

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltbrw.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltbrw-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	brw_area=dzbrw*DXTbrw
	brw_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*brw_area
	brw_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*brw_area
	brw_liq_flux_vol_cont=km3yr*brw_vol_cont.sum(dim=('lev','nlon'))
	brw_liq_flux_salt_cont=km3yr*brw_salt_cont.sum(dim=('lev','nlon'))
	#Split for salinities above and below s_ref
	brw_liq_flux_vol_cont_below_sref=km3yr*brw_vol_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','nlon'))
	brw_liq_flux_vol_cont_above_sref=km3yr*brw_vol_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','nlon'))
	brw_liq_flux_salt_cont_below_sref=km3yr*brw_salt_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','nlon'))
	brw_liq_flux_salt_cont_above_sref=km3yr*brw_salt_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','nlon'))

	
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
	#Can only do this if u,v arrays are same size and values are colocated
	fram_vel=(umid+vmid)*vel_scale
	fram_area=dzfram*np.sqrt(HTEfram*HTEfram+HTNfram*HTNfram)
	#Compute mean and anomalies for salinity and velocity
	vel_mn=fram_vel.mean('time')
	salt_mn=saltfram.mean('time')
	vel_anom=fram_vel-vel_mn
	salt_anom=saltfram-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	fram_vol_cont=((s_ref-salt_mn)/s_ref)*vel_anom*fram_area
	fram_salt_cont=(-1.0*salt_anom/s_ref)*vel_mn*fram_area
	fram_liq_flux_vol_cont=km3yr*fram_vol_cont.sum(dim=('lev','nlon','nlat'))
	fram_liq_flux_salt_cont=km3yr*fram_salt_cont.sum(dim=('lev','nlon','nlat'))
	#Split for salinities above and below s_ref
	fram_liq_flux_vol_cont_below_sref=km3yr*fram_vol_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	fram_liq_flux_vol_cont_above_sref=km3yr*fram_vol_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	fram_liq_flux_salt_cont_below_sref=km3yr*fram_salt_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	fram_liq_flux_salt_cont_above_sref=km3yr*fram_salt_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','nlon','nlat'))


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
	#Can only do this if u,v arrays are same size and values are colocated
	bso_vel=(umid+vmid)*vel_scale
	bso_area=dzbso*np.sqrt(HTEbso*HTEbso+HTNbso*HTNbso)
	#Compute mean and anomalies for salinity and velocity
	vel_mn=bso_vel.mean('time')
	salt_mn=saltbso.mean('time')
	vel_anom=bso_vel-vel_mn
	salt_anom=saltbso-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	bso_vol_cont=((s_ref-salt_mn)/s_ref)*vel_anom*bso_area
	bso_salt_cont=(-1.0*salt_anom/s_ref)*vel_mn*bso_area
	bso_liq_flux_vol_cont=km3yr*bso_vol_cont.sum(dim=('lev','nlon','nlat'))
	bso_liq_flux_salt_cont=km3yr*bso_salt_cont.sum(dim=('lev','nlon','nlat'))
	#Split for salinities above and below s_ref
	bso_liq_flux_vol_cont_below_sref=km3yr*bso_vol_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	bso_liq_flux_vol_cont_above_sref=km3yr*bso_vol_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	bso_liq_flux_salt_cont_below_sref=km3yr*bso_salt_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','nlon','nlat'))
	bso_liq_flux_salt_cont_above_sref=km3yr*bso_salt_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','nlon','nlat'))

	##--------Davis Strait
	#Liquid flux through Davis Strait: x=[293,304], y=364 #Indices on full world
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	j=364
	istart=293
	iend=304
	kmax=36 #maximum number of depth levels across the transect
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

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltdav.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltdav-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	dav_area=dzdav*DXTdav
	dav_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*dav_area
	dav_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*dav_area
	dav_liq_flux_vol_cont=km3yr*dav_vol_cont.sum(dim=('lev','nlon'))
	dav_liq_flux_salt_cont=km3yr*dav_salt_cont.sum(dim=('lev','nlon'))
	#Split for salinities above and below s_ref
	dav_liq_flux_vol_cont_below_sref=km3yr*dav_vol_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','nlon'))
	dav_liq_flux_vol_cont_above_sref=km3yr*dav_vol_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','nlon'))
	dav_liq_flux_salt_cont_below_sref=km3yr*dav_salt_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','nlon'))
	dav_liq_flux_salt_cont_above_sref=km3yr*dav_salt_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','nlon'))

	print('Liquid flux decomposition done!')

	#########2. Save the output
	sv_dims=['time'] #['ensemble','time']
	dvs={'ber_liq_flux_vol_cont':(sv_dims,ber_liq_flux_vol_cont),'ber_liq_flux_salt_cont':(sv_dims,ber_liq_flux_salt_cont),
	     'ber_liq_flux_vol_cont_above_sref':(sv_dims,ber_liq_flux_vol_cont_above_sref),
	     'ber_liq_flux_vol_cont_below_sref':(sv_dims,ber_liq_flux_vol_cont_below_sref),
	     'ber_liq_flux_salt_cont_above_sref':(sv_dims,ber_liq_flux_salt_cont_above_sref),
	     'ber_liq_flux_salt_cont_below_sref':(sv_dims,ber_liq_flux_salt_cont_below_sref),
	     'nares_liq_flux_vol_cont':(sv_dims,nar_liq_flux_vol_cont),'nares_liq_flux_salt_cont':(sv_dims,nar_liq_flux_salt_cont),
	     'nares_liq_flux_vol_cont_above_sref':(sv_dims,nar_liq_flux_vol_cont_above_sref),
	     'nares_liq_flux_vol_cont_below_sref':(sv_dims,nar_liq_flux_vol_cont_below_sref),
	     'nares_liq_flux_salt_cont_above_sref':(sv_dims,nar_liq_flux_salt_cont_above_sref),
	     'nares_liq_flux_salt_cont_below_sref':(sv_dims,nar_liq_flux_salt_cont_below_sref),
	     'brw_liq_flux_vol_cont':(sv_dims,brw_liq_flux_vol_cont),'brw_liq_flux_salt_cont':(sv_dims,brw_liq_flux_salt_cont),
	     'brw_liq_flux_vol_cont_above_sref':(sv_dims,brw_liq_flux_vol_cont_above_sref),
	     'brw_liq_flux_vol_cont_below_sref':(sv_dims,brw_liq_flux_vol_cont_below_sref),
	     'brw_liq_flux_salt_cont_above_sref':(sv_dims,brw_liq_flux_salt_cont_above_sref),
	     'brw_liq_flux_salt_cont_below_sref':(sv_dims,brw_liq_flux_salt_cont_below_sref),
	     'fram_liq_flux_vol_cont':(sv_dims,fram_liq_flux_vol_cont),'fram_liq_flux_salt_cont':(sv_dims,fram_liq_flux_salt_cont),
	     'fram_liq_flux_vol_cont_above_sref':(sv_dims,fram_liq_flux_vol_cont_above_sref),
	     'fram_liq_flux_vol_cont_below_sref':(sv_dims,fram_liq_flux_vol_cont_below_sref),
	     'fram_liq_flux_salt_cont_above_sref':(sv_dims,fram_liq_flux_salt_cont_above_sref),
	     'fram_liq_flux_salt_cont_below_sref':(sv_dims,fram_liq_flux_salt_cont_below_sref),
	     'bso_liq_flux_vol_cont':(sv_dims,bso_liq_flux_vol_cont),'bso_liq_flux_salt_cont':(sv_dims,bso_liq_flux_salt_cont),
	     'bso_liq_flux_vol_cont_above_sref':(sv_dims,bso_liq_flux_vol_cont_above_sref),
	     'bso_liq_flux_vol_cont_below_sref':(sv_dims,bso_liq_flux_vol_cont_below_sref),
	     'bso_liq_flux_salt_cont_above_sref':(sv_dims,bso_liq_flux_salt_cont_above_sref),
	     'bso_liq_flux_salt_cont_below_sref':(sv_dims,bso_liq_flux_salt_cont_below_sref),
	     'davis_liq_flux_vol_cont':(sv_dims,dav_liq_flux_vol_cont),'davis_liq_flux_salt_cont':(sv_dims,dav_liq_flux_salt_cont),
	     'davis_liq_flux_vol_cont_above_sref':(sv_dims,dav_liq_flux_vol_cont_above_sref),
	     'davis_liq_flux_vol_cont_below_sref':(sv_dims,dav_liq_flux_vol_cont_below_sref),
	     'davis_liq_flux_salt_cont_above_sref':(sv_dims,dav_liq_flux_salt_cont_above_sref),
	     'davis_liq_flux_salt_cont_below_sref':(sv_dims,dav_liq_flux_salt_cont_below_sref)}
	ds=xr.Dataset(data_vars=dvs, coords={'time':u.coords['time']}) #you can use any variable here that has a time coord
	#Change units attributes to be km3/yr, check the encoding params and the attributes
	for a in [v for v in ds.variables if 'flux' in v]:
		ds[a].attrs['units']='km3 yr-1'
	#Save it as a netcdf
	svdirec='/glade/u/home/zanowski/ArcticFW/'
	#Opening in 'a' mode overwrites exsiting variables and 'w' overwrites the whole file
	ds.to_netcdf(path=svdirec+'Arctic_fw_flux_salt_and_vol_decomp_ts_'+model+'_'+experiment.lower()+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')

