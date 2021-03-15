import xarray as xr
import numpy as np
import os
import get_models as gm

#####Arctic_liq_fw_salt_and_vol_decomp_MIROC6.py#######
#Decomposes the anomalous salt and volume contributions to the ocean
#fw flux through the Arctic gateways following Holland et al. 2006
#s=s_mean + s', u=u_mean+u' where u_mean,s_mean are the 2000-2100
#mean salinity and velocity in each grid cell and u',s' are temporal
#departures from that mean. The relevant terms in the integrals for the
#fluxes are then s'*u_mean (salinity contribution) and s_mean*u' (velocity contribution)
institution='MIROC'
model='MIROC6'
experiment='ssp585' #any of the ssps as well
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,11),'ssp126':[1,2,3],
                'ssp585':[1,2,3]}
nmax=int(len(ens_num_dict[experiment]))
ens_nums=ens_num_dict[experiment]

#Constants
s_ref=34.8 #Reference salinity [ppt]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Ocean grid variables
#Grid distance file
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+model+'_grid_distances.nc')
#get x and y distances: tracer x and y, y distances at v points and x distances at u points
DXV=xdf['dxv'] #x-spacing, u,v points
DYV=xdf['dyv'] #y-spacing, u,v points
DXT=xdf['dxt'] #x-spacing, tracer pts
DYT=xdf['dyt'] #y-spacing, tracer pts
HTN=xdf['htn'] #T-cell widths on North sides
HTE=xdf['hte'] #T-cell widths on East sides
HUS=xdf['hus'] #U-cell widths on South sides
HUW=xdf['huw'] #U-cell widths on West sides

#Loop over ensemble members for everything, including data read in
#NOTE:version numbers change with each ensemble member,
#so to account for this the code below sorts
#the contents of os.listdir() at the version level
#and chooses the last element, which should be the latest version date.
for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	print('Getting ocean variables')
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	#For getting the historical and ssp files
	simhist,miphist,tp_hist=gm.choose_model_and_sim_cmip6(model,'historical',variant)
	simssp, mipssp, tp_ssp=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean variables
	#Ocean directory for MIROC6 for CMIP6
	direc_hist='/glade/collections/cmip/CMIP6/'+miphist+'/'+institution+'/'+model+'/historical/'+variant+'/Omon/'
	direc_ssp='/glade/collections/cmip/CMIP6/'+mipssp+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	#Ocean variables
	salt_paths=[direc_hist+'so/gn/'+sorted(os.listdir(direc_hist+'so/gn/'))[-1]+'/so/so_Omon_'+simhist+'*.nc',
		direc_ssp+'so/gn/'+sorted(os.listdir(direc_ssp+'so/gn/'))[-1]+'/so/so_Omon_'+simssp+'*.nc']
	u_paths=[direc_hist+'uo/gn/'+sorted(os.listdir(direc_hist+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simhist+'*.nc',
                direc_ssp+'uo/gn/'+sorted(os.listdir(direc_ssp+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simssp+'*.nc']
	v_paths=[direc_hist+'vo/gn/'+sorted(os.listdir(direc_hist+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simhist+'*.nc',
                direc_ssp+'vo/gn/'+sorted(os.listdir(direc_ssp+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simssp+'*.nc']
	#Get salinity and velocities as well as the variables for computing layer thickness
	#historical
	xdo=xr.open_mfdataset(salt_paths[0],combine='by_coords',data_vars='minimal')
	salthist=xdo['so'].sel({'time':slice('2000-01-01',None)}) #ocean salinity in ppt
	etahist=xdo['eta'].sel({'time':slice('2000-01-01',None)})
	sigma_bnds=xdo['sigma_bnds']
	zlev_bnds=xdo['zlev_bnds']
	depth_c=xdo['depth_c']
	depth=xdo['depth']
	nsigma=xdo['nsigma']
	xdo=xr.open_mfdataset(u_paths[0],combine='by_coords',data_vars='minimal')
	uhist=xdo['uo'].sel({'time':slice('2000-01-01',None)}) #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(v_paths[0],combine='by_coords',data_vars='minimal')
	vhist=xdo['vo'].sel({'time':slice('2000-01-01',None)}) #ocean meridional velocity in [m/s]
	#ssp
	xdo=xr.open_mfdataset(salt_paths[1],combine='by_coords',data_vars='minimal')
	saltssp=xdo['so'] #.sel({'time':slice(None,'2100-01-01')}) #ocean salinity in ppt
	etassp=xdo['eta'] #.sel({'time':slice(None,'2100-01-01')})
	xdo=xr.open_mfdataset(u_paths[1],combine='by_coords',data_vars='minimal')
	ussp=xdo['uo'] #.sel({'time':slice(None,'2100-01-01')}) #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(v_paths[1],combine='by_coords',data_vars='minimal')
	vssp=xdo['vo'] #.sel({'time':slice(None,'2100-01-01')}) #ocean meridional velocity in [m/s]

	#Concatenate
	salt=xr.concat([salthist,saltssp],dim='time')
	u=xr.concat([uhist,ussp],dim='time')
	v=xr.concat([vhist,vssp],dim='time')
	eta=xr.concat([etahist,etassp],dim='time')	
	#Compute dz based on the function provided in the netcdf files:
	#k <= nsigma: z(n,k,j,i) = eta(n,j,i) + sigma(k)*(min(depth_c,depth(j,i))+eta(n,j,i))
	#k > nsigma: z(n,k,j,i) = zlev(k)
	dz=xr.DataArray(data=np.zeros(salt.shape),coords=salt.coords,dims=salt.dims)
	dmin=xr.ufuncs.minimum(eta+depth,depth_c) #no depth dependence
	for k in range(0,len(salt[0])):
		if k < nsigma: 
			#zu=eta+sigma_bnds[k,0]*(xr.ufuncs.minimum(eta+depth,depth_c))
			#zl=eta+sigma_bnds[k,1]*(xr.ufuncs.minimum(eta+depth,depth_c))
			#dz[:,k,:,:]=zu-zl
			#Simplified version of the above
			dz[:,k,:,:]=(sigma_bnds[k,0]-sigma_bnds[k,1])*dmin
		else:
			dz[:,k,:,:]=zlev_bnds[k,0]-zlev_bnds[k,1]

	print('Model output read in and ready for computation')


	#########1. Liquid flux decomposition through Arctic Gateways [km3/yr]

	##--------Bering Strait
	#Liquid flux through Bering Strait: x=[129:132], y=212 #Indices on the full world map
	#Compute the fw volume relative to 34.8. 
	#Flow out of Bering is southward so the sign convention is correct
	j=212
	istart=129
	iend=132
	kmax=10 #maximum number of depth levels across the transect
	saltber=salt[:,0:kmax,j,istart:iend+1].load()
	vber=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vber=vber.where(xr.ufuncs.isfinite(vber),other=0).compute()
	DYVber=DYV[j-1:j+1,istart-1:iend+1].load()
	DXVber=HUS[j,istart-1:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	dzber=dz[:,0:kmax,j,istart:iend+1].load()
	#Average meridional velocities in y first
	vber=(vber[:,:,0]*DYVber[0]+vber[:,:,1]*DYVber[1])/(DYVber[0]+DYVber[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltber),coords=saltber.coords,dims=saltber.dims)
	vmid[:,:,0]=(vber[:,:,0]*DXVber[0]+vber[:,:,1]*DXVber[1])/sum(DXVber[0:2])
	vmid[:,:,1]=(vber[:,:,1]*DXVber[1]+vber[:,:,2]*DXVber[2])/sum(DXVber[1:])

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltber.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltber-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	ber_area=dzber*DXTber
	ber_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*ber_area
	ber_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*ber_area
	ber_liq_flux_vol_cont=km3yr*ber_vol_cont.sum(dim=('lev','x'))
	ber_liq_flux_salt_cont=km3yr*ber_salt_cont.sum(dim=('lev','x'))
	#Split for salinities above and below s_ref
	ber_liq_flux_vol_cont_below_sref=km3yr*ber_vol_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','x'))
	ber_liq_flux_vol_cont_above_sref=km3yr*ber_vol_cont.where(saltber>s_ref,other=0).sum(dim=('lev','x'))
	ber_liq_flux_salt_cont_below_sref=km3yr*ber_salt_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','x'))
	ber_liq_flux_salt_cont_above_sref=km3yr*ber_salt_cont.where(saltber>s_ref,other=0).sum(dim=('lev','x'))


	##--------Nares Strait
	#Liquid flux through Nares: x=[241,243], y=[234,236]. Indices on the full world map
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is west and south (both negative) so the sign convention is correct
	jstart=234
	jend=236
	istart=241
	iend=243
	kmax=22  #Maximum number of depth levels across the transect
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)

	saltnar=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	unar=u[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	unar=unar.where(xr.ufuncs.isfinite(unar),other=0).compute()
	vnar=v[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	vnar=vnar.where(xr.ufuncs.isfinite(vnar),other=0).compute()
	DYVnar=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNnar=HTN[jstart:jend+1,istart:iend+1].load()
	DXVnar=DXV[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEnar=HTE[jstart:jend+1,istart:iend+1].load()
	dznar=dz[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(saltnar),coords=saltnar.coords,dims=saltnar.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltnar),coords=saltnar.coords,dims=saltnar.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
		#Avg in x first
		umidx=(unar[:,:,j:j+2,i]*DXVnar[j:j+2,i]+unar[:,:,j:j+2,i+1]*DXVnar[j:j+2,i+1])/DXVnar[j:j+2,i:i+2].sum(dim='x')
		vmidx=(vnar[:,:,j:j+2,i]*DXVnar[j:j+2,i]+vnar[:,:,j:j+2,i+1]*DXVnar[j:j+2,i+1])/DXVnar[j:j+2,i:i+2].sum(dim='x')
		#Then in y
		umid[:,:,j,i]=(umidx[:,:,0]*DYVnar[j,i]+umidx[:,:,1]*DYVnar[j+1,i])/DYVnar[j:j+2,i].sum(dim='y')
		vmid[:,:,j,i]=(vmidx[:,:,0]*DYVnar[j,i]+vmidx[:,:,1]*DYVnar[j+1,i])/DYVnar[j:j+2,i].sum(dim='y')
	#Can only do this if u,v arrays are same sizes and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Nares line makes a 45˚ angle with u and v
	nar_vel=(umid+vmid)*vel_scale
	nar_area=dznar*np.sqrt(HTEnar*HTEnar+HTNnar*HTNnar)
	#Compute mean and anomalies for salinity and velocity
	vel_mn=nar_vel.mean('time')
	salt_mn=saltnar.mean('time')
	vel_anom=nar_vel-vel_mn
	salt_anom=saltnar-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	nar_vol_cont=((s_ref-salt_mn)/s_ref)*vel_anom*nar_area
	nar_salt_cont=(-1.0*salt_anom/s_ref)*vel_mn*nar_area
	nar_liq_flux_vol_cont=km3yr*nar_vol_cont.sum(dim=('lev','x','y'))
	nar_liq_flux_salt_cont=km3yr*nar_salt_cont.sum(dim=('lev','x','y'))
	#Split for salinities above and below s_ref
	nar_liq_flux_vol_cont_below_sref=km3yr*nar_vol_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','x','y'))
	nar_liq_flux_vol_cont_above_sref=km3yr*nar_vol_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','x','y'))
	nar_liq_flux_salt_cont_below_sref=km3yr*nar_salt_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','x','y'))
	nar_liq_flux_salt_cont_above_sref=km3yr*nar_salt_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','x','y'))


	##--------Barrow Strait
	#Liquid fw flux through Barrow Strait: x=[212,218], y=242 #Indices on full world map
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	j=242
	istart=212
	iend=218
	kmax=25 #maximum number of depth levels across the transect
	irange=np.arange(istart,iend+1)
	saltbrw=salt[:,0:kmax,j,istart:iend+1].load()
	vbrw=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vbrw=vbrw.where(xr.ufuncs.isfinite(vbrw),other=0).compute()
	DYVbrw=DYV[j-1:j+1,istart-1:iend+1].load()
	DXVbrw=HUS[j,istart-1:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	dzbrw=dz[:,0:kmax,j,istart:iend+1].load()
	#Average meridional velocities in y first
	vbrw=(vbrw[:,:,0]*DYVbrw[0]+vbrw[:,:,1]*DYVbrw[1])/(DYVbrw[0]+DYVbrw[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltbrw),coords=saltbrw.coords,dims=saltbrw.dims)
	for i in range(0,len(irange)):
		vmid[:,:,i]=(vbrw[:,:,i]*DXVbrw[i]+vbrw[:,:,i+1]*DXVbrw[i+1])/sum(DXVbrw[i:i+2])

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltbrw.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltbrw-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	brw_area=dzbrw*DXTbrw
	brw_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*brw_area
	brw_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*brw_area
	brw_liq_flux_vol_cont=km3yr*brw_vol_cont.sum(dim=('lev','x'))
	brw_liq_flux_salt_cont=km3yr*brw_salt_cont.sum(dim=('lev','x'))
	#Split for salinities above and below s_ref
	brw_liq_flux_vol_cont_below_sref=km3yr*brw_vol_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','x'))
	brw_liq_flux_vol_cont_above_sref=km3yr*brw_vol_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','x'))
	brw_liq_flux_salt_cont_below_sref=km3yr*brw_salt_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','x'))
	brw_liq_flux_salt_cont_above_sref=km3yr*brw_salt_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','x'))


	##--------Fram Strait
	#Liquid flux through Fram Strait, x=[280,295], y=238. Indices on full world.
	#Fram is oriented such that northward (positive)
	#velocities are into the Arcitc, so the sign convention is correct
	#Compute the fw volume relative to 34.8
	istart=280
	iend=295
	j=238
	kmax=53 #maximum number of depth levels across the transect
	irange=np.arange(istart,iend+1)
	saltfram=salt[:,0:kmax,j,istart:iend+1].load()
	vfram=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vfram=vfram.where(xr.ufuncs.isfinite(vfram),other=0).compute()
	DYVfram=DYV[j-1:j+1,istart-1:iend+1].load()
	DXVfram=HUS[j,istart-1:iend+1].load()
	DXTfram=DXT[j,istart:iend+1].load()
	dzfram=dz[:,0:kmax,j,istart:iend+1].load()
	#Average meridional velocities in y first
	vfram=(vfram[:,:,0]*DYVfram[0]+vfram[:,:,1]*DYVfram[1])/(DYVfram[0]+DYVfram[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltfram),coords=saltfram.coords,dims=saltfram.dims)
	for i in range(0,len(irange)):
		vmid[:,:,i]=(vfram[:,:,i]*DXVfram[i]+vfram[:,:,i+1]*DXVfram[i+1])/sum(DXVfram[i:i+2])

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltfram.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltfram-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	fram_area=dzfram*DXTfram
	fram_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*fram_area
	fram_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*fram_area
	fram_liq_flux_vol_cont=km3yr*fram_vol_cont.sum(dim=('lev','x'))
	fram_liq_flux_salt_cont=km3yr*fram_salt_cont.sum(dim=('lev','x'))
	#Split for salinities above and below s_ref
	fram_liq_flux_vol_cont_below_sref=km3yr*fram_vol_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','x'))
	fram_liq_flux_vol_cont_above_sref=km3yr*fram_vol_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','x'))
	fram_liq_flux_salt_cont_below_sref=km3yr*fram_salt_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','x'))
	fram_liq_flux_salt_cont_above_sref=km3yr*fram_salt_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','x'))

	##--------Barents Sea Opening (BSO)
	#Liquid flux through the BSO, x=[304,317], y=[223,236]. Indices on full world.
	#This is is more or less a 45˚ line so it doesn't have a complex pattern
	#for the grid cells it passes through
	#The BSO is oriented such that eastward and northward (positive) velocities are
	#into the Arctic
	#Compute the fw volume relative to 34.8
	istart=304
	iend=317
	jstart=223
	jend=236
	kmax=29 #maximum number of depth levels across the transect
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)

	saltbso=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	ubso=u[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	ubso=ubso.where(xr.ufuncs.isfinite(ubso),other=0).compute()
	vbso=v[:,0:kmax,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	vbso=vbso.where(xr.ufuncs.isfinite(vbso),other=0).compute()
	DYVbso=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNbso=HTN[jstart:jend+1,istart:iend+1].load()
	DXVbso=DXV[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEbso=HTE[jstart:jend+1,istart:iend+1].load()
	dzbso=dz[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
		#Avg in x first
		umidx=(ubso[:,:,j:j+2,i]*DXVbso[j:j+2,i]+ubso[:,:,j:j+2,i+1]*DXVbso[j:j+2,i+1])/DXVbso[j:j+2,i:i+2].sum(dim='x')
		vmidx=(vbso[:,:,j:j+2,i]*DXVbso[j:j+2,i]+vbso[:,:,j:j+2,i+1]*DXVbso[j:j+2,i+1])/DXVbso[j:j+2,i:i+2].sum(dim='x')
		#Then in y
		umid[:,:,j,i]=(umidx[:,:,0]*DYVbso[j,i]+umidx[:,:,1]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')
		vmid[:,:,j,i]=(vmidx[:,:,0]*DYVbso[j,i]+vmidx[:,:,1]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')
	#Can only do this if u,v arrays are same size and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Fram line makes a 45˚ angle with u and v
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
	bso_liq_flux_vol_cont=km3yr*bso_vol_cont.sum(dim=('lev','x','y'))
	bso_liq_flux_salt_cont=km3yr*bso_salt_cont.sum(dim=('lev','x','y'))
	#Split for salinities above and below s_ref
	bso_liq_flux_vol_cont_below_sref=km3yr*bso_vol_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','x','y'))
	bso_liq_flux_vol_cont_above_sref=km3yr*bso_vol_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','x','y'))
	bso_liq_flux_salt_cont_below_sref=km3yr*bso_salt_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','x','y'))
	bso_liq_flux_salt_cont_above_sref=km3yr*bso_salt_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','x','y'))	

	##--------Davis Strait
	#Liquid flux through Davis Strait: x=[239:245], y=213 #Indices on the full world map
	#Compute the fw volume relative to 34.8.
	#Flow out of Davis is southward so the sign convention is correct
	j=213
	istart=239
	iend=245
	kmax=41 #maximum number of depth levels across the transect
	irange=np.arange(istart,iend+1)
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	vdav=v[:,0:kmax,j-1:j+1,istart-1:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	DYVdav=DYV[j-1:j+1,istart-1:iend+1].load()
	DXVdav=HUS[j,istart-1:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[:,0:kmax,j,istart:iend+1].load()
	#Average meridional velocities in y first
	vdav=(vdav[:,:,0]*DYVdav[0]+vdav[:,:,1]*DYVdav[1])/(DYVdav[0]+DYVdav[1])
	#Average in x
	vmid=xr.DataArray(data=np.zeros_like(saltdav),coords=saltdav.coords,dims=saltdav.dims)
	for i in range(0,len(irange)):
		vmid[:,:,i]=(vdav[:,:,i]*DXVdav[i]+vdav[:,:,i+1]*DXVdav[i+1])/sum(DXVdav[i:i+2])

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltdav.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltdav-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	dav_area=dzdav*DXTdav
	dav_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*dav_area
	dav_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*dav_area
	dav_liq_flux_vol_cont=km3yr*dav_vol_cont.sum(dim=('lev','x'))
	dav_liq_flux_salt_cont=km3yr*dav_salt_cont.sum(dim=('lev','x'))
	#Split for salinities above and below s_ref
	dav_liq_flux_vol_cont_below_sref=km3yr*dav_vol_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','x'))
	dav_liq_flux_vol_cont_above_sref=km3yr*dav_vol_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','x'))
	dav_liq_flux_salt_cont_below_sref=km3yr*dav_salt_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','x'))
	dav_liq_flux_salt_cont_above_sref=km3yr*dav_salt_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','x'))

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
	ds.to_netcdf(path=svdirec+'Arctic_fw_flux_salt_and_vol_decomp_ts_'+model+'_'+experiment+'_'+str('%02i' %ens_num)+'.nc')
	print('Output for ensemble member '+str('%02i' %ens_num)+' saved!')	

