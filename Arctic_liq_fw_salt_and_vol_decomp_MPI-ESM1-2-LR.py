import xarray as xr
import numpy as np
import os
import get_models as gm

####Arctic_liq_fw_salt_and_vol_decomp_MPI-ESM1-2-LR.py#######
#Decomposes the anomalous salt and volume contributions to the ocean
#fw flux through the Arctic gateways following Holland et al. 2006
#s=s_mean + s', u=u_mean+u' where u_mean,s_mean are the 2000-2100
#mean salinity and velocity in each grid cell and u',s' are temporal
#departures from that mean. The relevant terms in the integrals for the
#fluxes are then s'*u_mean (salinity contribution) and s_mean*u' (velocity contribution)
institution='MPI-M'
model='MPI-ESM1-2-LR'
experiment='ssp585' #ssp126 or ssp585
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,11),'ssp126':np.arange(1,11),
		'ssp585':np.arange(1,11)}
nmax=int(len(ens_num_dict[experiment]))
ens_nums=ens_num_dict[experiment]

#Constants.
s_ref=34.8 #Reference salinity [PSU]
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

####IMPORTANT:
#MPI's grid is upside down s.t. index zero corresponds to the
#farthest point in the NH and index -1 corresponds to the
#farthest point in the SH. Everything plots upside down, etc.
#The mask, grid distances, and indices for this code are all
#in the frame of this original indexing convention--I have
#NOT flipped the lats so that the NH points are last instead
#of first

for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)) #+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f1'

	#For getting the historical and ssp files
	simhist,miphist,tp_hist=gm.choose_model_and_sim_cmip6(model,'historical',variant)
	simssp, mipssp, tp_ssp=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean variables
	#Ocean directory for MPI-ESM1-2-LR for CMIP6
	direc_hist='/glade/collections/cmip/CMIP6/'+miphist+'/'+institution+'/'+model+'/historical/'+variant+'/Omon/'
	direc_ssp='/glade/collections/cmip/CMIP6/'+mipssp+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'

	#Ocean variables
	print('Getting ocean variables')
	#Historical
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_hist+'so/gn/'+sorted(os.listdir(direc_hist+'so/gn/'))[-1]+'/so/so_Omon_'+simhist+'*.nc',combine='by_coords',data_vars='minimal')
	salthist=xdo['so'].sel({'time':slice('2000-01-01',None)}) #ocean salinity in ppt
	xdo=xr.open_mfdataset(direc_hist+'uo/gn/'+sorted(os.listdir(direc_hist+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simhist+'*.nc',combine='by_coords',data_vars='minimal')
	uhist=xdo['uo'].sel({'time':slice('2000-01-01',None)}) #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_hist+'vo/gn/'+sorted(os.listdir(direc_hist+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simhist+'*.nc',combine='by_coords',data_vars='minimal')
	vhist=xdo['vo'].sel({'time':slice('2000-01-01',None)}) #ocean meridional velocity in [m/s]
	xdo=xr.open_mfdataset(direc_hist+'thkcello/gn/'+sorted(os.listdir(direc_hist+'thkcello/gn/'))[-1]+'/thkcello/thkcello_Omon_'+simhist+'*.nc',combine='by_coords',data_vars='minimal')
	dzhist=xdo['thkcello'].sel({'time':slice('2000-01-01',None)}) #ocean cell thicknesses in [m]
	#ssp
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_ssp+'so/gn/'+sorted(os.listdir(direc_ssp+'so/gn/'))[-1]+'/so/so_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	saltssp=xdo['so'] #ocean salinity in ppt
	xdo=xr.open_mfdataset(direc_ssp+'uo/gn/'+sorted(os.listdir(direc_ssp+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	ussp=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_ssp+'vo/gn/'+sorted(os.listdir(direc_ssp+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	vssp=xdo['vo'] #ocean meridional velocity in [m/s]
	xdo=xr.open_mfdataset(direc_ssp+'thkcello/gn/'+sorted(os.listdir(direc_ssp+'thkcello/gn/'))[-1]+'/thkcello/thkcello_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	dzssp=xdo['thkcello'] #ocean cell thicknesses in [m]

	#Concatenate
	salt=xr.concat([salthist,saltssp],dim='time')
	u=xr.concat([uhist,ussp],dim='time')
	v=xr.concat([vhist,vssp],dim='time')
	dz=xr.concat([dzhist,dzssp],dim='time')

	print('Model output read in and ready for computation')
	

	####-------------Compute the liquid flux decomposition through the gateways-------------####

	#-------Bering Strait-------#
	#Liquid flux through Bering Strait: x=[6,8], y=[85] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=85
	istart=6
	iend=8
	kmax=4 #maximum number of depth levels across the transect
	vber=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vber=vber.where(xr.ufuncs.isfinite(vber),other=0).compute()
	saltber=salt[:,0:kmax,j,istart:iend+1].load()
	DYVber=DYV[j-1:j+1,istart:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	dzber=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(saltber),coords=saltber.coords,dims=saltber.dims)
	vmid=(vber[:,:,0]*DYVber[0,:]+vber[:,:,1]*DYVber[1,:])/DYVber.sum(dim='j')

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltber.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltber-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	ber_area=dzber*DXTber
	ber_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*ber_area
	ber_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*ber_area
	ber_liq_flux_vol_cont=km3yr*ber_vol_cont.sum(dim=('lev','i'))
	ber_liq_flux_salt_cont=km3yr*ber_salt_cont.sum(dim=('lev','i'))
	#Split for salinities above and below s_ref
	ber_liq_flux_vol_cont_below_sref=km3yr*ber_vol_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','i'))
	ber_liq_flux_vol_cont_above_sref=km3yr*ber_vol_cont.where(saltber>s_ref,other=0).sum(dim=('lev','i'))
	ber_liq_flux_salt_cont_below_sref=km3yr*ber_salt_cont.where(saltber<=s_ref,other=0).sum(dim=('lev','i'))
	ber_liq_flux_salt_cont_above_sref=km3yr*ber_salt_cont.where(saltber>s_ref,other=0).sum(dim=('lev','i'))

	
	#-------Nares Strait-------#
	#Liquid flux through Nares: x=4, y=[31,32] full world indices
	#u is east of tracer for a given tracer index
	#Compute the fw volume relative to 34.8
	#Flow into Nares is west so multiply u by -1
	i=4
	jstart=31
	jend=32
	kmax=13  #Maximum number of deptth levels across the transect
	unar=u[:,0:kmax,jstart:jend+1,i-1:i+1]
	unar=unar.where(xr.ufuncs.isfinite(unar),other=0).compute()
	saltnar=salt[:,0:kmax,jstart:jend+1,i].load()
	DXUnar=DXU[jstart:jend+1,i-1:i+1].load()
	DYTnar=DYT[jstart:jend+1,i].load()
	dznar=dz[:,0:kmax,jstart:jend+1,i].load()
	umid=xr.DataArray(data=np.zeros_like(saltnar),coords=saltnar.coords,dims=saltnar.dims)
	umid=-1.0*(unar[:,:,:,0]*DXUnar[:,0]+unar[:,:,:,1]*DXUnar[:,1])/DXUnar.sum(dim='i')

	#Compute mean and anomalies for salinity and velocity
	umid_mn=umid.mean('time')
	salt_mn=saltnar.mean('time')
	umid_anom=umid-umid_mn
	salt_anom=saltnar-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	nar_area=dznar*DYTnar
	nar_vol_cont=((s_ref-salt_mn)/s_ref)*umid_anom*nar_area
	nar_salt_cont=(-1.0*salt_anom/s_ref)*umid_mn*nar_area
	nar_liq_flux_vol_cont=km3yr*nar_vol_cont.sum(dim=('lev','j'))
	nar_liq_flux_salt_cont=km3yr*nar_salt_cont.sum(dim=('lev','j'))
	#Split for salinities above and below s_ref
	nar_liq_flux_vol_cont_below_sref=km3yr*nar_vol_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','j'))
	nar_liq_flux_vol_cont_above_sref=km3yr*nar_vol_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','j'))
	nar_liq_flux_salt_cont_below_sref=km3yr*nar_salt_cont.where(saltnar<=s_ref,other=0).sum(dim=('lev','j'))
	nar_liq_flux_salt_cont_above_sref=km3yr*nar_salt_cont.where(saltnar>s_ref,other=0).sum(dim=('lev','j'))

	
	#-------Barrow Strait-------#
	#Liquid fw flux through Barrow Strait: x=[26,27],y=42, full map indices
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are out of the Arctic so multiply v by -1
	#v is north of tracer for a given tracer index
	j=42
	istart=26
	iend=27
	kmax=16
	vbrw=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vbrw=vbrw.where(xr.ufuncs.isfinite(vbrw),other=0).compute()
	saltbrw=salt[:,0:kmax,j,istart:iend+1].load()
	DYVbrw=DYV[j-1:j+1,istart:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	dzbrw=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(saltbrw),coords=saltbrw.coords,dims=saltbrw.dims)
	vmid=-1.0*(vbrw[:,:,0]*DYVbrw[0,:]+vbrw[:,:,1]*DYVbrw[1,:])/DYVbrw.sum(dim='j')

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltbrw.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltbrw-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	brw_area=dzbrw*DXTbrw
	brw_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*brw_area
	brw_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*brw_area
	brw_liq_flux_vol_cont=km3yr*brw_vol_cont.sum(dim=('lev','i'))
	brw_liq_flux_salt_cont=km3yr*brw_salt_cont.sum(dim=('lev','i'))
	#Split for salinities above and below s_ref
	brw_liq_flux_vol_cont_below_sref=km3yr*brw_vol_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','i'))
	brw_liq_flux_vol_cont_above_sref=km3yr*brw_vol_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','i'))
	brw_liq_flux_salt_cont_below_sref=km3yr*brw_salt_cont.where(saltbrw<=s_ref,other=0).sum(dim=('lev','i'))
	brw_liq_flux_salt_cont_above_sref=km3yr*brw_salt_cont.where(saltbrw>s_ref,other=0).sum(dim=('lev','i'))


	#-------Fram Strait-------#
	#Liquid flux through Fram Strait
	#Compute the fw volume relative to 34.8
	#Eastward velocities are into the Arctic so the sign convention is correct
	#Northward velocities are out of the Arctic so multiply v by -1
	#v is north of tracer for a given tracer index
	istart=216
	iend=224
	jstart=35
	jend=43
	kmax=35
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)
	ufram=u[:,0:kmax,jstart:jend+1,istart-1:iend+1] #one point larger in u for midpoint averaging
	ufram=ufram.where(xr.ufuncs.isfinite(ufram),other=0).compute()
	vfram=v[:,0:kmax,jstart-1:jend+1,istart:iend+1] #one point larger in v for midpoint averaging
	vfram=vfram.where(xr.ufuncs.isfinite(vfram),other=0).compute()
	saltfram=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	DYVfram=DYV[jstart-1:jend+1,istart:iend+1].load() #one point larger in v for midpoint averaging
	DXVfram=DXV[jstart:jend+1,istart:iend+1].load()
	DXUfram=DXU[jstart:jend+1,istart-1:iend+1].load() #one point larger in u for midpoint averaging
	DYUfram=DYU[jstart:jend+1,istart:iend+1].load()
	dzfram=dz[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(saltfram),coords=saltfram.coords,dims=saltfram.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltfram),coords=saltfram.coords,dims=saltfram.dims)

	#count backward for y
	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
	    umid[:,:,j,i]=(ufram[:,:,j,i]*DXUfram[j,i]+ufram[:,:,j,i+1]*DXUfram[j,i+1])/DXUfram[j,i:i+2].sum(dim='i')
	    vmid[:,:,j,i]=-1.0*(vfram[:,:,j,i]*DYVfram[j,i]+vfram[:,:,j+1,i]*DYVfram[j+1,i])/DYVfram[j:j+2,i].sum(dim='j')

	#Can only do this if u,v arrays are same size and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Fram line makes a 45˚ angle with u and v
	fram_vel=(umid+vmid)*vel_scale
	fram_area=dzfram*np.sqrt(DYUfram*DYUfram+DXVfram*DXVfram)
	#Compute mean and anomalies for salinity and velocity
	vel_mn=fram_vel.mean('time')
	salt_mn=saltfram.mean('time')
	vel_anom=fram_vel-vel_mn
	salt_anom=saltfram-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	fram_vol_cont=((s_ref-salt_mn)/s_ref)*vel_anom*fram_area
	fram_salt_cont=(-1.0*salt_anom/s_ref)*vel_mn*fram_area
	fram_liq_flux_vol_cont=km3yr*fram_vol_cont.sum(dim=('lev','i','j'))
	fram_liq_flux_salt_cont=km3yr*fram_salt_cont.sum(dim=('lev','i','j'))
	#Split for salinities above and below s_ref
	fram_liq_flux_vol_cont_below_sref=km3yr*fram_vol_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','i','j'))
	fram_liq_flux_vol_cont_above_sref=km3yr*fram_vol_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','i','j'))
	fram_liq_flux_salt_cont_below_sref=km3yr*fram_salt_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','i','j'))
	fram_liq_flux_salt_cont_above_sref=km3yr*fram_salt_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','i','j'))


	#-------Barents Sea Opening-------#
	#Liquid flux through Barents Sea Opening: indices in the full world
	#Compute the fw volume relative to 34.8
	#Eastward into the Arctic so the sign convention is correct 
	#Northward velocities are out of the Arctic so multiply v by -1
	#v is north of tracer for a given tracer index
	#u is east of tracer for a given tracer index
	istart=197
	iend=206
	jstart=48
	jend=57
	kmax=17
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)
	ubso=u[:,0:kmax,jstart:jend+1,istart-1:iend+1] #one point larger in u for midpoint averaging
	ubso=ubso.where(xr.ufuncs.isfinite(ubso),other=0).compute()
	vbso=v[:,0:kmax,jstart-1:jend+1,istart:iend+1] #one point larger in v for midpoint averaging
	vbso=vbso.where(xr.ufuncs.isfinite(vbso),other=0).compute()
	saltbso=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	DYVbso=DYV[jstart-1:jend+1,istart:iend+1].load() #one point larger in v for midpoint averaging
	DXVbso=DXV[jstart:jend+1,istart:iend+1].load()
	DXUbso=DXU[jstart:jend+1,istart-1:iend+1].load() #one point larger in u for midpoint averaging
	DYUbso=DYU[jstart:jend+1,istart:iend+1].load()
	dzbso=dz[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)
	vmid=xr.DataArray(data=np.zeros_like(saltbso),coords=saltbso.coords,dims=saltbso.dims)

	#count backward for y
	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
	    umid[:,:,j,i]=(ubso[:,:,j,i]*DXUbso[j,i]+ubso[:,:,j,i+1]*DXUbso[j,i+1])/DXUbso[j,i:i+2].sum(dim='i')
	    vmid[:,:,j,i]=-1.0*(vbso[:,:,j,i]*DYVbso[j,i]+vbso[:,:,j+1,i]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='j')

	#Can only do this if u,v arrays are same size and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes a 45˚ angle with u and v
	#Can only do this if u,v arrays are same size and values are colocated
	bso_vel=(umid+vmid)*vel_scale
	bso_area=dzbso*np.sqrt(DYUbso*DYUbso+DXVbso*DXVbso)
	#Compute mean and anomalies for salinity and velocity
	vel_mn=bso_vel.mean('time')
	salt_mn=saltbso.mean('time')
	vel_anom=bso_vel-vel_mn
	salt_anom=saltbso-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	bso_vol_cont=((s_ref-salt_mn)/s_ref)*vel_anom*bso_area
	bso_salt_cont=(-1.0*salt_anom/s_ref)*vel_mn*bso_area
	bso_liq_flux_vol_cont=km3yr*bso_vol_cont.sum(dim=('lev','i','j'))
	bso_liq_flux_salt_cont=km3yr*bso_salt_cont.sum(dim=('lev','i','j'))
	#Split for salinities above and below s_ref
	bso_liq_flux_vol_cont_below_sref=km3yr*bso_vol_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','i','j'))
	bso_liq_flux_vol_cont_above_sref=km3yr*bso_vol_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','i','j'))
	bso_liq_flux_salt_cont_below_sref=km3yr*bso_salt_cont.where(saltbso<=s_ref,other=0).sum(dim=('lev','i','j'))
	bso_liq_flux_salt_cont_above_sref=km3yr*bso_salt_cont.where(saltbso>s_ref,other=0).sum(dim=('lev','i','j'))

	#-------Davis Strait-------#
	#Liquid flux through Davis Strait: x=[67,88], y=[28] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. 
	#Flow out of Davis is southward so the sign convention is correct.
	j=28
	istart=67
	iend=88
	kmax=22 #maximum number of depth levels across the transect
	vdav=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(saltdav),coords=saltdav.coords,dims=saltdav.dims)
	vmid=(vdav[:,:,0]*DYVdav[0,:]+vdav[:,:,1]*DYVdav[1,:])/DYVdav.sum(dim='j')

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltdav.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltdav-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	dav_area=dzdav*DXTdav
	dav_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*dav_area
	dav_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*dav_area
	dav_liq_flux_vol_cont=km3yr*dav_vol_cont.sum(dim=('lev','i'))
	dav_liq_flux_salt_cont=km3yr*dav_salt_cont.sum(dim=('lev','i'))
	#Split for salinities above and below s_ref
	dav_liq_flux_vol_cont_below_sref=km3yr*dav_vol_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','i'))
	dav_liq_flux_vol_cont_above_sref=km3yr*dav_vol_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','i'))
	dav_liq_flux_salt_cont_below_sref=km3yr*dav_salt_cont.where(saltdav<=s_ref,other=0).sum(dim=('lev','i'))
	dav_liq_flux_salt_cont_above_sref=km3yr*dav_salt_cont.where(saltdav>s_ref,other=0).sum(dim=('lev','i'))

	print('Liquid flux decomposition done!')


	####-------------Save everything-------------####
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
