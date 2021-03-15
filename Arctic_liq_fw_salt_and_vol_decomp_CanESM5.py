import xarray as xr
import numpy as np
import get_models as gm

####Arctic_liq_fw_salt_and_vol_decomp_CanESM5.py#######
#Decomposes the anomalous salt and volume contributions to the ocean
#fw flux through the Arctic gateways following Holland et al. 2006
#s=s_mean + s', u=u_mean+u' where u_mean,s_mean are the 2000-2100
#mean salinity and velocity in each grid cell and u',s' are temporal
#departures from that mean. The relevant terms in the integrals for the
#fluxes are then s'*u_mean (salinity contribution) and s_mean*u' (velocity contribution)
institution='CCCma'
model='CanESM5'
experiment='ssp585' #ssp126 or ssp585
nmax=10 #Current number of ensemble members for eaach experiment
n_start=1
n_end=nmax #nmax,normally

#Constants
s_ref=34.8 #Reference salinity [PSU]
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

#Cell thicknesses
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p2f1/Ofx/'
xdo=xr.open_dataset(direc_static+'thkcello/gn/v20190429/thkcello/thkcello_Ofx_'+model+'_piControl_r1i1p2f1_gn.nc')
dz=xdo['thkcello'] #Ocean layer thickness [m] (this is 3D--z,y,x)

for ens_num in range(n_start,n_end+1):
	####-------------Read in the output-------------####
	#Loops over ensemble members. Technicaally I can skip the for loop for ice variables
	#by providing a list of the paths to the sea ice file to open_mfdataset instead
	#but easier to just do one member at a time

	print('Getting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p2f1'
	#For getting the historical and ssp files
	simhist,miphist,tp_hist=gm.choose_model_and_sim_cmip6(model,'historical',variant)
	simssp, mipssp, tp_ssp=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean directory for CanESM5 for CMIP6
	direc_hist='/glade/collections/cmip/CMIP6/'+miphist+'/'+institution+'/'+model+'/historical/'+variant+'/Omon/'
	direc_ssp='/glade/collections/cmip/CMIP6/'+mipssp+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'

	#Easiest to do this for historical paths
	salt_paths=[direc_hist+'so/gn/v20190429/so/so_Omon_'+simhist+'199101-200012.nc',
		direc_hist+'so/gn/v20190429/so/so_Omon_'+simhist+'200101-201012.nc',
		direc_hist+'so/gn/v20190429/so/so_Omon_'+simhist+'201101-201412.nc']
	u_paths=[direc_hist+'uo/gn/v20190429/uo/uo_Omon_'+simhist+'199101-200012.nc',
		direc_hist+'uo/gn/v20190429/uo/uo_Omon_'+simhist+'200101-201012.nc',
		direc_hist+'uo/gn/v20190429/uo/uo_Omon_'+simhist+'201101-201412.nc']
	v_paths=[direc_hist+'vo/gn/v20190429/vo/vo_Omon_'+simhist+'199101-200012.nc',
			direc_hist+'vo/gn/v20190429/vo/vo_Omon_'+simhist+'200101-201012.nc',
			direc_hist+'vo/gn/v20190429/vo/vo_Omon_'+simhist+'201101-201412.nc']
	#Ocean variables	
	print('Getting ocean variables')
	#Historical
	#Get salinity and velocities
	xdo=xr.open_mfdataset(salt_paths,combine='by_coords',data_vars='minimal')
	salthist=xdo['so'].sel({'time':slice('2000-01-01',None)}) #ocean salinity in PSU
	xdo=xr.open_mfdataset(u_paths,combine='by_coords',data_vars='minimal')
	uhist=xdo['uo'].sel({'time':slice('2000-01-01',None)}) #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(v_paths,combine='by_coords',data_vars='minimal')
	vhist=xdo['vo'].sel({'time':slice('2000-01-01',None)}) #ocean meridional velocity in [m/s]		
	#ssp
	#Get salinity and velocities
	xdo=xr.open_mfdataset(direc_ssp+'so/gn/v20190429/so/so_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	saltssp=xdo['so'] #ocean salinity in PSU
	xdo=xr.open_mfdataset(direc_ssp+'uo/gn/v20190429/uo/uo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	ussp=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_ssp+'vo/gn/v20190429/vo/vo_Omon_'+simssp+'*.nc',combine='by_coords',data_vars='minimal')
	vssp=xdo['vo'] #ocean meridional velocity in [m/s]

	#Concatenate
	salt=xr.concat([salthist,saltssp],dim='time')
	u=xr.concat([uhist,ussp],dim='time')
	v=xr.concat([vhist,vssp],dim='time')
		
	print('Model output read in and ready for computation')


	####-------------Compute the liquid flux decomposition through the gateways-------------####

	#-------Bering Strait-------#
	#Liquid flux through Bering Strait: x=[113,115], y=[245] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=245
	istart=113
	iend=114
	kmax=3 #maximum number of depth levels across the transect
	vber=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vber=vber.where(xr.ufuncs.isfinite(vber),other=0).compute()
	saltber=salt[:,0:kmax,j,istart:iend+1].load()
	DYVber=DYV[j-1:j+1,istart:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	dzber=dz[0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vber[:,:,1]),coords=vber[:,:,1].coords,dims=vber[:,:,1].dims)
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
	#Liquid flux through Nares: x=248, y=[276,277] full world indices
	#u is east of tracer for a given tracer index
	#Compute the fw volume relative to 34.8
	#Flow into Nares is east so sign convention is correct
	i=248
	jstart=276
	jend=277
	kmax=18  #Maximum number of deptth levels across the transect
	unar=u[:,0:kmax,jstart:jend+1,i-1:i+1]
	unar=unar.where(xr.ufuncs.isfinite(unar),other=0).compute()
	saltnar=salt[:,0:kmax,jstart:jend+1,i].load()
	DXUnar=DXU[jstart:jend+1,i-1:i+1].load()
	DYTnar=DYT[jstart:jend+1,i].load()
	dznar=dz[0:kmax,jstart:jend+1,i].load()
	umid=xr.DataArray(data=np.zeros_like(unar[:,:,:,1]),coords=unar[:,:,:,1].coords,dims=unar[:,:,:,1].dims)
	umid=(unar[:,:,:,0]*DXUnar[:,0]+unar[:,:,:,1]*DXUnar[:,1])/DXUnar.sum(dim='i')

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
	#Liquid fw flux through Barrow Strait: x=[227,229],y=278, full map indices
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	#v is north of tracer for a given tracer index
	j=278
	istart=227
	iend=229
	kmax=16
	vbrw=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vbrw=vbrw.where(xr.ufuncs.isfinite(vbrw),other=0).compute()
	saltbrw=salt[:,0:kmax,j,istart:iend+1].load()
	DYVbrw=DYV[j-1:j+1,istart:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	dzbrw=dz[0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vbrw[:,:,1]),coords=vbrw[:,:,1].coords,dims=vbrw[:,:,1].dims)
	vmid=(vbrw[:,:,0]*DYVbrw[0,:]+vbrw[:,:,1]*DYVbrw[1,:])/DYVbrw.sum(dim='j')

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
	#Liquid flux through Fram Strait, indices: x=[267,278],y=272, full world indices
	#Compute the fw volume relative to 34.8
	#Northward velocities are into the Arctic so the sign convention is correct
	#v is north of tracer for a given tracer index
	j=272
	istart=267
	iend=278
	kmax=35
	vfram=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vfram=vfram.where(xr.ufuncs.isfinite(vfram),other=0).compute()
	saltfram=salt[:,0:kmax,j,istart:iend+1].load()
	DYVfram=DYV[j-1:j+1,istart:iend+1].load()
	DXTfram=DXT[j,istart:iend+1].load()
	dzfram=dz[0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vfram[:,:,1]),coords=vfram[:,:,1].coords,dims=vfram[:,:,1].dims)
	vmid=(vfram[:,:,0]*DYVfram[0,:]+vfram[:,:,1]*DYVfram[1,:])/DYVfram.sum(dim='j')

	#Compute mean and anomalies for salinity and velocity
	vmid_mn=vmid.mean('time')
	salt_mn=saltfram.mean('time')
	vmid_anom=vmid-vmid_mn
	salt_anom=saltfram-salt_mn
	#Compute anomalous velocity (volume) and salinity contributions
	fram_area=dzfram*DXTfram
	fram_vol_cont=((s_ref-salt_mn)/s_ref)*vmid_anom*fram_area
	fram_salt_cont=(-1.0*salt_anom/s_ref)*vmid_mn*fram_area
	fram_liq_flux_vol_cont=km3yr*fram_vol_cont.sum(dim=('lev','i'))
	fram_liq_flux_salt_cont=km3yr*fram_salt_cont.sum(dim=('lev','i'))
	#Split for salinities above and below s_ref
	fram_liq_flux_vol_cont_below_sref=km3yr*fram_vol_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','i'))
	fram_liq_flux_vol_cont_above_sref=km3yr*fram_vol_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','i'))
	fram_liq_flux_salt_cont_below_sref=km3yr*fram_salt_cont.where(saltfram<=s_ref,other=0).sum(dim=('lev','i'))
	fram_liq_flux_salt_cont_above_sref=km3yr*fram_salt_cont.where(saltfram>s_ref,other=0).sum(dim=('lev','i'))


	#-------Barents Sea Opening-------#
	#Liquid flux through Barents Sea Opening: indices in the full world
	#Compute the fw volume relative to 34.8
	#Eastward and northward velocities are into the Arctic so the sign convention is correct for both
	#v is north of tracer for a given tracer index
	#u is east of tracer for a given tracer index
	istart=280
	iend=291
	jstart=258
	jend=269
	kmax=19
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)
	ubso=u[:,0:kmax,jstart:jend+1,istart:iend+1]
	ubso=ubso.where(xr.ufuncs.isfinite(ubso),other=0).compute()
	vbso=v[:,0:kmax,jstart:jend+1,istart:iend+1]
	vbso=vbso.where(xr.ufuncs.isfinite(vbso),other=0).compute()
	saltbso=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	DYVbso=DYV[jstart:jend+1,istart:iend+1].load()
	DXVbso=DXV[jstart:jend+1,istart:iend+1].load()
	DXUbso=DXU[jstart:jend+1,istart:iend+1].load()
	DYUbso=DYU[jstart:jend+1,istart:iend+1].load()
	dzbso=dz[0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(ubso),coords=ubso.coords,dims=ubso.dims)
	vmid=xr.DataArray(data=np.zeros_like(vbso),coords=vbso.coords,dims=vbso.dims)

	#count backward for y
	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
	    umid[:,:,j,i]=(ubso[:,:,j,i-1]*DXUbso[j,i-1]+ubso[:,:,j,i]*DXUbso[j,i])/DXUbso[j,i-1:i+1].sum(dim='i')
	    vmid[:,:,j,i]=(vbso[:,:,j-1,i]*DYVbso[j-1,i]+vbso[:,:,j,i]*DYVbso[j,i])/DYVbso[j-1:j+1,i].sum(dim='j')

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
	#Liquid flux through Davis Strait: x=[232,239], y=[249] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Davis is southward so the sign convention is correct.
	j=249
	istart=232
	iend=239
	kmax=21 #maximum number of depth levels across the transect
	vdav=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vdav=vdav.where(xr.ufuncs.isfinite(vdav),other=0).compute()
	saltdav=salt[:,0:kmax,j,istart:iend+1].load()
	DYVdav=DYV[j-1:j+1,istart:iend+1].load()
	DXTdav=DXT[j,istart:iend+1].load()
	dzdav=dz[0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vdav[:,:,1]),coords=vdav[:,:,1].coords,dims=vdav[:,:,1].dims)
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
