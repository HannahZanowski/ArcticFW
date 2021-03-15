import xarray as xr
import numpy as np
import os
import get_models as gm

##############Arctic_fw_gateways_CNRM-CM6-1.py#################
#Compute liquid and solid fw fluxes through the Arctic gateways
#for the ensemble members of CNRM-CM6-1
institution='CNRM-CERFACS'
model='CNRM-CM6-1'
experiment='ssp585' #any of the ssps as well
#Get the right number for the end value for the for loop because the number of ensemble members
#differs across the models and simulations
ens_num_dict={'historical':np.arange(1,20),'ssp126':np.arange(1,7),
                'ssp585':np.arange(1,6)}
nmax=int(len(ens_num_dict[experiment]))
ens_nums=ens_num_dict[experiment]

#Constants
s_ref=34.8 #Reference salinity [PSU]
s_si=4 #Sea ice salinity [PSU]. #NO CLUE IF THIS IS RIGHT. NOT IN ESDOCS
s_sn=0 #Snow salinity [PSU], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3] 
rho_si=917.0 #Sea ice density [kg/m3] #NO CLUE IF THIS IS RIGHT. NOT IN ESDOCS
rho_sn=330.0 #Snow denisty [kg/m3] #NO CLUE IF THIS IS RIGHT. NOT IN ESDOCS
rho_si_over_rho_fw=rho_si/rho_fw #Scale factor for converting ice volume to liquid freshwater volume [unitless]
rho_sn_over_rho_fw=rho_sn/rho_fw #Scale factor for converting snow volume to liquid freshwater volume [unitless]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Ocean and sea ice grid variables
#Grid distance file
direc_f2='/glade/work/zanowski/ArcticFW/'
xdf=xr.open_dataset(direc_f2+'CNRM_CM6_grid_distances.nc')
#get x and y distances: tracer x and y, y distances at v points and x distances at u points
DYV=xdf['dyv'] #y-spacing, v points
DXU=xdf['dxu'] #x-spacing, u points
DXT=xdf['dxt'] #x-spacing, tracer pts
DYT=xdf['dyt'] #y-spacing, tracer pts
DXV=xdf['dxv'] #x-spacing, v points
DYU=xdf['dyu'] #y-spacing, u points
#Arctic mask, cell areas and ocean depth for storage code
#Read in areacello
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f2/Ofx/'
xdo=xr.open_dataset(direc_static+'areacello/gn/v20180814/areacello/areacello_Ofx_'+model+'_piControl_r1i1p1f2_gn.nc')
#This [:-1,1:-1] gets rid of masked edge points
oarea=xdo['areacello'][:-1,1:-1] #Ocean tracer grid cell area [m^2]
xdo=xr.open_dataset(direc_static+'deptho/gn/v20180814/deptho/deptho_Ofx_'+model+'_piControl_r1i1p1f2_gn.nc')
odepth=xdo['deptho'][:,:-1,1:-1] #Ocean depth (tracer) [m]
#This is a 3D field so get the max level
odepth=odepth.max(dim='lev')
xda=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xda['arctic_mask']
oarea=oarea.where(amask==1,drop=True).load()
odepth=odepth.where(amask==1,drop=True).load()
arctic_vol=m3tokm3*(oarea*odepth).sum(dim=('x','y'))

#Loops over ensemble members. Technically I can skip the for loop for ice variables
#by providing a list of the paths to the sea ice file to open_mfdataset instead
#but I am optting not to do this

#NOTE:version numbers change with each ensemble member,
#so to account for this the code below sorts
#the contents of os.listdir() at the version level
#and chooses the last element, which should be the latest version date.

for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	variant='r'+str('%i' %ens_num)+'i1p1f2'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)
	#Ice directories for CNRM-CM6-1 for CMIP6
	direc_i='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/SImon/'
	sit_path=direc_i+'sithick/gn/'+sorted(os.listdir(direc_i+'sithick/gn/'))[-1]+'/sithick/sithick_SImon_'+simulation+time_period+'.nc'
	snt_path=direc_i+'sisnthick/gn/'+sorted(os.listdir(direc_i+'sisnthick/gn/'))[-1]+'/sisnthick/sisnthick_SImon_'+simulation+time_period+'.nc'
	siu_path=direc_i+'siu/gn/'+sorted(os.listdir(direc_i+'siu/gn/'))[-1]+'/siu/siu_SImon_'+simulation+time_period+'.nc'
	siv_path=direc_i+'siv/gn/'+sorted(os.listdir(direc_i+'siv/gn/'))[-1]+'/siv/siv_SImon_'+simulation+time_period+'.nc'
	sic_path=direc_i+'siconc/gn/'+sorted(os.listdir(direc_i+'siconc/gn/'))[-1]+'/siconc/siconc_SImon_'+simulation+time_period+'.nc'
	#Ocean directories for CNRM-CM6-1 for CMIP6	
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	salt_path=direc_o+'so/gn/'+sorted(os.listdir(direc_o+'so/gn/'))[-1]+'/so/so_Omon_'+simulation+'*.nc'
	u_path=direc_o+'uo/gn/'+sorted(os.listdir(direc_o+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simulation+'*.nc'
	v_path=direc_o+'vo/gn/'+sorted(os.listdir(direc_o+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simulation+'*.nc'
	dz_path=direc_o+'thkcello/gn/'+sorted(os.listdir(direc_o+'thkcello/gn/'))[-1]+'/thkcello/thkcello_Omon_'+simulation+'*.nc'

	#This model leaves extra edge points for the cyclic boundaries and the northern boundary 
	#so get rid of these by slicing like this: [:,:,:-1,1:-1]
	#Ocean variables
	print('Getting ocean variables')	
	#Get salinity, velocities, and cell thicknesses	
	xdo=xr.open_mfdataset(salt_path,combine='by_coords',data_vars='minimal')
	salt=xdo['so'][:,:,:-1,1:-1] #ocean salinity in PSU
	xdo=xr.open_mfdataset(u_path,combine='by_coords',data_vars='minimal')
	u=xdo['uo'][:,:,:-1,1:-1] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(v_path,combine='by_coords',data_vars='minimal')
	v=xdo['vo'][:,:,:-1,1:-1] #ocean meridional velocity in [m/s]
	xdo=xr.open_mfdataset(dz_path,combine='by_coords',data_vars='minimal')
	dz=xdo['thkcello'][:,:,:-1,1:-1] #ocean cell thicknesses in [m]
	
	#Ice variables
	print('Getting sea ice variables')
	#Get sea ice thickness, snow thickness, and sea ice velocities
	xdi=xr.open_dataset(sit_path,chunks={})
	sit=xdi['sithick'][:,:-1,1:-1] #sea ice thickness [m]
	xdi=xr.open_dataset(sic_path,chunks={})
	sic=xdi['siconc'][:,:-1,1:-1]*0.01 #sea ice concentration [%], converted to fraction
	xdi=xr.open_dataset(snt_path,chunks={})
	snt=xdi['sisnthick'][:,:-1,1:-1] #sea ice snow thickness [m]
	xdi=xr.open_dataset(siu_path,chunks={})
	siu=xdi['siu'][:,:-1,1:-1] #sea ice zonal velocity [m/s]
	xdi=xr.open_dataset(siv_path,chunks={})
	siv=xdi['siv'][:,:-1,1:-1] #sea ice meridional velocity [m/s]

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs in the for loop. Will have to do this here for
	#ice variables anyway because they are in one file
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
	liq_fw_storage=m3tokm3*liq_fw.sum(dim=('lev','y','x'))
	#keep track of the water with salinities above/below sref too
	liq_fw_storage_below_sref=m3tokm3*liq_fw.where(salt_arctic<=s_ref,other=0).sum(dim=('lev','y','x'))
	liq_fw_storage_above_sref=m3tokm3*liq_fw.where(salt_arctic>s_ref,other=0).sum(dim=('lev','y','x'))

	#Solid fw storage--no need for snow concentration correction
	solid_storage=(sit_arctic+snt_arctic)*sic_arctic*oarea #s_sn is zero so the salinity scale factor for snow is 1
	solid_fw_storage=m3tokm3*solid_storage.sum(dim=('x','y'))

	print('Storage done!')

	####-------------Compute the solid fluxes through the gateways-------------####

	#-------Bering Strait-------#
	# x=[112,114], y=[247] Indices on full world map
	# v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=247
	istart=112
	iend=114
	sivber=siv[:,j-1:j+1,istart:iend+1]
	sivber=sivber.where(xr.ufuncs.isfinite(sivber),other=0).compute()
	sicber=sic[:,j,istart:iend+1]
	sicber=sicber.where(xr.ufuncs.isfinite(sicber),other=0).compute()
	sitber=sit[:,j,istart:iend+1]
	sitber=sitber.where(xr.ufuncs.isfinite(sitber),other=0).compute()
	sntber=snt[:,j,istart:iend+1]
	sntber=sntber.where(xr.ufuncs.isfinite(sntber),other=0).compute()
	DYVber=DYV[j-1:j+1,istart:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sivber[:,1]),coords=sivber[:,1].coords,dims=sivber[:,1].dims)
	vmid=(sivber[:,0]*DYVber[0,:]+sivber[:,1]*DYVber[1,:])/DYVber.sum(dim='y')
	#Compute the fluxes
	ber_vel=vmid*DXTber
	ber_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitber
	#snow is pure fw (s_sn=0) so the salinity scale factor is 1
	ber_snow=rho_sn_over_rho_fw*sntber
	ber_solid=(ber_ice+ber_snow)*sicber*ber_vel
	ber_solid_flux=km3yr*ber_solid.sum(dim=('x'))

	#-------Nares Strait-------#
	#x=249, y=[278,279] full world indices
	#u is east of tracer for a given tracaer index
	#Compute the fw volume relative to 34.8
	#Flow into Nares is east so sign convention is correct
	i=249
	jstart=278
	jend=279
	siunar=siu[:,jstart:jend+1,i-1:i+1]
	siunar=siunar.where(xr.ufuncs.isfinite(siunar),other=0).compute()
	sicnar=sic[:,jstart:jend+1,i]
	sicnar=sicnar.where(xr.ufuncs.isfinite(sicnar),other=0).compute()
	sitnar=sit[:,jstart:jend+1,i]
	sitnar=sitnar.where(xr.ufuncs.isfinite(sitnar),other=0).compute()
	sntnar=snt[:,jstart:jend+1,i]
	sntnar=sntnar.where(xr.ufuncs.isfinite(sntnar),other=0).compute()
	DXUnar=DXU[jstart:jend+1,i-1:i+1].load()
	DYTnar=DYT[jstart:jend+1,i].load()
	umid=xr.DataArray(data=np.zeros_like(siunar[:,:,1]),coords=siunar[:,:,1].coords,dims=siunar[:,:,1].dims)
	umid=(siunar[:,:,0]*DXUnar[:,0]+siunar[:,:,1]*DXUnar[:,1])/DXUnar.sum(dim='x')
	#Compute the fluxes
	nar_vel=umid*DYTnar
	nares_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitnar
	#Snow is pure fw so the salinity scale factor is 1
	nares_snow=rho_sn_over_rho_fw*sntnar
	nares_solid=(nares_ice+nares_snow)*sicnar*nar_vel
	nares_solid_flux=km3yr*nares_solid.sum(dim=('y')) #sum over latitude

	#-------Barrow Strait-------#
	#x=[227,229],y=281, full map indices
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign conventtion is correct
	#v is north of tracer for a given tracer index
	j=281
	istart=227
	iend=229
	sivbrw=siv[:,j-1:j+1,istart:iend+1]
	sivbrw=sivbrw.where(xr.ufuncs.isfinite(sivbrw),other=0).compute()
	sicbrw=sic[:,j,istart:iend+1]
	sicbrw=sicbrw.where(xr.ufuncs.isfinite(sicbrw),other=0).compute()
	sitbrw=sit[:,j,istart:iend+1]
	sitbrw=sitbrw.where(xr.ufuncs.isfinite(sitbrw),other=0).compute()
	sntbrw=snt[:,j,istart:iend+1]
	sntbrw=sntbrw.where(xr.ufuncs.isfinite(sntbrw),other=0).compute()
	DYVbrw=DYV[j-1:j+1,istart:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sivbrw[:,1]),coords=sivbrw[:,1].coords,dims=sivbrw[:,1].dims)
	#This is only 3 points wide
	vmid=(sivbrw[:,0]*DYVbrw[0,:]+sivbrw[:,1]*DYVbrw[1,:])/DYVbrw.sum(dim='y')
	#Compute the fluxes
	brw_vel=vmid*DXTbrw
	brw_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbrw
	#Snow is pure fw so the salinity scale factor is 1
	brw_snow=rho_sn_over_rho_fw*sntbrw
	brw_solid=(brw_ice+brw_snow)*sicbrw*brw_vel
	brw_solid_flux=km3yr*brw_solid.sum(dim=('x')) #sum over longitude

	#-------Fram Strait-------#
	#Solid flux through Fram Strait, indices: x=[267,276],y=272, full world indices
	#Compute the fw volume relative to 34.8
	#Northward velocities are into the Arctic so the sign convention is correct
	#v is north of tracer for a given tracer index
	j=273
	istart=267
	iend=276
	sivfram=siv[:,j-1:j+1,istart:iend+1]
	sivfram=sivfram.where(xr.ufuncs.isfinite(sivfram),other=0).compute()
	sicfram=sic[:,j,istart:iend+1]
	sicfram=sicfram.where(xr.ufuncs.isfinite(sicfram),other=0).compute()
	sitfram=sit[:,j,istart:iend+1]
	sitfram=sitfram.where(xr.ufuncs.isfinite(sitfram),other=0).compute()
	sntfram=snt[:,j,istart:iend+1]
	sntfram=sntfram.where(xr.ufuncs.isfinite(sntfram),other=0).compute()
	DYVfram=DYV[j-1:j+1,istart:iend+1].load()
	DXTfram=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sivfram[:,1]),coords=sivfram[:,1].coords,dims=sivfram[:,1].dims)
	vmid=(sivfram[:,0]*DYVfram[0,:]+sivfram[:,1]*DYVfram[1,:])/DYVfram.sum(dim='y')
	#Compute the fluxes
	fram_vel=vmid*DXTfram
	fram_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitfram
	#Snow is pure fw so the salinity scale factor is 1
	fram_snow=rho_sn_over_rho_fw*sntfram
	fram_solid=(fram_ice+fram_snow)*sicfram*fram_vel
	fram_solid_flux=km3yr*fram_solid.sum(dim=('x'))

	#-------Barents Sea Opening-------#
	#Solid flux through Barents Sea Opening: indices in the full world
	#Compute the fw volume relative to 34.8
	#Eastward and northward velocities are into the Arctic so the sign convention is correct for both
	#v is north of tracer points for a given tracer index
	#u is east of tracer points for a given tracer index
	istart=281
	iend=291
	jstart=261
	jend=271
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)

	siubso=siu[:,jstart:jend+1,istart-1:iend+1] #one point larger in u for midpoint averaging
	siubso=siubso.where(xr.ufuncs.isfinite(siubso),other=0).compute()
	sivbso=siv[:,jstart-1:jend+1,istart:iend+1] #one point larger in v for midpoint averaging
	sivbso=sivbso.where(xr.ufuncs.isfinite(sivbso),other=0).compute()
	sicbso=sic[:,jstart:jend+1,istart:iend+1]
	sicbso=sicbso.where(xr.ufuncs.isfinite(sicbso),other=0).compute()
	sitbso=sit[:,jstart:jend+1,istart:iend+1]
	sitbso=sitbso.where(xr.ufuncs.isfinite(sitbso),other=0).compute()
	sntbso=snt[:,jstart:jend+1,istart:iend+1]
	sntbso=sntbso.where(xr.ufuncs.isfinite(sntbso),other=0).compute()
	DYVbso=DYV[jstart-1:jend+1,istart:iend+1].load() #one point larger in v for midpoint averaging
	DXVbso=DXV[jstart:jend+1,istart:iend+1].load()
	DXUbso=DXU[jstart:jend+1,istart-1:iend+1].load() #one point larger in u for midpoint averaging
	DYUbso=DYU[jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(siubso[:,:,1:]),coords=siubso[:,:,1:].coords,dims=siubso[:,:,1:].dims)
	vmid=xr.DataArray(data=np.zeros_like(sivbso[:,1:,:]),coords=sivbso[:,1:,:].coords,dims=sivbso[:,1:,:].dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
	    umid[:,j,i]=(siubso[:,j,i]*DXUbso[j,i]+siubso[:,j,i+1]*DXUbso[j,i+1])/DXUbso[j,i:i+2].sum(dim='x')
	    vmid[:,j,i]=(sivbso[:,j,i]*DYVbso[j,i]+sivbso[:,j+1,i]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')
	#Compute the fluxes
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes 45˚ angle with u and v
	#Can only do this if u,v are same size and values are colocated
	bso_vel=(umid+vmid)*vel_scale*np.sqrt(DXVbso*DXVbso+DYUbso*DYUbso)
	bso_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbso
	#snow is pure fw so salinity factor is 1
	bso_snow=rho_sn_over_rho_fw*sntbso
	bso_solid=(bso_ice+bso_snow)*sicbso*bso_vel
	bso_solid_flux=km3yr*bso_solid.sum(dim=('x','y')) #sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...

	print('Solid Fluxes Done!')


	####-------------Compute the liquid fluxes through the gateways-------------####

	#-------Bering Strait-------#
	#Liquid flux through Bering Strait: x=[112,114], y=[247] Indices on full world map
	#v is north of tracer for a given tracer index
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct.
	j=247
	istart=112
	iend=114
	kmax=18 #maximum number of depth levels across the transect
	vber=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vber=vber.where(xr.ufuncs.isfinite(vber),other=0).compute()
	saltber=salt[:,0:kmax,j,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again 
	#although it takes basically no time
	DYVber=DYV[j-1:j+1,istart:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	dzber=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vber[:,:,1]),coords=vber[:,:,1].coords,dims=vber[:,:,1].dims)
	vmid=(vber[:,:,0]*DYVber[0,:]+vber[:,:,1]*DYVber[1,:])/DYVber.sum(dim='y')
	#Compute the fluxes
	#Volume flux
	ber_vol=vmid*dzber*DXTber
	ber_vol_flux=km3yr*ber_vol.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#FW
	ber_liq=((s_ref-saltber)/s_ref)*ber_vol
	ber_liq_flux=km3yr*ber_liq.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	ber_liq_flux_below_sref=km3yr*ber_liq.where(saltber<=s_ref,other=0).sum(dim=('lev','x'))
	ber_liq_flux_above_sref=km3yr*ber_liq.where(saltber>s_ref,other=0).sum(dim=('lev','x'))	


	#-------Nares Strait-------#
	#Liquid flux through Nares: x=249, y=[278,279] full world indices
	#u is east of tracer for a given tracer index
	#Compute the fw volume relative to 34.8
	#Flow into Nares is east so sign convention is correct
	i=249
	jstart=278
	jend=279
	kmax=37  #Maximum number of deptth levels across the transect
	unar=u[:,0:kmax,jstart:jend+1,i-1:i+1]
	unar=unar.where(xr.ufuncs.isfinite(unar),other=0).compute()
	saltnar=salt[:,0:kmax,jstart:jend+1,i].load()
	DXUnar=DXU[jstart:jend+1,i-1:i+1].load()
	DYTnar=DYT[jstart:jend+1,i].load()
	dznar=dz[:,0:kmax,jstart:jend+1,i].load()
	umid=xr.DataArray(data=np.zeros_like(unar[:,:,:,1]),coords=unar[:,:,:,1].coords,dims=unar[:,:,:,1].dims)
	umid=(unar[:,:,:,0]*DXUnar[:,0]+unar[:,:,:,1]*DXUnar[:,1])/DXUnar.sum(dim='x')
	#Compute the fluxes
	#Volume flux
	nares_vol=umid*dznar*DYTnar
	nares_vol_flux=km3yr*nares_vol.sum(dim=('lev','y')) #only one lon pt so no sum over nlon
	#FW
	nares_liq=((s_ref-saltnar)/s_ref)*nares_vol
	nares_liq_flux=km3yr*nares_liq.sum(dim=('lev','y')) #only one lon pt so no sum over nlon
	#Fluxes above and below sref
	nares_liq_flux_below_sref=km3yr*nares_liq.where(saltnar<=s_ref,other=0).sum(dim=('lev','y'))
	nares_liq_flux_above_sref=km3yr*nares_liq.where(saltnar>s_ref,other=0).sum(dim=('lev','y'))	


	#-------Barrow Strait-------#
	#Liquid fw flux through Barrow Strait: x=[227,229],y=281, full map indices
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	#v is north of tracer for a given tracer index
	j=281
	istart=227
	iend=229
	kmax=31
	vbrw=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vbrw=vbrw.where(xr.ufuncs.isfinite(vbrw),other=0).compute()
	saltbrw=salt[:,0:kmax,j,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again
	#although it takes basically no time
	DYVbrw=DYV[j-1:j+1,istart:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	dzbrw=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vbrw[:,:,1]),coords=vbrw[:,:,1].coords,dims=vbrw[:,:,1].dims)
	vmid=(vbrw[:,:,0]*DYVbrw[0,:]+vbrw[:,:,1]*DYVbrw[1,:])/DYVbrw.sum(dim='y')
	#Compute the fluxes
	#Volume flux
	brw_vol=vmid*dzbrw*DXTbrw
	brw_vol_flux=km3yr*brw_vol.sum(dim=('lev','x')) #only one lat pt so no sum over nlat
	#FW
	brw_liq=((s_ref-saltbrw)/s_ref)*brw_vol
	brw_liq_flux=km3yr*brw_liq.sum(dim=('lev','x')) #only one lat pt so no sum over nlat
	#Fluxes above and below sref
	brw_liq_flux_below_sref=km3yr*brw_liq.where(saltbrw<=s_ref,other=0).sum(dim=('lev','x'))
	brw_liq_flux_above_sref=km3yr*brw_liq.where(saltbrw>s_ref,other=0).sum(dim=('lev','x'))	


	#-------Fram Strait-------#
	#Liquid flux through Fram Strait, indices: x=[267,276],y=273, full world indices
	#Compute the fw volume relative to 34.8
	#Northward velocities are into the Arctic so the sign convention is correct
	#v is north of tracer for a given tracer index
	j=273
	istart=267
	iend=276
	kmax=62
	vfram=v[:,0:kmax,j-1:j+1,istart:iend+1]
	vfram=vfram.where(xr.ufuncs.isfinite(vfram),other=0).compute()
	saltfram=salt[:,0:kmax,j,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again
	#although it takes basically no time
	DYVfram=DYV[j-1:j+1,istart:iend+1].load()
	DXTfram=DXT[j,istart:iend+1].load()
	dzfram=dz[:,0:kmax,j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(vfram[:,:,1]),coords=vfram[:,:,1].coords,dims=vfram[:,:,1].dims)
	vmid=(vfram[:,:,0]*DYVfram[0,:]+vfram[:,:,1]*DYVfram[1,:])/DYVfram.sum(dim='y')
	#Compute the fluxes
	#Volume flux
	fram_vol=vmid*dzfram*DXTfram
	fram_vol_flux=km3yr*fram_vol.sum(dim=('lev','x'))
	#FW
	fram_liq=((s_ref-saltfram)/s_ref)*fram_vol
	fram_liq_flux=km3yr*fram_liq.sum(dim=('lev','x'))
	#Fluxes above and below sref
	fram_liq_flux_below_sref=km3yr*fram_liq.where(saltfram<=s_ref,other=0).sum(dim=('lev','x'))
	fram_liq_flux_above_sref=km3yr*fram_liq.where(saltfram>s_ref,other=0).sum(dim=('lev','x'))	


	#-------Barents Sea Opening-------#
	#Liquid flux through Barents Sea Opening: indices in the full world
	#Compute the fw volume relative to 34.8
	#Eastward and northward velocities are into the Arctic so the sign convention is correct for both
	#v is north of tracer for a given tracer index
	#u is east of tracer for a given tracer index
	istart=281
	iend=291
	jstart=261
	jend=271
	kmax=39
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)
	ubso=u[:,0:kmax,jstart:jend+1,istart-1:iend+1] #one point larger in u for midpoint averaging
	ubso=ubso.where(xr.ufuncs.isfinite(ubso),other=0).compute()
	vbso=v[:,0:kmax,jstart-1:jend+1,istart:iend+1] #one point larger in v for midpoint averaging
	vbso=vbso.where(xr.ufuncs.isfinite(vbso),other=0).compute()
	saltbso=salt[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	#These are computed in the sea ice code so technically no need to read in again
	#although it takes basically no time
	DYVbso=DYV[jstart-1:jend+1,istart:iend+1].load() #one point larger in v for midpoint averaging
	DXVbso=DXV[jstart:jend+1,istart:iend+1].load()
	DXUbso=DXU[jstart:jend+1,istart-1:iend+1].load() #one point larger in u for midpoint averaging
	DYUbso=DYU[jstart:jend+1,istart:iend+1].load()
	dzbso=dz[:,0:kmax,jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(ubso[:,:,:,1:]),coords=ubso[:,:,:,1:].coords,dims=ubso[:,:,:,1:].dims)
	vmid=xr.DataArray(data=np.zeros_like(vbso[:,:,1:,:]),coords=vbso[:,:,1:,:].coords,dims=vbso[:,:,1:,:].dims)

	#count backward for y
	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
	    umid[:,:,j,i]=(ubso[:,:,j,i]*DXUbso[j,i]+ubso[:,:,j,i+1]*DXUbso[j,i+1])/DXUbso[j,i:i+2].sum(dim='x')
	    vmid[:,:,j,i]=(vbso[:,:,j,i]*DYVbso[j,i]+vbso[:,:,j+1,i]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')

	#Compute the fluxes
	#Can only do this if u,v are same size and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes a 45˚ angle with u and v
	#Volume flux
	bso_vol=dzbso*(umid+vmid)*vel_scale*np.sqrt(DXVbso*DXVbso+DYUbso*DYUbso)
	bso_vol_flux=km3yr*bso_vol.sum(dim=('lev','y','x'))#sum over latitude and lon.
	#FW flux
	bso_liq=((s_ref-saltbso)/s_ref)*bso_vol
	bso_liq_flux=km3yr*bso_liq.sum(dim=('lev','y','x'))#sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
	#Fluxes above and below sref
	bso_liq_flux_below_sref=km3yr*bso_liq.where(saltbso<=s_ref,other=0).sum(dim=('lev','y','x'))
	bso_liq_flux_above_sref=km3yr*bso_liq.where(saltbso>s_ref,other=0).sum(dim=('lev','y','x'))	

	print('Liquid fluxes done!')

	####-------------Compute the total fluxes and save everything-------------####
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
print('Ensemble members '+str('%02i' %n_start)+' to '+str('%02i' %n_end)+' done!')
