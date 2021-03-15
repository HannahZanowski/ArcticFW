import xarray as xr
import numpy as np
import os
import get_models as gm

#####Arctic_fw_gateways_MIROC6.py#######
# Compute liquid and solid fw fluxes through the Arctic gateways
#for the ensemble members of MIROC6
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
s_si=5 #Sea ice salinity [ppt]. From ESDOCS
s_sn=0 #Snow salinity [ppt], zero because it's pure freshwater
rho_fw=1000.0 #Density of freshwster [kg/m3]
rho_si=900.0 #Sea ice density [kg/m3]. From Yoshiki Komuro
rho_sn=330.0 #Snow denisty [kg/m3]. From Yoshiki Komuro
rho_si_over_rho_fw=rho_si/rho_fw #Scale factor for converting ice volume to liquid freshwater volume [unitless]
rho_sn_over_rho_fw=rho_sn/rho_fw #Scale factor for converting snow volume to liquid freshwater volume [unitless]
rho_fw_inv=0.001 #1/rho_fw [m3/kg] for converting freshwater mass to freshwater volume
km3yr=(1e-9)*(365*86400) #Convert from [m3/s] to [km3/yr]
m3sv=1e-6 #Convert from [m3/s] to [Sv] [1 Sv = 10^6 m3/s]
m3tokm3=1e-9 #Convert from [m3] to [km3]

#Ocean and sea ice grid variables
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

#Read in areacello, deptho, Arctic mask for storage
direc_static='/glade/collections/cmip/CMIP6/CMIP/'+institution+'/'+model+'/piControl/r1i1p1f1/Ofx/'
xdo=xr.open_dataset(direc_static+'areacello/gn/v20190311/areacello/areacello_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
oarea=xdo['areacello'] #Ocean tracer grid cell area [m^2]
xdo=xr.open_dataset(direc_static+'deptho/gn/v20190311/deptho/deptho_Ofx_'+model+'_piControl_r1i1p1f1_gn.nc')
odepth=xdo['deptho'] #Ocean depth (tracer) [m]
xda=xr.open_dataset(direc_f2+model+'_Arctic_Mask.nc')
amask=xda['arctic_mask']
oarea=oarea.where(amask==1,drop=True).load()
odepth=odepth.where(amask==1,drop=True).load()
arctic_vol=m3tokm3*(oarea*odepth).sum(dim=('y','x'))

#Loop over ensemble members for everything, including data read in
#NOTE:version numbers change with each ensemble member,
#so to account for this the code below sorts
#the contents of os.listdir() at the version level
#and chooses the last element, which should be the latest version date.
for ens_num in ens_nums:
	####-------------Read in the output-------------####
	print('Starting Ensemble Member '+str('%i' %ens_num)+'/'+str('%i' %nmax))
	print('Getting ocean and sea ice variables')
	variant='r'+str('%i' %ens_num)+'i1p1f1'
	simulation, mip, time_period=gm.choose_model_and_sim_cmip6(model,experiment,variant)

	#Ocean variables
	#Ocean directory for UKESM for CMIP6
	direc_o='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/Omon/'
	#Ocean variables
	#Get salinity and velocities as well as the variables for computing layer thickness
	xdo=xr.open_mfdataset(direc_o+'so/gn/'+sorted(os.listdir(direc_o+'so/gn/'))[-1]+'/so/so_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	salt=xdo['so'] #ocean salinity in ppt
	eta=xdo['eta']
	sigma_bnds=xdo['sigma_bnds']
	zlev_bnds=xdo['zlev_bnds']
	depth_c=xdo['depth_c']
	depth=xdo['depth']        
	nsigma=xdo['nsigma']
	xdo=xr.open_mfdataset(direc_o+'uo/gn/'+sorted(os.listdir(direc_o+'uo/gn/'))[-1]+'/uo/uo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	u=xdo['uo'] #ocean zonal velocity in [m/s]
	xdo=xr.open_mfdataset(direc_o+'vo/gn/'+sorted(os.listdir(direc_o+'vo/gn/'))[-1]+'/vo/vo_Omon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	v=xdo['vo'] #ocean meridional velocity in [m/s]

	#Ice variables
	#Ice directory for MIROC6 for CMIP6
	direc_i='/glade/collections/cmip/CMIP6/'+mip+'/'+institution+'/'+model+'/'+experiment+'/'+variant+'/SImon/'
	#Get sea ice concentration, thickness, snow thickness, and sea ice velocities
	xdi=xr.open_mfdataset(direc_i+'siconc/gn/'+sorted(os.listdir(direc_i+'siconc/gn/'))[-1]+'/siconc/siconc_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sic=xdi['siconc']*0.01 #sea ice concentration, converted from % to fraction
	xdi=xr.open_mfdataset(direc_i+'sithick/gn/'+sorted(os.listdir(direc_i+'sithick/gn/'))[-1]+'/sithick/sithick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	sit=xdi['sithick'] #sea ice thickness [m]
	xdi=xr.open_mfdataset(direc_i+'sisnthick/gn/'+sorted(os.listdir(direc_i+'sisnthick/gn/'))[-1]+'/sisnthick/sisnthick_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	snt=xdi['sisnthick'] #sea ice snow thickness [m]; already 'weighted' by sisnconc
	xdi=xr.open_mfdataset(direc_i+'siu/gn/'+sorted(os.listdir(direc_i+'siu/gn/'))[-1]+'/siu/siu_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	siu=xdi['siu'] #sea ice zonal velocity [m/s]
	xdi=xr.open_mfdataset(direc_i+'siv/gn/'+sorted(os.listdir(direc_i+'siv/gn/'))[-1]+'/siv/siv_SImon_'+simulation+'*.nc',combine='by_coords',data_vars='minimal')
	siv=xdi['siv'] #sea ice meridional velocity [m/s]

	#Easier to just cut the simulation lengths here rather than doing
	#it in the globs.
	if experiment=='historical':
		salt=salt.sel({'time':slice('1950-01-01',None)})
		eta=eta.sel({'time':slice('1950-01-01',None)})
		u=u.sel({'time':slice('1950-01-01',None)})
		v=v.sel({'time':slice('1950-01-01',None)})
		sic=sic.sel({'time':slice('1950-01-01',None)})
		sit=sit.sel({'time':slice('1950-01-01',None)})
		snt=snt.sel({'time':slice('1950-01-01',None)})
		siu=siu.sel({'time':slice('1950-01-01',None)})
		siv=siv.sel({'time':slice('1950-01-01',None)})
	

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

	#dz for Arctic storage calculations
	dz_arctic=dz.where(amask==1,drop=True).load()

	print('Model output read in and ready for computation')


	#########1. Arctic fw storage [km3]
	salt_arctic=salt.where(amask==1,drop=True).load()
	sit_arctic=rho_si_over_rho_fw*(1-(s_si/s_ref))*sit.where(amask==1,drop=True).load()
	sic_arctic=sic.where(amask==1,drop=True).load()
	snt_arctic=rho_sn_over_rho_fw*snt.where(amask==1,drop=True).load()

	#Compute fw storage (total, anywhere s is >, <= s_ref
	#ensures positive volumes only)
	liq_fw=((s_ref-salt_arctic)/s_ref)*dz_arctic*oarea
	liq_fw_storage=m3tokm3*liq_fw.sum(dim=('lev','y','x'))                    
	#keep track of the water with salinities above/below sref 
	liq_fw_storage_below_sref=m3tokm3*liq_fw.where(salt_arctic<=s_ref,other=0).sum(dim=('lev','y','x'))
	liq_fw_storage_above_sref=m3tokm3*liq_fw.where(salt_arctic>s_ref,other=0).sum(dim=('lev','y','x'))

	#Solid fw storage--no need for snow concentration correction
	solid_storage=(sit_arctic+snt_arctic)*sic_arctic*oarea
	solid_fw_storage=m3tokm3*solid_storage.sum(dim=('x','y'))

	print('Storage Done')

	#########2. Solid fluxes through Arctic Gateways [km3/yr]
	#u,v are not on same grid as tracers, etc so you need to compute velocities at tracer points, and you need to average
	#twice because velocities are catty-corner to tracer points on a B grid. This undoubtedly introduces errors.
	#For a given tracer index, u[index] is at the northeast corner

	##--------Bering Strait
	#Solid flux through Bering Strait: x=[129:132], y=212 #Indices on the full world map
	#Compute the fw volume relative to 34.8. Flow out of Bering is southward so the sign convention is correct
	j=212
	istart=129
	iend=132
	sicber=sic[:,j,istart:iend+1]
	sicber=sicber.where(xr.ufuncs.isfinite(sicber),other=0).compute()
	sitber=sit[:,j,istart:iend+1]
	sitber=sitber.where(xr.ufuncs.isfinite(sitber),other=0).compute()
	sntber=snt[:,j,istart:iend+1]
	sntber=sntber.where(xr.ufuncs.isfinite(sntber),other=0).compute()
	#Average meridional velocities in y first
	siv2=siv[:,j-1:j+1,istart-1:iend+1]
	siv2=siv2.where(xr.ufuncs.isfinite(siv2),other=0).compute()
	DYVber=DYV[j-1:j+1,istart-1:iend+1].load()
	sivber=(siv2[:,0]*DYVber[0]+siv2[:,1]*DYVber[1])/(DYVber[0]+DYVber[1])
	#Average in x
	DXVber=HUS[j,istart-1:iend+1].load()
	DXTber=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitber),coords=sitber.coords,dims=sitber.dims)
	vmid[:,0]=(sivber[:,0]*DXVber[0]+sivber[:,1]*DXVber[1])/sum(DXVber[0:2])
	vmid[:,1]=(sivber[:,1]*DXVber[1]+sivber[:,2]*DXVber[2])/sum(DXVber[1:])
	#Compute the fluxes
	ber_vel=vmid*DXTber
	ber_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitber
	ber_snow=rho_sn_over_rho_fw*sntber
	ber_solid=(ber_ice+ber_snow)*sicber*ber_vel
	ber_solid_flux=km3yr*ber_solid.sum(dim=('x'))


	##--------Nares Strait
	#Solid flux through Nares: x=[241,243], y=[234,236]. Indices on the full world map
	#This is a 45˚ line
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is west and south (negative) so the sign convention is correct
	jstart=234
	jend=236
	istart=241
	iend=243
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)
	
	siunar=siu[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	siunar=siunar.where(xr.ufuncs.isfinite(siunar),other=0).compute()
	sivnar=siv[:,jstart-1:jend+1,istart-1:iend+1] #one point larger in x,y for midpoint averaging
	sivnar=sivnar.where(xr.ufuncs.isfinite(sivnar),other=0).compute()
	sicnar=sic[:,jstart:jend+1,istart:iend+1]
	sicnar=sicnar.where(xr.ufuncs.isfinite(sicnar),other=0).compute()
	sitnar=sit[:,jstart:jend+1,istart:iend+1]
	sitnar=sitnar.where(xr.ufuncs.isfinite(sitnar),other=0).compute()
	sntnar=snt[:,jstart:jend+1,istart:iend+1]
	sntnar=sntnar.where(xr.ufuncs.isfinite(sntnar),other=0).compute()
	DYVnar=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNnar=HTN[jstart:jend+1,istart:iend+1].load()
	DXVnar=DXV[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEnar=HTE[jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(sitnar),coords=sitnar.coords,dims=sitnar.dims)
	vmid=xr.DataArray(data=np.zeros_like(sitnar),coords=sitnar.coords,dims=sitnar.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
		#Avg in x first
		umidx=(siunar[:,j:j+2,i]*DXVnar[j:j+2,i]+siunar[:,j:j+2,i+1]*DXVnar[j:j+2,i+1])/DXVnar[j:j+2,i:i+2].sum(dim='x')
		vmidx=(sivnar[:,j:j+2,i]*DXVnar[j:j+2,i]+sivnar[:,j:j+2,i+1]*DXVnar[j:j+2,i+1])/DXVnar[j:j+2,i:i+2].sum(dim='x')
		#Then in y
		umid[:,j,i]=(umidx[:,0]*DYVnar[j,i]+umidx[:,1]*DYVnar[j+1,i])/DYVnar[j:j+2,i].sum(dim='y')
		vmid[:,j,i]=(vmidx[:,0]*DYVnar[j,i]+vmidx[:,1]*DYVnar[j+1,i])/DYVnar[j:j+2,i].sum(dim='y')

	#Compute the fluxes
	#Can only do this if u,v (arrays) are the same size and colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Nares line makes a 45˚ angle with u and v
	nar_vel=(umid+vmid)*vel_scale*np.sqrt(HTEnar*HTEnar+HTNnar*HTNnar)
	nar_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitnar
	#snow is pure fw so salinity factor is 1
	nar_snow=rho_sn_over_rho_fw*sntnar
	#weight by sic in this step instead of in the separate snow and ice steps above
	nar_solid=(nar_ice+nar_snow)*sicnar*nar_vel
	nares_solid_flux=km3yr*nar_solid.sum(dim=('x','y')) #sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...


	##--------Barrow Strait
	#Solid flux through Barrow Strait: x=[212,218], y=242 #Indices on full world map
	#boundary out of Barrow
	#Compute the fw volume relative to 34.8
	#Northward (positive) velocities are into the Arctic so the sign convention is correct
	j=242
	istart=212
	iend=218
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
	DYVbrw=DYV[j-1:j+1,istart-1:iend+1].load()
	sivbrw=(siv2[:,0]*DYVbrw[0]+siv2[:,1]*DYVbrw[1])/(DYVbrw[0]+DYVbrw[1])
	#Average in x
	DXVbrw=HUS[j,istart-1:iend+1].load()
	DXTbrw=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitbrw),coords=sitbrw.coords,dims=sitbrw.dims)
	for i in range(0,len(irange)):
		vmid[:,i]=(sivbrw[:,i]*DXVbrw[i]+sivbrw[:,i+1]*DXVbrw[i+1])/sum(DXVbrw[i:i+2])
	#Compute the fluxes
	brw_vel=vmid*DXTbrw
	brw_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbrw
	#Snow is pure fw so the salinity scale factor is 1
	brw_snow=rho_sn_over_rho_fw*sntbrw
	brw_solid=(brw_ice+brw_snow)*sicbrw*brw_vel
	brw_solid_flux=km3yr*brw_solid.sum(dim=('x')) #sum over longitude


	##--------Fram Strait
	#Solid flux through Fram Strait, x=[280,295], y=238. Indices on full world. 
	#Fram is oriented such that northward (positive)
	#velocities are into the Arcitc, so the sign convention is correct
	#Compute the fw volume relative to 34.8
	istart=280
	iend=295
	j=238
	irange=np.arange(istart,iend+1)
	sicfram=sic[:,j,istart:iend+1]
	sicfram=sicfram.where(xr.ufuncs.isfinite(sicfram),other=0).compute()
	sitfram=sit[:,j,istart:iend+1]
	sitfram=sitfram.where(xr.ufuncs.isfinite(sitfram),other=0).compute()
	sntfram=snt[:,j,istart:iend+1]
	sntfram=sntfram.where(xr.ufuncs.isfinite(sntfram),other=0).compute()
	#Average meridional velocities in y first
	siv2=siv[:,j-1:j+1,istart-1:iend+1]
	siv2=siv2.where(xr.ufuncs.isfinite(siv2),other=0).compute()
	DYVfram=DYV[j-1:j+1,istart-1:iend+1].load()
	sivfram=(siv2[:,0]*DYVfram[0]+siv2[:,1]*DYVfram[1])/(DYVfram[0]+DYVfram[1])
	#Average in x
	DXVfram=HUS[j,istart-1:iend+1].load()
	DXTfram=DXT[j,istart:iend+1].load()
	vmid=xr.DataArray(data=np.zeros_like(sitfram),coords=sitfram.coords,dims=sitfram.dims)
	for i in range(0,len(irange)):
		vmid[:,i]=(sivfram[:,i]*DXVfram[i]+sivfram[:,i+1]*DXVfram[i+1])/sum(DXVfram[i:i+2])
	#Compute the fluxes
	fram_vel=vmid*DXTfram
	fram_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitfram
	#Snow is pure fw so the salinity scale factor is 1
	fram_snow=rho_sn_over_rho_fw*sntfram
	fram_solid=(fram_ice+fram_snow)*sicfram*fram_vel
	fram_solid_flux=km3yr*fram_solid.sum(dim=('x')) #sum over longitude


	##--------Barents Sea Opening (BSO)
	#Solid flux through the BSO, x=[304,317], y=[223,236]. Indices on full world. 
	#This is is more or less a 45˚ line so it doesn't have a complex pattern 
	#for the grid cells it passes through
	#The BSO is oriented such that eastward and northward (positive) velocities are 
	#into the Arctic
	#Compute the fw volume relative to 34.8
	istart=304
	iend=317
	jstart=223
	jend=236
	i_range=np.arange(istart,iend+1)
	j_range=np.arange(jend,jstart-1,-1)

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
	DYVbso=HUW[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTNbso=HTN[jstart:jend+1,istart:iend+1].load()
	DXVbso=DXV[jstart-1:jend+1,istart-1:iend+1].load() #one point larger in x,y for midpoint averaging
	HTEbso=HTE[jstart:jend+1,istart:iend+1].load()
	umid=xr.DataArray(data=np.zeros_like(sitbso),coords=sitbso.coords,dims=sitbso.dims)
	vmid=xr.DataArray(data=np.zeros_like(sitbso),coords=sitbso.coords,dims=sitbso.dims)

	i_range_enum=np.arange(len(i_range))
	for (i,j) in zip(i_range_enum,reversed(i_range_enum)):
		#Avg in x first
		umidx=(siubso[:,j:j+2,i]*DXVbso[j:j+2,i]+siubso[:,j:j+2,i+1]*DXVbso[j:j+2,i+1])/DXVbso[j:j+2,i:i+2].sum(dim='x')
		vmidx=(sivbso[:,j:j+2,i]*DXVbso[j:j+2,i]+sivbso[:,j:j+2,i+1]*DXVbso[j:j+2,i+1])/DXVbso[j:j+2,i:i+2].sum(dim='x')
		#Then in y
		umid[:,j,i]=(umidx[:,0]*DYVbso[j,i]+umidx[:,1]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')
		vmid[:,j,i]=(vmidx[:,0]*DYVbso[j,i]+vmidx[:,1]*DYVbso[j+1,i])/DYVbso[j:j+2,i].sum(dim='y')

	#Compute the fluxes
	#Can only do this if u,v arrays are same size and colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the BSO line makes a 45˚ angle with u and v
	bso_vel=(umid+vmid)*vel_scale*np.sqrt(HTEbso*HTEbso+HTNbso*HTNbso)
	bso_ice=rho_si_over_rho_fw*((s_ref-s_si)/s_ref)*sitbso
	#snow is pure fw so salinity factor is 1
	bso_snow=rho_sn_over_rho_fw*sntbso
	bso_solid=(bso_ice+bso_snow)*sicbso*bso_vel
	bso_solid_flux=km3yr*bso_solid.sum(dim=('x','y')) #sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
		

	print('Solid fluxes through gateways done!')

	#########3. Liquid fluxes through Arctic Gateways [km3/yr]

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

	#Compute the fluxes
	#Volume flux
	ber_vol=vmid*dzber*DXTber
	ber_vol_flux=km3yr*ber_vol.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#FW flux
	ber_liq=((s_ref-saltber)/s_ref)*ber_vol
	ber_liq_flux=km3yr*ber_liq.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	ber_liq_flux_below_sref=km3yr*ber_liq.where(saltber<=s_ref,other=0).sum(dim=('lev','x'))
	ber_liq_flux_above_sref=km3yr*ber_liq.where(saltber>s_ref,other=0).sum(dim=('lev','x')) 


	##--------Nares Strait
	#Liquid flux through Nares: x=[241,243], y=[234,236]. Indices on the full world map
	#Compute the fw volume relative to 34.8
	#Flow out of Nares is west and south (negative) so the sign convention is correct
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
	#Compute fluxes
	#Can only do this if u,v arrays are same sizes and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Nares line makes a 45˚ angle with u and v
	#Volume flux
	nar_vol=dznar*(umid+vmid)*vel_scale*np.sqrt(HTEnar*HTEnar+HTNnar*HTNnar)
	nares_vol_flux=km3yr*nar_vol.sum(dim=('lev','y','x'))
	#FW flux
	nar_liq=((s_ref-saltnar)/s_ref)*nar_vol
	nares_liq_flux=km3yr*nar_liq.sum(dim=('lev','y','x'))#sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
	#Fluxes above and below sref
	nares_liq_flux_below_sref=km3yr*nar_liq.where(saltnar<=s_ref,other=0).sum(dim=('lev','y','x'))
	nares_liq_flux_above_sref=km3yr*nar_liq.where(saltnar>s_ref,other=0).sum(dim=('lev','y','x')) 

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

	#Compute the fluxes
	#Volume flux
	brw_vol=vmid*dzbrw*DXTbrw
	brw_vol_flux=km3yr*brw_vol.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#FW flux
	brw_liq=((s_ref-saltbrw)/s_ref)*brw_vol
	brw_liq_flux=km3yr*brw_liq.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	brw_liq_flux_below_sref=km3yr*brw_liq.where(saltbrw<=s_ref,other=0).sum(dim=('lev','x'))
	brw_liq_flux_above_sref=km3yr*brw_liq.where(saltbrw>s_ref,other=0).sum(dim=('lev','x'))


	##--------Fram Strait
	#Liquid flux through Fram Strait, x=[280,295], y=238. Indices on full world.
	#Fram is oriented such that northward (positive)
	#velocities are into the Arctic, so the sign convention is correct
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

	#Compute the fluxes
	#Volume flux
	fram_vol=vmid*dzfram*DXTfram
	fram_vol_flux=km3yr*fram_vol.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#FW flux
	fram_liq=((s_ref-saltfram)/s_ref)*fram_vol
	fram_liq_flux=km3yr*fram_liq.sum(dim=('lev','x')) #only one lat pt so no sum over lat
	#Fluxes above and below sref and the volume flux
	fram_liq_flux_below_sref=km3yr*fram_liq.where(saltfram<=s_ref,other=0).sum(dim=('lev','x'))
	fram_liq_flux_above_sref=km3yr*fram_liq.where(saltfram>s_ref,other=0).sum(dim=('lev','x'))


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
	#Compute the fluxes
	#Can only do this if u,v arrays are same size and values are colocated
	vel_scale=0.5*np.sqrt(2) #The unit normal into the Arctic for the Fram line makes a 45˚ angle with u and v
	#Volume flux
	bso_vol=dzbso*(umid+vmid)*vel_scale*np.sqrt(HTEbso*HTEbso+HTNbso*HTNbso)
	bso_vol_flux=km3yr*bso_vol.sum(dim=('lev','y','x'))#sum over latitude and lon.
	#FW flux
	bso_liq=((s_ref-saltbso)/s_ref)*bso_vol
	bso_liq_flux=km3yr*bso_liq.sum(dim=('lev','y','x'))#sum over latitude and lon. This is fine as long as
	#umid and vmid are zero everywhere except where you computed it...
	#Fluxes above and below sref
	bso_liq_flux_below_sref=km3yr*bso_liq.where(saltbso<=s_ref,other=0).sum(dim=('lev','y','x'))
	bso_liq_flux_above_sref=km3yr*bso_liq.where(saltbso>s_ref,other=0).sum(dim=('lev','y','x'))

	print('Liquid fluxes through gateways done!')

	#########3.Compute flux totals [km3/yr]
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


	#########4. Save the output
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

