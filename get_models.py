
#######get_models.py########
#Houses a function for retrieving the correct strings
#for different models and filenames for CMIP6, mainly for the
#historical and SSPs 585 and 126

def choose_model_and_sim_cmip6(model_name='CESM2',simulation_name='historical', variant_label='r1i1p1f1'):
	if simulation_name=='historical':
		mip='CMIP'
		time_period='185001-201412'
		simulation=model_name+'_historical_'+variant_label+'_gn_'
	elif 'ssp' in simulation_name:
		mip='ScenarioMIP'
		time_period='201501-210012'
		simulation=model_name+'_'+simulation_name.lower()+'_'+variant_label+'_gn_'
	return simulation,mip,time_period


