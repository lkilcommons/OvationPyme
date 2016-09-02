import matplotlib.pyplot as pp
import ovation_prime
import ovation_utilities
import datetime
import numpy as np
from geospacepy import satplottools,special_datetime

def draw_interpolated_conductance(new_mlat_grid,new_mlt_grid,dt,startdt,enddt,hemi):
	"""
	Interpolate hall and pedersen conductance
	onto grid described by mlat_grid,mlt_grid
	""" 
	estimator = ovation_prime.ConductanceEstimator(startdt,enddt)
	
	mlatgrid,mltgrid,pedgrid,hallgrid = estimator.get_conductance(dt,hemi=hemi,auroral=True,solar=True)

	ped_interpolator = ovation_prime.LatLocaltimeInterpolator(mlatgrid,mltgrid,pedgrid)
	new_pedgrid = ped_interpolator.interpolate(new_mlat_grid,new_mlt_grid)

	hall_interpolator = ovation_prime.LatLocaltimeInterpolator(mlatgrid,mltgrid,hallgrid)
	new_hallgrid = hall_interpolator.interpolate(new_mlat_grid,new_mlt_grid)

	f = pp.figure(figsize=(11,5))
	aH = f.add_subplot(121)
	aP = f.add_subplot(122)

	X,Y = satplottools.latlt2cart(new_mlat_grid.flatten(),new_mlt_grid.flatten(),hemi)
	X = X.reshape(new_mlat_grid.shape)
	Y = Y.reshape(new_mlt_grid.shape)
	
	satplottools.draw_dialplot(aH)
	satplottools.draw_dialplot(aP)
	
	mappableH = aH.pcolormesh(X,Y,new_hallgrid,vmin=0.,vmax=20.)
	mappableP = aP.pcolormesh(X,Y,new_pedgrid,vmin=0.,vmax=15.)
	
	aH.set_title("Hall Conductance")
	aP.set_title("Pedersen Conductance")

	f.colorbar(mappableH,ax=aH)
	f.colorbar(mappableP,ax=aP)
	
	f.suptitle("OvationPyme Interpolated Conductance %s Hemisphere at %s" % (hemi,dt.strftime('%c')),
		fontweight='bold')
	return f

def draw_conductance(dt,startdt,enddt,hemi):
	"""
	Get the hall and pedersen conductance for one date and hemisphere
	"""
	estimator = ovation_prime.ConductanceEstimator(startdt,enddt)
	
	mlatgrid,mltgrid,pedgrid,hallgrid = estimator.get_conductance(dt,hemi=hemi,auroral=True,solar=True)

	f = pp.figure(figsize=(11,5))
	aH = f.add_subplot(121)
	aP = f.add_subplot(122)

	X,Y = satplottools.latlt2cart(mlatgrid.flatten(),mltgrid.flatten(),hemi)
	X = X.reshape(mlatgrid.shape)
	Y = Y.reshape(mltgrid.shape)
	
	satplottools.draw_dialplot(aH)
	satplottools.draw_dialplot(aP)
	
	mappableH = aH.pcolormesh(X,Y,hallgrid,vmin=0.,vmax=20.)
	mappableP = aP.pcolormesh(X,Y,pedgrid,vmin=0.,vmax=15.)
	
	aH.set_title("Hall Conductance")
	aP.set_title("Pedersen Conductance")

	f.colorbar(mappableH,ax=aH)
	f.colorbar(mappableP,ax=aP)
	
	f.suptitle("OvationPyme Conductance Output %s Hemisphere at %s \n" % (hemi,dt.strftime('%c')),
		fontweight='bold')

	return f


def draw_weighted_flux(dt,atype='diff',jtype='electron energy flux'):
	"""
	Test automatic generation of omni_intervals in ovation_utilities
	also by not specifying a start and end time for the FluxEstimator
	"""
	estimator = ovation_prime.FluxEstimator(atype,jtype)
	
	mlatgridN,mltgridN,fluxgridN = estimator.get_flux_for_time(dt,hemi='N')
	mlatgridS,mltgridS,fluxgridS = estimator.get_flux_for_time(dt,hemi='S')

	f = pp.figure(figsize=(11,5))
	aN = f.add_subplot(121)
	aS = f.add_subplot(122)

	XN,YN = satplottools.latlt2cart(mlatgridN.flatten(),mltgridN.flatten(),'N')
	XS,YS = satplottools.latlt2cart(mlatgridS.flatten(),mltgridS.flatten(),'S')
	XN = XN.reshape(mlatgridN.shape)
	YN = YN.reshape(mltgridN.shape)
	XS = XS.reshape(mlatgridS.shape)
	YS = YS.reshape(mltgridS.shape)
	
	satplottools.draw_dialplot(aN)
	satplottools.draw_dialplot(aS)
	
	mappableN = aN.pcolormesh(XN,YN,fluxgridN,vmin=0,vmax=2)
	mappableS = aS.pcolormesh(XS,YS,fluxgridS,vmin=0,vmax=2)
	
	#aN.set_title("Northern Hemisphere Flux")
	#aS.set_title("Southern Hemisphere Flux")

	f.colorbar(mappableN,ax=aN)
	f.colorbar(mappableS,ax=aS)
	
	f.suptitle("OvationPyme Auroral Model Flux Output at %s \n AuroralType:%s, FluxType:%s" % (dt.strftime('%c'),
																							atype,jtype),
																							fontweight='bold')
	return f

def draw_seasonal_flux(seasonN='summer',seasonS='winter',atype='diff',jtype='electron energy flux'):

	dF = 2134.17

	estimatorN = ovation_prime.SeasonalFluxEstimator(seasonN,atype,jtype)
	estimatorS = ovation_prime.SeasonalFluxEstimator(seasonS,atype,jtype)

	fluxtupleN = estimatorN.get_gridded_flux(dF,combined_N_and_S=False)
	(mlatgridN,mltgridN,fluxgridN) = fluxtupleN[:3]

	fluxtupleS = estimatorS.get_gridded_flux(dF,combined_N_and_S=False)
	(mlatgridS,mltgridS,fluxgridS) = fluxtupleS[3:]
	
	f = pp.figure(figsize=(11,5))
	aN = f.add_subplot(121)
	aS = f.add_subplot(122)

	f2 = pp.figure(figsize=(5,5))
	a2 = f2.add_subplot(111)

	XN,YN = satplottools.latlt2cart(mlatgridN.flatten(),mltgridN.flatten(),'N')
	XS,YS = satplottools.latlt2cart(mlatgridS.flatten(),mltgridS.flatten(),'S')
	XN = XN.reshape(mlatgridN.shape)
	YN = YN.reshape(mltgridN.shape)
	XS = XS.reshape(mlatgridS.shape)
	YS = YS.reshape(mltgridS.shape)
	
	satplottools.draw_dialplot(aN)
	satplottools.draw_dialplot(aS)
	satplottools.draw_dialplot(a2)
	
	mappableN = aN.pcolormesh(XN,YN,fluxgridN,vmin=0,vmax=2)
	mappableS = aS.pcolormesh(XS,YS,fluxgridS,vmin=0,vmax=2)
	mappableNS = a2.pcolormesh(XN,YN,(fluxgridS+fluxgridN)/2,vmin=0,vmax=2)
	
	#aN.set_title("Northern Hemisphere Flux")
	#aS.set_title("Southern Hemisphere Flux")

	f.colorbar(mappableN,ax=aN)
	f.colorbar(mappableS,ax=aS)
	f.colorbar(mappableNS,ax=a2)
	
	f.suptitle("OvationPyme Auroral Model Raw Flux Output \n Season:%s, AuroralType:%s, FluxType:%s, Newell Coupling:%.3f" % (seasonN,atype,jtype,dF),
		fontweight='bold')

	f2.suptitle("OvationPyme Combined Hemisphere Output \n Season:%s, AuroralType:%s, FluxType:%s, Newell Coupling:%.3f" % (seasonN,atype,jtype,dF),
		fontweight='bold')

	return f,f2


if __name__ == '__main__':

	seasonN='summer'
	seasonS='winter'
	atype='diff'
	jtype='electron energy flux'
	tfmt = '%Y%m%d'

	startdt = datetime.datetime(2011,11,29,0,0,0)
	enddt = datetime.datetime(2011,12,1,0,0,0)

	#Times given in figure 2 of Cousins et. al. 2015
	dt1 = datetime.datetime(2011,11,30,12,10,0)
	dt2 = datetime.datetime(2011,11,29,7,40,0)
	dt3 = datetime.datetime(2011,11,29,0,50,0)
	
	for dt in [dt1,dt2,dt3]:

		new_mlat,new_mlt = np.meshgrid(np.linspace(60.,80.,40.),np.linspace(2.,6.,30.))
		fiN = draw_interpolated_conductance(new_mlat,new_mlt,dt,startdt,enddt,'N')
		fiN.savefig('ovation_conductance_interp_N_%s.png' % (dt.strftime(tfmt)))
		
		fiS = draw_interpolated_conductance(new_mlat,new_mlt,dt,startdt,enddt,'N')
		fiS.savefig('ovation_conductance_interp_S_%s.png' % (dt.strftime(tfmt)))
		
		f2N = draw_conductance(dt,startdt,enddt,'N')
		f2N.savefig('ovation_conductance_N_%s.png' % (dt.strftime(tfmt)))
		
		f2S = draw_conductance(dt,startdt,enddt,'S')
		f2S.savefig('ovation_conductance_S_%s.png' % (dt.strftime(tfmt)))
		
		f1 = draw_weighted_flux(dt)
		f1.savefig('ovation_combflux_%s_%s_%s.png' % (atype,jtype.replace(' ','_'),dt.strftime(tfmt)))

		f3,f4 = draw_seasonal_flux(seasonN=seasonN,seasonS=seasonS,atype=atype,jtype=jtype)
		
		f3.savefig('ovation_rawflux_N%s_S%s_%s_%s.png' % (seasonN,seasonS,atype,jtype.replace(' ','_')))

		for f in [fiN,fiS,f2N,f2S,f1,f3,f4]:
			pp.close(f)

		