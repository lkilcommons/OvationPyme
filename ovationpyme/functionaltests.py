import matplotlib.pyplot as pp
import ovation_prime
import ovation_utilities
from geospacepy import satplottools

if __name__ == '__main__':

	seasonN,seasonS,atype,jtype = 'summer','winter','diff','electron energy flux'

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
	

	pp.show()
	f.savefig('ovation_rawflux_N%s_S%s_%s_%s.png' % (seasonN,seasonS,atype,jtype.replace(' ','_')))