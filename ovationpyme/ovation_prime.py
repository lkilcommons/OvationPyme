"""
This module contains the main model routines for 
Ovation Prime (historically called season_epoch.pro)
"""
import scipy.interpolate as interpolate
import numpy as np
from ovationpyme import ovation_utilities
import geospacepy
import os

#Determine where this module's source file is located
#to determine where to look for the tables
src_file_dir = os.path.dirname(os.path.realpath(__file__))
ovation_datadir = os.path.join(src_file_dir,'data')

class OvationConductance(object):
	"""
	Implements the 'Robinson Formula'
	for estimating Pedersen and Hall height integrated
	conductivity (conducance) from
	average electron energy and 
	total electron energy flux 
	(assumes a Maxwellian electron energy distribution)
	"""
	def __init__(self,dt):
		#Use diffuse aurora only
		self.eflux_estimator = FluxEstimator('diff','electron energy flux')
		self.eavg_estimator = FluxEstimator('diff','electron average energy')
		

class FluxEstimator(object):
	"""
	A class which estimates auroral flux
	based on the Ovation Prime regressions,
	at arbitrary locations and times.

	Locations are in magnetic latitude and local
	time, and are interpolated using a B-spline
	representation
	"""
	def __init__(self,atype,jtype,seasonal_estimators=None,start_dt=None,end_dt=None):
		"""

		doy - int
			day of year

		atype - str, ['diff','mono','wave','ions']
			type of aurora for which to load regression coeffients

		jtype - int or str
			1:"electron energy flux",
			2:"ion energy flux",
			3:"electron number flux",
			4:"ion number flux",
			5:"electron average energy",
			6:"ion average energy"

			Type of flux you want to estimate

		seasonal_estimators - dict, optional
			A dictionary of SeasonalFluxEstimators for seasons
			'spring','fall','summer','winter', if you 
			don't want to create them
			(for efficiency across multi-day calls)

		"""
		self.atype = atype #Type of aurora
		self.jtype = jtype #Type of flux
		if start_dt is not None and end_dt is not None:
			self.oi = geospacepy.omnireader.omni_interval(start_dt-datetime.timedelta(days=1),
															end_dt+datetime.timedelta(days=1),
															'5min',silent=True) #Give 1 day +- buffer because we need avg SW

			#Pre-create an omni interval (for speed if you are estimating auroral flux across many days)
		else:
			self.oi = None #omni_interval objects will be created on-the-fly (slow, but fine for single calls to get_flux)

		seasons = ['spring','summer','fall','winter']

		if seasonal_estimators is None:
			#Make a seasonal estimator for each season with nonzero weight
			self.seasonal_flux_estimators = {season:SeasonalFluxEstimator(seasons,atype,jtype) for season in seasons}				
		else:
			#Ensure the passed seasonal estimators are approriate for this atype and jtype
			for season,estimator in seasonal_flux_estimators.iteritems():
				jtype_atype_ok = jtype_atype_ok and (self.jtype == estimator.jtype and self.atype == estimator.atype)
			if not jtype_atype_ok:
				raise RuntimeError('Auroral and flux type of SeasonalFluxEstimators do not match %s and %s!' % (self.atype,self.jtype))

	def get_flux_for_time(dt,hemi='N'):
		"""
		doy must be single value
		mlats and mlts can be arbitary shape, but both must be same shape
		"""
		doy = dt.timetuple().tm_yday
		if hemi == 'S':
			doy = 365.-doy #Use opposite season coefficients to get southern hemisphere results

		weights = self.season_weights(doy)
		dF = ovation_utilities.calc_avg_solarwind(dt,oi=self.oi)
		seasonal_flux = {}
		n_mlat_bins = self.seasonal_estimators.items()[0].n_mlat_bins/2 #div by 2 because combined N&S hemispheres
		n_mlt_bins = self.seasonal_estimators.items()[0].n_mlt_bins
		gridflux = np.zeros(n_mlat_bins,n_mlt_bins)
		for season in weights:
			if weights[season] > 0.:
				flux_estimator = self.seasonal_flux_estimators[season]
				grid_mlats,grid_mlts,grid_seasonflux = flux_estimator.get_gridded_flux(dF)
				gridflux += grid_seasonflux*weights[season]

		if hemi == 'S':
			grid_mlats = -1.*grid_mlats #by default returns positive latitudes

		return grid_mlats,grid_mlts,gridflux

	def season_weights(self,doy):
		"""
		Determines the relative weighting of the 
		model coeffecients for the various seasons for a particular
		day of year (doy). Nominally, weights the seasons
		based on the difference between the doy and the peak
		of the season (solstice/equinox)

		Returns:
			a dictionary with a key for each season. 
			Each value in the dicionary is a float between 0 and 1
		"""
		weight = {'winter':0.,'spring':0.,'summer':0.,'fall':0.}
		winter_w,spring_w,summer_w,fall_w = 0.,0.,0.,0.

		if doy >= 79. and doy < 171:
		   weight['summer'] = 1. - (171.-doy)/92.
		   weight['spring'] = 1. - weight['summer']

		elif doy >= 171. and doy < 263.:
		   weight['fall'] = 1. - (263.-doy)/92.
		   weight['summer'] = 1. - weight['fall']
		   
		elif doy >= 263. and doy < 354.:
		   weight['winter'] = 1. - (354.-doy)/91.
		   weight['fall'] = 1. - weight['winter']
		
		elif doy >= 354 or doy < 79:
			#For days of year > 354, subtract 365 to get negative
			#day of year values for computation
			doy0 = doy- 365. if doy >= 354 else doy
			weight['spring'] = 1. - (79.-doy0)/90.
			weight['winter'] = 1. - weight['spring'] 

		return weight

class SeasonalFluxEstimator(object):
	"""
	A class to hold and caculate predictions from the regression coeffecients
	which are tabulated in the data/premodel/{season}_{atype}_*.txt
	files.

	Given a particular season, type of aurora ( one of ['diff','mono','wave','ions'])
	and type of flux, returns 
	"""
	def __init__(self,season,atype,jtype):
		"""

		season - str,['winter','spring','summer','fall']
			season for which to load regression coeffients

		atype - str, ['diff','mono','wave','ions']
			type of aurora for which to load regression coeffients

		jtype - int or str
			1:"electron energy flux",
			2:"ion energy flux",
			3:"electron number flux",
			4:"ion number flux",
			5:"electron average energy",
			6:"ion average energy"
		
		"""
		nmlt = 96                           #number of mag local times in arrays (resolution of 15 minutes)
		nmlat = 160                         #number of mag latitudes in arrays (resolution of 1/4 of a degree (.25))
		ndF = 12							#number of coupling strength bins
		self.jtype,self.atype = jtype,atype

		self.n_mlt_bins,self.n_mlat_bins,self.n_dF_bins = nmlt,nmlat,ndF

		self.mlats = np.concatenate([np.linspace(-90.,-50.,self.n_mlat_bins/2),
											np.linspace(50.,90.,self.n_mlat_bins/2)])

		self.mlts = np.linspace(0.,24.,self.n_mlt_bins)
		
		self.fluxtypes = {1:"electron energy flux",
							2:"ion energy flux",
							3:"electron number flux",
							4:"ion number flux",
							5:"electron average energy",
							6:"ion average energy"}

		#Determine file names
		file_suffix = '_n' if (jtype in [3,4] or 'number flux' in jtype) else ''
		self.afile = os.path.join(ovation_datadir,'premodel/%s_%s%s.txt' % (season,atype,file_suffix))
		self.pfile = os.path.join(ovation_datadir,'premodel/%s_prob_b_%s.txt' % (season,atype))
		#Defualt values of header (don't know why need yet)
		# b1 = 0. 
		# b2 = 0.
		# yend = 1900
		# dend = 1
		# y0 = 1900
		# d0 = 1
		# files_done = 0
		# sf0 = 0
		self.valid_atypes = ['diff','mono','wave','ions']

		with open(self.afile,'r') as f:
			aheader = f.readline() # y0,d0,yend,dend,files_done,sf0
			print "Read Auroral Flux Coefficient File %s,\n Header: %s" % (self.afile,aheader)
			# Don't know if it will read from where f pointer is after reading header line
			adata = np.genfromtxt(f,max_rows=nmlat*nmlt)
			print "First line was %s" % (str(adata[0,:]))

		self.b1a, self.b2a = np.zeros((nmlt,nmlat)), np.zeros((nmlt,nmlat))
		self.b1a.fill(np.nan)
		self.b2a.fill(np.nan)
		mlt_bin_inds,mlat_bin_inds = adata[:,0].astype(int),adata[:,1].astype(int) 
		self.b1a[mlt_bin_inds,mlat_bin_inds]=adata[:,2]
		self.b2a[mlt_bin_inds,mlat_bin_inds]=adata[:,3]

		self.b1p, self.b2p = np.zeros((nmlt,nmlat)), np.zeros((nmlt,nmlat))
		self.prob = np.zeros((nmlt,nmlat,ndF))
		self.b1p.fill(np.nan)
		self.b2p.fill(np.nan)
		self.prob.fill(np.nan)
		#pdata has 2 columns, b1, b2 for first 15361 rows
		#pdata has nmlat*nmlt rows (one for each positional bin)

		if atype in ['diff','mono','wave']:
			with open(self.pfile,'r') as f:
				pheader = f.readline() #y0,d0,yend,dend,files_done,sf0
				# Don't know if it will read from where f pointer is after reading header line
				pdata_b = np.genfromtxt(f,max_rows=nmlt*nmlat) # 2 columns, b1 and b2
				#print "Shape of b1p,b2p should be nmlt*nmlat=%d, is %s" % (nmlt*nmlat,len(pdata_b[:,0]))
				pdata_p = np.genfromtxt(f,max_rows=nmlt*nmlat*ndF) # 1 column, pval

			#in the file the probability is stored with coupling strength bin
			#varying fastest (this is Fortran indexing order)
			pdata_p_column_dFbin = pdata_p.reshape((-1,ndF),order='F')

			#I don't know why this is not used for atype 'ions'
			#mlt is first dimension
			self.b1p[mlt_bin_inds,mlat_bin_inds]=pdata_b[:,0]
			self.b2p[mlt_bin_inds,mlat_bin_inds]=pdata_b[:,1]
			for idF in range(ndF):
				self.prob[mlt_bin_inds,mlat_bin_inds,idF]=pdata_p_column_dFbin[:,idF]

			if season=='spring' and atype=='diff' and jtype=='electron energy flux':
				print self.b1p[22:26,138:142]
				print self.prob[22:26,138:142,5]

		#IDL original read
		#readf,20,i,j,b1,b2,rF
		#;;   b1a_all(atype, iseason,i,j) = b1
		#;;   b2a_all(atype, iseason,i,j) = b2
		#adata has 5 columns, mlt bin number, mlat bin number, b1, b2, rF
		#adata has nmlat*nmlt rows (one for each positional bin)

	def which_dF_bin(dF):
		"""
		
		Given a coupling strength value, finds the bin it falls into
		
		"""
		dFave = 4421. #Magic numbers!
		dFstep = dFave/8.
		i_dFbin = np.round(dF/dFstep)
		#Range check 0 <= i_dFbin <= n_dF_bins-1
		if i_dFbin < 0 or i_dFbin > self.n_dF_bins-1: 
			i_dFbin = 0 if i_dFbin < 0 else self.n_dF_bins-1
		return i_dFbin

	def prob_estimate(self,dF,i_mlt_bin,i_mlat_bin):
		"""

		Estimate probability of <something> by using tabulated
		linear regression coefficients ( from prob_b files ) 
		WRT coupling strength dF (which are different for each position bin)
		
		If p doesn't come out sensible by the initial regression, 
		(i.e both regression coefficients are zero)
		then tries loading from the probability array. If the value
		in the array is zero, then estimates a value using adjacent 
		coupling strength bins in the probability array.

		""" 
		#Look up the regression coefficients
		b1,b2 = self.b1p[i_mlt_bin,i_mlat_bin],self.b2p[i_mlt_bin,i_mlat_bin]

		p = b1 + b2*dF #What is this the probability of?
		
		#range check 0<=p<=1
		if p < 0. or p > 1.:
			p = 1. if p > 1. else 0.

		if b1 == 0. and b2 == 0.:
			i_dFbin = self.which_dF_bin(dF)
			#Get the tabulated probability
			p = self.prob[i_mlt_bin,i_mlat_bin,i_dFbin]
			
			if p == 0.:
				#If no tabulated probability we must estimate by interpolating
				#between adjacent coupling strength bins
				i_dFbin_1 = i_dFbin - 1 if i_dFbin > 0 else i_dFbin+2 #one dF bin before by preference, two after in extremis
				i_dFbin_2 = i_dFbin + 1 if i_dFbin < self.n_dF_bins-1 else i_dFbin-2 #one dF bin after by preference, two before in extremis
				p = (self.prob[i_mlt_bin,i_mlat_bin,i_dFbin_1] + self.prob[i_mlt_bin,i_mlat_bin,i_dFbin_2])/2.
		
		return p

	def estimate_auroral_flux(self,dF,i_mlt_bin,i_mlat_bin):
		"""
		Does what it says on the tin,
		estimates the flux using the regression coeffecients in the 'a' files

		fluxtype can use either numerical code or string description in self.fluxtypes
		"""
		b1,b2 = self.b1a[i_mlt_bin,i_mlat_bin],self.b2a[i_mlt_bin,i_mlat_bin]
		p = self.prob_estimate(dF,i_mlt_bin,i_mlat_bin)
		print p,b1,b2,dF
		flux = (b1+b2*dF)*p
		return self.correct_flux(flux)
		
	def correct_flux(self,flux):
		"""
		A series of magical (unexplained,unknown) corrections to flux given a particular
		type of flux
		"""
		fluxtype = self.jtype

		if flux < 0.:
			flux = 0.

		if self.atype is not 'ions':
			#Electron Energy Flux
			if fluxtype in [1,self.fluxtypes[1]]:
			  if flux > 10.:
				flux = 0.5
			  elif flux > 5.:
				flux = 5.
			
			#Electron Number Flux
			if fluxtype in [3,self.fluxtypes[3]]:
			  if flux > 2.0e9:
				flux = 1.0e9
			  elif flux > 2.0e10:
				flux = 0.
		else:
			#Ion Energy Flux
			if fluxtype in [2,self.fluxtypes[2]]:
			  if flux > 2.:
				flux = 2.
			  elif flux > 4.:
				flux = 0.25
			  
			#Ion Number Flux
			if fluxtype in [4,self.fluxtypes[4]]:
			  if flux > 1.0e8:
				flux = 1.0e8
			  elif flux > 5.0e8:
				flux = 0.
		return flux

	def get_gridded_flux(self,dF,combined_N_and_S=True):
		"""
		Return the flux interpolated onto arbitary locations
		in mlats and mlts
		
		combined_N_and_S, bool, optional

			Average the fluxes for northern and southern hemisphere
			and use them for both hemispheres (this is what standard
			ovation prime does always I think, so I've made it default)
			The original code says that this result is appropriate for 
			the northern hemisphere, and to use 365 - actual doy to
			get a combined result appropriate for the southern hemisphere

		"""

		fluxgridN = np.zeros((self.n_mlat_bins/2,self.n_mlt_bins))
		fluxgridN.fill(np.nan)
		#Make grid coordinates
		mlatgridN,mltgridN = np.meshgrid(self.mlats[self.n_mlat_bins/2:],self.mlts)
		
		fluxgridS = np.zeros((self.n_mlat_bins/2,self.n_mlt_bins))
		fluxgridS.fill(np.nan)
		#Make grid coordinates
		mlatgridS,mltgridS = np.meshgrid(self.mlats[:self.n_mlat_bins/2],self.mlts)
		
		for i_mlt in range(self.n_mlt_bins):
			for j_mlat in range(self.n_mlat_bins/2):
				#The mlat bins are orgainized like -90:dlat:-50,50:dlat:90
				north_grid[j_mlat,i_mlt] = self.estimate_auroral_flux(dF,i_mlt,self.n_mlat_bins/2+j_mlat)
				south_grid[j_mlat,i_mlt] = self.estimate_auroral_flux(dF,i_mlt,j_mlat)

		if not combined_N_and_S:
			return mlatgridN,mltgridN,fluxgridN,mlatgridS,mltgridS,fluxgridS
		else:
			return mlatgridN,mltgridN,(fluxgridN+fluxgridS)/2.










