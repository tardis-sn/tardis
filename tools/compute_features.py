#!/usr/bin/env python

############################  IMPORTS  #################################


import os																
import sys
import time
import warnings

path_S																	= os.environ["PATH_subtypes_main"]
sys.path.insert(1,path_S+'codes/modules/')

import numpy															as np
import pandas															as pd
import matplotlib.pyplot												as plt
import matplotlib														as mpl

import scipy

#from scipy.signal														import savgol_filter
from scipy.optimize														import curve_fit, OptimizeWarning
from matplotlib.ticker 													import MultipleLocator
from astropy															import constants as const
from scipy.interpolate													import interp1d

import read_data_v2														as re
import utility															as ut
import astropy
import math			
from math																import factorial	
								
mpl.rcParams['mathtext.fontset']										= 'stix'
mpl.rcParams['mathtext.fontset']										= 'stix'
mpl.rcParams['font.family']												= 'STIXGeneral'

warnings.simplefilter('error',OptimizeWarning)

#############################  CONSTANTS  ##############################

c																		= const.c.to('km/s').value

##########################  MASTER CODE  ###############################

class Utility(object):
	"""
	Utility functions
	"""
	
	def __init__(self, sep=20.):
		self.utility													= True
		self.sep														= sep
		self.keys														= ['6', '7']
		self.keys_to_fit												= ['6', '7']
		
	def rms(self,y_data,y_pred):
		rms_array														= math.sqrt(((y_data - y_pred)**2.).mean())
		rms_array														= np.asarray(rms_array).astype('float')
		rms_array[rms_array == 0]										= 1.e-5	
		return rms_array

	def make_parabola(self, x_ref):
		def parabola(x,a,b,c):
			return a*(x-x_ref)**2.+b*(x-x_ref)+c
		return parabola		
	
class analyse_spectra(Utility):

	"""
	Description
	----------
	TAKES A DATASET WITH COLUMNS ['wavelength_raw'] and ['flux_normalized']
	TO DETERMINE AND COMPUTE FEATURES. CAN BE CALLED BY STAND-ALONE
	PROGRAMS OR AS PART OF BSNIP_retrieve_observables.

	Parameters
	----------
	dataframe		:	Pandas dataframe containing a column of
						wavelengths and a column of fluxes - each row
						entry a spectrum.
	smoothing_mode	:	'savgol' for using savgol-golay filter.
						'wavelet' for using wavelets
						'None' for no smoothing - Not recommended.
					
	Returns
	-------
	run_SIM			:	Pandas dataframe where spectral properties are
						added to each row of the input frame.

	Notes
	-----
	*Parameters for the smoothing techniques are to be hard coded.
	
	Examples
	--------
	TBW.
	
	LOG history
	----------
	TBW.
	
	"""

	def __init__(self, dataframe, smoothing_mode='savgol', smoothing_window=21, deredshift_and_normalize=True,verbose=False):
		self.time														= time.time()
		
		'''
		CHECK SCIPY VERSION ACCEPTS 'param_bounds' (0.18 or higher)
		'''

		#if str(scipy.__version__)[0:4] != '0.18':
		#	sys.exit('\n\n\nError: Make sure scipy version is 0.18 or later, so that param_bounds is accepted to fit curves. Check whether the CONDA path is commented in bash file.\n\n\n')
		
		Utility.__init__(self)
			
		self.DF															= dataframe
		self.smoothing_mode												= smoothing_mode
		self.smoothing_window											= smoothing_window
		self.deredshift_and_normalize									= deredshift_and_normalize
		self.verbose													= verbose
		
		self.MD															= {}	#Table containing the boundaries of line regions. From Silverman+ 2012 (paper II).

		if self.verbose:
			print '\n*STARTING SPECTRAL ANALAYSIS.'
			
		"""
		Silicon lines can be found here:
		https://www.nist.gov/sites/default/files/documents/srd/jpcrd3720081501p.pdf
		"""

		#original
		self.MD['rest_f1'], self.MD['blue_lower_f1'], self.MD['blue_upper_f1'], self.MD['red_lower_f1'], self.MD['red_upper_f1']	= [3945.28], 3400., 3800., 3800., 4100.
		self.MD['rest_f2'], self.MD['blue_lower_f2'], self.MD['blue_upper_f2'], self.MD['red_lower_f2'], self.MD['red_upper_f2']	= [4129.73], 3850., 4000., 4000., 4150.
		self.MD['rest_f3'], self.MD['blue_lower_f3'], self.MD['blue_upper_f3'], self.MD['red_lower_f3'], self.MD['red_upper_f3']	= [4700.], 4000., 4150., 4350., 4700. #rest flux is the upper red bound for uniform selection criteria.  (was None before)
		self.MD['rest_f4'], self.MD['blue_lower_f4'], self.MD['blue_upper_f4'], self.MD['red_lower_f4'], self.MD['red_upper_f4']	= [5550.], 4350., 4700., 5050., 5550. #rest flux is the upper red bound for uniform selection criteria.
		self.MD['rest_f5'], self.MD['blue_lower_f5'], self.MD['blue_upper_f5'], self.MD['red_lower_f5'], self.MD['red_upper_f5']	= [5624.32], 5100., 5300., 5450., 5700.
		self.MD['rest_f6'], self.MD['blue_lower_f6'], self.MD['blue_upper_f6'], self.MD['red_lower_f6'], self.MD['red_upper_f6']	= [5971.85], 5400., 5700., 5750., 6000.
		self.MD['rest_f7'], self.MD['blue_lower_f7'], self.MD['blue_upper_f7'], self.MD['red_lower_f7'], self.MD['red_upper_f7']	= [6355.21], 5750., 6060., 6200., 6600.
		self.MD['rest_f8'], self.MD['blue_lower_f8'], self.MD['blue_upper_f8'], self.MD['red_lower_f8'], self.MD['red_upper_f8']	= [7773.37], 6800., 7450., 7600., 8000.
		self.MD['rest_f9'], self.MD['blue_lower_f9'], self.MD['blue_upper_f9'], self.MD['red_lower_f9'], self.MD['red_upper_f9']	= [8498., 8542., 8662.], 7500., 8100., 8200., 8900.
		
	def deredshift_spectrum(self):
		"""
		Remove cosmological redshift, which seems to not be removed in the downloaded data.
		BSNIP I paper says: "Plots of all of the fully reduced spectra as well as (for the objects with multiband SN and galaxy photometry) galaxy-subtracted spectra (as discussed in Section 3.3) presented in this work are available online ..."
		However it is not clear whether the wavelength is de-redshifted. Looking as objects such as 2006bu and 2006cj and comparing to the observed/rest wavelength of these objects in WISEREP, I strongly suspect that the data is *not* deredshifted.
		"""	
		start_time														= time.time()
		def remove_cosmological_redshift(wavelength,redshift):		
			try:
				wavelength												= np.asarray(wavelength).astype(np.float)
				redshift												= float(redshift)
				wavelength												= wavelength/(1.+redshift)
			except:
				wavelength												= np.full(len(wavelength), np.nan)
			return wavelength
		self.DF['wavelength_raw']										= self.DF.apply(lambda row: pd.Series([remove_cosmological_redshift(row['wavelength_raw'], row['host_redshift'])]), axis=1)	
		if self.verbose:
			print '  -RAN: De-redshifting the spectra. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None

	def normalize_flux_to_max(self):
		"""
		Normalize the flux to relative units so that methods such as
		wavelet smoothing can be applied if requested.
		"""	
		start_time														= time.time()

		def get_normalized_flux(flux):
			aux_flux													= np.asarray(flux).astype(np.float)
			normalization_factor										= max(aux_flux)
			aux_flux													= aux_flux/normalization_factor
			aux_flux													= list(aux_flux)
			return aux_flux
		self.DF['flux_normalized']										= self.DF.apply(lambda row: get_normalized_flux(row['flux_raw']), axis=1)	
		if self.verbose:
			print '  -RAN: Normalizing flux to maximum. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None

	def smooth_spectrum(self):
		"""
		Smooth the spectrum using either the savgol-golay filter or wavelets.
		"""									
		start_time														= time.time()

		def savitzky_golay(y, window_size, order, deriv=0, rate=1):		
			#FROM http://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
			try:
				window_size = np.abs(np.int(window_size))
				order = np.abs(np.int(order))
			except ValueError, msg:
				raise ValueError("window_size and order have to be of type int")
			if window_size % 2 != 1 or window_size < 1:
				raise TypeError("window_size size must be a positive odd number")
			if window_size < order + 2:
				raise TypeError("window_size is too small for the polynomials order")
			order_range = range(order+1)
			half_window = (window_size -1) // 2
			# precompute coefficients
			b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
			m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
			# pad the signal at the extremes with
			# values taken from the signal itself
			firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
			lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
			y = np.concatenate((firstvals, y, lastvals))
			return np.convolve( m[::-1], y, mode='valid')

		def do_smooth(wavelength,flux):
			wavelength, flux											= np.asarray(wavelength).astype(np.float), np.asarray(flux).astype(np.float)
			
			#f_smoothed													= savgol_filter(flux, self.smoothing_window, 3)
			f_smoothed													= savitzky_golay(flux, self.smoothing_window, 3)
			
			
			
			df															= np.asarray([f - f_next for f, f_next in zip(f_smoothed,f_smoothed[1:])])
			dw															= np.asarray([w - w_next for w, w_next in zip(wavelength,wavelength[1:])])						
			#derivative													= savgol_filter(np.asarray([np.nan]+list(df/dw)), self.smoothing_window, 3)			
			derivative													= savitzky_golay(np.asarray([np.nan]+list(df/dw)), self.smoothing_window, 3)			

			return f_smoothed, derivative

		flux_smoothed													= self.DF.apply(lambda row: pd.Series(do_smooth(row['wavelength_raw'],row['flux_normalized'])), axis=1)
		flux_smoothed.columns											= ['flux_smoothed', 'derivative']
		self.DF															= self.DF.join(flux_smoothed)
				
		self.DF['wavelength_smoothed']									= self.DF['wavelength_raw']								
		if self.verbose:
			print '  -RAN: Smoothing flux. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None

	def find_zeros_in_features(self):
		"""
		No feature boundary if no peak in the region.
		If only one peak, then maxima is trivially determined.
		If more than one peak in the line (blue or red) window, then
		usually choose the largest peak, but the proceedure is feature
		dependent. This needs the wavelength arrays to be ordered.
		"""			
		start_time														= time.time()
		def get_zeros(wavelength, flux, derivative, key):
			wavelength, flux, derivative								= np.asarray(wavelength).astype(np.float), np.asarray(flux).astype(np.float), np.asarray(derivative).astype(np.float)
			window_condition											= (wavelength >= self.MD['blue_lower_f'+key]) & (wavelength <= self.MD['red_upper_f'+key])			
			w_window, f_window, der_window								= wavelength[window_condition], flux[window_condition], derivative[window_condition]	
			
			idx_minima_window											= [idx for idx, (w,f,der,der_next) in enumerate(zip(w_window,f_window,der_window,der_window[1:])) if (np.sign(der) != np.sign(der_next)) and (der_next > 0.) and w < max(self.MD['rest_f'+key])]		
			idx_maxima_window											= [idx for idx, (w,f,der,der_next) in enumerate(zip(w_window,f_window,der_window,der_window[1:])) if (np.sign(der) != np.sign(der_next)) and (der_next < 0.)]			
	
			w_minima_window, f_minima_window							= np.asarray([w_window[idx] for idx in idx_minima_window]), np.asarray([f_window[idx] for idx in idx_minima_window])			
			w_maxima_window, f_maxima_window							= np.asarray([w_window[idx] for idx in idx_maxima_window]), np.asarray([f_window[idx] for idx in idx_maxima_window])

			'''
			Find where the true minimum of the feature is. Iterate over the minima to find the deepest deep
			that contains a maximum to the right and to the left.
			'''
			copy_w_minima_window, copy_f_minima_window					= w_minima_window[:], f_minima_window[:]

			for i in range(len(w_minima_window)):
				if len(copy_w_minima_window) > 0:
					
					'''
					Get the deepest minimum.
					'''
					w_min, f_min = copy_w_minima_window[copy_f_minima_window.argmin()], min(copy_f_minima_window)
					
					'''
					trimming minima and maxima in feature window: select only minima/maxima in the left (right) side of the true minimum for the blue (red) window.
					These are bounded by the pre-fixed limits for the window and the position of the true minimum. 
					'''
					min_blue_condition, min_red_condition					= (w_minima_window < w_min), (w_minima_window > w_min)
					max_blue_condition, max_red_condition					= (w_maxima_window < w_min), (w_maxima_window > w_min)

					minima_window_blue_condition							= min_blue_condition & (w_minima_window <= self.MD['blue_upper_f'+key]) & (w_minima_window >= self.MD['blue_lower_f'+key])				
					maxima_window_blue_condition							= max_blue_condition & (w_maxima_window <= self.MD['blue_upper_f'+key]) & (w_maxima_window >= self.MD['blue_lower_f'+key])				
					minima_window_red_condition								= min_red_condition & (w_minima_window <= self.MD['red_upper_f'+key]) & (w_minima_window >= self.MD['red_lower_f'+key])				
					maxima_window_red_condition								= max_red_condition & (w_maxima_window <= self.MD['red_upper_f'+key]) & (w_maxima_window >= self.MD['red_lower_f'+key])				
					
					w_minima_window_blue, f_minima_window_blue				= w_minima_window[minima_window_blue_condition], f_minima_window[minima_window_blue_condition]	
					w_maxima_window_blue, f_maxima_window_blue				= w_maxima_window[maxima_window_blue_condition], f_maxima_window[maxima_window_blue_condition]	
					w_minima_window_red, f_minima_window_red				= w_minima_window[minima_window_red_condition], f_minima_window[minima_window_red_condition]	
					w_maxima_window_red, f_maxima_window_red				= w_maxima_window[maxima_window_red_condition], f_maxima_window[maxima_window_red_condition]	

					'''
					Select the maxima to the right and to the left.
					'''	
					try:
						w_max_blue, f_max_blue								= w_maxima_window_blue[-1], f_maxima_window_blue[-1]
						w_max_red, f_max_red								= w_maxima_window_red[0], f_maxima_window_red[0]
					except:
						w_max_blue, f_max_blue								= np.nan, np.nan
						w_max_red, f_max_red								= np.nan, np.nan							

					'''
					If there is no maximum to either the left or to the right, remove the minimum from the list of minima and
					try the next deepest minimum.
					'''
					if np.isnan(w_max_blue) == False and np.isnan(w_max_red) == False:
						break
					else:
						copy_w_minima_window = np.asarray(filter(lambda x : x != w_min, copy_w_minima_window))
						copy_f_minima_window = np.asarray(filter(lambda x : x != f_min, copy_f_minima_window))	

			if len(copy_w_minima_window) == 0: 
				w_min, f_min, w_max_blue, f_max_blue, w_max_red, f_max_red = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan 	

			'''
			Once the true minimum is known, check whether the nearest maxima are spurious.
			For instance, if the maxima are just shoulders.
			'''

			if np.isnan(w_max_blue) == False and len(w_maxima_window_blue)>1:	
				#Compute distance from blue minima to the blue maximum.
				distance_w_minima_window_blue = w_minima_window_blue - w_max_blue	
				#Select only the minima which are bluer than the maximum and within the separation window.
				distance_w_minima_window_blue = distance_w_minima_window_blue[(distance_w_minima_window_blue < 0.) & (distance_w_minima_window_blue > -1.*self.sep)]
				#Take the bluest of the minima and check if there is another maximum bluer than that.
				if len(distance_w_minima_window_blue)>0:
					condition = (w_maxima_window_blue < w_max_blue + min(distance_w_minima_window_blue))
					w_maxima_window_blue = w_maxima_window_blue[condition]
					f_maxima_window_blue = f_maxima_window_blue[condition]
					if len(w_maxima_window_blue) >= 1:
						w_max_blue, f_max_blue								= w_maxima_window_blue[-1], f_maxima_window_blue[-1]

			if np.isnan(w_max_red) == False and len(w_maxima_window_red)>1:	
				
				#Compute distance from red minima to the red maximum.
				distance_w_minima_window_red = w_minima_window_red - w_max_red	
				#Select only the minima which are redr than the maximum and within the separation window.
				distance_w_minima_window_red = distance_w_minima_window_red[(distance_w_minima_window_red > 0.) & (distance_w_minima_window_red < 1.*self.sep)]
				#Take the redst of the minima and check if there is another maximum redr than that.
				if len(distance_w_minima_window_red)>0:
					condition = (w_maxima_window_red > w_max_red + max(distance_w_minima_window_red))
					w_maxima_window_red = w_maxima_window_red[condition]
					f_maxima_window_red = f_maxima_window_red[condition]
					if len(w_maxima_window_red) >= 1:
						w_max_red, f_max_red								= w_maxima_window_red[0], f_maxima_window_red[0]
										


			#print float(w_min), float(f_min), float(w_max_blue), float(f_max_blue), float(w_max_red), float(f_max_red)
			return float(w_min), float(f_min), float(w_max_blue), float(f_max_blue), float(w_max_red), float(f_max_red)
		
		for key in self.keys:
			#print 'key', key
			feature_zeros												= self.DF.apply(lambda row: pd.Series(get_zeros(row['wavelength_smoothed'],row['flux_smoothed'],row['derivative'],key)), axis=1)
			feature_zeros.columns										= ['wavelength_minima_f'+key, 'flux_minima_f'+key, 'wavelength_maxima_blue_f'+key, 'flux_maxima_blue_f'+key, 'wavelength_maxima_red_f'+key, 'flux_maxima_red_f'+key]
			self.DF														= self.DF.join(feature_zeros)							
		if self.verbose:
			print '  -RAN: Determining zeros that define the pseudo spectrum. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None
		
	def grab_feature_regions(self):
		"""
		Determined here whether to use raw or smoothed spectra for computing features.
		"""
		start_time														= time.time()
	
		def isolate_region(wavelength, flux_normalized, flux_smoothed, blue_boundary, red_boundary):				
			wavelength, flux_normalized, flux_smoothed							= np.asarray(wavelength).astype(np.float), np.asarray(flux_normalized).astype(np.float), np.asarray(flux_smoothed).astype(np.float)
			if not np.isnan(blue_boundary) and not np.isnan(red_boundary): 
				region_condition										= (wavelength >= blue_boundary) & (wavelength <= red_boundary)		
				wavelength_region, flux_normalized_region, flux_smoothed_region= wavelength[region_condition], flux_normalized[region_condition], flux_smoothed[region_condition]			
			else:
				wavelength_region, flux_normalized_region, flux_smoothed_region= [np.nan], [np.nan], [np.nan]	
			return wavelength_region, flux_normalized_region, flux_smoothed_region
		
		for key in self.keys:
			feature_region												= self.DF.apply(lambda row: pd.Series(isolate_region(row['wavelength_raw'],row['flux_normalized'],row['flux_smoothed'],row['wavelength_maxima_blue_f'+key],row['wavelength_maxima_red_f'+key])), axis=1)	
			feature_region.columns										= ['wavelength_region_f'+key, 'flux_normalized_region_f'+key, 'flux_smoothed_region_f'+key]
			self.DF														= self.DF.join(feature_region)		
		
		if self.verbose:
			print '  -RAN: Grabing the region of each feature. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None	

	def make_pseudo_continuum(self):
		"""
		The pseudo continuum slope is simply a line connecting the
		feature region boundaries. The pseudo continuum depends
		only on the wavelength array and boundary values, the latter
		coming from smoothed quantities, regardless of method chosen (raw or smoothed.)
		"""			
		start_time														= time.time()
		def get_psedo_continuum_flux(w,x1,y1,x2,y2,f_smoothed):
			try:	
				w, f_smoothed											= np.asarray(w), np.asarray(f_smoothed)
				slope													= (y2-y1)/(x2-x1)
				intercept												= y1 - slope*x1
				def pseudo_cont(x):
					return slope*x+intercept
				pseudo_flux												= pseudo_cont(w)			

				boolean_check											= [(f_s-f_c)>0.01 for (f_c,f_s) in zip(pseudo_flux, f_smoothed)]
				if True in boolean_check or len(pseudo_flux) < 1:
					pseudo_flux											= 'Failed'

			except:
				pseudo_flux												= 'Failed'
			return pseudo_flux											
		
		for key in self.keys:
			pseudo_cont_flux											= self.DF.apply(lambda row: pd.Series([get_psedo_continuum_flux(row['wavelength_region_f'+key],row['wavelength_maxima_blue_f'+key],row['flux_maxima_blue_f'+key],row['wavelength_maxima_red_f'+key],row['flux_maxima_red_f'+key],row['flux_smoothed_region_f'+key])]), axis=1)	
			pseudo_cont_flux.columns									= ['pseudo_cont_flux_f'+key]
			self.DF														= self.DF.join(pseudo_cont_flux)
		if self.verbose:
			print '  -RAN: Making the pseudo continuum. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None

	def compute_pEW(self):
		start_time														= time.time()
		def get_pEW(wavelength_region, flux_region, pseudo_flux):
			try:
				if len(wavelength_region) > 1:
					pEW													= 0.
					for i, (w, f, p) in enumerate(zip(wavelength_region[0:-1],flux_region[0:-1], pseudo_flux[0:-1])): #Not including last entry, so that differences can be computed.
						delta_lambda									= abs(wavelength_region[i+1] - wavelength_region[i])
						pEW												+= delta_lambda*(p - f)/p
				else:
					pEW													= np.nan									
			except:
				pEW														= np.nan		
			return pEW

		for key in self.keys:
			pEW_value													= self.DF.apply(lambda row: pd.Series(get_pEW(row['wavelength_region_f'+key],row['flux_normalized_region_f'+key],row['pseudo_cont_flux_f'+key])), axis=1)	
			pEW_value.columns											= ['pEW_f'+key]
			self.DF														= self.DF.join(pEW_value)
		if self.verbose:
			print '  -RAN: Computing the pEW. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None	

	def compute_smoothed_velocity_and_depth(self):
		"""
		Smoothed velocity is the velocity obtained from the wavelength
		of the minima in the feature region.
		"""			
		start_time														= time.time()
		def get_smoothed_velocity(wavelength_region,flux_region,pseudo_cont_flux,rest_wavelength):
			try:
				wavelength_region, flux_region, pseudo_cont_flux		= np.asarray(wavelength_region).astype(np.float), np.asarray(flux_region).astype(np.float), np.asarray(pseudo_cont_flux).astype(np.float)			
								
				flux_at_min, wavelength_at_min,  pseudo_cont_at_min		= min(flux_region), wavelength_region[flux_region.argmin()], pseudo_cont_flux[flux_region.argmin()] 
				
				wavelength_par_window									= wavelength_region[(wavelength_region >= wavelength_at_min - self.sep) & (wavelength_region <= wavelength_at_min + self.sep)]
				flux_par_window											= flux_region[(wavelength_region >= wavelength_at_min - self.sep) & (wavelength_region <= wavelength_at_min + self.sep)]
				par_to_fit												= self.make_parabola(wavelength_at_min)						
				popt, covt												= curve_fit(par_to_fit,wavelength_par_window,flux_par_window)
						
				rest_wavelength											= sum(rest_wavelength)/len(rest_wavelength)
				wavelength_par_min										= wavelength_at_min - popt[1]/(2*popt[0])
				flux_par_min											= par_to_fit(wavelength_par_min,popt[0],popt[1],popt[2])		
				velocity												= c*( (wavelength_par_min/rest_wavelength)**2. - 1.  )/( (wavelength_par_min/rest_wavelength)**2. + 1. )/1.e3
				depth													= 1. - flux_par_min/pseudo_cont_at_min
								
				if popt[0] < 0. or velocity > 0. or velocity < -30000.:			
					velocity											= np.nan					
			
			except:					
				wavelength_par_min, flux_par_min, velocity, depth		= np.nan, np.nan, np.nan, np.nan			
			
			return wavelength_par_min, flux_par_min, velocity, depth	
	
		for key in self.keys_to_fit:
			velocity_from_smoothing										= self.DF.apply(lambda row: pd.Series(get_smoothed_velocity(row['wavelength_region_f'+key],row['flux_normalized_region_f'+key],row['pseudo_cont_flux_f'+key],self.MD['rest_f'+key])), axis=1)	
			velocity_from_smoothing.columns								= ['wavelength_at_min_f'+key, 'flux_at_min_f'+key, 'velocity_f'+key, 'depth_f'+key]
			self.DF														= self.DF.join(velocity_from_smoothing)
		if self.verbose:
			print '  -RAN: Computing line velocities from minima in smoothed spectra. FINISHED IN ('+str(format(time.time()-start_time, '.1f'))+'s)'
		return None		
	
	def run_analysis(self):
		if self.deredshift_and_normalize:
			self.deredshift_spectrum()
			self.normalize_flux_to_max()	
		self.smooth_spectrum()	
		self.find_zeros_in_features()
		self.grab_feature_regions()
		self.make_pseudo_continuum()
		self.compute_pEW()
		self.compute_smoothed_velocity_and_depth()
		
		if self.verbose:
			print "    -TOTAL TIME IN SPECTRAL ANALYSIS: "+str(format(time.time()-self.time, '.1f'))+'s'			
			print '    *** RUN COMPLETED SUCCESSFULLY ***\n'						
		return self.DF	

class uncertainty(Utility):
	"""
	Description
	----------
	TAKES A DATASET WITH COLUMNS ['wavelength_raw'] and ['flux_normalized']
	TO DETERMINE AND COMPUTE FEATURES. CAN BE CALLED BY STAND-ALONE
	PROGRAMS OR AS PART OF BSNIP_retrieve_observables.

	Parameters
	----------
	dataframe		:	Pandas dataframe containing a column of
						wavelengths and a column of fluxes - each row
						entry a spectrum.
	smoothing_mode	:	'savgol' for using savgol-golay filter.
						'wavelet' for using wavelets
						'None' for no smoothing - Not recommended.
					
	Returns
	-------
	run_SIM			:	Pandas dataframe where spectral properties are
						added to each row of the input frame.

	Notes
	-----
	*Parameters for the smoothing techniques are to be hard coded.
	
	Examples
	--------
	TBW.
	
	LOG history
	----------
	TBW.
	
	"""
	def __init__(self, dataframe, smoothing_mode='savgol', smoothing_window=21, N_MC_runs=3000):
		self.time														= time.time()
		
		print '\n*STARTING CALCULATION OF UNCERTAINTIES.'
		
		Utility.__init__(self)
			
		self.df															= dataframe
		self.smoothing_mode												= smoothing_mode
		self.smoothing_window											= smoothing_window
		self.N_MC_runs													= N_MC_runs

		if smoothing_window == 21:
			self.smoothing_correction									= (1./0.93)
			#self.smoothing_correction									= (1./0.93)*1.47
		elif smoothing_window == 51:
			self.smoothing_correction									= 1./0.96	
		else:
			self.smoothing_correction									= 1.
			#raise ValueError("Smoothing correction not defined for this smoothing window.")

	def compute_flux_rms(self):
		"""
		Compute the flux uncertainty in each pixel using a simple rms
		in a bin defined by the self.sep parameter.
		"""									
		print '  -RUNNING: Computing flux rms pixel-wise.',		
		start_time														= time.time()
		
		def get_rms(wavelength, flux_normalized, flux_smoothed):
			try:
				wavelength, flux_normalized, flux_smoothed						= np.asarray(wavelength).astype(np.float), np.asarray(flux_normalized).astype(np.float), np.asarray(flux_smoothed).astype(np.float)		
				flux_normalized_bins											= [flux_normalized[(wavelength >= w - self.sep) & (wavelength <= w + self.sep)] for w in wavelength]
				flux_smoothed_bins										= [flux_smoothed[(wavelength >= w - self.sep) & (wavelength <= w + self.sep)] for w in wavelength]		
				rms														= [self.rms(f_raw_bins, f_smoothed_bins)*self.smoothing_correction for (f_raw_bins, f_smoothed_bins) in zip(flux_normalized_bins,flux_smoothed_bins)]
			except:
				rms														= np.nan
		
			print 'Median of the smoothing corrected rms = ', np.median(rms)
			print 'Medan of the smoothing corrected rms = ', np.mean(rms)
			return rms
	
		flux_rms														= self.df.apply(lambda row: pd.Series([get_rms(row['wavelength_raw'],row['flux_normalized'],row['flux_smoothed'])]), axis=1)	
		flux_rms.columns												= ['flux_rms']
		self.df															= self.df.join(flux_rms)		
	
		print 'DONE ('+str(format(time.time()-start_time, '.1f'))+'s)'			
		return None

	def compute_mock_spectra(self, index):
		"""
		Create N_runs mock spectra for each spectrum in the dataframe.
		The are to be used to compute uncertainties in the pEW and
		velocity measurements using a MC technique.
		"""		
		try:
			wavelength_raw, flux_normalized, flux_rms							= np.asarray(self.df.loc[index]['wavelength_raw']).astype(np.float), np.asarray(self.df.loc[index]['flux_normalized']).astype(np.float), np.asarray(self.df.loc[index]['flux_rms']).astype(np.float)		
			random_flux_draw											= [[np.random.normal(flux,noise) for (flux,noise) in zip(flux_normalized,flux_rms)] for i in range(self.N_MC_runs)]
			wavelength													= [wavelength_raw for i in range(self.N_MC_runs)]		
			mock_dict													= pd.DataFrame({'wavelength_raw': wavelength, 'flux_normalized': random_flux_draw})
		
		except:
			mock_dict													= {}		
		return mock_dict

	def compute_uncertainty(self,quantity,quantity_value,bin_size):
		def gaussian(x,A,mu,sigma):
			return A*np.exp(-(x-mu)**2./(2.*sigma**2.))

		if not np.isnan(quantity).all():
			flag														= False
			
			quantity													= np.asarray(quantity)
			quantity													= quantity[~np.isnan(quantity)]		
			quantity													= quantity[quantity != 0]		
			quantity_median												= np.median(quantity)	
			quantity_mean												= quantity.mean()

			#print quantity_value, quantity_mean, quantity_median

			bin_edges													= np.arange(math.floor(min(quantity)), math.ceil(max(quantity))+bin_size, bin_size)
			center_bins													= np.asarray([edge+(edge_next - edge)/2. for edge, edge_next in zip(bin_edges,bin_edges[1:])])	
			pEW_histogram, edges										= np.histogram(quantity,bins=bin_edges)
			#plt.plot(center_bins,pEW_histogram)
			
			try:				
				popt, pcov												= curve_fit(gaussian,center_bins,pEW_histogram,p0=[self.N_MC_runs/6.,quantity_median,abs(quantity_median/5.)])
				gaussian_mean, unc										= popt[1], abs(popt[2])
							
				#print quantity_value, quantity_mean, quantity_median, gaussian_mean
				#plt.plot(center_bins,gaussian(center_bins,popt[0], popt[1], popt[2]))


				if not np.isnan(quantity_value):
					#if abs(quantity_value - gaussian_mean) > 1.*unc: #But distribution is not necessarily symmetric
					if abs(quantity_value - quantity_median) > max([unc,bin_size]): #But distribution is not necessarily symmetric
					#if abs(quantity_value - quantity_mean) > 1.*unc: #But surious values comromise the mean
						#print 'Value far from gaussian mean.'
						flag											= True	
				elif np.isnan(quantity_value) and not np.isnan(quantity_median):  	
					flag												= True
					#print 'value from mock spectra.'
					quantity_value										= quantity_median
						
			except:
				#print 'Gaussian fit failed.'
				flag													= True
				unc														= 2.*abs(quantity_value - quantity_mean)	
			
							

			#if unc > quantity_value:
			#	flag													= True
			#	unc														= 2.*abs(quantity_value - quantity_mean)	
				


		else:
			unc, flag													= np.nan, True	

		#print 'Final', quantity_value, unc, flag
		#plt.show()
		
		return unc, flag, quantity_value

	def run_uncertainties(self):
		self.compute_flux_rms()			
				
		for idx in self.df.index.values:
			print '  -PROCESSING OBJECT WITH INDEX: '+str(idx)		
			mock_spectra_dict											= self.compute_mock_spectra(idx)

			if any(mock_spectra_dict):
				mock_spectra_dict										= analyse_spectra(mock_spectra_dict, smoothing_mode='savgol', smoothing_window=self.smoothing_window, deredshift_and_normalize=False, verbose=False).run_analysis()	
				
				
				for key in self.keys:
					pEW_unc, pEW_flag, pEW_value						= self.compute_uncertainty(mock_spectra_dict['pEW_f'+key].tolist(),self.df.loc[idx]['pEW_f'+key], 0.5)				
					#print 'pew', key, pEW_value, pEW_unc, pEW_flag
					self.df.loc[idx, 'pEW_unc_f'+key], self.df.loc[idx, 'pEW_flag_f'+key], self.df.loc[idx, 'pEW_f'+key] = pEW_unc, pEW_flag, pEW_value
					pass
				
				for key in self.keys_to_fit:				
					velocity_unc, velocity_flag, velocity_value			= self.compute_uncertainty(mock_spectra_dict['velocity_f'+key].tolist(),self.df.loc[idx]['velocity_f'+key], 0.1)				
					#print key, 'vel', key, velocity_value, velocity_unc, velocity_flag
					self.df.loc[idx, 'velocity_unc_f'+key], self.df.loc[idx, 'velocity_flag_f'+key], self.df.loc[idx, 'velocity_f'+key] = velocity_unc, velocity_flag, velocity_value
					
					depth_unc, depth_flag, depth_value					= self.compute_uncertainty(mock_spectra_dict['depth_f'+key].tolist(),self.df.loc[idx]['depth_f'+key], 0.01)				
					#print 'dep', key, depth_value, depth_unc, depth_flag
					self.df.loc[idx, 'depth_unc_f'+key], self.df.loc[idx, 'depth_flag_f'+key], self.df.loc[idx, 'depth_f'+key] = depth_unc, depth_flag, depth_value
					pass
			
			else:
				for key in self.keys:
					#print key
					self.df.loc[idx, 'pEW_unc_f'+key], self.df.loc[idx, 'pEW_flag_f'+key], self.df.loc[idx, 'pEW_f'+key] = np.nan, np.nan, np.nan
				for key in self.keys_to_fit:				
					self.df.loc[idx, 'velocity_unc_f'+key], self.df.loc[idx, 'velocity_flag_f'+key], self.df.loc[idx, 'velocity_f'+key] = np.nan, np.nan, np.nan
					self.df.loc[idx, 'depth_unc_f'+key], self.df.loc[idx, 'depth_flag_f'+key], self.df.loc[idx, 'depth_f'+key] = np.nan, np.nan, np.nan
			
			del mock_spectra_dict
							
		print "    -TOTAL TIME IN COMPUTING UNCERTAINTIES: "+str(format(time.time()-self.time, '.1f'))+'s'			
		print '    *** RUN COMPLETED SUCCESSFULLY ***\n'

		return self.df	



