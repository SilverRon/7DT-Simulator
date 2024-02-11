#============================================================
#	https://ps1images.stsci.edu/ps1image.html
#	By Sophia Kim 2019.01.22. based on code PS1 suggests on the link above
#	Pan-STARRS DR1 data query
#	from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
#	By CS Choi 
#	REVISED AND ORGANIZED BY GREGORY S.H. PAEK
#	UPDATE : 20.01.03
#============================================================
from __future__ import print_function

from astropy.io import fits
from astropy.io import ascii
import glob, os
# from imsng import tool
# from imsng import phot_tbd
import sys
sys.path.append('..')
from util import tool
from phot import gpphot
import numpy as np
from astropy.coordinates import Angle
import astropy.units as u
import astropy.coordinates as coord
from astropy.table import Table
import pandas as pd
#
from gaiaxpy import generate, PhotometricSystem
#
#-------------------------------------------------------------------------#
#
#	7DT
#
#-------------------------------------------------------------------------#
# def querybox(refcatname, obj, racent, decent, path_refcat, radius=0.5, refmagkey=''):
#-------------------------------------------------------------------------#
def calculate_snr(flux, flux_error):
    snr = flux / flux_error
    return snr
#-------------------------------------------------------------------------#
def calculate_magnitude_error(flux, flux_error):
    snr = calculate_snr(flux, flux_error)
    magnitude_error = 2.5*np.log10(1+(1/snr))    
    return magnitude_error
#-------------------------------------------------------------------------#
from gaiaxpy import generate, PhotometricSystem
from gaiaxpy import PhotometricSystem, load_additional_systems
from astroquery.gaia import GaiaClass
import pandas as pd
# path_to_filterset = '../config/filterset'
# PhotometricSystem = load_additional_systems(path_to_filterset)
# PhotometricSystem.get_available_systems().split(', ')[-2:]
#
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
#	Keys to Query
#-------------------------------------------------------------------------#
columns_to_query = [
	'gaia_source.solution_id',
	'gaia_source.source_id',
	'gaia_source.ra',
	'gaia_source.dec',
	'gaia_source.parallax',
	'gaia_source.l',
	'gaia_source.b',

	'gaia_source.phot_g_mean_mag',
	'gaia_source.bp_rp',
	'gaia_source.bp_g',
	'gaia_source.g_rp',

	'gaia_source.phot_variable_flag',
	'gaia_source.in_galaxy_candidates',
	'gaia_source.non_single_star',
	'gaia_source.has_xp_continuous',

	'gaia_source.has_rvs',
	# 'gaia_source.has_epoch_photometry',
	'gaia_source.ebpminrp_gspphot',
	'gaia_source.phot_bp_rp_excess_factor',
	'gaia_source.ruwe',
	'gaia_source.ipd_frac_multi_peak',
	'gaia_source.ipd_frac_odd_win',
	]

str_to_query = ', '.join(columns_to_query)
# print(f"{len(columns_to_query)} columns to query:")
# print(str_to_query)

#-------------------------------------------------------------------------#
def calculate_snr(flux, flux_error):
    snr = flux / flux_error
    return snr
#-------------------------------------------------------------------------#
def calculate_magnitude_error(flux, flux_error):
    snr = calculate_snr(flux, flux_error)
    magnitude_error = 2.5*np.log10(1+(1/snr))    
    return magnitude_error
#-------------------------------------------------------------------------#
path_to_filterset = '/home/gp/gppy/config/filterset'
PhotometricSystem = load_additional_systems(path_to_filterset)

path_save = '/large_data/factory/SPSS'
metadf = pd.read_csv(f'{path_save}/Gaia_SPSS.csv')
# fname = f'{path_save}/synphot.csv'
fname = f'{path_save}/XP_CONTINUOUS_RAW.csv'
phot_system_7dt = PhotometricSystem.USER_7DT_Edmund
synthetic_photometry_7dt = generate(fname, photometric_system=phot_system_7dt, save_file=False)


#============================================================
#	Synthetic Photometry
#------------------------------------------------------------
#	SDSS
#------------------------------------------------------------
phot_system_sdss = PhotometricSystem.SDSS
synthetic_photometry_sdss = generate(fname, photometric_system=phot_system_sdss, save_file=False)
synthetic_photometry_sdss
#------------------------------------------------------------
#	7DT
#------------------------------------------------------------
phot_system_7dt = PhotometricSystem.USER_7DT_Edmund
synthetic_photometry_7dt = generate(fname, photometric_system=phot_system_7dt, save_file=False)
synthetic_photometry_7dt
#============================================================
#	Merge the Table
#------------------------------------------------------------
#	Meta Table + SDSS + 7DT
#------------------------------------------------------------
synthetic_photometry = pd.merge(synthetic_photometry_sdss, synthetic_photometry_7dt, on='source_id')
merged_df = pd.merge(synthetic_photometry, metadf, on='source_id')
#------------------------------------------------------------
#	Filter List
#------------------------------------------------------------
filterlist_sdss = ['u', 'g', 'r', 'i', 'z',]
filterlist_7dt = []
for key in merged_df.keys():
	if 'USER_7DT_Edmund_mag_' in key:
		filte = key.replace('USER_7DT_Edmund_mag_', '')
		filterlist_7dt.append(filte)
#------------------------------------------------------------
#	Rename configurations
#------------------------------------------------------------
#	SDSS
prefix_mag_sdss = 'Sdss_mag_'
prefix_flux_sdss = 'Sdss_flux_'
prefix_fluxerr_sdss = 'Sdss_flux_error_'
#	7DT
prefix_mag_7dt = 'USER_7DT_Edmund_mag_'
prefix_flux_7dt = 'USER_7DT_Edmund_flux_'
prefix_fluxerr_7dt = 'USER_7DT_Edmund_flux_error_'
#------------------------------------------------------------
#	Rename the columns
#------------------------------------------------------------
key_to_rename_dict = {}
#	SDSS
for filte in filterlist_sdss:
	#	New Values
	# _mag = merged_df[f"{prefix_mag_sdss}{filte}"]
	_flux = merged_df[f"{prefix_flux_sdss}{filte}"]
	_fluxerr = merged_df[f"{prefix_fluxerr_sdss}{filte}"]
	_snr = calculate_snr(_flux, _fluxerr)
	_magerr = calculate_magnitude_error(_flux, _fluxerr)

	#	Former Keywords
	magkey = f"{prefix_mag_sdss}{filte}"
	fluxkey = f"{prefix_flux_sdss}{filte}"
	fluxerrkey = f"{prefix_fluxerr_sdss}{filte}"

	#	New Keywords
	newmagkey = f"{filte}_mag"
	newmagerrkey = f"{filte}_magerr"
	newfluxkey = f"{filte}_flux"
	newfluxerrkey = f"{filte}_fluxerr"
	newsnrkey = f"{filte}_snr"

	#	Keywords to Rename
	key_to_rename_dict[magkey] = newmagkey
	key_to_rename_dict[fluxkey] = newfluxkey
	key_to_rename_dict[fluxerrkey] = newfluxerrkey

	merged_df[newmagerrkey] = _magerr
	merged_df[newsnrkey] = _snr

new_columns = {}

# 7DT
for filte in filterlist_7dt:
	# New Values
	_flux = merged_df[f"{prefix_flux_7dt}{filte}"]
	_fluxerr = merged_df[f"{prefix_fluxerr_7dt}{filte}"]
	_snr = calculate_snr(_flux, _fluxerr)
	_magerr = calculate_magnitude_error(_flux, _fluxerr)

	# Former Keywords
	magkey = f"{prefix_mag_7dt}{filte}"
	fluxkey = f"{prefix_flux_7dt}{filte}"
	fluxerrkey = f"{prefix_fluxerr_7dt}{filte}"

	# New Keywords
	if '_50' in filte:
		replace_key = '_50'
		suffix = 'w'
	elif '_25' in filte:
		replace_key = '_25'
		suffix = ''
	else:
		replace_key = ''
		suffix = ''

	newmagkey = f"m{filte.replace(replace_key, suffix)}_mag"
	newmagerrkey = f"m{filte.replace(replace_key, suffix)}_magerr"
	newfluxkey = f"m{filte.replace(replace_key, suffix)}_flux"
	newfluxerrkey = f"m{filte.replace(replace_key, suffix)}_fluxerr"
	newsnrkey = f"m{filte.replace(replace_key, suffix)}_snr"

	# Keywords to Rename
	key_to_rename_dict[magkey] = newmagkey
	key_to_rename_dict[fluxkey] = newfluxkey
	key_to_rename_dict[fluxerrkey] = newfluxerrkey

	new_columns[newmagerrkey] = _magerr
	new_columns[newsnrkey] = _snr

# 새로운 열을 한 번에 추가
merged_df = pd.concat([merged_df, pd.DataFrame(new_columns)], axis=1)

#------------------------------------------------------------
#	Final Dataframe
#------------------------------------------------------------
merged_df.rename(columns=key_to_rename_dict, inplace=True)

merged_df.to_csv(f"{path_save}/gaiaxp_dr3_synphot_SPSS.csv")