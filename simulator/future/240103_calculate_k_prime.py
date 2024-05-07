#	Run at the proton server

#	Library
import os
import glob
import numpy as np
from ccdproc import ImageFileCollection
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.optimize import curve_fit
##	Plot presetting
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 20
plt.rcParams['savefig.dpi'] = 500
plt.rc('font', family='serif')
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# 에러바를 고려한 가중 선형 회귀 함수 정의
def linear_model(x, a, b):
	return a * x + b
#	Input
##	Target
# obj = "LTT9491"
# radec = '23 19 35.388 -17 05 28.47'
# date = '20231015'

obj = "LTT7987"
radec = '20 10 56.849 -30 13 06.63'
date = '20231015'
c = SkyCoord(radec, unit=(u.hourangle, u.deg))

##	Path & Pattern
path_processed_data = "/large_data/processed"
path_refcat = f"/large_data/factory/ref_cat"
suffix_image = 'com.fits'
suffix_photcat = f"com.phot.cat"

#	Data List
path_image = f"{path_processed_data}/{obj}/7DT??/*/calib*{date}*{suffix_image}"
path_photcat = f"{path_processed_data}/{obj}/7DT??/*/calib*{date}*{suffix_photcat}"
data_list = sorted(glob.glob(path_image))
cat_list = sorted(glob.glob(path_photcat))
print(f"Data List: {len(data_list)}")
print(f"Cat List: {len(cat_list)}")
refcat = f"{path_refcat}/gaiaxp_dr3_synphot_{obj}.csv"
if os.path.exists(refcat):
	reftbl = Table.read(refcat)
	c_ref = SkyCoord(reftbl['ra'], reftbl['dec'], unit="deg")
	print(f"Reference catalog for {obj} field found")
else:
	print(f"Reference catalog for {obj} field not found!")

#	Summarized Table
ic_all = ImageFileCollection(filenames=data_list, keywords='*',)
table = ic_all.summary
table['cat'] = cat_list

#	Collect Data
# nn = 0

def calc_k_prime(table, filte,):
	mag_obs_arr = []
	magerr_obs_arr = []
	mag_ref_arr = []
	airmass_arr = []

	_table = table[table['filter']==filte]
	for nn in range(len(_table)):
		incat = _table['cat'][nn]
		# filte = _table['filter'][nn]
		airmass = _table['airmass'][nn]
		intbl = Table.read(incat, format='ascii')
		c_cat = SkyCoord(intbl['ALPHA_J2000'], intbl['DELTA_J2000'], unit="deg")
		##	Matching
		indx_in, sep_in, _ = c.match_to_catalog_sky(c_cat)
		indx_ref, sep_ref, _ = c.match_to_catalog_sky(c_ref)

		mintbl = Table(intbl[indx_in])
		mreftbl = Table(reftbl[indx_ref])

		mag_obs = mintbl['MAG_AUTO'].item()
		magerr_obs = mintbl['MAGERR_AUTO'].item()
		mag_ref = mreftbl[f'{filte}_mag'].item()

		mag_obs_arr.append(mag_obs)
		magerr_obs_arr.append(magerr_obs)
		mag_ref_arr.append(mag_ref)
		airmass_arr.append(airmass)

	mag_obs_arr = np.array(mag_obs_arr)
	magerr_obs_arr = np.array(magerr_obs_arr)
	mag_ref_arr = np.array(mag_ref_arr)
	airmass_arr = np.array(airmass_arr)


	delta_mag = mag_obs_arr - mag_ref_arr

	# 가중치는 오차의 역수로 설정 (큰 오차를 가진 데이터 포인트의 가중치를 줄임)
	weights = 1 / magerr_obs_arr

	# 가중 선형 회귀 수행
	popt, pcov = curve_fit(linear_model, airmass_arr, delta_mag, sigma=magerr_obs_arr)

	# 회귀 결과
	slope_weighted = popt[0]
	intercept_weighted = popt[1]

	slope_error = np.sqrt(pcov[0,0])
	intercept_error = np.sqrt(pcov[1,1])

	outdict = {
		'mag_obs_arr': mag_obs_arr,
		'magerr_obs_arr': magerr_obs_arr,
		'mag_ref_arr': mag_ref_arr,
		'airmass_arr': airmass_arr,
		'k_p': slope_weighted,
		'k_p_err': slope_error,
		'Z': intercept_weighted, 
		'Z_err': intercept_error,
	}
	return outdict

def get_color(filte0, filte1):
	indx_ref, sep_ref, _ = c.match_to_catalog_sky(c_ref)
	mreftbl = Table(reftbl[indx_ref])
	mag_ref0 = mreftbl[f'{filte0}_mag'].item()
	mag_ref1 = mreftbl[f'{filte1}_mag'].item()

	color = mag_ref0 - mag_ref1

	return color



#

result_dict = {
	'u': calc_k_prime(table, filte='u'),
	'g': calc_k_prime(table, filte='g'),
	'r': calc_k_prime(table, filte='r'),
	'i': calc_k_prime(table, filte='i'),
	'z': calc_k_prime(table, filte='z'),
}

full_airmass_arr = np.arange(0.5, 2.5+0.1, 0.1)

# 플롯으로 결과 시각화
for ff, filte in enumerate(list(result_dict.keys())):
	_dict = result_dict[filte]
	airmass_arr = _dict['airmass_arr']
	delta_mag = _dict['mag_obs_arr'] - _dict['mag_ref_arr']
	magerr_obs_arr = _dict['magerr_obs_arr']
	Z = _dict['Z']
	k_p = _dict['k_p']

	plt.errorbar(airmass_arr, delta_mag, yerr=magerr_obs_arr, ls='none', color=default_colors[ff], ecolor='k', elinewidth=3, capsize=0)#, label='Data with error bars')
	plt.plot(airmass_arr, delta_mag, marker='s', mfc='none', mec=default_colors[ff], ls='none')
	plt.plot(full_airmass_arr, Z + k_p * full_airmass_arr, color=default_colors[ff], label=r"$\Delta$"+f'{filte}={Z:.3f}+{k_p:.3f}*X', alpha=0.5)

plt.xlabel('Airmass X')
plt.ylabel(r'$\Delta m$')
plt.xlim(0.5, 2.5)
# plt.ylim(np.mean(delta_mag)-0.1, np.mean(delta_mag)+0.1)
plt.title(f'{obj} ({radec})')
plt.legend(loc='upper center', ncol=2)
yl, yu = plt.ylim()
plt.ylim([yl, yu+1])
plt.tight_layout()
plt.show()



# color = get_color(filte0='g', filte1='r')

# filte = 'g'
# _dict = result_dict[filte]
# airmass_arr = _dict['airmass_arr']
# delta_mag = _dict['mag_obs_arr'] - _dict['mag_ref_arr']
# magerr_obs_arr = _dict['magerr_obs_arr']
# k_p = _dict['k_p']

# #	Fitting
# popt, pcov = curve_fit(linear_model, airmass_arr, delta_mag, sigma=magerr_obs_arr)

# # 회귀 결과
# slope_weighted = popt[0]
# intercept_weighted = popt[1]

# slope_error = np.sqrt(pcov[0,0])
# intercept_error = np.sqrt(pcov[1,1])

