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
from astropy.table import hstack
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
# for ff, filte in enumerate(list(result_dict.keys())):
# 	_dict = result_dict[filte]
# 	airmass_arr = _dict['airmass_arr']
# 	delta_mag = _dict['mag_obs_arr'] - _dict['mag_ref_arr']
# 	magerr_obs_arr = _dict['magerr_obs_arr']
# 	Z = _dict['Z']
# 	k_p = _dict['k_p']

# 	plt.errorbar(airmass_arr, delta_mag, yerr=magerr_obs_arr, ls='none', color=default_colors[ff], ecolor='k', elinewidth=3, capsize=0)#, label='Data with error bars')
# 	plt.plot(airmass_arr, delta_mag, marker='s', mfc='none', mec=default_colors[ff], ls='none')
# 	plt.plot(full_airmass_arr, Z + k_p * full_airmass_arr, color=default_colors[ff], label=r"$\Delta$"+f'{filte}={Z:.3f}+{k_p:.3f}*X', alpha=0.5)

# plt.xlabel('Airmass X')
# plt.ylabel(r'$\Delta m$')
# plt.xlim(0.5, 2.5)
# # plt.ylim(np.mean(delta_mag)-0.1, np.mean(delta_mag)+0.1)
# plt.title(f'{obj} ({radec})')
# plt.legend(loc='upper center', ncol=2)
# yl, yu = plt.ylim()
# plt.ylim([yl, yu+1])
# plt.tight_layout()
# plt.show()



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


#	Seperation Cut
sep_cut = 1.0

filte = 'g'
_table = table[table['filter']==filte]
nn = 0
# for nn in range(len(_table)):
	# incat = _table['cat'][nn]


matched_table_list = []
airmass_arr = []
for nn in range(len(_table)):
	incat = _table['cat'][nn]
	intbl = Table.read(incat, format='ascii')
	c_cat = SkyCoord(intbl['ALPHA_J2000'], intbl['DELTA_J2000'], unit="deg")
	indx_match, sep_match, _ = c_ref.match_to_catalog_sky(c_cat)
	airmass = _table['airmass'][nn]
	airmass_arr.append(airmass)


	mtbl = hstack([reftbl, intbl[indx_match]])[sep_match.arcsec < sep_cut]
	# delta_mag = mtbl[f"MAG_AUTO_{filte}"] - mtbl[f"{filte}_mag"]
	print(f"{len(mtbl)} sources matched")

	matched_table_list.append(mtbl)

airmass_arr = np.array(airmass_arr)

# 모든 테이블에서 겹치는 좌표의 인덱스를 찾기 위한 초기 설정
coords_list = [SkyCoord(ra=table['ALPHA_J2000'], dec=table['DELTA_J2000'], unit="deg") for table in matched_table_list]
base_coords = coords_list[0]  # 첫 번째 테이블을 기준으로 설정

# 기준 좌표와 나머지 테이블들의 좌표를 비교
common_indices = np.arange(len(base_coords))  # 초기에는 모든 인덱스를 포함
for coords in coords_list[1:]:
    idx, sep, _ = base_coords.match_to_catalog_sky(coords)
    # 일정 거리 이내에서 매칭되는 인덱스만 유지
    matched_idx = idx[sep.arcsec < sep_cut]
    # 교집합을 찾음
    common_indices = np.intersect1d(common_indices, matched_idx)

# 교집합에 해당하는 좌표만을 각 테이블에서 선택
filtered_tables = []
for mtable in matched_table_list:
	filtered_table = mtable[common_indices]
	filtered_tables.append(filtered_table)
	n_matched_source = len(filtered_table)
	print(f"{n_matched_source} sources filtered")

filtered_table0 = filtered_tables[0]

std_star_table = filtered_table0
# std_star_table['n'] = np.arange(n_matched_source)
std_star_table['g-r'] = [filtered_table0['g_mag'][ii] - filtered_table0['r_mag'][ii] for ii in range(len(filtered_table0))]
std_star_table['k_p'] = 0.0
std_star_table['k_p_err'] = 0.0
# std_star_table['FLAGS'] = filtered_table0['FLAGS']
# std_star_table['CLASS_STAR'] = filtered_table0['CLASS_STAR']
# std_star_table['C_term'] = filtered_table0['C_term']

for nn in range(n_matched_source):
	print(f"[{nn+1}/{n_matched_source} ({(nn+1)/n_matched_source:.1%})]     ", end='\r')

	mag_obs_arr = []
	magerr_obs_arr = []
	mag_ref_arr = []

	_table = table[table['filter']==filte]
	for intbl in filtered_tables:

		mag_obs = intbl['MAG_AUTO'][nn]
		magerr_obs = intbl['MAGERR_AUTO'][nn]
		mag_ref = intbl[f'{filte}_mag'][nn]

		mag_obs_arr.append(mag_obs)
		magerr_obs_arr.append(magerr_obs)
		mag_ref_arr.append(mag_ref)

	mag_obs_arr = np.array(mag_obs_arr)
	magerr_obs_arr = np.array(magerr_obs_arr)
	mag_ref_arr = np.array(mag_ref_arr)

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

	std_star_table['k_p'][nn] = slope_weighted
	std_star_table['k_p_err'][nn] = slope_error


#	
	


nn = 0
# print(f"[{nn+1}/{n_matched_source} ({(nn+1)/n_matched_source:.1%})]     ", end='\r')


_table = table[table['filter']==filte]

delta_mag = filtered_table0['MAG_AUTO'] - filtered_table0[f'{filte}_mag']
magerr_obs_arr = filtered_table0['MAGERR_AUTO']
plt.close()
plt.errorbar(std_star_table['g-r'], delta_mag, yerr=magerr_obs_arr, alpha=0.5, ls='none', marker='.')
plt.errorbar(std_star_table['g-r'][std_star_table['FLAGS']!=0], delta_mag[std_star_table['FLAGS']!=0], yerr=magerr_obs_arr[std_star_table['FLAGS']!=0], alpha=0.5, ls='none', color='r', marker='x')
plt.errorbar(std_star_table['g-r'][std_star_table['CLASS_STAR']<0.9], delta_mag[std_star_table['CLASS_STAR']<0.9], yerr=magerr_obs_arr[std_star_table['CLASS_STAR']<0.9], alpha=0.5, ls='none', color='g', marker='+')
plt.errorbar(std_star_table['g-r'][std_star_table['C_term']>2], delta_mag[std_star_table['C_term']>2], yerr=magerr_obs_arr[std_star_table['C_term']>2], alpha=0.5, ls='none', color='purple', marker='*')
plt.xlabel('(g-r) (mag)')
plt.ylabel(r'$M_{obs}-M_{std}$ (mag)')
plt.xlim(-0.5, +1.75)
plt.tight_layout()
plt.show()
# 가중치는 오차의 역수로 설정 (큰 오차를 가진 데이터 포인트의 가중치를 줄임)
# weights = 1 / magerr_obs_arr

# filtered_table0['FLAGS']
#
# plt.plot(std_star_table['g-r'], std_star_table[''])

indx_filter = np.where(
	(std_star_table['FLAGS']==0) &
	(std_star_table['CLASS_STAR']>0.95) &
	(std_star_table['C_term']<2) &
	(std_star_table['g_mag']>10) &
	(~np.isnan(std_star_table['k_p_err']))
)


# 회귀 분석을 위한 데이터 준비
X = airmass_arr[nn]  # 대기 질량 데이터
k_prime = std_star_table['k_p'][indx_filter]  # k' 데이터
color = std_star_table['g-r'][indx_filter]  # 색상 데이터
delta_m = delta_mag[indx_filter]  # 관측된 별의 크기 차이


# 수정된 회귀 모델 함수 정의
# def func_color_term(color, color_term, Z, k_prime):
#     return color_term * color + k_prime * X + Z
def func_color_term(X_input, color_term, Z):
	color, k_prime = X_input
	return Z + color_term * color + k_prime * X

# def func_color_term_2nd(X_input, color_term, Z, k_2prime):
# 	color, k_prime = X_input
# 	return  Z + color_term * color + k_prime * X + k_2prime * color * X

def func_color_term_2nd(X_input, color_term, Z, k_2prime):
	color, k_prime = X_input
	# color = X_input
	return Z + color_term * color + k_prime * X + k_2prime * color * X

# # curve_fit을 사용하여 회귀 분석 수행
# popt, pcov = curve_fit(func_color_term_, color, delta_m,)

# # 최적화된 매개변수 추출
# color_term, Z, k_prime = popt
# print("color_term:", color_term)
# print("Z:", Z)
# print("k_prime:", k_prime)


# curve_fit을 사용하여 회귀 분석 수행
popt, pcov = curve_fit(func_color_term, (color, k_prime), delta_m,)

x_color = np.arange(-2.0, +2.0+0.01, 0.01)

# 최적화된 매개변수 추출
color_term, Z = popt
print("color_term:", color_term)
print("Z:", Z)


# curve_fit을 사용하여 회귀 분석 수행
# popt, pcov = curve_fit(func_color_term_2nd, (color, k_prime), delta_m,)

# x_color = np.arange(-2.0, +2.0+0.01, 0.01)

# # 최적화된 매개변수 추출
# color_term, Z, k_2prime = popt
# print("color_term:", color_term)
# print("Z:", Z)
# print("k_2prime:", k_2prime)






# plt.close()
# plt.errorbar(std_star_table['g-r'][indx_filter], delta_mag[indx_filter], yerr=magerr_obs_arr[indx_filter], c='tomato', alpha=0.5, ls='none', marker='.')
# # plt.plot(np.sort(color), func_color_term((np.sort(color), k_prime), *popt), '-')
# plt.xlabel('(g-r) (mag)')
# plt.ylabel(r'$M_{obs}-M_{std}$ (mag)')
# plt.xlim(-0.5, +1.75)
# plt.ylim([np.median(delta_mag)-0.5, np.median(delta_mag)+0.5])
# plt.tight_layout()
# plt.show()




# color와 k_prime의 평균값을 사용하여 회귀선 계산
# 이 때, color는 x_color 범위를 사용하고, k_prime은 평균값을 사용합니다.
k_prime_median = np.median(k_prime)
regression_line = func_color_term((x_color, k_prime_median), color_term, Z)

# 오차 막대 그래프 그리기
plt.close()
plt.errorbar(std_star_table['g-r'][indx_filter], delta_mag[indx_filter], yerr=magerr_obs_arr[indx_filter], c='tomato', alpha=0.5, ls='none', marker='.')

# 회귀선 그리기
plt.plot(x_color, regression_line, color='blue', label=f"C={color_term:.3f}\nZ'={Z:.3f}")

# 그래프 레이블 및 제한 설정
plt.xlabel('(g-r) (mag)')
plt.ylabel(r'$M_{obs}-M_{std}$ (mag)')
plt.xlim(-0.5, +1.75)
plt.ylim([np.median(delta_mag)-0.5, np.median(delta_mag)+0.5])

# 범례 추가
plt.legend()

# 그래프 보여주기
plt.tight_layout()
plt.show()
