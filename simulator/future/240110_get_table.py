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
obj = "LTT9491"
# obj = "LTT7987"
date = '20231015'

##	Path & Pattern
path_processed_data = "/large_data/processed"
path_save = f"{path_processed_data}/{obj}"
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

table.write(f"{path_save}/summary.csv", format='csv')