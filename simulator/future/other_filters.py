import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.table import Table
from . import utils
from pathlib import Path
SCRIPT_DIR = str(Path(__file__).parent.absolute())

class Filter:

	_refdata = Path(f"{SCRIPT_DIR}/refdata")

	def __init__(self, filterset=None, verbose=False, **kwargs):
		if filterset == "sdss":
			if verbose: 
				print("Adopt the SDSS filterset")
			self.sdss_filterset()
		elif filterset == "tophat":
			if verbose: 
				print("Adopt the tophat filterset")
			self.tophat_filterset(**kwargs)
		else:
			if verbose: 
				print("Adopt the default 7DT filterset")
			self.default_filterset()
			
	def edmund_25nm_filterset(self):
		path = str(self._refdata)+f'/edmund_med_band/*_25.csv'
		infilterlist = sorted(glob.glob(self._refdata/path_filter))
		print(f"{len(infilterlist)} Filters found in {path_filter}")

		#	Create filter_set definition
		filter_set = {
			'wavelength': 0
			}

		cwl = []
		filterNameList = []
		bandwidth = []
		################################################################
		ff = 0
		filte = infilterlist[ff]
		for ff, filte in enumerate(infilterlist):
			part = os.path.basename(filte).replace('.csv', '').split('_')
			center_wavelength = float(part[0])*1e1
			bw = int(part[1])*1e1
			# print(filte, center_wavelength, bw)

			_intbl = Table.read(filte)
			_intbl = _intbl[_intbl['wavelength']<1000]
			intbl = Table()
			#	[nm] --> [Angstrom]
			intbl['lam'] = _intbl['wavelength']*1e1
			#	[%] --> ratio
			intbl['col2'] = _intbl['transmission']*1e-2
			#	Wavelength resolution
			lamres = np.mean(intbl['lam'][1:] - intbl['lam'][:-1])
			#	Wavelength min & max
			lammin = intbl['lam'].min()
			lammax = intbl['lam'].max()

			# cwl.append(np.sum(intbl['lam']*intbl['col2']*lamres)/np.sum(intbl['col2']*lamres))
			cwl.append(center_wavelength)
			filterNameList.append(f"m{center_wavelength:g}")
			filter_set.update({f"m{center_wavelength:g}": intbl['col2']})

		cwl = np.array(cwl)

		step = 500
		ticks = np.arange(round(intbl['lam'].min(), -3)-step, round(intbl['lam'].max(), -3)+step, step)

		
		filter_set['wavelength'] = np.array(intbl['lam'].data)

		self.filterset = Table(filter_set)
		self.filterNameList = filterNameList
		self.lam = intbl['lam'].data
		self.lammin = lammin
		self.lammax = lammax
		# self.lammin = 3500
		# self.lammax = 9250
		self.lamres = lamres
		self.bandstep = 12.5
		self.bandwidth = bw
		self.color = utils.makeSpecColors(len(cwl))
		self.ticks = ticks
		return filter_set


	def edmund_50nm_filterset(self):
		path = str(self._refdata)+f'/edmund_med_band/*_50.csv'
		infilterlist = sorted(glob.glob(path))
		print(f"{len(infilterlist)} Filters found in {path_filter}")

		#	Create filter_set definition
		filter_set = {
			'wavelength': 0
			}

		cwl = []
		filterNameList = []
		bandwidth = []
		################################################################
		ff = 0
		filte = infilterlist[ff]
		for ff, filte in enumerate(infilterlist):
			part = os.path.basename(filte).replace('.csv', '').split('_')
			center_wavelength = float(part[0])*1e1
			bw = int(part[1])*1e1
			# print(filte, center_wavelength, bw)

			_intbl = Table.read(filte)
			_intbl = _intbl[_intbl['wavelength']<1000]
			intbl = Table()
			#	[nm] --> [Angstrom]
			intbl['lam'] = _intbl['wavelength']*1e1
			#	[%] --> ratio
			intbl['col2'] = _intbl['transmission']*1e-2
			#	Wavelength resolution
			lamres = np.mean(intbl['lam'][1:] - intbl['lam'][:-1])
			#	Wavelength min & max
			lammin = intbl['lam'].min()
			lammax = intbl['lam'].max()

			# cwl.append(np.sum(intbl['lam']*intbl['col2']*lamres)/np.sum(intbl['col2']*lamres))
			cwl.append(center_wavelength)
			filterNameList.append(f"m{center_wavelength:g}")
			filter_set.update({f"m{center_wavelength:g}": intbl['col2']})

		cwl = np.array(cwl)

		step = 500
		ticks = np.arange(round(intbl['lam'].min(), -3)-step, round(intbl['lam'].max(), -3)+step, step)

		filter_set['wavelength'] = np.array(intbl['lam'].data)

		self.filterset = Table(filter_set)
		self.filterNameList = filterNameList
		self.lam = intbl['lam'].data
		self.lammin = lammin
		self.lammax = lammax
		# self.lammin = 3500
		# self.lammax = 9250
		self.lamres = lamres
		self.bandcenter = cwl
		self.bandstep = 25
		self.bandwidth = bw
		self.color = utils.makeSpecColors(len(cwl))
		self.ticks = ticks
		return filter_set