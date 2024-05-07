# 7DT-Simulator

## Overview
`7DT-Simulator` is to help user to calculate the expected photometric results for a specific set of optical spectra, considering the realistic sensitifity of below components:
1. Optics
2. Filter Transmission Curve
3. Sky Transmission
4. Camera Quantum Efficiency

- *7-Dimensional Telescope (`7DT`) is one of the most biggest multiple telescope system, 20 0.5-m telescopes with 20-40 medium-band filters in Chile. The main goal of `7DT` is is to rapidly identify the optical counterpart of gravitaitonal-wave counteraprt, kilonova.
- **7-Dimensional Sky Survey (`7DS`) is survey project with `7DT` to spectral mapping of all-sky.

If you are interested in our project, please visit Google Document about details 
- https://docs.google.com/document/d/1WZ4sfvgd8PFF6FsxNkKlroQqSgf5JS7jZd158eS7RLQ/edit?usp=sharing

## Installation
The required version of Python is `>=3.1X`. For the lower version of Python, it needs to be checked in near future. If you install this package, **I higly recommend to install it with a new conda environment.**

Note: *The installation using `setup.py` is under construction.*

1. Clone this repository.
	```
	git clone https://github.com/SilverRon/7DT-Simulator
	```

2. Install the required libraries and modules.
	```
	cd 7DT-Simulator/
	python setyp.py install
	```

Here are requirements in `requirements.txt`:
```
astropy==6.0.0
astroquery==0.4.6
ccdproc==2.4.1
GaiaXPy==2.1.0
matplotlib==3.7.1
numpy==1.24.2
pandas==2.0.1
phot==0.2.12
scikit_learn==1.2.2
scipy==1.12.0
seaborn==0.13.2
setuptools==65.5.0
speclite==0.14
```
If you want to install dependencies, run the following command:
```
pip install -r requirements.txt
```

## Usage
This version of `7DT-Simulator` has limited functionality, but there is no problem to generate a simulated photometry. Below example makes you start your work quickly:
```
test/Example_7DT_Synthetic_Photometry.ipynb
```
If you have a spectrum data that you want to convert to simulated `7DT` photometry containing `lam` and `flam`, it is easy! 

**Note that this version works only in relative path like `7DT-Simulator/test` or `7DT-Simulator/code`**.

## Authors and Acknowledgment

If you are interested in `7DT`/`7DS`, or need some supports, please let us know.

| Name                | Role                   | Institution | Email                     |
|---------------------|------------------------|-------------|---------------------------|
| Mr. Gregory S.H. Paek | Project SW Main Developer | SNU         | gregorypaek94_at_g_mail |
| Prof. Myungshin Im  | Principal Investigator (PI) | SNU         |                           |
| Dr. Ji Hoon Kim     | Project Manager (PM)   | SNU         |                           |

This code is based on the modified and advanced code from the Dr. Yujin Yang (KASI) in 2022. We thank to the Dr. Yujin Yang.

## Changelog
### v0.1.0 - 2024-02-11
- Initial Commit
### v0.1.1 - 2024-02-13
- Modify `README.md`
### v0.1.2 - 2024-05-07
- Simplify codes

## License

This project is provided under the MIT license. For more details, please refer to the `LICENSE` file.

## Features
TBD

## Documentation
TBD

## Contributing
TBD

