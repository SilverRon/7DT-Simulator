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

## Usage
This version of `7DT-Simulator` has limited functionality, but there is no problem to generate a simulated photometry. Below example makes you start your work quickly:
```
test/Example_7DT_Synthetic_Photometry.ipynb
```
If you have a spectrum data that you want to convert to simulated `7DT` photometry containing `lam` and `flam`, it is easy! 

**Note that this version works only in relative path like `7DT-Simulator/test` or `7DT-Simulator/code`**.

## Features
TBD

## Documentation
TBD

## Contributing
TBD

## License

This project is provided under the MIT license. For more details, please refer to the `LICENSE` file.

## Authors and Acknowledgment

- Mr. Gregory S.H. Paek (SNU): Project SW Main Developer
	- Email: gregorypaek94__at__g__mail__
- Prof. Myungshin Im (SNU): Principal Investigator (PI)
- Dr. Ji Hoon Kim (SNU): Project Manager (PM)

This code is based on the modified and advanced code from the Dr. Yujin Yang (KASI) in 2022.

## Changelog
### v0.1.0 - 2024-02-11
- Initial Commit
### v0.1.1 - 2024-02-13
- Modify `README.md`