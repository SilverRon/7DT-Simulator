def extract_param_kn_sim_cube(knsp):
	part = os.path.basename(knsp).split('_')

	if part[1] == 'TP':
		dshape = 'toroidal'
	elif part[1] == 'TS':
		dshape = 'spherical'
	else:
		dshape = ''

	#	Latitude
	if part[5] == 'wind1':
		lat = 'Axial'
	elif part[5] == 'wind2':
		lat = 'Edge'
	else:
		lat = ''

	#	Ejecta mass for low-Ye [solar mass]
	md = float(part[7].replace('md', ''))

	#	Ejecta velocity for low-Ye [N*c]
	vd = float(part[8].replace('vd', ''))

	#	Ejecta mass for high-Ye [solar mass]
	mw = float(part[9].replace('mw', ''))

	#	Ejecta velocity for high-Ye [N*c]
	vw = float(part[10].replace('vw', ''))

	#	Angle
	try:
		if 'angularbin' not in knsp:
			angle = float(part[11].replace('angle', ''))
		else:
			angle = int(part[11].replace('angularbin', ''))
	except:
		angle = None

	return (dshape, lat, md, vd, mw, vw, angle)
