def extract_info_pos( inlist ):
	ipos = {}

	with open( inlist, 'r' ) as intx:

		lineas = intx.readlines()

		for linea in lineas:
			l = linea.strip().split(' ')

			chr = l[1].split(':')[0]

			ipos[ chr ] = {}

		for linea in lineas:
			l = linea.strip().split(' ')
		
			chr = l[1].split(':')[0]
			pos = l[2]

			ipos[ chr ][ pos ] = {}

		for linea in lineas:
			l = linea.strip().split(' ')

			chr = l[1].split(':')[0]
			pos = l[2]
			if max( [ float(l[4]), float(l[5]) ] ) >= 0.75:
				if float(l[4]) > float(l[5]):
					allele = l[3]
				else:
					allele = l[6]
				ipos[ chr ][ pos ][ allele ] = []
			else:
				pass


		for linea in lineas:

			l = linea.strip().split(' ')
		
			chr = l[1].split(':')[0]
			pos = l[2]
			if max( [ float(l[4]), float(l[5]) ] ) >= 0.75:
				if float(l[4]) > float(l[5]):
					allele = l[3]
				else:
					allele = l[6]
				ipos[ chr ][ pos ][ allele ].append( l[0] )
			else:
				pass


	return( ipos )
