import sys
from os import listdir
from os.path import isfile, join

import info_pos_extractor
import bam_parser


indire = sys.argv[1]
infpos = sys.argv[2]
#refsp = sys.argv[3]

snps = info_pos_extractor.extract_info_pos( infpos )

onlyfiles = [ f for f in listdir( indire ) if isfile( join( indire,f ) ) and 'mkdups' in f.split('.') and  'bam' in f.split('.') and 'bai' not in f.split('.')]

#print(onlyfiles)

########## Function for finding the most common element in list (taken from 
########## https://stackoverflow.com/a/1520716)

import itertools
import operator

def most_common(L):
  # get an iterable of (item, iterable) pairs
	SL = sorted((x, i) for i, x in enumerate(L))
  # print 'SL:', SL
	groups = itertools.groupby(SL, key=operator.itemgetter(0))
  # auxiliary function to get "quality" for an item
	def _auxfun(g):
		item, iterable = g
		count = 0
		min_index = len(L)
		for _, where in iterable:
			count += 1
			min_index = min(min_index, where)
    # print 'item %r, count %r, minind %r' % (item, count, min_index)
		return count, -min_index
  # pick the highest-count/earliest item
	return max(groups, key=_auxfun)[0]
	
########## Parsing BAM file and classifying each mapped sequence



for file in onlyfiles:

	f_file = join( indire, file )
	lib = file.split('.')[0]
	bam = bam_parser.parse_bam( f_file )

	for b in bam:

		chr = b.split(',')[0]
		read = b.split(',')[1]
		seq = b.split(',')[2]
		offset = int(b.strip().split(',')[3])
		sp_list = []

		if chr in snps:
			#print( snps[ chr ] )
			poss = [pos for pos in snps[ chr ] ]
			for i in range(0, len( poss ) ):
				pos = int( poss[i] )
				
				if pos >= offset and pos <= offset + len( seq ) - 1:
					#try:
					if  seq[ pos - offset - 1 ].upper() in snps[ chr ][ poss[i] ] :

						base = seq[ pos - offset - 1].upper()
	
						if  base in snps[ chr ][ poss[i] ]:#len( set(snps[ chr ][ poss[i] ][ base ] ))  > 0 :
							sp_list.append( snps[ chr ][ poss[i]][ base ] )
						else:
							sp_list.append('NA')
					#except IndexError:
						#sp_list.append('NA')
					#	print('Position out of sequence BAM:%s READ:%s SEQ_LEN:%d POS:%d C_POS:%d OFF:%d' % (file, read, len(seq), pos, pos - offset -1, offset ) )
											
			if len( sp_list ) > 0:
				if len( most_common( sp_list) ) == 1:
					print(lib, chr, read,''.join( most_common( sp_list ) ) )
				else:
					print(lib, chr, read, 'NA' )

			else:
				print(lib, chr, read, 'NA' )
		else:
			print(lib, chr, read, 'NA' )
