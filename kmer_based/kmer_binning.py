from Bio import SeqIO
import sys
import re
from collections import Counter

diagk = sys.argv[1]
fastq = sys.argv[2]
k = 25
kmatch_th = 10

kmers = {}

with open( diagk, 'r' ) as infile:
	tx = infile.readlines()
	for t in tx[1:]:
		e = t.replace('"', '').split(',')
		key = tuple([e[2], e[3]])
		if key in kmers:
			kmers[key].append(e[1])
		else:
			kmers[key] = [e[1]]



for r in SeqIO.parse( fastq, 'fastq' ):			#### Looping through the fastq file

	sp = []
	klist = []
	kmatch_list = []

	for i in range( len(str(r.seq)) - k + 1):
		klist.append( str(r.seq)[i:i+k] )

	klist = set(klist)				#### klist: the set of k-mers in a read
	
	for splocus in kmers:				#### Looping through the sp/locus pairs, each of which has a list of diagnostic k-mers

		intsct = klist.intersection( kmers[splocus] )	#### intsct: the intersection between the list of read k-mers and the list of diagnostic k-mers of each sp/locus pairs		
		kmatch_list.append( len(intsct) )		#### kmatch_list: stores the length of all intsct

		if len(intsct) > kmatch_th:			#### Check that intsct is higher than threshold (Set to 10)
			kmatch = tuple( [ splocus, len(intsct) ] )	#### Append a tuple with the sp/locus pair and the length of intsct to the sp list
			sp.append( kmatch )

	assigned_sp = []

	for kmatch in sp:

		if kmatch[1] == max(kmatch_list):
			assigned_sp.append(kmatch)

#		#else:
#			pass
	

	if len( set( assigned_sp ) ) == 1:
		
		print( '%s,%s,%d' % (r.description, ','.join(list( assigned_sp[0][0] )), assigned_sp[0][1] ) ) 
	else:
		print( '%s,NA,NA,NA' % (r.description) )

