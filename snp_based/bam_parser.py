import pysam
import sys


def parse_bam( inbam ):
	entries = []
	samfile = pysam.AlignmentFile( inbam, "rb" )
	complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}

	for r in samfile:


		haplotypes = []

		nucs = ['A', 'T', 'G', 'C']
		seq = str(r.seq)
		read = r.qname
		cigar = r.cigarstring
		flag = r.flag
		cigar_n = cigar.replace('X', '=').replace('I', '=').replace('D', '=').replace('N', '=').replace('M', '=').replace('S', '=').replace('H', '=').split('=')
		offset = r.pos
		cigar_l = [i for i in cigar if not i.isdigit()]

		cigar_list = []

		for i in range(0, len(cigar_l)):
			cigar_element = tuple([cigar_n[i], cigar_l[i]])
			cigar_list.append(cigar_element)

		seq_tr = ''

		for c in cigar_list:

			if c[1] == 'M' or c[1] == '=' or c[1] == 'X':
				seq_tr += seq[0:int(c[0])]
				seq = seq[int(c[0]) :len(seq)  ]



			elif c[1] == 'D' or c[1] == 'N':
				seq_tr += int(c[0]) * '-'


			elif c[1] == 'I' or c[1] == 'S' or c[1] == 'H' :
				seq = seq[int(c[0]) :len(seq) ]						

		

			if flag == '16' or flag == '1040':

				try:

					seq_tr = "".join([complement[x] for x in seq_tr])

				except KeyError:

					pass

			else:

				seq_tr = seq_tr			
		entries.append( '%s,%s,%s,%s' % (r.reference_name, read, seq_tr, offset) )
	return(entries)
