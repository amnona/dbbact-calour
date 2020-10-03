#!/usr/bin/env python

# add database based whole sequence ids to dbbact sequences

'''Add SILVA/Greengenes IDs to each dbBact sequence
The IDs are added to the WholeSeqIDs table which contains the following columns:
dbID: SILVA/GreenGenes
dbbactid: the dbbact id of the sequence
wholeseqid: the whole sequence database id

NOTE: each dbbact sequence can match more than one silva/GG id! so can have multiple entries for same dbbactid

Adding is based on complete match of dbbact sequence to database (as a subsequence)
'''

import sys
from collections import defaultdict
import random

import argparse
import psycopg2
import psycopg2.extras

__version__ = "0.1"


def debug(level, msg):
	print(msg)


def hash_sequences(filename, short_len=100):
	'''hash all the sequences in a fasta file

	Parameters
	----------
	filename: str
		the fasta file

	Returns
	-------
	seq_hash: dict of {seq: seqid}
	seq_lens : list of int
		all the sequence lengths in the fasta file (so we can hash all the lengths in the queries)
	short_hash: dict of {short_seq: seq_hash dict}
	'''
	num_too_short = 0
	seq_hash = {}
	seq_lens = set()
	short_hash = defaultdict(dict)
	for cseq, chead in iter_fasta_seqs(filename):
		clen = len(cseq)
		if clen < short_len:
			num_too_short += 1
			continue
		short_seq = cseq[:short_len]
		short_hash[short_seq][cseq] = chead
		if clen not in seq_lens:
			seq_lens.add(clen)
		seq_hash[cseq] = chead
	print('processed %d sequences.' % len(seq_hash))
	print('lens: %s' % seq_lens)
	print('num too short: %d' % num_too_short)
	return seq_hash, seq_lens, short_hash


def iter_fasta_seqs(filename):
	"""
	iterate a fasta file and return header,sequence
	input:
	filename - the fasta file name

	output:
	seq - the sequence
	header - the header
	"""

	fl = open(filename, "rU")
	cseq = ''
	chead = ''
	for cline in fl:
		if cline[0] == '>':
			if chead:
				yield(cseq.lower(), chead)
			cseq = ''
			chead = cline[1:].rstrip()
		else:
			cline = cline.strip().lower()
			cline = cline.replace('u', 't')
			cseq += cline.strip()
	if cseq:
		yield(cseq, chead)
	fl.close()


def seqs_to_silva_ids(dbfile, inputfile, outputfile, short_len=100, max_match=0):
	seq_hash, seq_lens, short_hash = hash_sequences(filename=inputfile, short_len=short_len)
	idx = 0
	num_matches = 0
	matches = defaultdict(list)
	for cseq, chead in iter_fasta_seqs(dbfile):
		idx += 1
		if idx % 1000 == 0:
			print(idx)
		for cpos in range(len(cseq) - short_len):
				ccseq = cseq[cpos:cpos + short_len]
				if ccseq in short_hash:
					# print('found short with %d sequences' % len(short_hash[ccseq]))
					for k, v in short_hash[ccseq].items():
						if k in cseq:
							cid = chead.split(' ')[0]
							matches[v].append(cid)
							num_matches += 1
	print('found %d dbbact sequences in database' % num_matches)
	with open(outputfile, 'w') as ofl:
		for cseq, cids in matches.items():
			if max_match > 0:
				if len(cids) > max_match:
					random.shuffle(cids)
					cids = cids[:max_match]
			ofl.write('>%s\n%s\n' % (','.join(cids), cseq))
	print('done')


def main(argv):
	parser = argparse.ArgumentParser(description='seqs_to_silva_ids - create a fasta file with the silva ids of each sequence. version ' + __version__)
	parser.add_argument('--dbfile', help='name of whole sequence database fasta file (i.e. greenegenes or silva)')
	parser.add_argument('-i', '--input', help='name of the input fasta file')
	parser.add_argument('-m', '--max', help='maximal number of db matches to keep (otherwise select random subset of matches. 0 to keep all', type=int, default=0)
	parser.add_argument('-o', '--output', help='name of the output fasta file')
	args = parser.parse_args(argv)
	seqs_to_silva_ids(dbfile=args.dbfile, inputfile=args.input, outputfile=args.output, max_match=args.max)


if __name__ == "__main__":
	main(sys.argv[1:])
