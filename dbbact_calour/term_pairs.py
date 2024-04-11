from collections import defaultdict
import requests

from logging import getLogger, NOTSET, basicConfig
from logging.config import fileConfig
from pkg_resources import resource_filename

import numpy as np
import scipy as sp

# prepare the logging
logger = getLogger(__name__)

try:
	# get the logger config file location
	log = resource_filename(__name__, 'log.cfg')
	# log = path.join(path.dirname(path.abspath(__file__)), 'log.cfg')
	# set the logger output according to log.cfg
	# setting False allows other logger to print log.
	fileConfig(log, disable_existing_loggers=False)
	# set the log level to same value as calour log level
	clog = getLogger('calour')
	calour_log_level = clog.getEffectiveLevel()
	if calour_log_level != NOTSET:
		logger.setLevel(calour_log_level)
except:
	print('FAILED loading logging config file %s' % log)
	basicConfig(format='%(levelname)s:%(message)s')


def debug(level, msg):
	'''to make compatible with the scdb website
	'''
	if level < 4:
		logger.debug(msg)
	elif level < 7:
		logger.info(msg)
	else:
		logger.warning(msg)


def get_enrichment_score(annotations, seqannotations, ignore_exp=[], term_info=None, term_types=('single'), threshold=None, dbbact_server_url='https://api.dbbact.org', low_number_correction=0):
	'''Get f score, recall and precision for set of annotations

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	ingore_exp: list of str, optional
		list of experiment ids to ignore in the analysis
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info
	term_types: list of str
		the terms to calculate enrichment scores for. options are:
		'single': each single term (i.e. 'feces')
		'pairs': all term pairs (i.e. 'feces+homo sapiens')
	threshold: float or None, optional
		if not None, return only terms that are significantly enriched in the annotations compared to complete database null with p-val <= threshold
	dbbact_server_url: str, optional
		the dbbact server url for getting the term_info for recall (if term_info is None)
	low_number_correction: int, optional
		the constant to penalize low number of annotations in the precision. used as precision=obs/(total+low_number_correction)
	

	Returns
	-------
	fscore: dict of {term(str): fscore(float)}
	recall: dict of {term(str): recall(float)}
	precision: dict of {term(str): precision(float)}
	term_count: dict of {term(str): total experiments(float)}
		the number of experiments where annotations for each term appear
	reduced_f
	'''
	debug(2, 'getting enrichment scores from %d sequences' % len(seqannotations))
	debug(1, 'getting recall')
	recall = get_recall(annotations, seqannotations, ignore_exp=ignore_exp, term_info=term_info, term_types=term_types, dbbact_server_url=dbbact_server_url)
	debug(1, 'getting precision')
	precision = get_precision(annotations, seqannotations, ignore_exp=ignore_exp, term_types=term_types, low_number_correction=low_number_correction)
	debug(1, 'getting term count from get_enrichent_score()')
	term_count = get_term_total_counts(annotations, seqannotations, ignore_exp=ignore_exp, term_types=term_types)
	debug(1, 'calculating the enrichment scores')

	fscore = {}
	for cterm, crecall in recall.items():
		cprecision = precision[cterm]
		fscore[cterm] = 2 * (crecall * cprecision) / (crecall + cprecision)
		# fscore[cterm] = 2 * (crecall * cprecision) / (crecall + cprecision) * ((cprecision + 1.1) / (cprecision + 0.1))
		# fscore[cterm] = np.sqrt(crecall + cprecision)

	# create the reduced f-scores that contain each term only once (for term pairs)
	if 'pairs' in term_types:
		zz = sorted(fscore.items(), key=lambda x: x[1], reverse=True)
		found_items = set()
		reduced_f = {}
		for (cterm, cscore) in zz:
			isok = True
			for ccterm in cterm.split('+'):
				if ccterm in found_items:
					isok = False
					continue
			if isok:
				reduced_f[cterm] = cscore
				for ccterm in cterm.split('+'):
					found_items.add(ccterm)
	else:
		reduced_f = fscore.copy()

	if threshold is not None:
		# filter away non-significant terms
		pvals = get_term_pvals(annotations, seqannotations, term_info, threshold=threshold)
		for cterm, creject in pvals.items():
			if not creject:
				fscore.pop(cterm)
				recall.pop(cterm)
				precision.pop(cterm)
				term_count.pop(cterm)
				if cterm in reduced_f:
					reduced_f.pop(cterm)

	return fscore, recall, precision, term_count, reduced_f


def get_term_total_counts(annotations, seqannotations, ignore_exp=[], term_types=('single')):
	'''Get the number of experiments containing each term from our annotations

	Used to calculate the color for the wordcloud

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	ingore_exp: list of str, optional
		list of experiment ids to ignore in the analysis

	Returns
	-------
	dict of {term(str): number of experiments(int)}
	'''
	term_exps = defaultdict(set)
	for cseqid, cseq_annotations in seqannotations:
		for cannotationid in cseq_annotations:
			cannotation = annotations[str(cannotationid)]
			cexpid = cannotation['expid']
			if cexpid in ignore_exp:
				continue
			for cterm in get_terms(cannotation, term_types=term_types):
				term_exps[cterm].add(cannotation['expid'])

	term_exp_count = {}
	for cterm, cexps in term_exps.items():
		term_exp_count[cterm] = len(cexps)
	return term_exp_count


def get_recall(annotations, seqannotations, method='exp-mean', ignore_exp=[], term_info=None, term_types=('single'), low_num_correction=1, dbbact_server_url='https://api.dbbact.org'):
	'''Calculate the recall (what fraction of the database enteries for this term are covered by the query)

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	method: str, optional
		the method to calculate the recall. options are:
		'exp-mean': calculate the per-experiment mean for the term
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info
	low_num_correction: int, optional
		the constant to penalize low number of experiments in the recall. used as recall=obs/(total+low_num_correction)
	dbbact_server_url: str, optional
		the url of the dbbact server to use for getting the term_info if term_info is None

	Returns
	-------
	dict of {term (str): recall(float)}
	'''
	import itertools
	debug(1, 'calculating recall')
	recall = defaultdict(float)

	# get the term counts for all the terms and store in term_info
	termslist = []
	for cannotation in annotations.values():
		cterms = get_terms(cannotation, term_types=term_types)
		termslist.append(cterms)
	# we use this as much faster than multiple unions
	all_terms = set(itertools.chain.from_iterable(termslist))
	debug(1, 'total terms in all annotations: %d' % len(all_terms))

	if term_info is None:
		debug(2, 'term_info was None, getting from dbbact')
		term_info = get_term_info(list(all_terms), term_types=term_types, dbbact_server_url=dbbact_server_url)
		# term_info = get_term_info(all_terms_positive, term_types=term_types)
	else:
		debug(1, 'term_info already supplied')

	num_sequences = len(seqannotations)
	debug(1, 'total sequences: %d' % num_sequences)
	for cseq, cseq_annotations in seqannotations:
		debug(1, 'processing seq %s' % cseq)
		# for each term (in all the annotations), get all the experiments where it appears
		cseq_term_exps = defaultdict(set)
		for cannotationid in cseq_annotations:
			cannotation = annotations[str(cannotationid)]
			terms = get_terms(cannotation, term_types=term_types)
			cexp = cannotation['expid']
			if cexp in ignore_exp:
				continue
			for cterm in terms:
				cseq_term_exps[cterm].add(cexp)
		# and add the normalized count
		debug(1, 'going over exp list')
		for cterm, cexplist in cseq_term_exps.items():
			debug(1, 'processing term %s' % cterm)
			if cterm not in term_info:
					debug(2, 'term %s not in term_info' % cterm)
					continue
			try:
				observed = term_info[cterm]['total_experiments']
				crecall = len(cexplist) / (observed + low_num_correction)
			except:
				observed = 0
				debug(1, 'term %s does not have total_experiments ok' % cterm)
				# crecall = len(cexplist)
				crecall = 0

			recall[cterm] += crecall / num_sequences
			debug(1, 'term %s observed in %s, total in db %s, recall %s' % (cterm, cexplist, observed, crecall))
	debug(1, 'recall contains %d terms' % len(recall))
	return recall


def get_term_info(terms, term_types=('single'), dbbact_server_url='https://api.dbbact.org'):
	'''
	Get the statistics about each term in annotations

	Parameters
	----------
	terms: list of str
		the terms to get the info about
	term_types: list of str, optional
		Not implemented yet
	dbbact_server_url: str, optional
		the url of the dbbact server to use for getting the term_info

	Returns
	-------
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats)
		The statistics about each term
	'''
	debug(2, 'getting term_info for %d terms' % len(terms))
	res = requests.get(dbbact_server_url + '/ontology/get_term_stats', json={'terms': terms})
	if res.status_code != 200:
		debug(6, 'error encountered in get_term_stats: %s' % res.reason)
		return []
	ans = res.json()
	term_info = ans.get('term_info')
	return term_info


def get_precision(annotations, seqannotations, method='total-annotation', ignore_exp=[], term_types=('single'), low_number_correction=0):
	'''Calculate the precision (how many of the sequences contain the term) for each term in annotations.

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	method: str, optional
		the method to calculate the precision. options are:
		'per-sequence': what fraction of the sequences contain this term at least once
		'total-annotation': what fraction of all sequences annotations contain this term (annotation can be counted more than once since iterating over all seqannotations)
	ignore_exp: list of int, optional:
		the experimentIDs to ignore for the precision calculation (if empty use all experiments)
	term_types: list of str, optional
		types of terms to use. can include the following (including combinations):
			'single': use each term
			'pairs': use term pairs
	low_number_correction: int, optional
		the constant to penalize low number of annotations in the precision. used as precision=obs/(total+low_number_correction)


	Returns
	-------
	dict of {term (str): precision(float)}
	'''
	# get the sequences where each term appears (at least once in their annotations)
	if method == 'per-sequence':
		term_seqs = defaultdict(set)
		for cseqid, cseq_annotations in seqannotations:
			for cannotationid in cseq_annotations:
				cannotation = annotations[str(cannotationid)]
				if cannotation['expid'] in ignore_exp:
					continue
				for cterm in get_terms(cannotation, term_types=term_types):
					term_seqs[cterm].add(cseqid)
		# and calculate the precision (what fraction of the sequences have this term)
		precision = {}
		total_seqs = len(seqannotations)
		for cterm, cterm_seqs in term_seqs.items():
			precision[cterm] = len(cterm_seqs) / total_seqs

	elif method == 'total-annotation':
		term_counts = defaultdict(float)
		for cseqid, cseq_annotations in seqannotations:
			cseq_term_counts = defaultdict(float)
			cseq_total_annotations = 0
			for cannotationid in cseq_annotations:
				cannotation = annotations[str(cannotationid)]
				if cannotation['expid'] in ignore_exp:
					continue
				cseq_total_annotations += 1
				for cterm in get_terms(cannotation, term_types=term_types):
					# we weigh each annotation by the number of annotations for this sequence (since we want mean over all sequences)
					cseq_term_counts[cterm] += 1
					# if we use the annotation type score - must fix normalization!!!!! need to do
					# term_counts[cterm] += get_annotation_type_score(cannotation) / cseq_total_annotations
			if cseq_total_annotations == 0:
				continue
			for cterm in cseq_term_counts.keys():
				term_counts[cterm] += cseq_term_counts[cterm] / (cseq_total_annotations + low_number_correction)
		precision = {}
		total_seqs = len(seqannotations)
		for cterm, cterm_counts in term_counts.items():
			precision[cterm] = cterm_counts / total_seqs

	else:
		raise ValueError('method %s unknown' % method)

	return precision


def get_annotation_type_score(annotation):
	'''Get the score factor associated with an annotation type.
	Score is based on the annotation type (i.e. "common/dominant/diffexp/contam/other/positive correlation/negative correlation")

	Parameters
	----------
	annotation: dict
		as from dbbact rest-api annotations/get_annotation.
		should contain at least:
			"annotationtype"
	Returns
	-------
	float: the annotation score factor
	'''
	score = 1
	anntationtype = annotation['anntationtype']
	if anntationtype == 'highfreq':
		score = 2
	elif anntationtype == 'common':
		score = 1
	elif anntationtype == 'other':
		score = 1
	elif anntationtype == 'diffexp':
		score = 1
	elif anntationtype == 'positive correlation':
		score = 1
	elif anntationtype == 'negative correlation':
		score = 1
	elif anntationtype == 'contamination':
		score = 1
	else:
		debug(4, 'unknown annotation type %s' % anntationtype)
	return score


def tessa(source):
	'''get all pairs from a list
	'''
	result = []
	for p1 in range(len(source)):
			for p2 in range(p1 + 1, len(source)):
					result.append([source[p1], source[p2]])
	return result


def get_terms(annotation, term_types=('single')):
	'''Get a list of terms present in the annotation. terms that are "lower in" are preceded by a "-"

	Parameters
	----------
	annotation: dict
		as from dbbact rest-api annotations/get_annotation.
		should contain at least:
			"annotationid" (str)
			"annotationtype"
			"details" (list of [detail_type, term])
	term_types: list of str, optional
		types of terms to use. can include the following (including combinations):
			'single': use each term
			'pairs': use term pairs

	Returns
	-------
	list of str - the terms in the annotation
	'''
	terms = []
	if 'single' in term_types:
		details = annotation['details']
		for cdetail in details:
			cterm = cdetail[1]
			ctype = cdetail[0]
			if ctype == 'low':
				cterm = '-' + cterm
			terms.append(cterm)

		# handle the contamination annotation as well
		if annotation['annotationtype'] == 'contamination':
			terms.append('contamination')
		debug(1, 'found %d single terms for annotation %s' % (len(terms), annotation['annotationid']))

	if 'pairs' in term_types:
		# add pairs. of the form (term1+term2) where term1 is before term2 alphabetically
		if len(details) < 10:
			single_terms = []
			for cdetail in details:
				cterm = cdetail[1]
				ctype = cdetail[0]
				if ctype == 'low':
					cterm = '-' + cterm
				single_terms.append(cterm)
			pairs = tessa(single_terms)
			pairs = ['+'.join(sorted([x, y])) for (x, y) in pairs]
			debug(1, 'found %d term pairs for annotation %s' % (len(pairs), annotation['annotationid']))
			terms.extend(pairs)
		else:
			debug(1, 'too many terms (%d) for term pair calculation for annotation %s' % (len(details), annotation['annotationid']))

	debug(1, 'found %d total terms for annotation %s' % (len(terms), annotation['annotationid']))
	return terms


def get_term_pvals(annotations, seqannotations, term_info=None, ignore_exp=[], term_types=('single'), threshold=0.01):
	'''Calculate the precision (how many of the sequences contain the term) for each term in annotations.

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info
	ignore_exp: list of int, optional:
		the experimentIDs to ignore for the precision calculation (if empty use all experiments)
	term_types: list of str, optional
		types of terms to use. can include the following (including combinations):
			'single': use each term
			'pairs': use term pairs
	threshold: float, optional
		maximal p-value to filter

	Returns
	-------
	dict of {term (str): reject (bool)}
		if True, reject null hypothesis (not random)
	'''
	# get the sequences where each term appears (at least once in their annotations)
	term_counts = defaultdict(int)
	total_annotations = 0
	for cseqid, cseq_annotations in seqannotations:
		cseq_term_counts = defaultdict(float)
		cseq_total_annotations = 0
		for cannotationid in cseq_annotations:
			cannotation = annotations[str(cannotationid)]
			if cannotation['expid'] in ignore_exp:
				continue
			cseq_total_annotations += 1
			for cterm in get_terms(cannotation, term_types=term_types):
				cseq_term_counts[cterm] += 1
		if cseq_total_annotations == 0:
			continue
		for cterm in cseq_term_counts.keys():
			term_counts[cterm] += cseq_term_counts[cterm]
		total_annotations += cseq_total_annotations

	# TODO: fix this!!!!
	total_db_annotations = 3925
	debug(2, 'total db annotations: %d' % total_db_annotations)
	debug(2, 'total annotations: %d' % total_annotations)
	pvals = {}
	for cterm, cterm_counts in term_counts.items():
		if cterm not in term_info:
			continue
		p_null = term_info[cterm]['total_annotations'] / total_db_annotations
		pv = 1 - sp.stats.binom.cdf(cterm_counts - 1, total_annotations, p_null)
		if pv > threshold:
			pvals[cterm] = False
			# print('term %s, counts %d, observed p: %f, p_null: %f, pv:%f' % (cterm, cterm_counts, cterm_counts / total_annotations, p_null, pv))
		else:
			pvals[cterm] = True
	return pvals
