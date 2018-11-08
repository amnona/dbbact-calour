from collections import defaultdict
import requests

from logging import getLogger, NOTSET, basicConfig
from logging.config import fileConfig
from pkg_resources import resource_filename


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


def get_enrichment_score(annotations, seqannotations, ignore_exp=[], term_info=None, term_types=('single')):
	'''Get f score, recall and precision for set of annotations

	Parameters
	----------
	annotations: dict of {annotationid (str): annotation(dict)}
	seqannotations: list of (seqid, [annotation ids])
	ingore_exp: list of str, optional
		list of experiment ids to ignore in the analysis
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats) or None, optional
		The statistics about each term. if None, the function will contact dbbact to get the term_info

	Returns
	-------
	fscore: dict of {term(str): fscore(float)}
	recall: dict of {term(str): recall(float)}
	precision: dict of {term(str): precision(float)}
	term_count: dict of {term(str): total experiments(float)}
		the number of experiments where annotations for each term appear
	reduced_f
	'''
	# term_info=None
	logger.debug(2, 'getting enrichment scores from %d sequences' % len(seqannotations))
	debug(1, 'getting recall')
	recall = get_recall(annotations, seqannotations, ignore_exp=ignore_exp, term_info=term_info, term_types=term_types)
	debug(1, 'getting precision')
	precision = get_precision(annotations, seqannotations, ignore_exp=ignore_exp, term_types=term_types)
	debug(1, 'getting term count from get_enrichent_score()')
	term_count = get_term_total_counts(annotations, seqannotations, ignore_exp=ignore_exp, term_types=term_types)
	debug(1, 'calculating the enrichment scores')

	fscore = {}
	for cterm, crecall in recall.items():
		cprecision = precision[cterm]
		fscore[cterm] = 2 * (crecall * cprecision) / (crecall + cprecision)

	# create the reduced f-scores that contain each term only once (for term pairs)
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


def get_recall(annotations, seqannotations, method='exp-mean', ignore_exp=[], term_info=None, term_types=('single'), low_num_correction=1):
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

	Returns
	-------
	dict of {term (str): recall(float)}
	'''
	# get the term counts for all the terns
	debug(1, 'calculating recall')
	recall = defaultdict(float)
	all_terms = set()
	for cannotation in annotations.values():
		cterms = get_terms(cannotation, term_types=term_types)
		all_terms = all_terms.union(set(cterms))
	# all_terms_positive = [x[1:] if x[0] == '-' else x for x in all_terms]
	debug(1, 'total terms in all annotations: %d' % len(all_terms))

	if term_info is None:
		debug(2, 'term_info was None, getting from dbbact')
		term_info = get_term_info(list(all_terms), term_types=term_types)
		# term_info = get_term_info(all_terms_positive, term_types=term_types)

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
		# if 'whale blow' not in cseq_term_exps:
		# 	print('*whale blow not found for sequence')
		# 	print(cseq_term_exps)
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
				debug(4, 'term %s does not have total_experiments ok' % cterm)
				# crecall = len(cexplist)
				crecall = 0

			recall[cterm] += crecall / num_sequences
			debug(1, 'term %s observed in %s, total in db %s, recall %s' % (cterm, cexplist, observed, crecall))
	debug(1, 'recall contains %d terms' % len(recall))
	return recall


def get_term_info(terms, term_types=('single')):
	'''
	Get the statistics about each term in annotations

	Parameters
	----------
	terms: list of str
		the terms to get the info about

	Returns
	-------
	term_info: dict of {term (str): details {"total_annotations": float, "total_sequences": float}} (see dbbact rest-api /ontology/get_term_stats)
		The statistics about each term
	'''
	debug(2, 'getting term_info for %d terms' % len(terms))
	res = requests.get(get_db_address() + '/ontology/get_term_stats', json={'terms': terms})
	if res.status_code != 200:
		debug(6, 'error encountered in get_term_stats: %s' % res.reason)
		return []
	ans = res.json()
	term_info = ans.get('term_info')
	return term_info


def get_precision(annotations, seqannotations, method='total-annotation', ignore_exp=[], term_types=('single')):
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
				term_counts[cterm] += cseq_term_counts[cterm] / cseq_total_annotations
		precision = {}
		total_seqs = len(seqannotations)
		for cterm, cterm_counts in term_counts.items():
			precision[cterm] = cterm_counts / total_seqs

	else:
		raise ValueError('method %s unknown' % method)

	return precision


def get_annotation_type_score(annotation):
	'''Get the score factor associated with an annotation type.
	Score is based on the annotation type (i.e. "common/highfreq/diffexp/contamination/other")

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
