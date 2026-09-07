import importlib
import os
from logging.config import fileConfig
from logging import getLogger, NOTSET, basicConfig

__credits__ = "https://github.com/amnona/dbbact-calour/graphs/contributors"
__version__ = '2026.09.07'


def get_resource_filename(package, resource):
	'''Return the filesystem path for a package resource.'''
	module = importlib.import_module(package)
	if module.__file__ is None:
		raise ImportError('Package %s has no filesystem location' % package)
	return os.path.join(os.path.dirname(module.__file__), resource)

# from .dbbact import DBBact

# load the logger config
try:
	# get the logger config file location
	log_file = get_resource_filename(__package__, 'log.cfg')
	# log = path.join(path.dirname(path.abspath(__file__)), 'log.cfg')
	# set the logger output according to log.cfg
	# setting False allows other logger to print log.
	fileConfig(log_file, disable_existing_loggers=False)
except:
	print('failed to load logging config file')
	basicConfig(format='%(levelname)s:%(message)s')

print('creating logger')
logger = getLogger(__package__)
# set the log level to the same as calour module if present
try:
	clog = getLogger('calour')
	calour_log_level = clog.getEffectiveLevel()
	if calour_log_level != NOTSET:
		logger.setLevel(calour_log_level)
except:
	print('calour module not found for log level setting. Level not set')
