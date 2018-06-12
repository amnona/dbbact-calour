.. dbbact_calour documentation master file, created by
   sphinx-quickstart on Wed Jan 11 16:08:35 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calour dbBact database interface
================================

This is a plugin module for the microbiome interactive analysis tool `Calour <https://github.com/biocore/calour>`_.

This plugin enables Calour to interact with the `dbBact <http://dbbact.org>`_ sequence annotation database in order to display what is known about each bacteria.

Things you can do with Calour-dbBact interface
----------------------------------------------

* Get information about bacterial sequences from the Calour interactive heatmap

* Identify database terms enriched in subgroups of bacteria

* Get information about the distribution of a given dbbact term in the experiment bacteria

Installing Calour-dbBact
------------------------
* First you need to install Calour (see Calour installation instructions `here <https://github.com/biocore/calour/blob/master/INSTALL.md>`_)

* Install the [dbBact](http://www.dbbact.org) calour interface:

```
pip install git+git://github.com/amnona/dbbact-calour
```

When you run Calour (via python/jupyter notebook/EZCalour), the dbbact-calour interface will be available. 

When running the Calour.plot() command, you can use the databases=['dbbact'] parameter (set by default).

To get the dbbact interface class from calour, use:

```
import calour.database
```

```
db = calour.database._get_database_class('dbbact')
```


Tutorials
---------
.. toctree::
   :maxdepth: 2

   tutorials_microbiome


The Jupyter notebooks in this tutorial section can be downloaded `here <https://github.com/biocore/calour/tree/master/doc/source/notebooks>`_.


API Documentation
-----------------
.. toctree::
   :maxdepth: 2

   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
