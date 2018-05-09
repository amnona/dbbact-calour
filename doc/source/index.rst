.. calour documentation master file, created by
   sphinx-quickstart on Wed Jan 11 16:08:35 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calour dbBact database interface
================================

This is a plugin module for the microbiome interactive analysis tool `Calour <https://github.com/biocore/calour>`_.

This plugin enables Calour to interact with the `dbBact <http://dbbact.org>`_ sequence annotation database in order to display what is known about each bacteria.

Things you can do with Calour-dbBact interface
----------------------------------------------

* Normalize, filter, reorder and cluster your data.

* Permutation based differential abundance testing with powerful `dsFDR <http://msystems.asm.org/content/2/6/e00092-17>`_ correction.

* Interactive heatmap plotting withh convenient zoom, pan, multiple feature selection and information about selected feature/sequence

* Integration with databases (`dbBact <http://dbbact.org>`_, `phenoDB <http://msphere.asm.org/content/2/4/e00237-17>`_, `SpongeEMP <http://www.spongeemp.com/main>`_ for microbiome data, `GNPS <https://gnps.ucsd.edu/ProteoSAFe/static/gnps-splash.jsp>`_ for metabolomics) enables viewing and statistical analysis of database information about the experiment features.


Installing Calour-dbBact
------------------------
* First you need to install Calour (see Calour installation instructions `here <https://github.com/biocore/calour/blob/master/INSTALL.md>`_)

* Install the [dbBact](http://www.dbbact.org) calour interface:

```
pip install git+git://github.com/amnona/dbbact-calour
```

When you run Calour (via python/jupyter notebook/EZCalour), the dbbact-calour interface will be available.

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

Extending Calour
----------------
.. toctree::
   :maxdepth: 2

   database_interface

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

