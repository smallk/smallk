.. SmallK documentation master file, created by
   sphinx-quickstart on Thu Aug 17 13:41:06 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. It is likely that several requirements must be installed first.
   sudo -H pip install -r requirements.txt	#install the theme, and so on
   sudo apt-get install texlive-latex-base  	#latex for math
   sudo apt-get install texlive-latex-extra 	#fix a utf8x error
   sudo apt-get install dvipng		 	#for rendering math stuff
   sudo apt-get install latexm                  #for rendering pdf

.. sort of trying to follow style guide here: http://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html

######
SmallK
######

.. toctree::
   :maxdepth: 8
   :numbered:
   :hidden:

   self
   pages_about
   pages_introduction
   pages_quickstartInstall
   pages_quickstartSmallkAPI
   pages_installation
   pages_commandLineTools
   pages_smallkAPI
   pages_pysmallkAPI
   pages_tests
   pages_benchmarks_results
   pages_publications
   pages_software_repo


.. only:: html
   
   .. include:: index_splash.rst

Partners
--------

----------------------------------------------------------------------

.. raw:: html

   <center>

|gtri|_

.. raw:: html

   </center>


.. raw:: html

   <center>

|gt|_
|darpa|_

.. raw:: html

   </center>

----------------------------------------------------------------------

..
   Indices and tables
   ==================
   
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

.. |darpa| image:: _static/darpa.png
   :height: 588 px
   :width: 1096 px
   :scale: 15 %
.. _darpa: http://www.darpa.mil/

.. |darpa_orig| image:: _static/darpa_original.png
   :height: 300 px
   :width: 600 px
   :scale: 70 %
.. _darpa_orig: http://www.darpa.mil/

.. |gtri| image:: _static/gtri.png
   :height: 300 px
   :width: 1418 px
   :scale: 20 %
.. _gtri: https://gtri.gatech.edu/

.. |gt| image:: _static/georgiatech.png
   :height: 235 px
   :width: 532 px
   :scale: 40 %
.. _gt: http://www.gatech.edu/

