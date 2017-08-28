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
   pages_documentation
   pages_publications
   pages_software_repo


.. only:: html
   
   .. include:: index_splash.rst

..
   Indices and tables
   ==================
   
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

