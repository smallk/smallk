
main smallk repo     = https://github.com/smallk/smallk 
documentation repo = https://github.com/smallk/smallk.github.io

# Setting up the documentation build environment

Clone the main smallk repo and enter the sphinx directory.

    sudo -H pip install -r requirements.txt	#install the theme, and so on

The html documentation and the pdf documentation both rely on Latex -- the html for generating svg images of math equations, and the pdf for compiling auto generated tex files into a pdf. It is necessary to install latex if not already present. If the html or pdf build fails, one potential fix is to install the extended/full latex, which is over a Gigabyte:

    sudo apt-get install texlive-latex-base  	#latex for math
    sudo apt-get install texlive-latex-extra	#fix a utf8x error
    sudo apt-get install dvipng			#for rendering math stuff
    sudo apt-get install latexm			#for rendering pdf


# Editing style guide for source rst documents

We have been attempting to follow the style guide here:

http://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html

# Building html and pdf

Change the working directory to `<main smallk repo>/sphinx` and run three `make` commands:

    make html
    make latex
    make latexpdf

The generated html and pdf documentation will be found in `sphinx/build/html/*` and 'sphinx/build/latex/SmallK.pdf`. 

# Updating online documentation
	
After building the html and pdf, change the working directory to <documentation repo>. You may choose to delete all current html files:

**caution! this assumes all html and related files (js, css, images, etc) are obsolete; it does not clean up the `doc`, `papers` nor `code` directories**

    rm -rf _images _sources _static *.html *.js 

Then, copy the newly generated html and pdf documentation over: 

    cp -r <main smallk repo>/sphinx/build/html/* .
    cp <main smallk repo>/sphinx/build/latex/SmallK.pdf doc/smallk_readme.pdf 

And finally commit the changes to git and push it online:

    git add -A
    git commit -m 'updated documentation'
    git push
