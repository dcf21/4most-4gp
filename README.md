# 4most-4gp

This is the code for the IWG7 Galactic pipeline for the 4MOST multi-object
spectrograph. It comprises a collection of Python modules which use a common data format to
store and manipulate spectra and their associated metadata. It makes it easy to pass spectra
between a range of spectral synthesis and processing tools including Turbospectrum, the
4MOST Facility Simulator (4FS) and the Cannon, without the need for manual data format
conversion. It includes the ability to store spectra in libraries and search for them by arbitrary
metadata constraints, making it easy to create new tests on subgroups of stars filtered from
larger samples.

In addition, a simple web interface allows the contents of spectrum libraries to be searched
and viewed quickly for diagnostic purposes.

The framework is available in two repositories on GitHub, and includes step-by-step
installation instructions. The first repository contains the Python modules which provide
programmatic interfaces for creating and manipulating libraries of spectra, including wrappers
for passing them to various analysis tools:

[https://github.com/dcf21/4most-4gp]

The second repository contains python scripts which utilise these modules to synthesise the
spectra and perform the tests described in this report:

[https://github.com/dcf21/4most-4gp-scripts]

This code is under active development, and stable releases are periodically
made.

Visiting the GitHub URLs above will present you with the `master` branch of our
code, which should always correspond to the latest stable release. If you click
on the "branches" dropdown menu, you can select a different version of the code
to download.

Stable releases are given date stamps, for example, `release-2018-01-09-1`. The
master branch points to the most recent release. The `dev` branch is not stable
and should be used with extreme caution.

# Code structure

The code is organised into a collection of python modules, most of which do not
depend on each other. These can be found in the directory `src/pythonModules`.

The code is divided up in this way since each module has different dependencies,
and it allows functions to be used without installing every dependency.

It is strongly recommended that they be installed in the python virtual
environment as described below, and that you do not tamper with your
system-wide python installation.

The python modules are as follows:

**fourgp_speclib** - Defines core classes representing spectra and libraries to keep them in. Functionality to search for spectra within libraries is provided using both SQLite and MySQL. All other 4GP modules depend on this core module.

**fourgp_specsynth** -- Provides a wrapper for synthesising spectra with given stellar parameters using Turbospectrum.

**fourgp_fourfs** -- Provides a wrapper for degrading spectra using 4FS; specifically, using the 4MOST Exposure Time Calculator (ETC).

**fourgp_cannon** -- Provides a wrapper for passing arrays of spectra to the Cannon.

**fourgp_rv** -- Provides a simple MCMC code for estimating the radial velocities of spectra.

**fourgp_degrade** -- Provides very simple functions for convolving spectra with Gaussians and for interpolating them onto new wavelength rasters.

**fourgp_telescope_data** -- Provides some basic data about 4MOST, including its wavelength coverage and spectral resolution.

# Getting 4GP to do useful things

All of this code in this repository is simply a collection of python modules which provide a programmatic interface for handling spectra. There is no example code here to show how to use these modules in practice.

Such examples can be found separately in the respository [4most-4gp-scripts](https://github.com/dcf21/4most-4gp-scripts).

# Installing 4GP

For more complete instructions, visit the [Wiki](https://github.com/dcf21/4most-4gp/wiki/dependencies) on our GitHib page.

## Dependencies

Before you start, you need to have all the dependencies you need to build Turbospectrum, 4FS, and the other tools that 4GP wraps.

4GP currently only runs on python 2.7.

These packages and libraries are required:

* *git* (required to check the code out from GitHub)
* *SQLite3* (including the python-sqlite3 binding; you can test for this by typing `import sqlite3` into a python terminal)
* *python-tk* (required to use the 4GP Spectrum Browser and the Cannon; you can test for this by typing `import tkinter` into a python terminal)

The following packages are strongly recommended:

* *python-virtualenv* (required to set up a python virtual environment)

The following packages are needed to run certain parts of 4GP:

* *libchealpix-dev* (required to build and install 4FS)
* *python-healpy* (required to build and install 4FS)
* *libcfitsio3-dev* (required to build and install 4FS)
* *gfortran* (required to build and install Turbospectrum)
* *MySQL* (including the libmysqlclient and python-mysql bindings -- currently only required by the 4GP unit tests, so not very important)

## Ubuntu installation instructions

Under Ubuntu Linux, you can install all these packages with a single command, as follows:

Ubuntu 16.04:

```
apt-get install git python-sqlite3 mysql-server libmysqlclient-dev python-virtualenv libhealpix-cxx-dev libchealpix-dev libcfitsio3-dev python-healpy gfortran python-tk
```

Ubuntu 14.04:

```
apt-get install git python-sqlite mysql-server libmysqlclient-dev python-virtualenv libcfitsio3-dev libblas-dev liblapack-dev libblas3gf liblapack3gf gfortran python-tk
```

## Installing 4GP in a python virtual environment

4GP is distributed with a standard `setuptools` installation script `setup.py` which allows you to install its
constituent modules into your local python environment.

We strongly recommended that you install them in a python virtual
environment, rather than tampering with your
system-wide python installation.

Follow these steps in a Linux shell to do this:

```
# Check out code from GitHub
git clone https://github.com/dcf21/4most-4gp.git
 
# Set up a python virtual environment
virtualenv virtualenv
source virtualenv/bin/activate
pip install numpy scipy astropy MySQL-python flask
 
# Install 4GP code
cd src/pythonModules/fourgp_speclib
python setup.py install
cd ../fourgp_cannon
python setup.py install
cd ../fourgp_degrade
python setup.py install
cd ../fourgp_rv
python setup.py install
cd ../fourgp_specsynth
python setup.py install
cd ../fourgp_telescope_data
python setup.py install
cd ../fourgp_fourfs
python setup.py install
 
# Create API documentation using sphinx
pip install Sphinx
cd ../../../docs
make html
 
# View HTML documentation
# At this point you need to edit your Apache configuration and point a webserver at the directory
# docs/_build/html
```

## Installing the tools which 4GP wraps

You will probably also want to install Turbospectrum and 4FS.

Your 4GP installation will need to know the paths to these tools so that it can invoke them as required. You can configure search paths whenever you invoke the wrappers for these tools.

### Installing Turbospectrum

If you want to synthesize spectra using Turbospectrum, you need to do the following. Download the code, as follows:

```
wget http://www.pages-perso-bertrand-plez.univ-montp2.fr/DATA/Turbospectrum-v15.1.tar.gz
wget http://marcs.astro.uu.se/documents/auxiliary/interpol_marcs.tar.gz
```

Proceed to build the tools as follows:

```
tar xvfz Turbospectrum-v15.1.tar.gz
mv EXPORT-15.1 turbospectrum-15.1
cd turbospectrum-15.1/exec-gf-v15.1/
make
 
tar xvfz interpol_marcs.tar.gz
cd interpol_marcs
gfortran interpol_modeles.f -o interpol_modeles
```

For more information, read the README.md file in the fourgp_specsynth directory.

### Installing 4FS

If you want to use 4FS, build it as follows:

```
git checkout https://dcf21@gitlab.4most.eu/tdwelly/OpSys.git
sudo apt-get install libhealpix-cxx-dev libchealpix-dev libcfitsio3-dev python-healpy
cd OpSys
make
```

Note that the 4MOST GitLab account is password protected, so you will need to get your own account before you will be able to check out the code.

### Installing the Cannon

At the time of writing their are various branches of Annie's Lasso with different APIs. The following branch uses Andy Casey's latest development API, which is the one 4GP expects:

```
git checkout https://github.com/dcf21/AnniesLasso.git develop
cd AnniesLasso
python setup.py install
```

# Using the web-based spectrum browser

4GP comes with a web-based tool for browsing the contents of spectrum
libraries. The tool is based on flask, and requires additional Javascript
dependencies to be installed using bower. bower is a part of nodeJS and can be
installed with the node package manager (npm).

```
apt-get install nodejs
npm update
sudo npm install -g bower
cd src/spectrumBrowser
bower install
python spectrumBrowser.py --library-path ../../../workspace
```

Flask will then tell you what web address to point your web browser: the
address is usually `http://127.0.0.1:5000`.

# Further information

Once you've followed the steps above, you should have lots of HTML
documentation of the 4GP API autogenerated by Sphinx.

# Testing your installation

The 4GP code comes with a set of unit tests to validate your installation.

The tests include building spectrum libraries using both SQLite and MySQL databases. For all the tests to pass, you need to create a local MySQL database called ```fourgp_unittest```. A MySQL user account with username ```fourgp_unittest``` and password ```fourgp_unittest``` is needed with full access to this database:
 
```
CREATE USER 'fourgp_unittest'@'localhost' IDENTIFIED BY 'fourgp_unittest';
CREATE DATABASE fourgp_unittest;
GRANT ALL ON fourgp_unittest.* TO 'fourgp_unittest'@'localhost';
```

You can run the unit tests as follows:

```
source virtualenv/bin/activate
cd src/pythonModules/fourgp_speclib/fourgp_speclib/tests
python -m unittest discover
```

# Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden
