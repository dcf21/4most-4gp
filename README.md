**Please note that this is Dominic Ford's personal copy of the 4GP code. It is basically a frozen copy of how the code stood in April 2019, when Dominic stopped working for 4MOST. The up-to-date code is in the official 4MOST IWG7 repository here: [https://github.com/4most-iwg7/4most-4gp](https://github.com/4most-iwg7/4most-4gp)**

The text on this page is a brief introduction to `4most-4gp`. For a more
complete tutorial, please visit the
[Wiki](https://github.com/dcf21/4most-4gp/wiki).

# 4most-4gp

4GP is the Galactic pipeline for the 4MOST multi-object
spectrograph, developed by Infrastructure Working Group 7 (IWG7).

The pipeline comprises a collection of Python modules which use a common
data format to store and manipulate spectra and their associated metadata. It
makes it easy to pass spectra between a range of spectral synthesis, processing and
analysis tools including Turbospectrum, the 4MOST Facility Simulator (4FS),
and abundance analysis codes such as the Cannon, without the need for manual data format conversion.

It includes
the ability to store spectra in libraries and search for them by arbitrary
metadata constraints, making it easy to create new tests on subgroups of stars
filtered from larger samples.

In addition, 4GP includes a simple web interface which allows the contents of spectrum libraries
to be searched and viewed quickly for diagnostic purposes.

### Code structure

The 4GP framework is available in two repositories on GitHub, and these Wiki pages provide
step-by-step installation instructions. The first repository contains the
Python modules which provide programmatic interfaces for creating and
manipulating libraries of spectra. It includes wrappers for passing them to
various analysis tools:

<https://github.com/dcf21/4most-4gp>

The second repository contains python scripts which utilise these modules to
create command-line tools for
synthesising spectra, manipulating them, and then testing abundance analysis codes such as the
Cannon and the Payne on them:

<https://github.com/dcf21/4most-4gp-scripts>

### Getting started

This code is under active development, but stable releases are periodically
made.

Visiting the GitHub URLs above will present you with the `master` branch of our
code, which should always correspond to the latest stable release. If you click
on the "branches" dropdown menu, you can select a different version of the code
to download.

Stable releases are given date stamps, for example, `release-2019-03-01-1`. The
master branch points to the most recent release. The `dev` branch may contain experimental code
and should be used with extreme caution.

## Code structure

4GP is organised into a collection of python modules, most of which do not
depend on each other. These can be found in the directory `src/pythonModules`.

The code is divided up in this way since each module has different dependencies,
and it allows various parts of 4GP to be used without installing every dependency.

Before you start, you will need to install these python modules into your python environment.

It is recommended that they be installed in a python virtual
environment as described in the [installation instructions](installing_4gp),
and that you do not install them into your
system-wide python installation.
Alternatively, you can use a python distribution such as Conda to create your own python environment.

The python modules that make up 4GP are as follows:

**fourgp_speclib** - Core classes representing spectra and libraries to keep them in. These include functionality to search for spectra within libraries based on arbitrary metadata constraints. Behind the scenes, this metadata is stored in a SQL database, and both SQLite (recommended) and MySQL are supported. All other 4GP modules depend on this core module.

**fourgp_specsynth** - A wrapper for synthesising spectra with given stellar parameters using Turbospectrum.

**fourgp_fourfs** - A wrapper for create mock 4MOST observations of spectra using 4FS; specifically, using the 4MOST Exposure Time Calculator (ETC). This reduces the resolution of the input spectra to match 4MOST's instrumental profile, and also adds noise to the spectra.

**fourgp_cannon** - A wrapper for the Cannon, a machine learning framework for deriving stellar parameters and abundances from spectra.

**fourgp_payne** - A wrapper for the Payne, a machine learning framework for deriving stellar parameters and abundances from spectra.

**fourgp_degrade** - Some simple functions for convolving spectra with custom instrumental profiles, for resampling them onto new wavelength rasters, and for adding simple noise to spectra.

**fourgp_telescope_data** - A container for basic instrumental data about 4MOST, including its wavelength coverage and spectral resolution.

**fourgp_rv** - A module which implements a simple cross-correlation algorithm for determining the radial velocities of stars (based loosely on the GUESS code used by GALAH).

**fourgp_pipeline** - A module which implements a simple pipeline for analysing FGK stars. The example pipeline uses the `fourgp_rv` module to determine each star's RV, before using the Cannon to do abundance analysis. The pipeline is designed to be easily extensible, so that additional tasks can be added to the work flow, or existing tasks replaced with others.

## Getting 4GP to do useful things

The code in this repository (i.e. `4most-4gp`) is simply a collection of python modules
which provide a programmatic interface for handling spectra. There are no command-line
interfaces or example code here to show how to use these modules in practice. To use this
code, you need a python script to invoke these modules.

However, in the separate repository [4most-4gp-scripts](https://github.com/dcf21/4most-4gp-scripts) you will
find lots of command line interfaces to `4most-4gp` which probably already do most of the tasks
that you're likely to need to do.

## Installing 4GP

For complete installation instructions, please visit the
[Wiki](https://github.com/dcf21/4most-4gp/wiki/dependencies) on our GitHib
page.

### Dependencies

Before installing 4GP, you need to make sure that you have all the software packages which are required to build Turbospectrum, 4FS, and the other tools that 4GP wraps.

4GP has currently only been tested on python 3.5, and so you will need to have access to this version of Python.

In addition, the following external packages and libraries are required:

* **git** - required to check the code out from GitHub
* **SQLite3** - including the python-sqlite3 binding; you can test for this by typing `import sqlite3` into a python terminal
* **python-matplotlib** - required to use the 4GP Spectrum Browser and the Cannon; you can test for this by typing `import matplotlib` into a python terminal
* **python-tk** - required to use the 4GP Spectrum Browser and the Cannon; you can test for this by typing `import tkinter` into a python terminal
* **pyxplot** - required to produce plots of the Cannon's performance

The following packages are strongly recommended:

* **python-virtualenv** - required to set up a python virtual environment

The following packages are needed to run certain parts of 4GP:

* **libchealpix-dev** - required to build and install 4FS
* **python-healpy** - required to build and install 4FS
* **libcfitsio3-dev** - required to build and install 4FS
* **gfortran** - required to build and install Turbospectrum
* **pyphot** - required to do photometry on spectra
* **MySQL** - including the libmysqlclient and python-mysql bindings - currently only required by the 4GP unit tests, so not very important

## Ubuntu installation instructions

Under Ubuntu Linux, you can install all these packages with the in-built package manager, as follows:

#### Ubuntu 18.04

```
apt-get install git python-sqlite mysql-server libmysqlclient-dev python-virtualenv libhealpix-cxx-dev libchealpix-dev libcfitsio-dev python-healpy gfortran python-tk python-matplotlib sqlite3 python3-dev
```

The visualisation scripts use a plotting package called `Pyxplot`, which is unfortunately not packaged within Ubuntu 18.04. You can download it and build it from source here:

https://github.com/dcf21/pyxplot9

## Installing 4GP in a python virtual environment

4GP is distributed with a standard `setuptools` installation script `setup.py` which allows you to install its
constituent modules into your local python environment.

Note that 4GP is currently only tested against python 3.5.

We recommended that you install them in a python virtual
environment, rather than tampering with your
system-wide python installation.

Follow these steps in a Linux shell to do this:

```
# Check out code from GitHub
git clone https://github.com/dcf21/4most-4gp.git

# Sometimes this line is necessary, if your locale settings are broken
export LC_ALL=C
 
# Set up a python 3 virtual environment
virtualenv -p python3 virtualenv
source virtualenv/bin/activate

# Install some of the python packages we required
pip install numpy scipy astropy mysqlclient flask matplotlib tables

# Install 4GP code
cd 4most-4gp/src/pythonModules/fourgp_speclib
python setup.py develop
cd ../fourgp_cannon
python setup.py develop
cd ../fourgp_degrade
python setup.py develop
cd ../fourgp_rv
python setup.py develop
cd ../fourgp_specsynth
python setup.py develop
cd ../fourgp_telescope_data
python setup.py develop
cd ../fourgp_fourfs
python setup.py develop
cd ../fourgp_payne
python setup.py develop
cd ../fourgp_pipeline
python setup.py develop
```

### 4GP API documentation

4GP includes a set of scripts which produce auto-generated HTML documentation of its programmatic API.

These use a tool called Sphinx to automatically extract comments from the 4GP source code, and turn them into documentation of the methods which are available within each Python class.

After running the commands above, you can create the API documentation as follows:

```
# Create API documentation using sphinx
pip install Sphinx
cd ../../../docs
make html
```

After doing this, a directory `_build` will have appeared within the `docs` directory containing HTML documentation. To view its contents, you should open a web browser and enter the address:

```file:///path_to_your_4gp_installation/docs/_build/index.html```

### Installing the tools which 4GP wraps

Depending what tasks you envisage doing with 4GP, you will almost certainly need to install some of the software
packages it provides wrappers for - for example, the Cannon, Turbospectrum or 4FS.

Note that you do not need to install tools that you are not going to use.

### Installation paths

When you install these tools, 4GP needs to be able to find them.

By far the simplest approach is if you install each tool in the same directory where you keep your working copies of `4most-4gp` and `4most-4gp-scripts`. Thus, after installing various tools, this directory might look as follows:

```bash
dcf21@astrolabe:~/iwg7_pipeline$ ls
4most-4gp                  downloads         idl_packages    rvspecfit
4most-4gp-scripts          forwardModelling  interpol_marcs  sme
4most-iwg7-pipeline-tests  fromBengt         OpSys           TheCannon
4MOST_testspectra          fromKeith         pepsi           turbospectrum-15.1
AnniesLasso                hot_stars         pyphot          virtualenv
```

This is the default place where 4GP looks to find each tool, and if they are installed in this way then you will not need to explicitly tell it where to look.

You are permitted to install these tools in a different location, but then you will need to pass configuration
parameters to a number of 4GP's modules (usually called `binary_path`) to tell it where to find each tool.

### Installing the Cannon

At the time of writing, there are various branches of the Cannon (also known as Annie's Lasso), and each has its
own different API. To avoid confusion, we have our own 4GP [GitHub repository](https://github.com/dcf21/AnniesLasso) containing the version of the Cannon
that we use. This is a fork of a version released by Andy Casey in 2007.

Within this repository, there are branches named after each release of 4GP, e.g. `release-2019-03-01-1`. These
contain the versions of the Cannon expected by each release of 4GP with the same name. The `master` branch
always represents the most recent stable release.  

```
git clone https://github.com/dcf21/AnniesLasso.git
cd AnniesLasso
python setup.py install
```

Do not use the `dev` branch in this repository. It contains the latest release of the Cannon by Andy Casey,
which has some bugs which are unresolved at the time of writing. It produces *worse* results than before.

### Installing pyphot


The pyphot source code can be obtained from a GitHub repository. As with the Cannon,
we maintain our own fork of the repository with branches labelled `release-2018-01-12-1`,
etc, to indicate versions which are compatible with each release of 4GP.

```
git clone https://github.com/dcf21/pyphot.git
cd pyphot
python setup.py install
```

### Installing Turbospectrum

If you want to synthesize spectra using Turbospectrum, the following commands will download and install it for you.

Download the code, as follows:

```
# Fetch the Turbospectrum source code
wget http://www.pages-perso-bertrand-plez.univ-montp2.fr/DATA/Turbospectrum-v15.1.tar.gz

# Fetch the interpol source code, which is used to interpolate MARCS models
wget http://marcs.astro.uu.se/documents/auxiliary/interpol_marcs.tar.gz
```

Proceed to build the tools as follows:

```
# Unpack Turbospectrum and build it from source
tar xvfz Turbospectrum-v15.1.tar.gz
mv EXPORT-15.1 turbospectrum-15.1
cd turbospectrum-15.1/exec-gf-v15.1/
make

# Unpack interpol and build it from source
cd ../..
tar xvfz interpol_marcs.tar.gz
cd interpol_marcs
gfortran interpol_modeles.f -o interpol_modeles
```

For more information, read the README.md file in the fourgp_specsynth directory.

You will also need to obtain a copy of the MARCS grid of model atmospheres, and a suitable line list for the
4MOST spectral range. Your best bet is to ask us to supply you with this data, which is stored in the following
directory on a computer in Lund:

```rsync -av dcf21@astrolabe.astro.lu.se:iwg7_pipeline/fromBengt .```

### Installing 4FS

If you want to use 4FS, build it as follows:

```
git checkout https://dcf21@gitlab.4most.eu/tdwelly/OpSys.git
sudo apt-get install libhealpix-cxx-dev libchealpix-dev libcfitsio3-dev python-healpy
cd OpSys
make
```

Note that the 4MOST GitLab account is password protected, so you will need to get your own account before you will be able to check out the code.

## Using the web-based spectrum browser

4GP comes with a simple web-based tool for browsing the contents of spectrum
libraries.

This is very handy as it lets you very quickly search for spectra, view graphs of their flux vs wavelength, and also export spectra as text files for use in other tools.

### Installation

The tool is based on python / flask, which is a simple framework for hosting websites from within a short python script.

It requires some additional Javascript dependencies, which for simplicity we supply in a `.tar.gz` file. You can extract this as follows:

```
cd 4most-4gp/src/spectrumBrowser/static
tar xvfz vendor_code.tar.gz
```

### Running the browser

Once you have all the dependencies installed, you can start the browser as follows:

```
python spectrumBrowser.py --library-path ../../../workspace
```

Once you start the python script,
Flask will then tell you what web address to point your web browser: the
address is usually `http://127.0.0.1:5000`.

## Further information

Once you've followed the steps above, you should have lots of HTML
documentation of the 4GP API autogenerated by Sphinx.

## Testing your installation

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

## Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden

