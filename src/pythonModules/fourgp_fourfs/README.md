# fourgp_fourfs

This python package provides a wrapper for functions provided by the 4MOST Facility Simulator, 4FS.

It requires the 4FS code to have been downloaded and compiled, and provides functions for degrading high-resolution spectra to the resolution and signal-to-noise ratios expected from 4MOST.

# Installation instructions

Build 4FS as follows:

```
git checkout https://dcf21@gitlab.4most.eu/tdwelly/OpSys.git
sudo apt-get install libhealpix-cxx-dev libchealpix-dev libcfitsio3-dev python-healpy
cd OpSys
make
```

Note that the 4MOST GitLab account is password protected, so you will need to get your own account before you will be able to check out the code.

# Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden
