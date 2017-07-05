# fourgp_specsynth

This python package wraps spectral synthesis tools.

Currently TurboSpectrum is the only tool which is supported.

# Installation instructions

To use this wrapper, you first need to install Turbospectrum and interpol - a tool used to interpolate MARCS models. You can find more information about these tools on their respective websites:

```
http://www.pages-perso-bertrand-plez.univ-montp2.fr
http://marcs.astro.uu.se/software.php
```

First download the code, as follows:

```
wget http://www.pages-perso-bertrand-plez.univ-montp2.fr/DATA/Turbospectrum-v15.1.tar.gz
wget http://marcs.astro.uu.se/documents/auxiliary/interpol_marcs.tar.gz
```

Proceed to build the tools as follows:

```
sudo apt-get install gfortran
 
tar xvfz Turbospectrum-v15.1.tar.gz
mv EXPORT-15.1 turbospectrum-15.1
cd turbospectrum-15.1/exec-gf-v15.1/
make
 
tar xvfz interpol_marcs.tar.gz
cd interpol_marcs
gfortran interpol_modeles.f -o interpol_modeles
```

# Contact details
This code is maintained by:

Dominic Ford  
Lund Observatory  
Box 43  
SE-221 00 Lund  
Sweden
