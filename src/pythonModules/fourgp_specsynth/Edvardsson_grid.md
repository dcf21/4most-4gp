# The Edvarddson grid of MARCS models

The Edvardsson grid of MARCS models contains 51,994 model atmospheres:

```
dcf21@thebe:~/iwg7_pipeline/fromBengt/marcs_grid$ ls | wc
  51994   51994 4003538
```

The filename of each model enhances the following parameters: effective temperatures, surface gravity, mass,
microturbulence (km/s), metallicity class, metallicity, alpha enrichment, r-process elements, s-process elements.

The coverage of parameter space can be queried as follows:

```
python
>>> from fourgp_specsynth import TurboSpectrum
>>> x = TurboSpectrum()
>>> x.marcs_values
```

```
{'spherical': ['p', 's'],
 'metallicity': [-5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.05, 0.25, 0.5, 0.75, 1.0],
 'turbulence': [0.0, 1.0, 2.0, 5.0],
 'log_g': [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5],
 'a': [-0.4, 0.0, 0.1, 0.11, 0.2, 0.3, 0.4],
 'c': [-0.38, -0.13, 0.0, 0.08],
 'temperature': [2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0,
                 3700.0, 3800.0, 3900.0, 4000.0, 4250.0, 4500.0, 4750.0, 5000.0, 5250.0, 5500.0, 5750.0, 6000.0,
                 6250.0, 6500.0, 6750.0, 7000.0, 7250.0, 7500.0, 7750.0, 8000.0],
 'o': [-0.4, 0.0, 0.1, 0.12, 0.2, 0.3, 0.4],
 'n': [0.0, 0.09, 0.31, 0.53],
 's': [0.0],
 'r': [0.0],
 'mass': [0.0, 0.5, 1.0, 2.0, 5.0, 15.0],
 'model_type': ['ae', 'an', 'ap', 'gs', 'hc', 'mc', 'st']
 }
```

The model types are metallicity classes, e.g. alpha enhanced, alpha poor, carbon enhanced. Distribution of the number
of models per type:

```
ae:  7797
an: 11945
ap:  7516
gs:  1140
hc:  3985
mc:  3931
st: 15680
```

The distribution of number of models per mass is:

```
 0.0: 18209
 0.5:  3755
 1.0: 21773
 2.0:  3978
 5.0:  3999
15.0:   280
```

The distribution of number of models per turbulence setting is:

```
0.0:  4518
1.0:  7077
2.0: 23493
5.0: 16906
```