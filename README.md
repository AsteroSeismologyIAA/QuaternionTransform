# Asteroseismology / Quaternion

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/manuparra/asteroseismology/HEAD)

Asteroseismology, that is, the use of the frequency content of a time series caused by variations in brightness or radial velocity of a stellar object, is based on the hypothesis that such a series is harmonic and therefore can be described by a sum of sines and cosines. If this were not the case (e.g. the oscillations of an ellipsoid of revolution) it cannot be guaranteed that the Discrete Fourier transform is the least squares approximation to the time series.

This report study the effect of extending the Fourier kernel to a particular quaternion and exploring the impact when it is applied to the best time series that we have (GOLF/SoHO) from the closest star, our Sun. The results are consistent with a notable improvement in the signal-to-noise ratio in the low frequency range. This opens the possibility of detecting the elusive g modes of the Sun in future works.

![PowerSpec](./docs/powerspectrum_18_5.png)
