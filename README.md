TLE Estimator by GPS Measurements
=================================

Simple batch filter estimation algorithm via SVD to determine the orbit of a satellite by GPS single-point positioning measurements.
Code provides the following features:
* Read GPS receiver telemetry files.
* Batch filter to calculate estimate position and velocity of satellite at initial measurement epoch.
* Convert position and velocity in xyz to Kepler orbital parameters.
* Write and read satellite TLE files based on estimated ephemeris.

How to run
----------

This code is developed and tested on Python3. There is a number of class dependencies:
* numpy
* pandas
* pyatmos
* spiceypy
* numdifftools
* scipy
* pyshtools

These may be simply installed using:

`pip3 install numpy pandas pyatmos spiceypy numdifftools scipy pyshtools`

An example file is included  _main_example.py_  that can be run to demonstrate use of the BlockEstimator class.

The Spice library is used to calculate positions of the Moon and Sun, as well as the
orientation of the Earth. The  _kernel_  folder contains a number of required files, including:
* [de440.bsp](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/)
* [earth_assoc_itrf93.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/)
* [earth_latest_high_prec.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/)
* [naif0012.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/)

Each file must be downloaded and placed into the kernel folder before running. The location is set using  _batch_estimator.tm_ .

Files are usually contained in the  _data_  folder.