# Ambient Noise Generator & Evaluation
Last modified 23 April 2019

**1. Research Subject**

The aim of this software is to generate synthetic ambient noise in physical dimension and evaluate the cross-correlation function for the structural monitoring.

**2. Specifications**

* Generate noise waveform using representations of ambient noise.
* Incorporating scattering, attenuation, disspersion, variety of source wavelet, Gaussian noise.
* Automatic evaluation of cross-correlation (1st & 2nd order, full wave and coda wave) and its signal to noise ratio
* Flexible source type, source locations and receiver locations
* Selectable outputs (noise waveform, cc function, beamforming analysis, Signal to noise ratio)
* Available for virtual dv/v analysis

**3. Way to run the simulation**

Please go through the work-flow as following:

**Make model**

* run `source_and_receiver.m` to prescribe source and receiver location
* run `model_structure.m` to prescribe structural information such as
	- phase velocity
	- dispersion curve
	- attenuation
* run `make_config` to prescribe model parameters

**Run Simulation**

  
