# Architecture aware phase polynomial synthesis for NISQ devices.
This repository contains the source code to replicate the results from the paper 
"Architecture-aware phase polynomial synthesis for NISQ devices" by Arianne Meijer - van de Griend and Ross Duncan.

## Installation
To run our experiments, make a python 3.6 environment and install the packages from the requirements.txt.

### Installing Staq
Retrieve the Staq source code from <https://github.com/softwareQinc/staq> and build it as instructed.
Then copy the executable to this folder such that the python file can find it.

## Replicating results
Our results can be replicated using [the jupyter notebook](Architecture%20aware%20phase%20polynomial%20synthesis.ipynb). 
It also contains our script for generating phase polynomials and to analyse the raw metrics from the synthesised circuits.

## The proposed algorithm
The code of our proposed algorithm can be found in [phase_poly.py](phase_poly.py) the corresponding method in the `PhasePoly` class is called `Ariannes_synth()`.
