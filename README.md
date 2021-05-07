# ries: resonances integrated over energy and space

A python library for the estimation of photonuclear reaction rates in low-energy nuclear physics experiments.

## Introduction

> When I increase my sample ('target') size by a factor of two, my particle detectors will record twice as many reactions.

If you fully agree to this statement, you have a rough idea of the scaling laws in beam-on-target experiments, but you may have overlooked some effects that make photonuclear experiments both challenging and exciting.
The `ries` library is all about these finite-temperature, thick-target, isolated-resonance, multiple-scattering ... effects that call for a reconsideration of statements like the one above.
When utilized in a clever way, second-order effects can even be the basis for more efficient experiments.

In fact, now is the right time to get involved, because currently existing [1,2] and planned [3] intense, polarized quasi-monoenergetic photon sources ('gamma-ray lasers') in the mega-electronvolt energy range have led to a 'renaissance' of photonuclear reaction studies.

## Description

The object-orientated and modular library `ries` allows users to simulate the passage of a photon beam through a sample in a nuclear physics experiment.
It aims to bridge the gap between pen-and-paper calculations and full-scale Monte Carlo particle simulations like Geant4 [4] or MCNP [5]:
Precise enough for serious work, and fast enough to act as the cost function of an optimization problem or the data model in a Bayesian analysis.

`ries` uses `ℏc = 197.3269804 MeV` fm for the product of the reduced Planck constant and the speed of light.
This means that all energies are assumed to be given in Mega electron volts, and all lengths in femtometers.
Temperatures are assumed to be given in Kelvin.

## Prerequisites

`ries` is [python3](https://www.python.org/) code and requires the following packages:

* [numpy](https://numpy.org/)
* [setuptools](https://setuptools.readthedocs.io/)
* [scipy](https://www.scipy.org/)

The [tox](https://tox.readthedocs.io/) tool is used to run self tests, build the documentation, and check whether the code is in standard format.
This requires the following packages:

* [black](https://black.readthedocs.io)
* [matplotlib](https://matplotlib.org/)
* [pytest](https://docs.pytest.org/)
* [pytest-cov](https://pytest-cov.readthedocs.io/)
* [sphinx](https://www.sphinx-doc.org/)
* [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/)
* [tox](https://tox.readthedocs.io/) 

## Installation

Clone the repository:

```
git clone https://github.com/uga-uga/ries.git
```

Assuming that the `setup.py` file of `ries` is located in `RIES_DIR`, the code can be installed by executing:

```
cd $RIES_DIR
python setup.py install
```

To run self tests and build the documentation, execute

```
tox
```

in the same directory.
By default, the documentation will be generated in html format.
It can be found in `$RIES_DIR/build` and opened in a web browser.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

Copyright (C) 2020, 2021 Udo Friman-Gayer (ufg@email.unc.edu)

## References

[1] High-Intensity Gamma-Ray Source (HIγS) at TUNL [https://tunl.duke.edu/research/our-facilities](https://tunl.duke.edu/research/our-facilities), accessed on 12/02/2020

[2] NewSUBARU at University of Hyogo [https://www.lasti.u-hyogo.ac.jp/NS-en/newsubaru/](https://www.lasti.u-hyogo.ac.jp/NS-en/newsubaru/), accessed on 12/02/2020

[3] Variable Energy Gamma (VEGA) System at ELI-NP [https://www.eli-np.ro/rd2.php](https://www.eli-np.ro/rd2.php), accessed on 12/02/2020

[4] [https://geant4.web.cern.ch/](https://geant4.web.cern.ch/), accessed on 12/02/2020

[5] [https://mcnp.lanl.gov/](https://mcnp.lanl.gov/), accessed on 12/02/2020
