# Arbor simulation of the induction of early- and late-phase plasticity at a single synapse

This code simulates a current-based, plastic synapse that connects two leaky integrate-and-fire neurons. The synapse model reproduces the outcome of established plasticity-inducing protocols by featuring early-phase plasticity based on calcium dynamics and late-phase plasticity based on synaptic tagging and capture.
 
It has been tested with [Arbor](https://arbor-sim.org/) version 0.8.1.

 The following protocols are provided and can be run by the bash scripts of the same name:
 * "run\_arbor\_basic\_early": demonstrates the simulation of basic calcium and early-phase plasticity dynamics. The presynaptic neuron spikes at predefined times and thereby stimulates the postsynaptic site and induces plasticity.
 * "run\_arbor\_classical\_protocols": reproduces classical induction protocols for early- and late-phase long-term potentiation and depression. The presynaptic neuron spikes in a Poissonian manner according to one of the protocols \{STET, WTET, SLFS, WLFS\}, thereby stimulating the postsynaptic site and inducing plasticity.

The underlying synapse and neuron model has been described in detail in [Luboeinski \& Tetzlaff, 2021](https://doi.org/10.1038/s42003-021-01778-y).

## Basic early-phase dynamics
Noise in plasticity generated from seed 0:
![Resulting basic early-phase dynamics without noise in plasticity](data_predef/arbor_2N1S_traces_0.svg)

No noise in plasticity:
![Resulting basic early-phase dynamics without noise in plasticity](data_predef_no_pl_noise/arbor_2N1S_traces_0.svg)


## Classical plasticity-inducing protocols
For details see, for example, [Sajikumar et al., 2005](https://doi.org/10.1523/JNEUROSCI.1104-05.2005) or [Luboeinski \& Tetzlaff, 2021](https://doi.org/10.1038/s42003-021-01778-y).

STET:
![Strong tetanic stimulation inducing late-phase LTP](arbor_averaged_STET.png)

WTET:
![Weak tetanic stimulation inducing early-phase LTP](arbor_averaged_WTET.png)

SLFS:
![Strong low-frequency stimulation inducing late-phase LTD](arbor_averaged_SLFS.png)

WLFS:
![Weak low-frequency stimulation inducing early-phase LTD](arbor_averaged_WLFS.png)

### Arbor installation
The code has been tested with Arbor version v0.9.1-dev, state of commit [2f4c325](https://github.com/arbor-sim/arbor/commit/2f4c32598d37f9852978c76952b0a09aeb84385b).

To install this Arbor version from source code (with SIMD support), you can run the following:
```
git clone --recursive https://github.com/arbor-sim/arbor/ arbor_source_repo
mkdir arbor_source_repo/build && cd arbor_source_repo/build
git checkout 2f4c32598d37f9852978c76952b0a09aeb84385b -b plastic_arbor_v1
cmake -DARB_WITH_PYTHON=ON -DARB_USE_BUNDLED_LIBS=ON -DARB_VECTORIZE=ON -DCMAKE_INSTALL_PREFIX=$(readlink -f ~/arbor_v0.9.1-dev-plastic_arbor_v1-simd) -DPYTHON_EXECUTABLE:FILEPATH=`which python3.10` -S .. -B .
#make tests && ./bin/unit # optionally: testing
make install
```
