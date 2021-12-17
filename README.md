# Arbor simulation of the induction of early- and late-phase plasticity at a single synapse

This code simulates a current-based, plastic synapse that connects two leaky integrate-and-fire neurons.
 
 The following protocols are provided and can be run by the bash scripts of the same name:
 * "run\_arbor\_basic\_early": demonstrates the simulation of basic calcium and early-phase plasticity dynamics. The presynaptic neuron spikes at predefined times and thereby stimulates the postsynaptic site and induces plasticity.
 * "run\_arbor\_classical\_protocols": reproduces classical induction protocols for early- and late-phase long-term potentiation and depression. The presynaptic neuron spikes in a Poissonian manner according to one of the protocols \{STET, WTET, SLFS, WLFS\}, thereby stimulating the postsynaptic site and inducing plasticity.

The underlying synapse and neuron model has been described in detail in [Luboeinski and Tetzlaff, 2021](https://doi.org/10.1038/s42003-021-01778-y).

Please note that the early-phase plasticity dynamics in the model feature a noise term (also cf. [Graupner and Brunel, 2012](https://doi.org/10.1073/pnas.1109359109)).
Since Arbor does not yet support stochastic differential equations (see issues [\#1643](https://github.com/arbor-sim/arbor/issues/1643) and [\#1655](https://github.com/arbor-sim/arbor/issues/1655)), 
it is necessary to manually add the noise term to the C++ code that Arbor's modcc generates from the 'expsyn\_curr\_calcium\_plasticity.mod' mechanism.
To do this, the following lines can be added to 'expsyn\_curr\_calcium\_plasticity\_cpu.cpp':

```
#include <random>
#include <chrono>
std::minstd_rand0 rg; // random number generator
std::normal_distribution<double> norm_dist; // normal distribution to obtain Gaussian white noise

static unsigned int getClockSeed()
{
	static int last_seed;
	int seed;

	while ( (seed = std::chrono::system_clock::now().time_since_epoch().count()) == last_seed ) {}
	last_seed = seed;

	return seed;
}

...

// in init():
rg.seed(_s); // sets seed of the random generator, '_s' has to be replaced by the seed value
norm_dist = std::normal_distribution<double>(0.0,1.0); // initializes normal distribution with mean 0 and standard deviation 1 

...

// in advance_state():
xi *= norm_dist(rg) / pow(dt, 0.5);
```

The seed value `_s` should be chosen `0` for the 'basic\_early' simulation and `getClockSeed()` for the 'classical\_protocols' simulations.
