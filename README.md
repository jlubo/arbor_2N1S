# Arbor simulation of the induction of early- and late-phase plasticity at a single synapse

This code simulates a current-based, plastic synapse that connects two leaky integrate-and-fire neurons. The synapse model reproduces the outcome of established plasticity-inducing protocols by featuring early-phase plasticity based on calcium dynamics and late-phase plasticity based on synaptic tagging and capture.

It requires Arbor v0.6 or later with https://github.com/boeschf/arbor/tree/SDE (as of commit 1b66a602bb3efef955ea441f02807c049d8c63c4).
 
 The following protocols are provided and can be run by the bash scripts of the same name:
 * "run\_arbor\_basic\_early": demonstrates the simulation of basic calcium and early-phase plasticity dynamics. The presynaptic neuron spikes at predefined times and thereby stimulates the postsynaptic site and induces plasticity.
 * "run\_arbor\_classical\_protocols": reproduces classical induction protocols for early- and late-phase long-term potentiation and depression. The presynaptic neuron spikes in a Poissonian manner according to one of the protocols \{STET, WTET, SLFS, WLFS\}, thereby stimulating the postsynaptic site and inducing plasticity.

The underlying synapse and neuron model has been described in detail in [Luboeinski and Tetzlaff, 2021](https://doi.org/10.1038/s42003-021-01778-y).

The `SEED` value in "expsyn_curr_calcium_plasticity.mod" should be chosen `0` for the deterministic "basic\_early" simulation and be left out for the "classical\_protocols" and other simulations. Following this, the script "build\_cat\_only" can be used to recompile the catalogue from the already generated C++ code.
