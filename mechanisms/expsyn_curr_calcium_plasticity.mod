: Exponential, current-based, synapse with early-phase plasticity based on calcium dynamics
: and late-phase plasticity based on synaptic tagging and capture

NEURON {
	POINT_PROCESS expsyn_curr_calcium_plasticity
	RANGE h_0, tau_syn, Ca_pre, Ca_post, theta_p, theta_d, theta_pro, theta_tag, R_mem, sigma_pl
	NONSPECIFIC_CURRENT I
}

UNITS {
	(ms) = (milliseconds)
	(mV) = (millivolt)
	(MOhm) = (megaohm)
}

PARAMETER {
	R_mem = 10.0 (MOhm) : membrane resistance
	tau_syn = 5.0 (ms) : synaptic time constant
	Ca_pre = 1.0 : pre-synaptic calcium contribution
	Ca_post = 0.2758  : post-synaptic calcium contribution

	tau_h = 688400 (ms) : early-phase time constant
	tau_z = 3600000 (ms) : late-phase time constant
	tau_Ca = 48.8 (ms) : calcium time constant
	tau_p = 3600000 (ms) : protein time constant
	h_0 = 4.20075 (mV) : initial weight
	gamma_p = 1645.6 : potentiation rate
	gamma_d = 313.1 : depression rate
	t_Ca_delay = 18.8 (ms) : delay of postsynaptic calcium increase after presynaptic spike
	theta_p = 3.0 : calcium threshold for potentiation
	theta_d = 1.2 : calcium threshold for depression
	sigma_pl = 2.90436 (mV) : standard deviation for plasticity fluctuations
	alpha = 1 : protein synthesis rate
	theta_pro = 2.10037 (mV) : protein synthesis threshold
	theta_tag = 0.840149 (mV) : tagging threshold
}

STATE {
	g (mV) : instantaneous synaptic conductance
	w (mV) : total synaptic weight
	h (mV) : early-phase weight
	z : late-phase weight
	p : protein amount
	Ca : calcium amount
}

INITIAL {
	g = 0
	w = h_0
	h = h_0
	z = 0
	p = 0
	Ca = 0
}

BREAKPOINT {
	:SOLVE state METHOD cnexp : solver not compatible with late-phase equation
	:SOLVE state METHOD sparse
	SOLVE state METHOD stochastic
	
	w = h + z*h_0
	I = -g / R_mem
}

WHITE_NOISE {
    zeta
}

DERIVATIVE state {
	LOCAL xi
	
	: Exponential decay of postsynaptic potential
	g' = -g / tau_syn
	
	: Early-phase dynamics
	xi = ( tau_h * (heaviside(Ca - theta_p) + heaviside(Ca - theta_d)) )^(1/2) * sigma_pl * zeta
	h' = ( 0.1 * (h_0 - h) + gamma_p * (10-h) * heaviside(Ca - theta_p) - gamma_d * h * heaviside(Ca - theta_d) + xi ) / tau_h
	
	: Late-phase dynamics
	z' = ( p * (1 - z) * heaviside( (h-h_0) - theta_tag ) - p * (z + 0.5) * heaviside( (h_0-h) - theta_tag ) ) / tau_z
		
	: Protein dynamics
	p' = ( -p + alpha * heaviside( abs(h-h_0) - theta_pro ) ) / tau_p 
	
	: Exponential decay of calcium concentration
	Ca' = -Ca/tau_Ca
}

NET_RECEIVE(weight) {
	if (weight >= 0) {
		: Start of postsynaptic potential
		g = g + w
	}
	else {
		: Increase of calcium amount by presynaptic spike
		Ca = Ca + Ca_pre
	}
}

POST_EVENT(time) {
	: Increase of calcium amount by postsynaptic spike
	Ca = Ca + Ca_post
}

FUNCTION heaviside(x) { : the Heaviside Theta function
	if (x >= 0) {
		heaviside = 1
	}
	else {
		heaviside = 0
	}
}

FUNCTION abs(x) { : the absolute value of a number
	if (x >= 0) {
		abs = x
	}
	else {
		abs = -x
	}
}
