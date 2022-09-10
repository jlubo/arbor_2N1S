#!/bin/python3

# Arbor simulation of a single current-based synapse connecting two leaky integrate-and-fire neurons.
# Depending on the neuronal activities, the synapse undergoes calcium-based early-phase plasticity,
# or late-phase plasticity described by synaptic tagging and capture.

import arbor
import numpy as np
from plotResults import plotResults
import time
from datetime import datetime
import gc
import json
import argparse
import os

#####################################
# SingleRecipe
# Implementation of Arbor simulation recipe
class SingleRecipe(arbor.recipe):

	# constructor
	# config: dictionary containing configuration data
	# learn_prot: protocol for the learning stimulus
	# delta_t: duration of one timestep in ms
	def __init__(self, config, learn_prot, delta_t):

		# The base C++ class constructor must be called first, to ensure that
		# all memory in the C++ class is initialized correctly. (see https://github.com/tetzlab/FIPPA/blob/main/STDP/arbor_lif_stdp.py)
		arbor.recipe.__init__(self)
		
		self.ncells = 2 # employ two neurons

		self.the_props = arbor.neuron_cable_properties() # initialize the cell properties to match Neuron's defaults 
		                                                 # (cf. https://docs.arbor-sim.org/en/v0.5.2/tutorial single_cell_recipe.html)	
		
		the_cat = arbor.load_catalogue("./custom-catalogue.so") # load the catalogue of custom mechanisms
		the_cat.extend(arbor.default_catalogue(), "") # add the default catalogue
		self.the_props.catalogue = the_cat
		  
		self.the_probes = [arbor.cable_probe_membrane_voltage('"center"'), \
		                   arbor.cable_probe_total_ion_current_density('"center"'), \
		                   arbor.cable_probe_total_ion_current_cell(), \
		                   arbor.cable_probe_point_state_cell("expsyn_curr_calcium_plasticity", "Ca"), \
		                   arbor.cable_probe_point_state_cell("expsyn_curr_calcium_plasticity", "h"), \
		                   arbor.cable_probe_point_state_cell("expsyn_curr_calcium_plasticity", "z"), \
		                   arbor.cable_probe_point_state_cell("expsyn_curr_calcium_plasticity", "p")] # set the probes 

		self.neuron_config = config["neuron"]
		self.syn_config_exc = config["synapses"]["syn_exc_calcium_plasticity"]
		self.syn_exc_input0_config = config["synapses"]["syn_exc_input0"]
		self.syn_exc_input1_config = config["synapses"]["syn_exc_input1"]
		self.learn_prot = learn_prot
		self.dt = delta_t
	
	# cell_kind
	# gid: global identifier of the cell
	# return: type of the cell
	def cell_kind(self, gid):
		
		return arbor.cell_kind.cable # note: implementation of arbor.cell_kind.lif is not ready to use yet

	# cell_description
	# gid: global identifier of the cell
	# return: description of the cell
	def cell_description(self, gid):
		assert gid == 0 or gid == 1
		
		# cylinder morphology
		tree = arbor.segment_tree()
		radius = self.neuron_config["radius"] # radius of cylinder (in µm)
		height = 2*radius # height of cylinder (in µm)
		tree.append(arbor.mnpos,
		            arbor.mpoint(-height/2, 0, 0, radius),
		            arbor.mpoint(height/2, 0, 0, radius),
		            tag=1)
		labels = arbor.label_dict({'center': '(location 0 0.5)'})
		area_m2 = 2 * np.pi * (radius * 1e-6) * (height * 1e-6) # surface area of the cylinder in m^2 (excluding the circle-shaped ends, since Arbor does not consider current flux there)
		area_cm2 = 2 * np.pi * (radius * 1e-4) * (height * 1e-4) # surface area of the cylinder in cm^2 (excluding the circle-shaped ends, since Arbor does not consider current flux there)
		i_factor = (1e-9/1e-3) / area_cm2 # conversion factor from nA to mA/cm^2; for point neurons
		c_mem = self.neuron_config["C_mem"] / area_m2 # specific capacitance in F/m^2, computed from absolute capacitance of a point neuron
		
		# cell mechanism
		decor = arbor.decor()
		decor.set_property(Vm=self.neuron_config["V_init"], cm=c_mem)
		neuron = arbor.mechanism(self.neuron_config["mechanism"])
		V_th = self.neuron_config["V_th"]
		V_reset = self.neuron_config["V_reset"]
		R_leak = self.neuron_config["R_leak"]
		tau_mem = R_leak*10**9 * self.neuron_config["C_mem"] # membrane time constant in ms
		neuron.set("V_th", V_th)
		neuron.set("V_reset", V_reset)
		neuron.set("R_reset", self.neuron_config["R_reset"])
		neuron.set("R_leak", R_leak)
		neuron.set("I_0", self.neuron_config["I_0"])
		neuron.set("i_factor", i_factor)
		neuron.set("t_ref", self.neuron_config["t_ref"])
		decor.paint('(all)', arbor.density(neuron))
			
		if gid == 0:
			# output information
			#print("area =", area_m2, "m^2")
			#print("i_factor =", i_factor, "(mA/cm^2) / (nA)")
			#print("c_mem =", c_mem, "F/m^2")
			#print("tau_mem =", tau_mem, "ms")
		
			# plastic excitatory exponential synapse
			mech_expsyn_exc = arbor.mechanism('expsyn_curr_calcium_plasticity')
			mech_expsyn_exc.set('h_0', self.syn_config_exc["h_0"])
			mech_expsyn_exc.set('R_mem', R_leak)
			mech_expsyn_exc.set('tau_syn', self.syn_config_exc["tau_syn"])
			mech_expsyn_exc.set('Ca_pre', self.syn_config_exc["Ca_pre"])
			mech_expsyn_exc.set('Ca_post', self.syn_config_exc["Ca_post"])
			mech_expsyn_exc.set('theta_p', self.syn_config_exc["theta_p"])
			mech_expsyn_exc.set('theta_d', self.syn_config_exc["theta_d"])
			mech_expsyn_exc.set('theta_pro', self.syn_config_exc["theta_pro"])
			mech_expsyn_exc.set('theta_tag', self.syn_config_exc["theta_tag"])
			decor.place('"center"', arbor.synapse(mech_expsyn_exc), "syn_exc_calcium_plasticity")
			
			# additional excitatory delta synapse
			mech_deltasyn_exc = arbor.mechanism('deltasyn')
			mech_deltasyn_exc.set('g_spike', 100*(V_th-V_reset)*np.exp(self.dt/tau_mem)) # choose sufficently large increase in conductance
			decor.place('"center"', arbor.synapse(mech_deltasyn_exc), "syn_exc_input0")
			
			# place spike detector
			decor.place('"center"', arbor.spike_detector(V_th), "spike_detector0")
			
		elif gid == 1:
			# excitatory delta synapse
			mech_deltasyn_exc = arbor.mechanism('deltasyn')
			mech_deltasyn_exc.set('g_spike', 100*(V_th-V_reset)*np.exp(self.dt/tau_mem)) # choose sufficently large increase in conductance
			decor.place('"center"', arbor.synapse(mech_deltasyn_exc), "syn_exc_input1")
			
			# place spike detector
			decor.place('"center"', arbor.spike_detector(V_th), "spike_detector1")
			
		return arbor.cable_cell(tree, labels, decor)
		
	# connections_on
	# gid: global identifier of the cell
	# return: connections to the given neuron
	def connections_on(self, gid):
		
		if gid == 0: # add connection from neuron 1 to neuron 0
			src = 1 # gid of the source neuron
			d0 = self.syn_config_exc["t_ax_delay"] # delay time of the postsynaptic potential in ms
			d1 = self.syn_config_exc["t_Ca_delay"] # delay time of the calcium increase in ms
			return [arbor.connection((src,"spike_detector1"), "syn_exc_calcium_plasticity", 1, d0), \
			        arbor.connection((src,"spike_detector1"), "syn_exc_calcium_plasticity", -1, d1)] # the weight argument is used to indicate which dynamics is targeted (>=0: voltage, <0: calcium)
		else:
			return []
	
	# event_generators
	# gid: global identifier of the cell
	# return: events generated from Arbor schedule
	def event_generators(self, gid):
			
		if gid == 0:
			stimulus_times_exc = self.syn_exc_input0_config["stimulus_times"]
			
			# create event generator
			spike_exc = arbor.event_generator(
				"syn_exc_input0",
				0.,
				arbor.explicit_schedule(stimulus_times_exc))
				
			return [spike_exc]
				
		elif gid == 1:

			if self.learn_prot == "STET":
				t_start = np.linspace(3600, 4800, num=3, endpoint=True) # start times of the pulses (in s)
				t_end = np.linspace(3601, 4801, num=3, endpoint=True) # end times of the pulses (in s)
				freq = 100 # average spike frequency in Hz
				seed = int(datetime.now().timestamp() * 1e6)

				# average number of spikes (random number drawn for every timestep, then filtered with probability):
				stimulus_times_exc = np.array([])
				rng = np.random.default_rng(seed)
				num_timesteps = np.int_(np.round_((t_end[0]-t_start[0])*1000./self.dt))
				for i in range(len(t_start)):
					spike_mask = rng.random(size=num_timesteps) < freq*self.dt/1000.
					timestep_values = np.linspace(t_start[i]*1000, t_end[i]*1000, num=num_timesteps, endpoint=False)
					spikes = timestep_values[spike_mask]
					stimulus_times_exc = np.concatenate([stimulus_times_exc, spikes])
						      
				#np.savetxt('arbor_2N1S_stimulation_' + self.learn_prot + '.txt', stimulus_times_exc, fmt="%i") # to check the specific stimulation
				
			elif self.learn_prot == "WTET":
				t_start = 3600 # start time in s
				t_end = 3600.2 # end time in s
				freq = 100 # average spike frequency in Hz
				seed = int(datetime.now().timestamp() * 1e6)
				
				# average number of spikes (random number drawn for every timestep, then filtered with probability):
				rng = np.random.default_rng(seed)
				num_timesteps = np.int_(np.round_((t_end-t_start)*1000./self.dt))
				spike_mask = rng.random(size=num_timesteps) < freq*self.dt/1000.
				timestep_values = np.linspace(t_start*1000, t_end*1000, num=num_timesteps, endpoint=False)
				stimulus_times_exc = timestep_values[spike_mask]	
				
				#np.savetxt('arbor_2N1S_stimulation_' + self.learn_prot + '.txt', stimulus_times_exc, fmt="%i") # to check the specific stimulation
				
			elif self.learn_prot == "SLFS":
				t_start = np.linspace(3600, 4633.85, num=900, endpoint=True) # start times of the pulses (in s)
				t_end = np.linspace(3600.15, 4634, num=900, endpoint=True) # end times of the pulses (in s)
				freq = 20 # average spike frequency in Hz
				seed = int(datetime.now().timestamp() * 1e6)
						      
				# average number of spikes (random number drawn for every timestep, then filtered with probability):
				stimulus_times_exc = np.array([])
				rng = np.random.default_rng(seed)
				num_timesteps = np.int_(np.round_((t_end[0]-t_start[0])*1000./self.dt))
				for i in range(len(t_start)):
					spike_mask = rng.random(size=num_timesteps) < freq*self.dt/1000.
					timestep_values = np.linspace(t_start[i]*1000, t_end[i]*1000, num=num_timesteps, endpoint=False)
					spikes = timestep_values[spike_mask]
					stimulus_times_exc = np.concatenate([stimulus_times_exc, spikes])
				
				#np.savetxt('arbor_2N1S_stimulation_' + self.learn_prot + '.txt', stimulus_times_exc, fmt="%i") # to check the specific stimulation
				
			elif self.learn_prot == "WLFS":
				t_start = 3600 # start time in s
				t_end = 4500 # end time in s
				freq = 1 # average spike frequency in Hz
				seed = int(datetime.now().timestamp() * 1e6)
					                     
				# average number of spikes (random number drawn for every timestep, then filtered with probability):
				rng = np.random.default_rng(seed)
				num_timesteps = np.int_(np.round_((t_end-t_start)*1000./self.dt))
				spike_mask = rng.random(size=num_timesteps) < freq*self.dt/1000.
				timestep_values = np.linspace(t_start*1000, t_end*1000, num=num_timesteps, endpoint=False)
				stimulus_times_exc = timestep_values[spike_mask]			
				
				#np.savetxt('arbor_2N1S_stimulation_' + self.learn_prot + '.txt', stimulus_times_exc, fmt="%i") # to check the specific stimulation

			else: # use predefined
				stimulus_times_exc = self.syn_exc_input1_config["stimulus_times"]
				
			# create event generator
			spike_exc = arbor.event_generator(
				"syn_exc_input1",
				0.,
				arbor.explicit_schedule(stimulus_times_exc))

			return [spike_exc]
			
		else:
		
			return []
		
	# global_properties
	# Sets properties that will be applied to all cells of that kind in the model
	# gid: global identifier of the cell
	# return: the cell properties 
	def global_properties(self, kind): 

		assert kind == arbor.cell_kind.cable 

		return self.the_props
	
	# num_cells
	# return: the number of cells
	def num_cells(self):
		
		return self.ncells

	# probes
	# gid: global identifier of the cell
	# return: the probes on the given cell
	def probes(self, gid):

		return self.the_probes

#####################################
# arbor2N1S
# Runs simulation of two neurons with one synapse, stores results, and plots them
# trial: number of the current trial
# learn_prot: protocol for the learning stimulus
# runtime: runtime of the simulation in milliseconds of biological time
# config_file: configuration of simulation (JSON file)
# data_saving [optional]: specifies if fewer data (e.g., no membrane potential and current, fewer samples) shall be recorded
def arbor2N1S(trial, learn_prot, runtime, config_file, data_saving):
	
	#####################################
	# start taking the time
	t0 = time.time()

	#####################################
	# set up and run simulation
	config = json.load(open(config_file + ".json", "r"))	
	delta_t = config["simulation"]["dt"]	
	recipe = SingleRecipe(config, learn_prot, delta_t)

	clockseed = int(t0*10000)
	print("Random seed " + str(clockseed))

	context = arbor.context() # constructs a local context with one thread, no GPU, no MPI (cf. https://docs.arbor-sim.org/en/v0.5.2/python/hardware.html#arbor.context)
	domains = arbor.partition_load_balance(recipe, context) # constructs a domain_decomposition that distributes the cells in the model described by an arbor.recipe over the distributed and local hardware resources described by an arbor.context (cf. https://docs.arbor-sim.org/en/v0.5.2/python/domdec.html#arbor.partition_load_balance)
	sim = arbor.simulation(recipe, context, domains, seed = clockseed)

	sim.record(arbor.spike_recording.all)
	
	if not data_saving:
		sample_size = delta_t # sample size
	else:
		sample_size = 100*delta_t # sample size

	reg_sched = arbor.regular_schedule(0, sample_size, runtime)
	
	if not data_saving:
		handle_mem = sim.sample((0, 0), reg_sched)
		handle_curr_dens = sim.sample((0, 1), reg_sched) # current density
		handle_curr = sim.sample((0, 2), reg_sched) # current
	handle_Ca = sim.sample((0, 3), reg_sched) # calcium amount
	handle_h = sim.sample((0, 4), reg_sched) # early-phase weight
	handle_z = sim.sample((0, 5), reg_sched) # late-phase weight
	handle_p = sim.sample((0, 6), reg_sched) # protein amount
	
	sim.run(tfinal=runtime, dt=delta_t)
	
	print("Simulation completed.")

	#####################################
	# get traces and spikes from simulator
	if not data_saving:
		data_mem, _ = sim.samples(handle_mem)[0]
		data_curr_dens, _ = sim.samples(handle_curr_dens)[0]
		data_curr, _ = sim.samples(handle_curr)[0]

	data_Ca, _ = sim.samples(handle_Ca)[0]
	data_h, _ = sim.samples(handle_h)[0]
	data_z, _ = sim.samples(handle_z)[0]
	data_p, _ = sim.samples(handle_p)[0]

	#####################################
	# assemble and store data
	if not data_saving:
		data_stacked = np.column_stack(
		      [data_mem[:, 0], data_mem[:, 1], np.negative(data_curr[:, 1]), data_h[:, 1], data_z[:, 1], data_Ca[:, 1], data_p[:, 1]])
	else:
		data_stacked = np.column_stack(
	              [data_h[:, 0], data_h[:, 1], data_z[:, 1], data_Ca[:, 1], data_p[:, 1]])
	
	spike_times = np.array([])
	for s in sim.spikes():
		spike_times = np.append(spike_times, np.array([s[0][0], s[1]]), axis=0)
	spike_times = np.reshape(spike_times, (-1,2))

	del sim # free the 'sim' object
	gc.collect() # collect garbage

	out_path = "data_" + learn_prot
	if not os.path.isdir(out_path): # if the directory does not exist yet
		os.mkdir(out_path)
		
	np.savetxt(os.path.join(out_path, 'arbor_2N1S_' + str(trial) + '_traces.txt'), data_stacked, fmt="%.4f")
	np.savetxt(os.path.join(out_path, 'arbor_2N1S_' + str(trial) + '_spikes.txt'), spike_times, fmt="%.0f %.4f") # integer formatting for neuron number

	#####################################
	# do the plotting
	if not data_saving:
		plotResults(config, data_stacked, trial, mem_dyn_data = True, store_path=out_path, figure_fmt = 'svg')
	else:
		plotResults(config, data_stacked, trial, mem_dyn_data = False, store_path=out_path, figure_fmt = 'png')
	
	print("Data output completed.")

	#####################################
    	# format and print elapsed time
	time_el = round(time.time()-t0) # total elapsed time in seconds
	
	if time_el < 60:
		time_el_str = str(time_el) + " s"
	elif time_el < 3600:
		time_el_str = str(time_el // 60) + " m " + str(time_el % 60) + " s"
	else:
		time_el_str = str(time_el // 3600) + " h " + str((time_el % 3600) // 60) + " m " + str((time_el % 3600) % 60) + " s"
		
	print("Total elapsed time: " + time_el_str + ".")

#####################################
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-trial', type=int, default=0, help="number of the current trial")
	parser.add_argument('-learn', default="predefined", help="protocol for the learning stimulus")
	parser.add_argument('-runtime', type=int, default=100, help="runtime of the simulation in ms")
	parser.add_argument('-config', default="config_arbor_2N1S", help="configuration of simulation (JSON file)")
	parser.add_argument('-data_saving', type=int, choices=[0, 1], default=0, help="switch for recording fewer data")
	args = parser.parse_args()
	
	print("Running trial #" + str(args.trial) + " with protocol " + str(args.learn))
	
	arbor2N1S(args.trial, args.learn, args.runtime, args.config, args.data_saving)
	
