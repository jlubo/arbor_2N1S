#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import json
import os
	
#####################################
# plotResults
# Plots the results of a simulation of synaptic weights changing based on calcium dynamics
# config: configuration parameters in JSON format
# data_stacked: two-dimensional array containing the values of the membrane potential, weights, calcium amount, etc. over time
# trial [optional]: number of the current trial
# mem_dyn_data [optional]: specifies if membrane potential and current shall be plotted
# figure_fmt [optional]: format of resulting graphics file
# store_path [optional]: path to store resulting graphics file
def plotResults(config, data_stacked, trial = 0, mem_dyn_data = False, store_path = ".", figure_fmt = 'png'):

	h_0 = config["synapses"]["syn_exc_calcium_plasticity"]["h_0"] # get reference value of synaptic weight
	
	# adapt depending on whether membrane dynamics shall be plotted or not
	if mem_dyn_data:
		num_rows = 3 # number of rows in figure array
		data_ptr = 3 # pointer to the first column of synapse data
	else:
		num_rows = 2 # number of rows in figure array
		data_ptr = 1
		
	fig, axes = plt.subplots(nrows=num_rows, ncols=1, sharex=False, figsize=(10, 10)) # create figure with 'num_rows' subfigures

	# set axis labels for plasticity dynamics
	axes[0].set_xlabel("Time (ms)")
	axes[0].set_ylabel("Synaptic weight (%)")
	
	# plot data for plasticity dynamics
	axes[0].plot(data_stacked[:,0], data_stacked[:,data_ptr]/h_0*100, color="#800000", label='h', marker='None', zorder=10)
	axes[0].plot(data_stacked[:,0], (data_stacked[:,data_ptr+1]+1)*100, color="#1f77b4", label='z', marker='None', zorder=9)
	axes[0].axhline(y=(config["synapses"]["syn_exc_calcium_plasticity"]["theta_pro"]/h_0+1)*100, label='Protein thresh.', linestyle='-.', color="#dddddd", zorder=5)
	axes[0].axhline(y=(config["synapses"]["syn_exc_calcium_plasticity"]["theta_tag"]/h_0+1)*100, label='Tag thresh.', linestyle='dashed', color="#dddddd", zorder=4)
	
	# create legend for plasticity dynamics
	axes[0].legend(loc='upper right')

	if mem_dyn_data:
		# set axis labels for membrane potential and current plots
		axes[1].set_xlabel("Time (ms)")
		axes[1].set_ylabel("Membrane potential (mV)")
		ax1twin = axes[1].twinx() # create twin axis for axes[1]
		ax1twin.set_ylabel("Current (nA)")
		
		# plot data for membrane potential and current data
		ax1g1 = axes[1].plot(data_stacked[:,0], data_stacked[:,1], color="#ff0000", label='Membrane pot.', marker='None', zorder=10)
		ax1g2 = ax1twin.plot(data_stacked[:,0], data_stacked[:,2], color="#ffee00", label='Membrane curr.', marker='None', zorder=9)
		
		# create common legend for axes[1] and ax1twin
		handles, labels = axes[1].get_legend_handles_labels()
		handles_twin, labels_twin = ax1twin.get_legend_handles_labels()
		axes[1].legend(handles + handles_twin, labels + labels_twin, loc='upper right')

	# set axis labels for calcium and protein plots
	axes[num_rows-1].set_xlabel("Time (ms)")
	axes[num_rows-1].set_ylabel("Calcium or protein amount")
	
	# plot data for calcium and protein dynamics
	axes[num_rows-1].plot(data_stacked[:,0], data_stacked[:,data_ptr+2], color="#c8c896", label='Ca', marker='None', zorder=10)
	axes[num_rows-1].plot(data_stacked[:,0], data_stacked[:,data_ptr+3], color="#008000", label='p', marker='None', zorder=9)
	axes[num_rows-1].axhline(y=config["synapses"]["syn_exc_calcium_plasticity"]["theta_p"], label='LTP thresh.', linestyle='dashed', color="#969664", zorder=8)
	axes[num_rows-1].axhline(y=config["synapses"]["syn_exc_calcium_plasticity"]["theta_d"], label='LTD thresh.', linestyle='dashed', color="#969696", zorder=7)
	
	# create legend for calcium and protein plots
	axes[num_rows-1].legend(loc='upper right')
	
	# save figure in given format (e.g., 'png' or 'svg')
	fig.savefig(os.path.join(store_path, 'arbor_2N1S_traces_' + str(trial) + '.' + figure_fmt))	
	
#####################################
if __name__ == '__main__':
	config = json.load(open("config_arbor_2N1S.json", "r"))
	
	os.chdir("traces_predef/")
	data_stacked = np.loadtxt(f'arbor_2N1S_0_traces.txt')
	
	plotResults(config, data_stacked, mem_dyn_data = True, figure_fmt = 'svg')
