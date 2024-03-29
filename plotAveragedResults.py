#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import json
import argparse
	
#####################################
# plotAveragedResults
# Plots the averaged results of simulations of synaptic weights changing based on calcium dynamics
# config: configuration parameters in JSON format
# data_stacked: two-dimensional array containing the values of the early-phase weight, late-phase weight, and calcium amount over time
# protocol: name of the protocol that was used to obtain the data
# out_format: format of the graphics output
# sampling [optional]: sampling (1: best resolution, >1: samples are left out)
def plotAveragedResults(config, data_stacked, protocol, out_format, sampling=1):

	h_0 = config["synapses"]["syn_exc_calcium_plasticity"]["h_0"]
	fig, axes = plt.subplots(nrows=3, ncols=1, sharex=False, figsize=(10, 10))

	# prepare data
	data_stacked[:,0] /= (1000*3600) # change time to hours
	data_stacked[:,0] -= 1 # shift time by 1 hour
	data_stacked = data_stacked[data_stacked[:,0] > -0.5] # strip some of the idle initial time span
	data_stacked = data_stacked [np.arange(len(data_stacked[:,0])) % sampling == 0] # sample only every x-th value
	time = data_stacked[:,0]
	h_mean = data_stacked[:,1]/h_0*100 # normalize
	h_stdev = data_stacked[:,2]/h_0*100
	z_mean = (data_stacked[:,3]/h_0+1)*100
	z_stdev = (data_stacked[:,4]/h_0)*100
	Ca_mean = data_stacked[:,5] 
	Ca_stdev = data_stacked[:,6]
	p_mean = data_stacked[:,7]
	p_stdev = data_stacked[:,8]
	w_mean = (data_stacked[:,9]/h_0)*100
	w_stdev = (data_stacked[:,10]/h_0)*100
	
	# set axis labels for plasticity dynamics
	axes[0].set_xlabel("Time (h)")
	axes[0].set_ylabel("Synaptic weight (%)")
	axes[0].set_xlim(-0.02, 7.01)
	
	# plot data for plasticity dynamics
	axes[0].plot(time, h_mean, color="#800000", label='h', marker='None', zorder=10)
	axes[0].fill_between(time, h_mean+h_stdev, h_mean-h_stdev, color="#800000", alpha=0.5)
	axes[0].plot(time, z_mean, color="#1f77b4", label='z', marker='None', zorder=9)
	axes[0].fill_between(time, z_mean+z_stdev, z_mean-z_stdev, color="#1f77b4", alpha=0.5)
	axes[0].plot(time, w_mean, color="#ff7f0e", label='w', marker='None', zorder=8)
	axes[0].fill_between(time, w_mean+w_stdev, w_mean-w_stdev, color="#ff7f0e", alpha=0.5)
	
	# adjust frame of axes[0] and add lines for tag and protein synthesis thresholds
	if np.sum(h_mean-100) > 0.: # if LTP dominates
		axes[0].axhline(y=(1+config["synapses"]["syn_exc_calcium_plasticity"]["theta_pro"]/h_0)*100, label='Protein thresh.', linestyle='-.', color="#dddddd", zorder=5)
		axes[0].axhline(y=(1+config["synapses"]["syn_exc_calcium_plasticity"]["theta_tag"]/h_0)*100, label='Tag thresh.', linestyle='dashed', color="#dddddd", zorder=4)
		axes[0].set_ylim(95,205)
		axes[0].legend(loc=(0.79,0.24)) # create legend for axes[0]
	else: # if LTD dominates
		axes[0].axhline(y=(1-config["synapses"]["syn_exc_calcium_plasticity"]["theta_pro"]/h_0)*100, label='Protein thresh.', linestyle='-.', color="#dddddd", zorder=5)
		axes[0].axhline(y=(1-config["synapses"]["syn_exc_calcium_plasticity"]["theta_tag"]/h_0)*100, label='Tag thresh.', linestyle='dashed', color="#dddddd", zorder=4)
		axes[0].set_ylim(-5,105)
		axes[0].legend(loc=(0.79,0.05)) # create legend for axes[0]

	# set axis labels for calcium dynamics
	axes[1].set_xlabel("Time (h)")
	axes[1].set_ylabel("Calcium concentration")
	axes[1].set_xlim(-0.01, 0.31)
	
	# plot data for calcium and protein dynamics
	axes[1].plot(time, Ca_mean, color="#c8c896", label='Ca', marker='None', zorder=10)
	axes[1].fill_between(time, Ca_mean+Ca_stdev, Ca_mean-Ca_stdev, color="#c8c896", alpha=0.5)
	axes[1].axhline(y=config["synapses"]["syn_exc_calcium_plasticity"]["theta_p"], label='LTP thresh.', linestyle='dashed', color="#969664", zorder=9)
	axes[1].axhline(y=config["synapses"]["syn_exc_calcium_plasticity"]["theta_d"], label='LTD thresh.', linestyle='dashed', color="#969696", zorder=8)
	axes[1].legend(loc=(0.78,0.65)) # create legend for axes[1]

	# set axis labels for protein dynamics
	axes[2].set_xlabel("Time (h)")
	axes[2].set_ylabel("Protein concentration")
	axes[2].set_xlim(-0.02, 7.01)

	# plot data for protein dynamics
	axes[2].plot(time, p_mean, color="#008000", label='p', marker='None', zorder=10)
	axes[2].fill_between(time, p_mean+p_stdev, p_mean-p_stdev, color="#008000", alpha=0.5)
	axes[2].legend(loc=(0.78,0.65)) # create legend for axes[1]
	
	# save figure as vector graphics
	fig.savefig('arbor_averaged_' + protocol + '.' + out_format)
	
	
#####################################
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-protocol', type=str, required=True, help='name of the protocol being used as a label for the data')
	parser.add_argument('-out_format', type=str, default='png', help='format of the graphics output')
	parser.add_argument('-sampling', type=int, default=1, help='sampling (1: best resolution, >1: samples are left out)')
	parser.add_argument('-config', type=str, default="config_arbor_2N1S", help="configuration of simulation (JSON file)")
	args = parser.parse_args()
	
	data_stacked = np.loadtxt('averaged_' + str(args.protocol) + '.txt')
	config = json.load(open(args.config + ".json", "r"))	
	
	plotAveragedResults(config, data_stacked, str(args.protocol), str(args.out_format), int(args.sampling))
