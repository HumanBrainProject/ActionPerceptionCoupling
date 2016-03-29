# Author : Dr. Espen Hagen, Juelich, Germany
# Modified by : Mohammad Hovaidi Ardestani, Tuebingen, Germany
# To be evolved :
# This code will be used as the main core of the model implemented in NEST. 
# Due to some restrictions of NEST the whole design implemented already in MATLAB
# has not been yet transfered to NEST.


# -*- coding: utf-8 -*-
#!/usr/bin/env python


import sys
sys.path.append('/opt/nest/lib/python2.7/site-packages')
import numpy as np
import matplotlib.pyplot as plt
import nest
from nest import topology as tp
import pylab as pl
import math

#simulation time and resolution
dt = 0.1
T = 1000
start = 0.0
stop = 800.0

#model neuron
model = 'aeif_cond_alpha'

#size of excitatory and inhibitory populations within cell assembly
NE = 80
NI = 20

#Input current to the neuron models
#amplitude = 100.0    #pA
I_e_ex = 800.        #input current in pA to excitatory neurons
I_e_in = 0.          #inh. neurons

#default synapse weight times 10, as we're running with weight = 0.1
syn_g = 1E-1

#delay and "weight" of connections
delay = 1.
weight = 0.1

#number of Neural Ensembles
n_ensembles = 30

#noise parameters
noise_mean = 0.
noise_std = I_e_ex*0.1
noise_dt = dt


# defining the weight function
def  funW(x) :

    gama = 5
    alpha = .7
    C = .7
    fac = 20e-9 # Kernel strength

    ''' Weight Kernel '''
    y = ((1+math.cos(alpha*x))/2)**gama; 
    '''               '''
#    y = fac * (y-C)
    y = fac * y
    return y
vec_funW = pl.vectorize(funW)
s = pl.linspace(-pl.pi,pl.pi,n_ensembles)
knl = vec_funW(s)    # weight matrix  

#==============================================================================
# ''' define the stimulation  '''
# 
# def funS(x):
#     c =0.6
#     mean = 0.0
#     I = 10
#     y = I/(pl.sqrt(2*pl.pi)*c) * pl.exp(-(x-mean)**2/(2*c**2))
#     return y
#     
# vec_funS = pl.vectorize(funS)
# s = pl.linspace(-pl.pi,pl.pi,n_ensembles)
# stim = amplitude*vec_funS(s) # stimulation matrix 1 by pop
# params_dc =[{'amplitude':stim[i]} for i in range(0,n_ensembles)] # parameter of multiple nodes is a list of dictionary
# nest.SetDefaults('dc_generator',{"start":start, "stop":stop})
# dc = nest.Create('dc_generator',n_ensembles,params_dc)
#==============================================================================


nest.ResetKernel()
nest.SetKernelStatus({'resolution' : dt,
                      'total_num_virtual_procs' : 4})

#update defaults for neuron model
nest.SetDefaults(model, dict(
    E_ex = 0.,
    g_ex = syn_g,
    V_peak = 20.0,
    tau_syn_ex = 2.0,
    E_in = -70.,
    g_in = syn_g,
    tau_syn_in = 10.0,
    gsl_error_tol = 1.E-6,
    ))


#the input currents differ between excitatory and inhibitory cells, so we make
#copies of the base model
nest.CopyModel(model, 'exc', {'I_e' : I_e_ex})
nest.CopyModel(model, 'inh', {'I_e' : I_e_in})

#make sure we're recording to memory
nest.SetDefaults('spike_detector',{"start":start, "stop":stop},{'to_memory' : True})
#nest.SetDefaults('spike_detector', {
#   'to_memory' : True,
#})


#create 1D layer of n_ensembles
layer = tp.CreateLayer({
    		'rows' : n_ensembles,
    		'columns' : 1,
			'edge_wrap': True,			#ring topology(to reduce the effect of boundary conditions)
    		'elements' : ['exc', NE,    #excitatory population
			'inh', NI,    				#inhibitory population
            'parrot_neuron', 1,  		#one parrot neuron per ensemble
                  		 ]
})

# shows the network topology
# tp.PlotLayer(layer) 					

#record to memory, create one spike detector for each main population
nest.SetDefaults('spike_detector',{"start":start, "stop":stop},{'to_memory' : True})
#nest.SetDefaults('spike_detector', {
#    'to_memory' : True,
#})
spikes_ex = nest.Create('spike_detector')
spikes_in = nest.Create('spike_detector')
spikes_parrot = nest.Create('spike_detector')

#create a voltmeter to record potentials in neurons
#nest.SetDefaults('Voltmeter',{"start":start, "stop":stop})
voltmeter = nest.Create('voltmeter')

#white noise generator, set defaults and create object
nest.SetDefaults('noise_generator', {
    'mean' : noise_mean,
    'std' : noise_std,
    'dt' : noise_dt,    
},{"start":start, "stop":stop})
noise = nest.Create('noise_generator')

#poisson generator object
nest.SetDefaults('poisson_generator',{"start":start, "stop":stop})
poisson = nest.Create('poisson_generator')
nest.SetStatus(poisson, dict(rate=20.))


#expose all nodes of network, including parrot neuron
nodes = nest.GetLeaves(layer)[0]

#Create a list of all the parrot neurons
i_parrot_list = []

for i in range(n_ensembles):
    #first and last indices for excitatory and inhibitory units in subnet,
    #and for parrot neuron in this subnet
    i_ex_0 = i*NE
    i_ex_1 = (i+1)*NE
    i_in_0 = i*NI + NE*n_ensembles
    i_in_1 = (i+1)*NI + NE*n_ensembles
    i_parrot_list.append(n_ensembles*(NE+NI) + i)
#==============================================================================
#     ''' In case of different arrengment of Exc and Inh neurons''' 
#	  i_parrot = n_ensembles*(NE+NI) + i 
#     i_ex_0 = i*(NE+NI)
#     i_ex_1 = (i*(NE + NI)) + NE
#     i_in_0 = (i*(NE + NI)) + NE
#     i_in_1 = (i+1)*(NE+NI)  
#     
#     
#     '''connect stimulation to neuorns'''
#     nest.Connect([dc[i]],nodes[i_ex_0:i_ex_1])
#==============================================================================
    
    print i_ex_0, i_ex_1, i_in_0, i_in_1, i_parrot_list[i]
    
    #connect exc to exc
    nest.Connect(nodes[i_ex_0:i_ex_1], nodes[i_ex_0:i_ex_1],
                 conn_spec='all_to_all',
                 syn_spec={'model' : 'static_synapse',
                           'weight' : weight,
                           'delay' : delay})
    
    #connect exc to inh
    nest.Connect(nodes[i_ex_0:i_ex_1], nodes[i_in_0:i_in_1],
                 conn_spec='all_to_all',
                 syn_spec={'model' : 'static_synapse',
                           'weight' : weight,
                           'delay' : delay})                 

    #connect inh to inh
    nest.Connect(nodes[i_in_0:i_in_1], nodes[i_in_0:i_in_1], 
                 conn_spec='all_to_all',
                 syn_spec={'model' : 'static_synapse',
                           'weight' : weight,
                           'delay' : delay})

    #connect inh to exc
    nest.Connect(nodes[i_in_0:i_in_1], nodes[i_ex_0:i_ex_1], 
                 conn_spec='all_to_all',
                 syn_spec={'model' : 'static_synapse',
                           'weight' : -weight,
                           'delay' : delay})
    
    #connect noise to all neurons
#    nest.Connect(noise, nodes[i_ex_0:i_in_1])
    nest.Connect(noise, nodes[i_ex_0:i_ex_1])
    nest.Connect(noise, nodes[i_in_0:i_in_1])
    
    #connect all neurons to parrot neuron
#    nest.Connect(nodes[i_ex_0:i_in_1], (nodes[i_parrot],))
#    nest.Connect(nodes[i_ex_0:i_ex_1], (nodes[i_parrot_list[i]],))
#	 nest.Connect(nodes[i_in_0:i_in_1], (nodes[i_parrot_list[i]],))

    #connect Exc neurons to parrot neuron
    nest.Connect(nodes[i_ex_0:i_ex_1], (nodes[i_parrot_list[i]],))


    #connect spike detector excitatory and inhibitory neurons
    nest.Connect(nodes[i_ex_0:i_ex_1], spikes_ex)
    nest.Connect(nodes[i_in_0:i_in_1], spikes_in)

    #connect parrot neuron to corresponding spike detector
    #nest.Connect((nodes[i_parrot], ), spikes_parrot)
    nest.Connect((nodes[i_parrot_list[i]], ), spikes_parrot)

    #connect one excitatory and inhibitory neuron to voltmeter
    nest.Connect(voltmeter, (nodes[i_ex_0], nodes[i_in_0]))
 
#==============================================================================
# #Connecting parrot neurons using Gaussuian Kernel
# for i in range (n_ensembles):
#     for j in range(n_ensembles):
#         para = {'weight': knl[i-j]}
#         conn_para = {'rule':'one_to_one'}
#         #nest.Connect([nodes[i]],[nodes[j]],conn_para, para) 
#         nest.Connect([nodes[i_parrot_list[i]]],[nodes[i_parrot_list[j]]],conn_para, para) 
#==============================================================================

#print out some info
nest.PrintNetwork(depth=3)

#run simulation
nest.Simulate(T)

#extract events and do a simple plot of activity
#events_ex = nest.GetStatus(spikes_ex)[0]['events']
events_ex = nest.GetStatus(spikes_ex,'events')[0]
events_in = nest.GetStatus(spikes_in)[0]['events']
events_parrot = nest.GetStatus(spikes_parrot)[0]['events']
events_voltmeter = nest.GetStatus(voltmeter)[0]['events']




#==============================================================================
# 
# for i in range(n_ensembles):
#     #first and last indices for excitatory and inhibitory units in subnet,
#     #and for parrot neuron in this subnet
#     i_ex_0 = i*NE
#     i_ex_1 = (i+1)*NE
#     i_in_0 = i*NI + NE*n_ensembles
#     i_in_1 = (i+1)*NI + NE*n_ensembles
# #==============================================================================
# #    i_parrot = n_ensembles*(NE+NI) + i
# #    i_ex_0 = i*(NE+NI)
# #    i_ex_1 = (i*(NE + NI)) + NE
# #    i_in_0 = (i*(NE + NI)) + NE
# #    i_in_1 = (i+1)*(NE+NI)  
# #==============================================================================
# 
#     
#     fig, ax = plt.subplots(4)
#     fig.suptitle('spike raster, rate, voltages for Neural Ensemble %i' % (i+1))
#     
#     #raster plots
#     inds_ex = (events_ex['senders'] >= nodes[i_ex_0]) & (events_ex['senders'] < nodes[i_ex_1])
#     ax[0].plot(events_ex['times'][inds_ex], events_ex['senders'][inds_ex], 'b.')
#     ax[0].set_xlim(0, T)
#     ax[0].set_ylim(nodes[i_ex_0], nodes[i_ex_1])
#     ax[0].set_ylabel('exc. units')
# 
#     #raster plots
#     inds_in = (events_in['senders'] >= nodes[i_in_0]) & (events_in['senders'] < nodes[i_in_1])
#     ax[1].plot(events_in['times'][inds_in], events_in['senders'][inds_in], 'r.')
#     ax[1].set_xlim(0, T)
#     ax[1].set_ylim(nodes[i_in_0], nodes[i_in_1])
#     ax[1].set_ylabel('inh. units')
#     
#     #rate plot
#     #inds_parrot = events_parrot['senders'] == nodes[i_parrot]
#     inds_parrot = events_parrot['senders'] == nodes[i_parrot_list[i]]
#     ax[2].hist(events_parrot['times'][inds_parrot], bins=np.arange(T+dt))
#     ax[2].set_xlim(0, T)
#     ax[2].set_ylabel('spike #')
# 
#     #membrane voltage traces
# 
#     inds_ex = events_voltmeter['senders'] == nodes[i_ex_0]
#     inds_in = events_voltmeter['senders'] == nodes[i_in_0]
#     ax[3].plot(events_voltmeter['times'][inds_ex], events_voltmeter['V_m'][inds_ex], 'b')
#     ax[3].plot(events_voltmeter['times'][inds_in], events_voltmeter['V_m'][inds_in], 'r')
#     ax[3].set_ylabel('Vm')
#==============================================================================
i = 2

i_ex_0 = i*NE
i_ex_1 = (i+1)*NE
i_in_0 = i*NI + NE*n_ensembles
i_in_1 = (i+1)*NI + NE*n_ensembles

fig, ax = plt.subplots(4)
fig.suptitle('spike raster, rate, voltages for Neural Ensemble %i' % (i+1))
 
#raster plots
inds_ex = (events_ex['senders'] >= nodes[i_ex_0]) & (events_ex['senders'] < nodes[i_ex_1])
ax[0].plot(events_ex['times'][inds_ex], events_ex['senders'][inds_ex], 'b.')
ax[0].set_xlim(0, T)
ax[0].set_ylim(nodes[i_ex_0], nodes[i_ex_1])
ax[0].set_ylabel('exc. units')

#raster plots
inds_in = (events_in['senders'] >= nodes[i_in_0]) & (events_in['senders'] < nodes[i_in_1])
ax[1].plot(events_in['times'][inds_in], events_in['senders'][inds_in], 'r.')
ax[1].set_xlim(0, T)
ax[1].set_ylim(nodes[i_in_0], nodes[i_in_1])
ax[1].set_ylabel('inh. units')

#rate plot
#inds_parrot = events_parrot['senders'] == nodes[i_parrot]
inds_parrot = events_parrot['senders'] == nodes[i_parrot_list[i]]
ax[2].hist(events_parrot['times'][inds_parrot], bins=np.arange(T+dt))
ax[2].set_xlim(0, T)
ax[2].set_ylabel('spike #')

#membrane voltage traces
 
inds_ex = events_voltmeter['senders'] == nodes[i_ex_0]
inds_in = events_voltmeter['senders'] == nodes[i_in_0]
ax[3].plot(events_voltmeter['times'][inds_ex], events_voltmeter['V_m'][inds_ex], 'b')
ax[3].plot(events_voltmeter['times'][inds_in], events_voltmeter['V_m'][inds_in], 'r')
ax[3].set_ylabel('Vm')    

plt.show()
