# ActionPerceptionCoupling
﻿
Main routine implementing two coupled neural fields, one corresponding to the visual action representation and to a motor representation of actions. Both fields are implemented by neural ensembles of spiking neurons, modeled by adaptive exponential integrate and fire units (Jolivet et al. 2008). The total number of ensembles for each field is set to thirty (NN_Mtr= 30, NN_Vsn = 30). Each ensemble consists of 100 neurons with 80% excitatory and 20% inhibitory neurons. Connection strength between the neurons is chosen randomly using uniform distributions centered around appropriate mean values. The interactions between the ensembles is defined by the lateral interaction kernels (WW_Mtr, WW_Vsn), which are asymmetric in order to support traveling solutions. 

#Note:
	Implementation uses the Forward Euler Method. 
	Output neural neural activity of each ensemble is defined based on the excitatory population activity, calculating the fraction of       neurons that are activated within a short time interval [t, t+∆t] divided by ∆t.
	Parameters with ending “_Vsn” are related to the Visual neural field and the ones  with “_Mtr” to the Motor field.
	This script captures the neural activities in an AVI format video. The default path for storing the video can be changed.
