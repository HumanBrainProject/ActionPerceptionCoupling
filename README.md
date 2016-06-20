# ActionPerceptionCoupling
﻿
1. ActionPerceptionCycle();

Main routine implementing two coupled neural fields, one corresponding to the visual action representation and to a motor representation of actions. Both fields are implemented by neural ensembles of spiking neurons, modeled by adaptive exponential integrate and fire units (Jolivet et al. 2008). The total number of ensembles for each field is set to thirty (NN_Mtr= 30, NN_Vsn = 30). Each ensemble consists of 100 neurons with 80% excitatory and 20% inhibitory neurons. Connection strength between the neurons is chosen randomly using uniform distributions centered around appropriate mean values. The interactions between the ensembles is defined by the lateral interaction kernels (WW_Mtr, WW_Vsn), which are asymmetric in order to support traveling solutions. 

#Note:
  1. Implementation uses the Forward Euler Method. 
  2. Output neural neural activity of each ensemble is defined based on the excitatory population activity, calculating the fraction of   neurons that are activated within a short time interval [t, t+∆t] divided by ∆t.
  3. Parameters with ending “_Vsn” are related to the Visual neural field and the ones  with “_Mtr” to the Motor field.
  4. This script captures the neural activities in an AVI format video. The default path for storing the video can be changed.


2.	AsymmetricCouplingMotor

Syntax 

WW = AsymmetricCouplingMotor(NN, gama, alpha, C)

Description
This function defines the lateral interaction kernel of the motor neural field, as depicted in figure 1.  
Note
The parameter “fac” in the routine specifies the amplitude of the lateral interaction kernel. It has to be chosen large enough to stabilize a propagating solution peak.

Input Arguments


NN: number of neural ensembles placed in the motor field
gama: width of the kernel
alpha: amplitude of the Gaussian function
C: determines the DC offset and adjusts the amplitude of the positive maximum. 

3. AsymmetricCouplingVision

Syntax 

WW = AsymmetricCouplingVision(NN, gama, alpha, C)

Description
This function defines the lateral interaction kernel for the vision field, as illustrated in figure 1. 
Input Arguments


NN: number of neural ensembles placed in the motor field
gama:  width of the kernel
alpha: amplitude of the Gaussian function
C: determines the DC offset and adjusts the amplitude of the positive maximum. 

4.	MtrInput
Syntax 

GaussianInput = MtrInput(M_size,sigma1,ST_v)

Description

This function produces a transient Gaussian bump (GO signal).

Input Arguments
M_size: number of neural ensembles
sigma1: standard deviation of Gaussian bump
ST_v: Sets the Go signal onset time

5.	MtrToVsnKernel

Syntax 
	
WW = MtrToVsnKernel(NN, gama, alpha, C)

Description
This function defines a unidirectional coupler between motor and vision fields.
 Note
Since the coupling strength from motor to field direction is different from the reverse direction, two distinct functions of “MtrToVsnKernel” and “VsnToMtrKernel” are defined.
Input Arguments


NN: number of neural ensembles placed in the motor field
gama: width of the kernel
alpha: amplitude of the Gaussian function
C: determines the DC offset and adjusts the amplitude of the positive maximum. 

6.	VsnToMtrKernel

Syntax 

WW = VsnToMtrKernel(NN, gama, alpha, C)

Description
This function defines a unidirectional coupler between vision and motor fields.
Input Arguments


NN: number of neural ensembles placed in the motor field
gama: width of the kernel
alpha: amplitude of the Gaussian function
C: determines the DC offset and adjusts the amplitude of the positive maximum. 

7.	GaussianFilter
Syntax 
GaussianFilter = GaussianFilter(sigma,G_size)

Description

This function implements a Gaussian filter to smooth the excitatory population activities.

Input Arguments
G_size: Number of neurons in each ensemble
sigma: standard deviation of Gaussian filter

8.	DetectionRate
Syntax 

[Detections, DetectionRate] 
= DetectionRate(VsnOutput, WinSize, Threshold,SNR)

Description
This function calculates a detection rate based on the percentage of time that the activity of the visual field is above a defined threshold.
Note
The routine adds noise to the external visual stimulus to create a randomly fluctuating activity. The noise strength is set by the parameter “NoiseSD” and the detection threshold is set by the parameter “Threshold” in the main script.  In the current version, we define the perceived signal as the vision field activity that is above the “Threshold”. Since the activity fluctuates, a window detection parameter,”WinSize”, is set to avoid taking immediate activity above threshold into consideration. 
Input Arguments
VsnOutput: timeseris of activity in the visual field 
Winsize: size of the detection windows 
Threshold: detection threshold
SNR: Signal to Noise ratio (Parameter only used for DetectioRate plot.)

	NEST implementation:
This implementation is not yet completed.  Due to some restrictions of the NEST framework the whole design implemented already in MATLAB is not yet transferable to NEST. Currently there are two scripts provided:

1.	Amari_field

Syntax 
Amari_field.py

Description

This script is a test version for implementing a neural field in NEST. A simple neural field consisting of thirty one identical neurons, which are connected through the interaction kernel (funW), is defined. The input signal of the neural field is defined by a traveling input peak (funS). This input signal determines the currents injected in the corresponding neurons. 
The neural activity u(x,t) is interpreted as membrane potential, and the threshold function of the neural field follows from the multistable characteristics of the individual neurons.
 
2.	ActionPerception

Syntax 
ActionPerception.py

Description
This script defines the skeleton of the full architecture.  It specifies the required neural ensembles and their interactions, and computes the relevant population activities as well.
The neural field topology is created using “topology” module which defines a ring topology with composite nodes. Each node in this topology-defined 1D layer is a neural ensemble composed of measurement devices, 80 excitatory (E) and 20 inhibitory (I) aEIF neurons that are fully connected. All E and I neurons send spikes to their correspondents  “parrot” neurons (P) that represent the activity of related neural ensembles. Since a parrot neuron simply emits one spike for every incoming spike, we can interpret the neural activity of each ensemble by the number of incoming spikes to its correspondent “parrot” neuron. The current step is to accomplish the goal of neural field design by connecting “parrot” neurons to other neurons based on the lateral interaction kernel (funW) which is already defined in the Amari_field script

