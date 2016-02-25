%{
Author : Mohammad Hovaidi Ardestani
Date and Place: 22.08.2015, CIN, Tuebingen, Germany
the last modification date : 02.01.2016
.
.
.
This is the core part of the model in which the model architecture has been 
designed. 
%}


tic
clear all
close all

%% Parameters setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Simulation Time Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1e-3;                                 % time step
Time = 1;                                  % simulation time
T = dt:dt:Time;                            % time step vector   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%% Synaptic current related parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% They are extracted from (Rolls and Deco 2014)

tau_syn_e_Mtr = 2e-3;                      % tua_GABA = 2ms
tau_syn_e_Vsn = 2e-3;                      % tua_GABA = 2ms
tau_syn_i_Mtr = 10e-3;                     % tua_GABA = 10ms
tau_syn_i_Vsn = 10e-3;                     % tua_GABA = 10ms

E_e_Mtr = 0;                              
E_e_Vsn = 0;                             
E_i_Mtr = -70e-3;   
E_i_Vsn = -70e-3;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Neuron Model Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Default 80 
N_e_Mtr = 80;                          N_i_Mtr = N_e_Mtr/4;
N_e_Vsn = 80;                          N_i_Vsn = N_e_Mtr/4;

%%% Default 0.45e-9;  
C_m_Mtr = 281e-12;                     % membrane capacitance(F)
C_m_Vsn = 281e-12;                     % membrane capacitance(F)

%%% Default 25e-9;
g_L_Mtr = 30e-9;                       % leak conductance(S)
g_L_Vsn = 30e-9;                       % leak conductance(S)
 
E_L_Mtr = -70.6e-3;                    % leak reversal potential(V)
E_L_Vsn = -70.6e-3;                    % leak reversal potential(V)

V_T_Mtr = -50.4e-3;                    % Spike threshold(V)
V_T_Vsn = -50.4e-3;                    % Spike threshold(V)

Delta_T_Mtr = 2e-3;                    % slope factor(V)
Delta_T_Vsn = 2e-3;                    % slope factor(V)

tau_w_Mtr = 144e-3;                    % adaptation time constant(S)
tau_w_Vsn = 144e-3;                    % adaptation time constant(S)

a_s_Mtr = 4e-9;                        % The level of subthreshold adaptation(S)
a_s_Vsn = 4e-9;                        % The level of subthreshold adaptation(S)

a_Mtr = [a_s_Mtr*ones(N_e_Mtr,1);      a_s_Mtr*ones(N_i_Mtr,1)];
a_Vsn = [a_s_Vsn*ones(N_e_Vsn,1);      a_s_Vsn*ones(N_i_Vsn,1)];

b_s_Mtr = 0.0805e-9;                   % spike-triggered adaptation(A)
b_s_Vsn = 0.0805e-9;                   % spike-triggered adaptation(A)

b_Mtr = [b_s_Mtr*ones(N_e_Mtr,1);      b_s_Mtr*ones(N_i_Mtr,1)];
b_Vsn = [b_s_Vsn*ones(N_e_Vsn,1);      b_s_Vsn*ones(N_i_Vsn,1)];

Vr_Mtr = [E_L_Mtr*ones(N_e_Mtr,1);     E_L_Mtr*ones(N_i_Mtr,1)]; 
Vr_Vsn = [E_L_Vsn*ones(N_e_Vsn,1);     E_L_Vsn*ones(N_i_Vsn,1)]; 

V_peak_Mtr = 20e-3;                    % peak potential(V)
V_peak_Vsn = 20e-3;                    % peak potential(V)

N_Mtr = N_e_Mtr + N_i_Mtr;             % Number of all neurons
N_Vsn = N_e_Vsn + N_i_Vsn;             % Number of all neurons 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% Coupling Matrix parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN_Mtr = 30;                           % Number of ensembles in Motor
gama_Mtr = 10;     alpha_Mtr = .2;      C_Mtr =.7;


NN_Vsn = 30;                           % Number of ensembles in Vision
gama_Vsn = 10;     alpha_Vsn = .2;      C_Vsn =.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% SNR = .4;                            % Signal to Noise Ratio
% NoiseSD = .8;
%% Coupling neurons inter and intra fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% W : weight matrix of synaptic connections strength within each ensemble
%%% WW : lateral connection kernel used in each neural field
%%% WWW : lateral connection kernel used to interconnect both neural fields 

%%% Wfac :specifies the strength of the connections among neurons
Wfac  =  1e-11;                         % Weight factor 

II_Mtr = Wfac*rand(N_i_Mtr,N_i_Mtr);
EE_Mtr = Wfac*rand(N_e_Mtr,N_e_Mtr);
EI_Mtr = Wfac*rand(N_i_Mtr,N_e_Mtr);
IE_Mtr = -Wfac*rand(N_e_Mtr,N_i_Mtr);

Inhibitory_Mtr = [IE_Mtr;II_Mtr];
Excitatory_Mtr = [EE_Mtr;EI_Mtr];
W_Mtr = horzcat(Excitatory_Mtr,Inhibitory_Mtr);   
W_Mtr(1:N_Mtr+1:N_Mtr*N_Mtr) = 0;       % No neuron can have effect on itself

II_Vsn = Wfac*rand(N_i_Vsn,N_i_Vsn);
EE_Vsn = Wfac*rand(N_e_Vsn,N_e_Vsn);
EI_Vsn = Wfac*rand(N_i_Vsn,N_e_Vsn);
IE_Vsn = -Wfac*rand(N_e_Vsn,N_i_Vsn);

Inhibitory_Vsn = [IE_Vsn;II_Vsn];
Excitatory_Vsn = [EE_Vsn;EI_Vsn];
W_Vsn = horzcat(Excitatory_Vsn,Inhibitory_Vsn); 
W_Vsn(1:N_Vsn+1:N_Vsn*N_Vsn) = 0;      % No neuron can have effect on itself

WW_Mtr = AsymmetricCouplingMotor(NN_Mtr, gama_Mtr, alpha_Mtr, C_Mtr);
WW_Vsn = AsymmetricCouplingVision(NN_Vsn, gama_Vsn, alpha_Vsn, C_Vsn);

WWW_Mtr = MtrToVsnKernel(NN_Mtr, gama_Mtr, alpha_Mtr, C_Mtr);
WWW_Vsn = VsnToMtrKernel(NN_Vsn, gama_Vsn, alpha_Vsn, C_Vsn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% memory allocation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Population Activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_e_Mtr = zeros(NN_Mtr,length(T));     % excitatory population activity of each ensemble
A_i_Mtr = zeros(NN_Mtr,length(T));     % inhibitory population activity of each ensemble
A_Mtr = zeros(NN_Mtr,length(T));       % total population activity of each ensemble
A_fnl_Mtr = zeros(NN_Mtr,length(T));   % guassian filtered population activity

A_e_Vsn = zeros(NN_Vsn,length(T));     % excitatory population activity of each ensemble         
A_i_Vsn = zeros(NN_Vsn,length(T));     % inhibitory population activity of each ensemble        
A_Vsn = zeros(NN_Vsn,length(T));       % total population activity of each ensemble            
A_fnl_Vsn = zeros(NN_Mtr,length(T));   % guassian filtered population activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%% membrane potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_Mtr = zeros(NN_Mtr,N_Mtr,length(T));  
v_Mtr(:,:,1) = E_L_Mtr;                % initialize membrane potential  
v_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); 
v_Vsn(:,:,1) = E_L_Vsn;                % initialize membrane potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_Mtr = zeros(NN_Mtr,N_Mtr,length(T)); % adaptaion variable (~Ca-K current)
w_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); % adaptaion variable (~Ca-K current)

%%% Currents %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_ext_Mtr = zeros(NN_Mtr,N_Mtr,length(T));
I_ext_Vsn = zeros(NN_Vsn,N_Vsn,length(T));

I_syn_Mtr = zeros(NN_Mtr,N_Mtr,length(T));
I_syn_e_Mtr = zeros(NN_Mtr,N_Mtr,length(T)); 
I_syn_i_Mtr = zeros(NN_Mtr,N_Mtr,length(T)); 

I_syn_Vsn = zeros(NN_Vsn,N_Vsn,length(T));
I_syn_e_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); 
I_syn_i_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); 

I_MtrToVsn = zeros(NN_Vsn,N_Vsn,length(T));
I_VsnToMtr = zeros(NN_Vsn,N_Vsn,length(T));

WE_Mtr = zeros(NN_Mtr,length(T));      % intrafield current among Ensembles
WF_Mtr = zeros(NN_Mtr,length(T));      % interFields current from Motor to Vision

WE_Vsn = zeros(NN_Vsn,length(T));
WF_Vsn = zeros(NN_Vsn,length(T));

I_Field_Mtr = zeros(NN_Mtr,N_Mtr,length(T));
I_Field_Vsn = zeros(NN_Vsn,N_Vsn,length(T));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Synaptic Cunductance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_syn_e_Mtr = zeros(NN_Mtr,N_Mtr,length(T)); g_syn_e_Mtr(:,:,1) = 0;
g_syn_i_Mtr = zeros(NN_Mtr,N_Mtr,length(T)); g_syn_i_Mtr(:,:,1) = 0;

g_syn_e_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); g_syn_e_Vsn(:,:,1) = 0;
g_syn_i_Vsn = zeros(NN_Vsn,N_Vsn,length(T)); g_syn_i_Vsn(:,:,1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SpikingNeurons_Mtr = cell(NN_Mtr,length(T));
FiringNeurons_e_Mtr = zeros(NN_Mtr,length(T));
FiringNeurons_i_Mtr = zeros(NN_Mtr,length(T));
  
SpikingNeurons_Vsn = cell(NN_Vsn,length(T));
FiringNeurons_e_Vsn = zeros(NN_Vsn,length(T));
FiringNeurons_i_Vsn = zeros(NN_Vsn,length(T));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%% Stimuli settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VM_delay = 0;               % Delay that defines the Ascunchronisity between vision and  motor simuli
VM_delay = 280;
% VM_delay = 560;
% delay = 0;                  % Delay that works for the other experiments.
ST_v = 1;                   % Visiual Stimulus starting time
Sf = 30e-9;                 % Stimulus_factor of strength
%%%%%%%%%%%%%% Waiting bar settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0,'Please be patient!');

for t = 2:length(T)    
    for EnId = 1:NN_Mtr     % Ensemble ID   
      
%%% If different connection weights needed in each ensemble %%%%%%%%%%%%%%% 
%         II_Mtr = 1e-11*rand(N_i_Mtr,N_i_Mtr);
%         EE_Mtr = 1e-11*rand(N_e_Mtr,N_e_Mtr);
%         EI_Mtr = 1e-11*rand(N_i_Mtr,N_e_Mtr);
%         IE_Mtr = -1e-11*rand(N_e_Mtr,N_i_Mtr);
% 
%         Inhibitory_Mtr = [IE_Mtr;II_Mtr];
%         Excitatory_Mtr = [EE_Mtr;EI_Mtr];
%         W_Mtr = horzcat(Excitatory_Mtr,Inhibitory_Mtr);   
% 
%         II_Vsn = 1e-11*rand(N_i_Vsn,N_i_Vsn);
%         EE_Vsn = 1e-11*rand(N_e_Vsn,N_e_Vsn);
%         EI_Vsn = 1e-11*rand(N_i_Vsn,N_e_Vsn);
%         IE_Vsn = -1e-11*rand(N_e_Vsn,N_i_Vsn);
% 
%         Inhibitory_Vsn = [IE_Vsn;II_Vsn];
%         Excitatory_Vsn = [EE_Vsn;EI_Vsn];
%         W_Vsn = horzcat(Excitatory_Vsn,Inhibitory_Vsn); 
        
        
        %%% External Stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        Stimulus_Mtr = Sf*MtrInput(NN_Mtr,1,ST_v);
%         Stimulus_Vsn = Sf*TravelingBumpVsn(NN_Vsn,1,ST_v);  
        %%% add noise to external Stimulus
%         Stimulus_Vsn = sf*NoisyInputVsn(NN_Vsn,1,NoiseSD);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_ext_Mtr(EnId,1:N_e_Mtr,t) = repmat(Stimulus_Mtr(EnId,t),N_e_Mtr,1);
%         I_ext_Mtr(EnId,1:N_e_Vsn,t) = 0; %Baseline  
%         I_ext_Vsn(EnId,1:N_e_Vsn,t) = repmat(Stimulus_Vsn(EnId,t),N_e_Vsn,1);
          %%% In case of seperated visual stimulus needed.
        %%% add noise to external Stimulus
%         I_ext_Mtr(EnId,1:N_e_Mtr,t) = awgn(I_ext_Mtr(EnId,1:N_e_Mtr,t),200);
%         I_ext_Vsn(EnId,1:N_e_Vsn,t) = awgn(I_ext_Vsn(EnId,1:N_e_Vsn,t),200);  
     
        %%% Current added from other populations within a Neural Field %%%%%%%
        I_Field_Mtr(EnId,:,t) = repmat(WE_Mtr(EnId,t-1),N_Mtr,1);   
        I_Field_Vsn(EnId,:,t) = repmat(WE_Vsn(EnId,t-1),N_Vsn,1);  
     
        %%% Current added from other Neural Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_MtrToVsn(EnId,:,t) = repmat(WF_Mtr(EnId,t-1),N_Mtr,1); 
        I_VsnToMtr(EnId,:,t) = repmat(WF_Vsn(EnId,t-1),N_Vsn,1); 
        %%% Subjects control the visual stimulus by their own hand movements
%         I_ext_Vsn(EnId,1:N_e_Vsn,t) = I_ext_Vsn(EnId,1:N_e_Vsn,t) + I_MtrToVsn(EnId,1:N_e_Vsn,t);
%         shift = 5;
%         PoId = EnId + shift;
%         if (PoId > NN_Mtr)
%             PoId = PoId - NN_Mtr;
%         end
        if (t > VM_delay)
           I_ext_Vsn(EnId,1:N_e_Vsn,t) = 1e-3*I_MtrToVsn(EnId,1:N_e_Vsn,t - VM_delay);
          %%% add noise to external Stimulus(This part should be evolved)
%           I_ext_Vsn(EnId,1:N_e_Vsn,t) = I_ext_Vsn(EnId,1:N_e_Vsn,t) + (mean(I_ext_Vsn(:))*randn(1)); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Synaptic Current calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %%% satisfy Exponential decay of Excitatory Conductance %%%%%%%%%%%%%%
        g_syn_e_Mtr(EnId,:,t) = g_syn_e_Mtr(EnId,:,t-1) - (dt/tau_syn_e_Mtr)*g_syn_e_Mtr(EnId,:,t-1);
        g_syn_e_Vsn(EnId,:,t) = g_syn_e_Vsn(EnId,:,t-1) - (dt/tau_syn_e_Vsn)*g_syn_e_Vsn(EnId,:,t-1);
        %%% Calculating Excitatory Synaptic Current%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        I_syn_e_Mtr(EnId,:,t) = sum(g_syn_e_Mtr(EnId,:,t-1) * (v_Mtr(EnId,:,t-1)-E_e_Mtr)');
        I_syn_e_Vsn(EnId,:,t) = sum(g_syn_e_Vsn(EnId,:,t-1) * (v_Vsn(EnId,:,t-1)-E_e_Vsn)');
        %%% satisfy Exponential decay of inhibitory Conductance %%%%%%%%%%%%%%
        g_syn_i_Mtr(EnId,:,t) = g_syn_i_Mtr(EnId,:,t-1) - (dt/tau_syn_i_Mtr)*g_syn_i_Mtr(EnId,:,t-1);
        g_syn_i_Vsn(EnId,:,t) = g_syn_i_Vsn(EnId,:,t-1) - (dt/tau_syn_i_Vsn)*g_syn_i_Vsn(EnId,:,t-1);
        %%% Calculating Inhibitory Synaptic Current %%%%%%%%%%%%%%%%%%%%%%%%%% 
        I_syn_i_Mtr(EnId,:,t) = sum(g_syn_i_Mtr(EnId,:,t-1) * (v_Mtr(EnId,:,t-1)-E_i_Mtr)');
        I_syn_i_Vsn(EnId,:,t) = sum(g_syn_i_Vsn(EnId,:,t-1) * (v_Vsn(EnId,:,t-1)-E_i_Vsn)');
      
     
     
        %%% aEIF differential equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Motor Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v_Mtr(EnId,:,t) = v_Mtr(EnId,:,t-1) + (dt/C_m_Mtr)*(g_L_Mtr*(E_L_Mtr-v_Mtr(EnId,:,t-1)) + g_L_Mtr*Delta_T_Mtr*exp((v_Mtr(EnId,:,t-1) - V_T_Mtr)/Delta_T_Mtr) - w_Mtr(EnId,:,t-1) ...
        + I_ext_Mtr(EnId,:,t) + I_syn_e_Mtr(EnId,:,t-1) + I_syn_i_Mtr(EnId,:,t-1) + I_Field_Mtr(EnId,:,t)) + I_VsnToMtr(EnId,:,t);                            
        w_Mtr(EnId,:,t) = w_Mtr(EnId,:,t-1) + (dt/tau_w_Mtr)*(a_s_Mtr.*(v_Mtr(EnId,:,t)-E_L_Mtr) - w_Mtr(EnId,:,t)); 
        %%% Vision Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        v_Vsn(EnId,:,t) = v_Vsn(EnId,:,t-1) + (dt/C_m_Vsn)*(g_L_Vsn*(E_L_Vsn-v_Vsn(EnId,:,t-1)) + g_L_Vsn*Delta_T_Vsn*exp((v_Vsn(EnId,:,t-1) - V_T_Vsn)/Delta_T_Vsn) - w_Vsn(EnId,:,t-1) ...
        + I_ext_Vsn(EnId,:,t) + I_syn_e_Vsn(EnId,:,t-1) + I_syn_i_Vsn(EnId,:,t-1) + I_Field_Vsn(EnId,:,t)) + I_MtrToVsn(EnId,:,t);
        %%%  
        w_Vsn(EnId,:,t) = w_Vsn(EnId,:,t-1) + (dt/tau_w_Vsn)*(a_s_Vsn.*(v_Vsn(EnId,:,t)-E_L_Vsn) - w_Vsn(EnId,:,t-1));
        %%%%%%%%%%%%%%%% When a neuron spikes!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fired_Mtr = find(v_Mtr(EnId,:,t) > V_peak_Mtr);
%       SpikingNeurons_Mtr{EnId,t} = [SpikingNeurons_Mtr{EnId,t}; t+0*Fired_Mtr,Fired_Mtr];
        SpikingNeurons_Mtr{EnId,t} = Fired_Mtr;
        Fired_Mtr = 0; 
        Fired_Vsn = find(v_Vsn(EnId,:,t) > V_peak_Vsn);
        SpikingNeurons_Vsn{EnId,t} = Fired_Vsn;
        Fired_Vsn = 0;
        %%% Excitatory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Motor Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fired_e_Mtr = find(v_Mtr(EnId,1:N_e_Mtr,t) > V_peak_Mtr);
        %%% Reset Conditions! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v_Mtr(EnId,Fired_e_Mtr,t) = E_L_Mtr; 
        w_Mtr(EnId,Fired_e_Mtr,t) = w_Mtr(EnId,Fired_e_Mtr,t-1)+ b_s_Mtr; 
        g_syn_e_Mtr(EnId,Fired_e_Mtr,t) = g_syn_e_Mtr(EnId,Fired_e_Mtr,t-1) + sum(W_Mtr(:,Fired_e_Mtr));
        FiringNeurons_e_Mtr(EnId,t) = length(Fired_e_Mtr);
        %%% Vision Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fired_e_Vsn = find(v_Vsn(EnId,1:N_e_Vsn,t) > V_peak_Vsn);
        %%% Reset Conditions! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v_Vsn(EnId,Fired_e_Vsn,t) = E_L_Vsn;
        w_Vsn(EnId,Fired_e_Vsn,t) = w_Vsn(EnId,Fired_e_Vsn,t-1)+ b_s_Vsn;
        g_syn_e_Vsn(EnId,Fired_e_Vsn,t) = g_syn_e_Vsn(EnId,Fired_e_Vsn,t-1) + sum(W_Vsn(:,Fired_e_Vsn));
        FiringNeurons_e_Vsn(EnId,t) = length(Fired_e_Vsn);
        %%% Inhibitory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Motor Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fired_i_Mtr = find(v_Mtr(EnId,N_e_Mtr+1:end,t) > V_peak_Mtr);
        %%% Reset Conditions! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        v_Mtr(EnId,Fired_i_Mtr+N_e_Mtr,t) = E_L_Mtr;
        w_Mtr(EnId,Fired_i_Mtr+N_e_Mtr,t) = w_Mtr(EnId,Fired_i_Mtr+N_e_Mtr,t-1)+ b_s_Mtr;
        g_syn_i_Mtr(EnId,N_e_Mtr+Fired_i_Mtr,t) = g_syn_i_Mtr(EnId,N_e_Mtr+Fired_i_Mtr,t-1) + sum(W_Mtr(:,N_e_Mtr+Fired_i_Mtr));
        FiringNeurons_i_Mtr(EnId,t) = length(Fired_i_Mtr); 
        %%% Vision Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fired_i_Vsn = find(v_Vsn(EnId,N_e_Vsn+1:end,t) > V_peak_Vsn);
        %%% Reset Conditions! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        v_Vsn(EnId,Fired_i_Vsn+N_e_Vsn,t) = E_L_Vsn;
        w_Vsn(EnId,Fired_i_Vsn+N_e_Vsn,t) = w_Vsn(EnId,Fired_i_Vsn+N_e_Vsn,t-1)+ b_s_Vsn;
        g_syn_i_Vsn(EnId,N_e_Vsn+Fired_i_Vsn,t) = g_syn_i_Vsn(EnId,N_e_Vsn+Fired_i_Vsn,t-1) + sum(W_Vsn(:,N_e_Vsn+Fired_i_Vsn));
        FiringNeurons_i_Vsn(EnId,t) = length(Fired_i_Vsn);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Population Activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Motor Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_e_Mtr(EnId,t) = FiringNeurons_e_Mtr(EnId,t)/(N_Mtr);        
        A_i_Mtr(EnId,t) = FiringNeurons_i_Mtr(EnId,t)/(N_Mtr);
        A_Mtr(EnId,t) = A_e_Mtr(EnId,t);  % Ensemble activity   
        %%% Vision Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_e_Vsn(EnId,t) = FiringNeurons_e_Vsn(EnId,t)/(N_Vsn);        
        A_i_Vsn(EnId,t) = FiringNeurons_i_Vsn(EnId,t)/(N_Vsn);
        A_Vsn(EnId,t) = A_e_Vsn(EnId,t);  % Ensemble activity 
      
    end
    
    %%% Applying connections kernel within a Neural Field %%%%%
    WE_Mtr(:,t) = WW_Mtr*A_Mtr(:,t);
    WE_Vsn(:,t) = WW_Vsn*A_Vsn(:,t);
    %%% Applying connections kernel between Neural Fields %%%%%
    WF_Mtr(:,t) = WWW_Mtr*A_Mtr(:,t);
    WF_Vsn(:,t) = WWW_Vsn*A_Vsn(:,t);
    
    waitbar(t / length(T))
    
end    
%%   Gaussian Filter
sigma = 10;
G_size = 100;
mu = 0;
for EnId = 1:NN_Mtr
    
%     A_fnl_Mtr(EnId,:) = conv (A_Mtr(EnId,:), GaussianFilter(mu,sigma,size), 'same');
A_fnl_Mtr(EnId,:) = conv (A_Mtr(EnId,:), GaussianFilter(sigma,G_size), 'same');
%     A_fnl_Vsn(EnId,:) = conv (A_Vsn(EnId,:), GaussianFilter(mu,sigma,size), 'same');
A_fnl_Vsn(EnId,:) = conv (A_Vsn(EnId,:), GaussianFilter(sigma,G_size), 'same');
       
end

close(h)
%% Capture the video here :
MaxActivity = 5*N_i_Mtr/100;
nframe=t;
mov(1:nframe)= struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren')
for k = 1:nframe

    plot(A_fnl_Vsn(:,k),'b','LineWidth',2.5)
    legend('Vision','Location','northeast')
    ylim([0,MaxActivity]);
    hold on
    plot(A_fnl_Mtr(:,k),'r','LineWidth',2.5)
    ylim([0,MaxActivity]);
    hold off
    legend('Visual activity','Motor activity','Location','northeast')
    ylim([0,MaxActivity]);
    title(['Time = ',num2str(k),'  (ms)']);
    ylabel('Neural Field Activity')
    xlabel('Ensemble Number')
    mov(k)= getframe(gcf);
  
end
%  movie2avi(mov, 'C:\Users\Mohammad\Desktop\Sync.avi', 'compression', 'None');
% movie2avi(mov, 'C:\Users\Mohammad\Desktop\test.avi', 'compression', 'None');
%  movie2avi(mov, 'C:\Users\Mohammad\Desktop\Delay1.avi', 'compression', 'None');
%% Population codes :
% PHI_Vsn = zeros(1,length(T));
% PHI_Mtr = zeros(1,length(T));
% x_n = linspace(0,.5,NN_Mtr/2);
% X_n = horzcat(x_n, fliplr(x_n));  %%% Due to ring topology
% % A_t_Vsn = A_fnl_Vsn.*floor(0.5*(1+sign(A_e_Vsn)));  %%% Threshold_Linear
% % A_t_Mtr = A_fnl_Mtr.*floor(0.5*(1+sign(A_e_Mtr)));  %%% Threshold_Linear
% % A_t_Vsn = sigmoid(A_fnl_Vsn);  %%% Threshold_Sigmoid
% % A_t_Mtr = sigmoid(A_fnl_Mtr);  %%% Threshold_Sigmoid
% 
% for tt = 1:length(T)
%     PHI_Vsn(1,tt) = X_n*A_fnl_Vsn(:,tt) / sum(A_fnl_Vsn(:,tt));  %%% Linear
%     PHI_Mtr(1,tt) = X_n*A_fnl_Mtr(:,tt) / sum(A_fnl_Mtr(:,tt));  %%% Linear
% %     PHI_Vsn(1,tt) = X_n*A_t_Vsn(:,tt) / sum(A_t_Vsn(:,tt));
% %     PHI_Mtr(1,tt) = X_n*A_t_Mtr(:,tt) / sum(A_t_Mtr(:,tt));
% end
% figure
% plot(T,PHI_Vsn,'-b')
% hold on 
% plot(T,PHI_Mtr,'.r')
% legend('Vision','Motor','Location','northwest')
% title('Population Code');
% ylabel('Population Activity')
% xlabel('Time [ms]')
% hold off
%% Action Detection :
% % Noise should be added at the output, also calculation analytically. 
% Threshold = 0.35;
% % Threshold = 0.45;
% WinSize = 10;                                       %%% Detection window size
% [sumSignal,detectionrate] = DetectionRate(A_fnl_Vsn,WinSize,Threshold,NoiseSD,VM_delay);
% detectionrate = sum(A_fnl_Vsn(:)>Threshold);
% % plot(sum(A_fnl_Vsn>Threshold))
%% figures are here:
% figure
% % imagesc(A_fnl_Mtr)
% imagesc(A_fnl_Vsn)
% % surf(A_fnl_Vsn
% % surf(A_fnl_Mtr)
% colorbar
% % title('Motor Field Activity');
% title('Visual Field Activity');
% ylabel('Neural Field Activity')
% xlabel('Time [ms]')
%%
toc;
