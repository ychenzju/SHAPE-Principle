%% Simulation of single bubbles
% This script can also be used to simulate microbubble populations when properly configured.
% This script has been tailed to single bubble simulations.
% This script is shared for manuscript "Explaining the correlation between Subharmonic Amplitude ..."
% This script can be configured to simulate multiple cases, mainly through
% -------------------------------------------------------------------------
% bubble_model            = 'EEM';
% rads = [1:0.5:3,4,5] * 1e-6;  %for multiple single bubbles [1:0.5:3,4,5]  
% povs = [0] * 1e3;             %0-25KPa [0,10,25]
% frqs = [2.5] * 1e6;           %1-10MHz [1:0.5:6]
% pacs = [350] * 1e3;           %10-800kPa [10:10:50,100:50:800]
% pulse_negative = 1;           %only used for non-modulated pulses
% to_solve_bubble_resonance:    %1: use small-signal simulation to estimate
%                                   the resonance frequencies
% -------------------------------------------------------------------------
% Please contact ychen.zju@outlook.com or yao.chen1@ge.com for discussion.
% Copyright reserved by Yao Chen. 2024-07-01

%% environmental setting
format short g
format compact
% close all;
clearvars;
global plot_on
plot_on = 0;

%% variables to be stored for simulation
Treceived = {};
Psimulated = {};
Preceived = {};
Tincident = {};
Rincident = {};
Tnodes = {};
Pnodes = {};
PreceivedPI = {};
PreceivedAM = {};

global ATTEN TotalBubbles Xnodes Ynodes Znodes TotalNodes R_rand Rs Xs Ys Zs DX DP
global incitime incident
global Rinit Finit Fresn R00 Pov Pac Frq bubble_model
global Pia Psi T R Rdot Rdot2 P R00 R0 tt r00
global recv_time recv_sigl
global tx_window_beg tx_window_end ix_window_beg ix_window_end
global rx_window_beg rx_window_end sg_window_beg sg_window_end


global probe_fs probe_fc probe_bandpassFreq1 probe_bandpassFreq2 probe_bandpassfilter
global probe_lead_cycles probe_cease_cycles probe_pulse_repeats probe_pulse_type probe_pulse_profiling probe_pulse_polarity probe_pulse_length
global probe_tran_enabled probe_recv_enabled
global atten_tran_enabled atten_recv_enabled atten_water_alpha atten_water_gamma atten_blood_beta
global eem_sigma_0 eem_E_s_0 eem_alpha_s eem_Kappa_s
global marmottant_chi marmottant_sigma_R00 marmottant_Kappa_s
global churchhoff_G_s churchhoff_mu_s churchhoff_d_sh_0 churchhoff_Kappa_s
global sound_speed Patm kappa mu_liquid rho_liquid sigma_water

%% Bubble Model Parameters
probe_fs                = 100e6;    % sampling rate
probe_fc                = 2.5e6;    % driving frequency
probe_bandpassFreq1     = 1e6;
probe_bandpassFreq2     = 10e6;
probe_bandpassfilter    = designfilt('bandpassiir','FilterOrder',10,...
    'HalfPowerFrequency1',probe_bandpassFreq1,...
    'HalfPowerFrequency2',probe_bandpassFreq2,...
    'SampleRate',probe_fs);
probe_lead_cycles       = 0;
probe_cease_cycles      = 12;
probe_pulse_repeats     = 1;
probe_pulse_type        = 'sin';
probe_pulse_profiling   = 'step';
probe_pulse_polarity    = 'bipolar';
probe_pulse_length      = 3;
probe_tran_enabled      = 1;
probe_recv_enabled      = 1;
atten_tran_enabled      = 1;
atten_recv_enabled      = 1;
atten_water_alpha       = 0.002;
atten_water_gamma       = 2;
atten_blood_beta        = 0.18;

eem_sigma_0             = 0.019;    % N/m; constant interfacial tension
eem_E_s_0               = 0.55;     % N/m; elastic modulus of shell
eem_alpha_s             = 1.5;      % the coeffient to describe the progressively loss of elasticity with increasing buble bubble area fraction
                                    % 2.5 or 3 may be used for an early growing stage in calibration curv
eem_Kappa_s             = 1.2e-8;   % N*s/m; surface dilatational viscosity of shell
marmottant_chi          = 0.53;     % N/m; elastic modulus of shell
marmottant_sigma_R00    = 0.02;     % N/m; initial interfacial tension at Patm
marmottant_Kappa_s      = 1.2e-8;   % Ns/m; change with radius, will be real-time calculated.
churchhoff_G_s          = 52e6;     % Pa; shear modulus of shell
churchhoff_mu_s         = 0.99;     % Ns/m2 = Pa*s; shear viscosity of shell
churchhoff_d_sh_0       = 4e-9;     % m; shell width
churchhoff_Kappa_s      = 1.2e-8;   % Ns/m; change with radius, will be real-time calculated.

sound_speed             = 1485;     % m/s;
Patm                    = 1.01e5;   % Pascal = N/m^2
kappa                   = 1.07;     % polytropic exponential index
mu_liquid               = 0.001;    % N*s/m^2; liquid viscosity
rho_liquid              = 1000;     % kg/m^3; liquid density
sigma_water             = 0.072;    % N/m;

%% bubble size relevant parameters
BubbleConc              = 0.78e9;   % number of bubbles per mL
DilutionRatio           = 1/3000;   % dilution ratio 1:3000, concentration = 2.6e5 ppmL
lognorm_mean            = 0.39;
lognorm_std             = 0.50;

%% =====================================================================%%
%% Parameters to be adjusted for different simulation cases
%% bubble_model:    'EEM', 'Marmottant', 'Church-Hoff'
%% rads:            bubble radii
%% povs:            ambient overpressures
%% frqs:            excitation frequencies
%% pacs:            pulse magnitudes
%% pulse_negative:  selection of positive (0) and negative pulse (1)
%% to_solve_bubble_resonance:   1: use small-signal simulation to estimate
%%                                 the resonance frequencies
%% =====================================================================%%
bubble_model            = 'EEM';
rads = [1:0.5:3,4,5] * 1e-6; %for multiple single bubbles [1:0.5:3,4,5]  
povs = [0] * 1e3; %0-25KPa [0,10,25]
frqs = [2.5] * 1e6; %1-10MHz [1:0.5:6]
pacs = [350] * 1e3; %10-800kPa [10:10:50,100:50:800]

pulse_negative        = 1; %only used for non-modulated pulses
pulse_modulated       = 0; %when using modulated pulse, the pulse_negative should be set to 0
pulse_pi              = 0; %pulse inversion for 2nd-harmonic imaging
pulse_am              = 0; %amplitude modulation for 2nd-harmonic imaging
pulse_cps             = 0; %contrast pulse sequencing for 2nd-harmonic imaging
modulating_pulse = [...
    1,      1,      1,      1,      1,      1;...   %magnitude
    1,      1,      1,      1,      1,      1;...   %frequency
    0,      1,      0,      1,      0,     1];     %phase: 1: 180度

% 通过小摄动激励仿真获得微泡的共振频率
to_solve_bubble_resonance = 0; %enable or disable the perturbation-simulation


%% Simulation Settings
%% Combined settings with microbubble populations

% seed = 1; rng(seed); GeneratedSeeds = randi(1000,2,50); %Monte-Carlo simulation
% GeneratedSeeds = GeneratedSeeds(:,21:30);
% GeneratedSeeds = [1999 * ones(1,10); (1122:100:2022)]; %seeds for random locations
% GeneratedSeeds = [(1999:100:2899); 1122 * ones(1,10)]; %seeds for random radius
GeneratedSeeds = [1999; 1122]; %fixed seed for bubble location and radius

for testcase = 1:size(GeneratedSeeds,2)
    Rseed = GeneratedSeeds(1,testcase); %seed for random bubble radius
    Lseed = GeneratedSeeds(2,testcase); %seed for random bubble locations
    
    %% single bubble definition
    single_bubble_simu      = 1;
    if single_bubble_simu
        %% following parameters for single bubble simulation
        Probe2Vessel_Depth      = 4e-2;     % m;
        Vessel2ROI_Depth        = 0e-3;     % m;
        ROI2Bottom_Depth        = 3e-3;     %m
        ROI_Xwidth              = 0e-3;     % m;
        ROI_Ywidth              = 0e-3;     % m;
        ROI_Zwidth              = 3e-3;     % m;
        TotalBubbles = 1;
        Xnodes = 1; Ynodes = 1; Znodes = 1;
        R_rand = 3;
        Rs = rads;
        TotalNodes = length(rads);
        Xs = 0; Ys = 0; Zs = -Probe2Vessel_Depth;
        DX = 1; DP = Probe2Vessel_Depth * ones(1,TotalNodes);
        [~,water_attn] = nonlinear_attenuation(atten_water_alpha,...
            atten_water_gamma,...
            Probe2Vessel_Depth*1e2,...
            1/probe_fs);
        ATTEN = ones(TotalNodes,1) * water_attn(1:100)';
    else
        %% following parameters for bubble populations
        %% ROI relevant parameters
        Probe2Vessel_Depth      = 4e-2;     % m;
        Vessel2ROI_Depth        = 3e-3;     % m;
        ROI2Bottom_Depth        = 3e-3;     % m;
        ROI_Xwidth              = 0e-3;     % m;
        ROI_Ywidth              = 0e-3;     % m;
        ROI_Zwidth              = 3e-3;     % m;
        gridInterval            = 1e-2 / ((BubbleConc * DilutionRatio)^(1/3));  % m;
        
        config_str = [num2str(ROI_Xwidth*1e3),num2str(ROI_Ywidth*1e3),num2str(ROI_Zwidth*1e3)];
        montecarlo_str = ['R3MC_R',num2str(Rseed),'_L',num2str(Lseed)];
        if size(GeneratedSeeds,2)>=10
            population_bubble_file = ['population_',config_str,'_',montecarlo_str];
        elseif size(GeneratedSeeds,2)>1
            if length(unique(GeneratedSeeds(1,:)))==1
                population_bubble_file = ['population_',config_str,'_L',num2str(Lseed)];
            elseif length(unique(GeneratedSeeds(2,:)))==1
                population_bubble_file = ['population_',config_str,'_R',num2str(Rseed)];
            end
        else
            population_bubble_file= ['population_',config_str];
        end
        
        if length(povs)>1 && length(pacs)==1 && length(frqs)==1
            population_bubble_file = [population_bubble_file,'_Povs','.mat'];
        elseif length(povs)>1 && length(frqs)>1 && length(pacs)==1
            population_bubble_file = [population_bubble_file,'_Povs_Frqs','.mat'];
        elseif length(povs)>1 && length(pacs)>1 && length(frqs)==1
            population_bubble_file = [population_bubble_file,'_Povs_Pacs','.mat'];
        elseif length(povs)>1 && length(pacs)>1 && length(frqs)>1
            population_bubble_file = [population_bubble_file,'_Povs_Frqs_Pacs','.mat'];
        else
            population_bubble_file = [population_bubble_file,'.mat'];
        end
        disp(population_bubble_file);

        
        if ~exist(population_bubble_file,'file')
            
            %% population relevant parameters
            TotalBubbles = BubbleConc * DilutionRatio * (ROI_Xwidth * ROI_Ywidth * (ROI_Zwidth + Vessel2ROI_Depth + ROI2Bottom_Depth)) * 100^3;
            Xnodes = floor(ROI_Xwidth / gridInterval) + 1;
            Ynodes = floor(ROI_Ywidth / gridInterval) + 1;
            Znodes = floor((ROI_Zwidth + Vessel2ROI_Depth + ROI2Bottom_Depth) / gridInterval) + 1;
            TotalNodes = round(Xnodes * Ynodes * Znodes);
            
            %% radius randomly generated
            rng(Rseed);
            R_rand = lognrnd(lognorm_mean,lognorm_std,Xnodes,Ynodes,Znodes);
            Rs = reshape(R_rand,Xnodes*Ynodes*Znodes,1) * 1e-6; %m; bubble radius,
            
            %% bubble location randomly generated
            rng(Lseed);
            Xs = (rand(TotalNodes,1) - 0.5) .* ROI_Xwidth; %X-coordinates for all bubbles
            Ys = (rand(TotalNodes,1) - 0.5) .* ROI_Ywidth; %Y-coordinates for all bubbles
            Zs = -((rand(TotalNodes,1) * (Vessel2ROI_Depth + ROI_Zwidth + ROI2Bottom_Depth)) + Probe2Vessel_Depth); %Z-coordinates for all bubbles
            Zs = sort(Zs,'descend');
                       
            overlap_remained = 1;
            while overlap_remained
                %% Remove the bubbles that are too close
                overlapped_ndx = 0;
                jnode = 1;
                while jnode<=length(Rs)
                    for inode = 1:length(Rs)
                        if jnode~=inode
                            distance = sqrt((Xs(jnode)-Xs(inode))^2 + (Ys(jnode)-Ys(inode))^2 + (Zs(jnode)-Zs(inode))^2);
                        else
                            distance = 1000;
                        end
                        if distance<( 10e-6 + eps)
                            if Rs(jnode)<Rs(inode)
                                overlapped_ndx = jnode;
                            else
                                overlapped_ndx = inode;
                            end
                            break;
                        end
                    end
                    jnode = jnode + 1;
                end
                if overlapped_ndx==1
                    Rs = Rs(2:end);
                    Xs = Xs(2:end);
                    Ys = Ys(2:end);
                    Zs = Zs(2:end);
                elseif overlapped_ndx==length(Rs)
                    Rs = Rs(1:end-1);
                    Xs = Xs(1:end-1);
                    Ys = Ys(1:end-1);
                    Zs = Zs(1:end-1);
                elseif overlapped_ndx>0
                    Rs = [Rs(1:overlapped_ndx-1);Rs(overlapped_ndx+1:end)];
                    Xs = [Xs(1:overlapped_ndx-1);Xs(overlapped_ndx+1:end)];
                    Ys = [Ys(1:overlapped_ndx-1);Ys(overlapped_ndx+1:end)];
                    Zs = [Zs(1:overlapped_ndx-1);Zs(overlapped_ndx+1:end)];
                else
                    overlap_remained = 0;
                end
                disp(['removing overlapped nodes = [',num2str(overlapped_ndx),']']);
                
            end
            
            %% Calculate the inter-bubble distance and probe-to-bubble distance
            DX = zeros(length(Rs),length(Rs)); %inter-bubble distances
            DP = zeros(length(Rs),1); %probe-to-bubble distances, the coordinate of probe is (0,0,0);
            for jnode = 1:length(Rs)
                for inode = 1:length(Rs)
                    DX(jnode,inode) = sqrt((Xs(jnode)-Xs(inode))^2 + (Ys(jnode)-Ys(inode))^2 + (Zs(jnode)-Zs(inode))^2);
                end
                DP(jnode) = sqrt((Xs(jnode))^2 + (Ys(jnode))^2 + (Zs(jnode))^2);
                DX(jnode,jnode) = 1;
            end
            TotalNodes = length(Rs);
            
            %% Calculate the attenuation model for each bubble
            ATTEN = zeros(TotalNodes,100);
            for jnode = 1:TotalNodes
                disp([jnode, TotalNodes]);
                [~,water_attn] = nonlinear_attenuation(atten_water_alpha,...
                    atten_water_gamma,...
                    DP(jnode)*1e2 * abs(Probe2Vessel_Depth/Zs(jnode)),...
                    1/probe_fs);
                water_attn = water_attn(1:100);
                [~,blood_attn] = linear_attenuation(atten_blood_beta,...
                    DP(jnode)*1e2 * (1 - abs(Probe2Vessel_Depth/Zs(jnode))),...
                    1/probe_fs);
                blood_attn = blood_attn(1:100);
                attn = conv(water_attn,blood_attn);
                ATTEN(jnode,:) = attn(1:100);
            end
            % save population_bubbles ATTEN TotalBubbles Xnodes Ynodes Znodes TotalNodes R_rand Rs Xs Ys Zs DX DP ATTEN
            save(population_bubble_file, 'ATTEN', 'TotalBubbles', ...
                'Xnodes', 'Ynodes', 'Znodes', 'TotalNodes', 'R_rand', 'Rs', ...
                'Xs', 'Ys', 'Zs', ...
                'DX', 'DP', 'ATTEN');
        else
            % load population_bubbles
            load(population_bubble_file);
        end
    end
    
    %% Simulating monodisperse populations
    %% assume all bubbles to be 3um
    % Rs = ones(size(Rs)) * 3e-6;
    
    %% Calculate the initial radius under hydrostatic pressure
    disp('Initiation of Bubbles')
    R00 = Rs; %Rs has unit of meter
    solver_options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);
    Rinit = zeros(length(R00),length(povs));
    Finit = Rinit;
    
    % 通过小摄动激励仿真获得微泡的共振频率
    if to_solve_bubble_resonance==1
        Fresn = Finit;
    end
    
    % 计算静息态微泡半径及其共振频率
    for rad_ndx = 1:length(R00)
        for pov_ndx = 1:length(povs)
            
            Pov = povs(pov_ndx);
            Rinit(rad_ndx,pov_ndx) = fsolve(@(x) gaslaw(x,R00(rad_ndx),povs(pov_ndx),bubble_model), R00(rad_ndx), solver_options);
            Finit(rad_ndx,pov_ndx) = resonanceFreq(Rinit(rad_ndx,pov_ndx),bubble_model);
            
            % 通过小摄动激励获得微泡的共振频率
            if to_solve_bubble_resonance==1
                disp([R00(rad_ndx)*1e6 povs(pov_ndx)]);
                Festm = Finit(rad_ndx,pov_ndx);
                if pov_ndx==1
                    Festm = solve_resonance(R00(rad_ndx),Rinit(rad_ndx,pov_ndx),1e6:0.5e6:10e6);
                    freqspan = Festm-0.5e6:0.01e6:Festm+0.5e6;
                else
                    Festm = Fresn(rad_ndx,pov_ndx-1);
                    freqspan = Festm-0.1e6:0.01e6:Festm+0.4e6;
                end
                Fresn(rad_ndx,pov_ndx) = solve_resonance(R00(rad_ndx),Rinit(rad_ndx,pov_ndx),freqspan);
            end
        end
    end
    % 单独保存小摄动激励仿真获得的微泡共振频率
    resn_filename = ['population1_bubble_resonance_',bubble_model,'.mat'];
    if to_solve_bubble_resonance
        save(resn_filename, 'rads', 'povs', 'R00', 'Rinit', 'Finit', 'Fresn');
        return
    end
    
    %% simulate the population model
    disp('Simulation of bubble models')

    % for all excitation frequency
    recv = cell(length(povs),length(frqs),length(pacs));
    for frq_ndx = 1:length(frqs)
        Frq = frqs(frq_ndx);
        
        % calculate the time windows
        tx_window_beg = probe_lead_cycles / Frq;
        tx_window_end = tx_window_beg + 2 * ROI_Zwidth / sound_speed;
        ix_window_beg = probe_lead_cycles / Frq + (Probe2Vessel_Depth) / sound_speed;
        ix_window_end = ix_window_beg + (Vessel2ROI_Depth + ROI_Zwidth + ROI2Bottom_Depth) / sound_speed; % + Param.probe.pulse_length / Param.probe.fc;
        rx_window_beg = probe_lead_cycles / Frq + 2 * Probe2Vessel_Depth / sound_speed; %the begin time of the time window
        rx_window_end = rx_window_beg + 2 * (Vessel2ROI_Depth + ROI_Zwidth + ROI2Bottom_Depth + 3e-3) / sound_speed; % + Param.probe.pulse_length / Param.probe.fc + ...
        sg_window_beg = probe_lead_cycles / Frq + 2 * (Probe2Vessel_Depth + Vessel2ROI_Depth) / sound_speed; %the begin time of the time window
        sg_window_end = sg_window_beg + 2 * ROI_Zwidth / sound_speed; % + Param.probe.pulse_length / Param.probe.fc + ...
        
        % generate the input signals
        if ~pulse_modulated
            pulse_option = probe_pulse_length;
        else
            pulse_option = modulating_pulse;
        end
        [tex, uex] = gnrtPulseSignal_x(Frq,pulse_option,...
            probe_pulse_repeats,...
            probe_fs,...
            probe_pulse_type,...
            probe_pulse_polarity,...
            probe_pulse_profiling,...
            probe_lead_cycles,probe_cease_cycles); %pressure signal
        
        %transmitting acoustic pressure
        if pulse_negative
            uex = -uex;
        end
        
        xex0 = uex;
        if probe_tran_enabled
            % simulate the probe
            % xs = lsim(Param.probe.sys_probe,us,ts);
            xex0 = filter(probe_bandpassfilter,xex0);
        end
        
        transmit_wv = [tex(:), uex(:), xex0(:)];
        
        % for all excitation magnitudes
        for pac_ndx = 1:length(pacs)
            Pac = pacs(pac_ndx);
            
            % generate transmitting acoustic signal
            xex = xex0 * Pac;
            
            % for all ambient overpressures
            for pov_ndx = 1:length(povs)
                Pov = povs(pov_ndx);
                disp([Frq/1e6, Pac/1e3, Pov/1e3]);
                
                %% initialization
                simulation_tbeg = ix_window_beg;
                simulation_tend = ix_window_end + 6e-6; %64 microseconds
                
                tt = (simulation_tbeg:1/probe_fs:simulation_tend)';
                xt = interp1((tex-tex(1)+simulation_tbeg),xex,tt,'linear',0);
                
                incitime0 = zeros(length(tt),TotalNodes);
                incident0 = zeros(length(tt),TotalNodes);
                for inode = 1:TotalNodes
                    ts = (tt - tt(1)) + DP(inode) / sound_speed;
                    if atten_tran_enabled
                        xs = conv(xt,ATTEN(inode,:)');
                        xs = xs(1:length(ts));
                    else
                        xs = xt;
                    end
                    
                    incitime0(:,inode) = ts(:);
                    incident0(:,inode) = xs(:);
                end
                
                incitime = incitime0;
                incident = incident0;
                
                R0 = Rinit(:,pov_ndx);
                solve_model();
                
                %% store the simulation
                Treceived{pov_ndx,frq_ndx,pac_ndx} = recv_time;  %receiving time series
                Psimulated{pov_ndx,frq_ndx,pac_ndx} = recv_sigl; %bubble-wise receiving pressure
                Preceived{pov_ndx,frq_ndx,pac_ndx} = sum(recv_sigl,2); %summed receiving pressure
                
                Tincident{pov_ndx,frq_ndx,pac_ndx} = tt; %incident time series
                Rincident{pov_ndx,frq_ndx,pac_ndx} = R;  %bubble-wise radius dynamics at incident time
                
                Tnodes{pov_ndx,frq_ndx,pac_ndx} = T; %bubble-wise simulation time series
                Pnodes{pov_ndx,frq_ndx,pac_ndx} = P; %bubble-wise scattering pressure
                
                if pulse_pi==1
                    incident = -1 * incident0;
                    R0 = Rinit(:,pov_ndx);
                    solve_model();
                    PreceivedPI{pov_ndx,frq_ndx,pac_ndx} = Preceived{pov_ndx,frq_ndx,pac_ndx} + sum(recv_sigl,2);
                end
                
                if pulse_am==1
                    incident = 0.5 * incident0;
                    R0 = Rinit(:,pov_ndx);
                    solve_model();
                    PreceivedAM{pov_ndx,frq_ndx,pac_ndx} = Preceived{pov_ndx,frq_ndx,pac_ndx} - 2.0 * sum(recv_sigl,2);
                end
                
                if pulse_cps==1
                    incident = -0.5 * incident0;
                    R0 = Rinit(:,pov_ndx);
                    solve_model();
                    PreceivedCPS{pov_ndx,frq_ndx,pac_ndx} = Preceived{pov_ndx,frq_ndx,pac_ndx} + 2.0 * sum(recv_sigl,2);
                end
                
            end %of povs
        end %of pacs
    end %of frqs
    
    if exist('population_bubble_file','var') && exist(population_bubble_file,'file')
        save(population_bubble_file,'*');
    end
    
end %end of Monte-Carlo Simulation

if size(GeneratedSeeds,2)>1
    disp('This is a Monte-Carlo simulation. Return without plotting!');
    return
end
if length(povs)>1 && length(frqs)>1 && length(pacs)>1
    disp('This is an overall simulation for data collection!')
    save population1 *
    return
end

%% plotting the results
plot_population1;

%% function of solving resonance frequency through small-perturbation simulation
function [resonance_freq] = solve_resonance(aR00,aRinit,freqs)
global bubble_model
global Patm Pov Pac Frq Pia Psi kappa mu_liquid rho_liquid sound_speed
global r00

magds = zeros(size(freqs));
for ndx=1:length(freqs)
    disp([ndx, freqs(ndx)/1e6]);
    afreq = freqs(ndx);
    sampling_rate = afreq * 20;
    % 通过仿真来获得微泡在特定平衡状态下的共振频率
    [atex, auex] = gnrtPulseSignal_x(afreq,...
        32,...              %pulse length
        1,...               %probe_pulse_repeats,...
        sampling_rate,...   %sampling rate
        'sin',...           %probe_pulse_type,...
        'bipolar',...       %probe_pulse_polarity,...
        'step',...          %probe_pulse_profiling,...
        0,...               %probe_lead_cycles,...
        0);                 %probe_cease_cycles
    auex = auex * 1;
    
    solver_options = odeset('RelTol',1e-4,'AbsTol',1e-10,'InitialStep',0.1/sampling_rate,'MaxStep',0.5/sampling_rate);
    radiidyn = zeros(length(atex),2);
    for step=1:length(atex)
        t = atex(step);
        r00 = aR00; %used in odeFcn
        Pia = interp1(atex,auex,t,'linear',0); %used in odeFcn
        Psi = 0; %used in odeFcn
        if step==1
            y0 = [aRinit;0];
        else
            y0 = radiidyn(step-1,:);
        end
        [~,yt] = ode45(@(t,y) odeFcn(t,y), [t-1/sampling_rate,t], y0, solver_options);
        radiidyn(step,:) = yt(end,:);
    end
    
%     %画图调试
%     figure(900); plot(atex(:)*1e6,radiidyn(:,1)*1e6);
%     title('Bubble Radius Dyanmics');
%     xlabel('Time (\mus)'); ylabel('Radii (\mum)');
%     pause(0.5);
%     
    segment_beg = round(length(atex)/2);
    magds(ndx) = max(radiidyn(segment_beg:end,1)) - min(radiidyn(segment_beg:end,1));
    
    if magds(ndx)<max(magds(1:ndx))*0.8
        break;
    end
end

figure(900); plot(freqs,magds);
[~,maxidx] = max(magds);
resonance_freq = freqs(maxidx(1));
end

%% function of solving bubble population
function solve_model()
% prepared variables
global tt incitime incident
global TotalNodes DX DP ATTEN probe_bandpassfilter
global Frq Pac Pov Rad Patm
global probe_fs rho_liquid sound_speed
global rx_window_beg rx_window_end
global T R Rdot Rdot2 P Pia Psi R00 R0 r00
global recv_time recv_sigl
global atten_recv_enabled probe_recv_enabled

% predefine the variables
T     = tt(:) * ones(1,TotalNodes);
R     = ones(length(tt),1) * R0';
Rdot  = zeros(length(tt),TotalNodes);
Rdot2 = zeros(length(tt),TotalNodes);
P     = zeros(length(tt),TotalNodes);
R(1,:) = R0';

% solve the equation
solver_options = odeset('RelTol',1e-3,'AbsTol',1e-8,'InitialStep',0.1/probe_fs,'MaxStep',0.5/probe_fs);
for step=1:length(tt)
    t = tt(step);
    % disp(t);
    
    for inode=1:TotalNodes

        r00 = R00(inode); %used in odeFcn
        Pia = interp1(incitime(:,inode),incident(:,inode),t,'linear',0); %used in odeFcn
        Psi = 0; %used in odeFcn
%         for jnode=1:TotalNodes
%             if jnode~=inode
%                 tao = DX(inode,jnode)/sound_speed;
%                 if (t-tao) < incitime(1,jnode)
%                     pssum = 0;
%                 else
%                     R_jMinusTao = interp1(tt,R(:,jnode),t-tao,'linear',R0(jnode));
%                     Rdot_jMinusTao = interp1(tt,Rdot(:,jnode),t-tao,'linear',0);
%                     Rdot2_jMinusTao = interp1(tt,Rdot2(:,jnode),t-tao,'linear',0);
%                     pssum = rho_liquid * (R_jMinusTao^2 * Rdot2_jMinusTao + ...
%                         2.0 * R_jMinusTao * Rdot_jMinusTao^2) / DX(inode,jnode);
%                 end
%                 Psi = Psi + pssum;
%             end
%             
%         end
        % Psi = 0; %% assume no bubble interactions
        
        if step==1
            y0 = [R0(inode);0];
        else
            y0 = [R(step-1,inode); Rdot(step-1,inode)];
        end
        [~,yt] = ode45(@(t,y) odeFcn(t,y), [t-1/probe_fs,t], y0, solver_options);
        dydt = odeFcn(t,yt(end,:));
        
        if isnan(yt(end,1))
            disp('DEBUG: Incorrect ODE solving!');
            keyboard;
        end
        
        R(step,inode)     = yt(end,1);
        Rdot(step,inode)  = dydt(1);
        Rdot2(step,inode) = dydt(2);
        P(step,inode)     = rho_liquid * (R(step,inode) / DP(inode)) *...
            (2.0 * Rdot(step,inode)^2 + R(step,inode) * Rdot2(step,inode));
        T(step,inode)     = t + DP(inode) / sound_speed;
    end %end of Other Nodes
end %end of time steps

% reconstruct the receiving signal
recv_time = (rx_window_beg:1/probe_fs:rx_window_end)';
recv_sigl = zeros(length(recv_time),TotalNodes);
for inode = 1:TotalNodes
    ps = P(:,inode);
    if atten_recv_enabled
        ps = conv(ps,ATTEN(inode,:)');
        ps = ps(1:length(P(:,inode)));
    end
    
    if probe_recv_enabled
        ps = filter(probe_bandpassfilter,ps);
    end
    
    recv_sigl(:,inode) = interp1(T(:,inode),ps,recv_time,'linear',0);
end

end


%% function of [rdot,rdot2] = odeFcn(t,[r,rdot])
function dydt = odeFcn(t,y)
global bubble_model
global Patm Pov Pac Frq Pia Psi kappa mu_liquid rho_liquid sound_speed
global r00
% global R Rdot Rdot2 DX incitime incident inode tt TotalNodes

r = y(1);
rdot = y(2);
dydt = [0;0];

% Pia = interp1(incitime(:,inode),incident(:,inode),t,'linear',0);
% Psi = 0;
% for jnode=1:TotalNodes
%     if jnode~=inode
%         tao = DX(inode,jnode)/sound_speed;
%         R_jMinusTao = interp1(tt,R(:,jnode),t-tao,'linear',R(1,jnode));
%         Rdot_jMinusTao = interp1(tt,Rdot(:,jnode),t-tao,'linear',0);
%         Rdot2_jMinusTao = interp1(tt,Rdot2(:,jnode),t-tao,'linear',0);
%         pssum = rho_liquid * (R_jMinusTao^2 * Rdot2_jMinusTao + ...
%             2.0 * R_jMinusTao * Rdot_jMinusTao^2) / DX(inode,jnode);
%         Psi = Psi + pssum;
%     end
% end

switch bubble_model
    case 'EEM'
        [sigmaOfR0]             = rheologyFunc_eem(r00,r00); %at resting state
        [sigmaOfR,Kappa_s]      = rheologyFunc_eem(r,r00);
    case 'Marmottant'
        [sigmaOfR0]             = rheologyFunc_marmottant(r00,r00); %at resting state
        [sigmaOfR,Kappa_s]      = rheologyFunc_marmottant(r,r00);
    case 'Church-Hoff'
        [sigmaOfR0]             = rheologyFunc_churchhoff(r00,r00); %at resting state
        [sigmaOfR,Kappa_s]      = rheologyFunc_churchhoff(r,r00);
    otherwise
        disp('Param.model is not well defined!');
        return
end

P_G0 = (Patm + 2.0 * sigmaOfR0 / r00);

dydt(1) = rdot;
dydt(2) = ...
    ( ((P_G0 * (r00 / r)^(3.0 * kappa) * (1.0 - 3.0* kappa * rdot / sound_speed)) ...
    - 2.0 * sigmaOfR / r - 4.0 * Kappa_s * rdot / r^2 - 4 * mu_liquid * rdot / r - (Patm + Pov + Pia)) ...
    - Psi ...
    - (3.0/2.0 * rho_liquid * rdot^2) ) / (rho_liquid * r); %Rdotdot

end

%% Function: gaslaw
function eq = gaslaw(R,R00,Pov,bubble_model)
% function of ideal gas condition law to calculate the static bubble radius
% under different quadrostatic pressure
global Patm kappa

switch bubble_model
    case 'EEM'
        sigmaOfR = rheologyFunc_eem(R,R00); %model specific calculation
        sigmaOfR00 = rheologyFunc_eem(R00,R00);
    case 'Church-Hoff'
        sigmaOfR = rheologyFunc_churchhoff(R,R00); %model specific calculation
        sigmaOfR00 = rheologyFunc_churchhoff(R00,R00);
    case 'Marmottant'
        sigmaOfR = rheologyFunc_marmottant(R,R00); %model specific calculation
        sigmaOfR00 = rheologyFunc_marmottant(R00,R00);
    otherwise
        disp('Param.model is not well defined!');
        return;
end
eq = (Patm + Pov + 2.0 * sigmaOfR ./ R) - ...
    (Patm + 2.0 * sigmaOfR00 ./ R00) .* ((R00 ./ R) .^ (3.0 * kappa));
end
% R0 = fsolve(@(x), gaslaw_eem(x,R00), R00);

%% Function: sigmaFunc
% function of interfacial tension based on radius of microbubble
function [sigmaOfR, Kappa_s] = rheologyFunc_eem(R,R00)
global eem_sigma_0 eem_E_s_0 eem_alpha_s eem_Kappa_s

sigma_0     = eem_sigma_0;
E_s_0       = eem_E_s_0;
alpha_s     = eem_alpha_s;

N = length(R);
if N>1
    sigmaOfR = zeros(N,1);
    Kappa_s = zeros(N,1);
end

for ndx=1:N
    k1 = sqrt( 1.0 + ( 1.0 - sqrt(1.0 + 4.0 * sigma_0 * alpha_s / E_s_0) ) / (2 * alpha_s) );
    R_E = R00(ndx) / k1; %equilibrium radius
    beta = (R(ndx) / R_E)^2 - 1; %area fraction change
    E_s = E_s_0 * exp(-alpha_s * beta);
    
    sigmaOfR(ndx) = sigma_0 + E_s * beta;
    Kappa_s(ndx)  = eem_Kappa_s;
end
end

function [sigmaOfR, Kappa_s] = rheologyFunc_churchhoff(R,R00)

global churchhoff_G_s churchhoff_mu_s churchhoff_d_sh_0

G_s         = churchhoff_G_s;
mu_s        = churchhoff_mu_s;
d_sh_0      = churchhoff_d_sh_0;

N = length(R);
if N>1
    sigmaOfR = zeros(N,1);
    Kappa_s = zeros(N,1);
end

for ndx=1:N
    sigmaOfR(ndx) = (6.0 * G_s * d_sh_0) * (R00(ndx) / R(ndx))^2 * ( 1 - R00(ndx) / R(ndx));
    Kappa_s(ndx)  = (3.0 * mu_s * d_sh_0) * (R00(ndx) / R(ndx))^2;
end
end

function [sigmaOfR, Kappa_s] = rheologyFunc_marmottant(R0,R00)

global marmottant_chi marmottant_sigma_R00 marmottant_Kappa_s
global sigma_water

chi = marmottant_chi;
sigma_R00 = marmottant_sigma_R00;

N = length(R0);
if N>1
    sigmaOfR = zeros(N,1);
    Kappa_s = zeros(N,1);
end

for ndx=1:N
    
    %% Raw Marmottant model
    % R_buckling = R00 / sqrt(1 + sigma_R00 / chi);
    % R_rupture = R_buckling * sqrt(1 + sigma_water / chi);
    %% Paul Parameters
    % R_buckling = R00(ndx);
    % R_rupture = 1.5 * R_buckling;
    %% CY modified Marmottant model
    R_buckling = R00(ndx) * sqrt(1 - sigma_R00 / chi);
    R_rupture  = R00(ndx) * sqrt(1 + (sigma_water - sigma_R00) / chi);
    
    if (R0(ndx) <= R_buckling)
        sigmaOfR(ndx) = 0;
    elseif (R0(ndx) < R_rupture) && (R0(ndx) >= R_buckling)
        %% Raw Marmottant model
        % sigmaOfR(ndx) = chi * ((R0(ndx) / R00(ndx))^2 - 1);
        %% CY modified Marmottant model
        sigmaOfR(ndx) = sigma_R00 + chi * ((R0(ndx) / R00(ndx))^2 - 1);
    elseif (R0(ndx) >= R_rupture)
        sigmaOfR(ndx) = sigma_water;
    end
    
    Kappa_s(ndx) = marmottant_Kappa_s;
end
end

%% Function：resonanceFreq
function rsnFrq = resonanceFreq(R00,bubble_model)
switch bubble_model
    case 'EEM'
        rsnFrq = resonanceFreq_eem(R00);
    case 'Church-Hoff'
        rsnFrq = resonanceFreq_churchhoff(R00);
    case 'Marmottant'
        rsnFrq = resonanceFreq_marmottant(R00);
    otherwise
        disp('bubble_model is not well defined!');
        return;
end
end

function rsnFreq = resonanceFreq_eem(R00)
global eem_sigma_0 eem_E_s_0 eem_alpha_s
global Patm kappa rho_liquid

v1 = sqrt(1.0 + 4.0 * eem_sigma_0 * eem_alpha_s / eem_E_s_0);
v2 = 1.0 + 2.0 * eem_alpha_s - v1;
v3 = (2.0 * eem_E_s_0) ./ R00 * ( v1 / eem_alpha_s);
v4 = 3.0 * kappa * Patm + v3 * v2;
rsnFreq = sqrt(v4 / rho_liquid) / (2.0 * pi * R00);
end

function rsnFrq = resonanceFreq_churchhoff(R00)
global churchhoff_G_s churchhoff_d_sh_0
global Patm kappa rho_liquid

v1 = 12.0 * churchhoff_G_s * churchhoff_d_sh_0 ./ R00;
v2 = 3.0 * kappa * Patm + v1;
rsnFrq = sqrt(v2 / rho_liquid) / (2.0 * pi * R00);
end

function rsnFrq = resonanceFreq_marmottant(R00)
global marmottant_chi marmottant_sigma_R00
global Patm kappa rho_liquid

Kv = kappa * Patm + 4.0 / 3.0 * marmottant_chi ./ R00 + ...
    (3.0 * kappa - 1) / 3.0 * 2.0 * marmottant_sigma_R00 ./ R00;
v1 = 3.0 * Kv / rho_liquid;
rsnFrq = sqrt(v1) / (2.0 * pi * R00);
end


