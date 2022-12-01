%runRateChangeHarvestScenarios.m
%E.W. Tekwa Apr 26, 2022

%Define parameters for set of stochastic harvest runs, plotting biomass outcomes
%Each set contains reps for each of two treatments. Trajectories start at upper
%stable equilibrium biomass for initial bifurcation value (initVs). They then experience
%harvest rate noise (D), bifurcation parameter noise (DDV), bifurcation parameter
%directional change (amplitude: DV, rates of change: rates), and bifurcation parameter
%cycle (amplitude CV,frequency cycleFreq). Once directional change is
%completed, cycle and noise continues for an additional time (endTimes).
set(0,'defaulttextinterpreter','tex');
set(0, 'defaultAxesTickLabelInterpreter','tex');
set(0, 'defaultLegendInterpreter','tex');
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');

figs=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)]);


rng(1); %set random number generator seed

% C=1; %process marginal cost (marginal cost with catch)
% C2=0.02; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.02
% C3=0; %multiplier for diminishing marginal cost with catch: 0.5
%V=1; %reference price (multiplier for diminishing marginal benefit) - this
%is being varied in the plot
%V2=0; %constant marginal benefit with catch: 0.1
Ds=[0.005 0.05]; %gaussian process noise variance in F: 2
%Vsteps=200; %number of steps in V to plot

%run trajectories with parameter change
reps=100; %replicates per rate of change in V
%plotReps=3; %number of trajectories to plot
rates=-[0.01, 1]; %rates of change in V
cycleFreq=[0.01 1]; %frequency of cycles in V
initV=0.9; %1.5
DV=0; %total directional V change
CV=[0 0]; %cyclical amplitude in V
%CV=[0.2 0.2]; %cyclical amplitude in V: 0.25
%DDVs=[0.001 1]; %temporal variance in V
DDVs=[0 0]; %temporal variance in V
%dVdts=[rates(1)*ones(1,reps) rates(2)*ones(1,reps)]; %directional rate of change in V
%cycVs=[CV(1)*sin(2*pi*t/cycleFreq(1))*ones(1,reps)  CV(2)*sin(2*pi*t/cycleFreq(2))*ones(1,reps)]; %cyclical change in V
%endTimes=abs(DV./rates([2 1])); %constant environment time after change
endTimes=[10 10];
%simTimeStep=0.01; %time step in Euler–Maruyama stochastic simulation

%define scenarios in rows (parameters as elements or in columns):
% initVs=[];
% rates
% cycleFreq
% endTimes
% D
% DDV
% DV
% CV

% p=1; %rho, (-1 to Inf, with 1 being logistic growth and <1 moving Smsy closer to R=0) *note values other than 0.5 and 1 take longer to compute
% Vs=0; %linearly diminishing cost (increasing benefit) with stock size ("stock ecosystem service")
% Cf=0; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.01
% bet=-1; %risk aversion controlling diminishing returns
p=[1 1 1 1 0.5 0.5 0.5 0.5]; %rho, (-1 to Inf, with 1 being logistic growth and <1 moving Smsy closer to R=0) *note values other than 0.5 and 1 take longer to compute
Vs=[0 0 0 0 0 0.01 0 0.01]; %stock ecosystem service
Cf=[0 0.01 0 0 0 0 0.01 0.01]; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.01
bet=[1 1 5 -1 1 1 1 5]; %risk aversion
alph=[1 1 1 1 1 1 0.5 0.5]; %effort or labour elasticity

for pl=1:8
    params=[p(pl) Vs(pl) Cf(pl) bet(pl) alph(pl)];
    subplot(4,2,pl)
    harvestEvolPlotMixedCostBenefit_solutions(initV,rates,cycleFreq,endTimes,Ds,DDVs,DV,CV,params);
end
harvestEvolPlotUs_solutions(initVs,rates,cycleFreq,endTimes,Ds,DDVs,DV,CV,params);
harvestEvolPlotMixedCostBenefit_rateChange_2Dstochastic(initV,rates,cycleFreq,endTimes,Ds,DDVs,DV,CV);