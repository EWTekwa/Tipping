%runRateChangeHarvestScenarios.m
%E.W. Tekwa Apr 26, 2022

%Define parameters for set of stochastic harvest runs, plotting biomass outcomes
%Each set contains reps for each of two treatments. Trajectories start at upper
%stable equilibrium biomass for initial bifurcation value (initu_betaFSs). They then experience
%harvest rate noise (D), bifurcation parameter noise (DDu_betaFS), bifurcation parameter
%directional change (amplitude: Du_betaFS, rates of change: rates), and bifurcation parameter
%cycle (amplitude Cu_betaFS,frequency cycleFreq). Once directional change is
%completed, cycle and noise continues for an additional time (endTimes).
set(0,'defaulttextinterpreter','tex');
set(0, 'defaultAxesTickLabelInterpreter','tex');
set(0, 'defaultLegendInterpreter','tex');
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');

%figs=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)]);


rng(1); %set random number generator seed

% C=1; %process marginal cost (marginal cost with catch)
% C2=0.02; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.02
% C3=0; %multiplier for diminishing marginal cost with catch: 0.5
%u_betaFS=1; %reference price (multiplier for diminishing marginal benefit) - this
%is being varied in the plot
%u_betaFS2=0; %constant marginal benefit with catch: 0.1
var_F=[0.005 0.05; 0.001 0.001;0.005 0.005; 0.005 0.005]; %gaussian process noise variance in F: 2
%u_betaFSsteps=200; %number of steps in u_betaFS to plot

%run trajectories with parameter change
R=[0 0; -0.001 -1; 0 0; 0 0]; %rates of change in u_betaFS
f=[0 0; 0 0; 0 0; 0.1 1]; %frequency of cycles in u_betaFS
initLambda_c=[0.9 1.5 1 1.5]; %1.5
DLambda_c=[0 -1 0 0]; %total directional u_betaFS change
A=[0 0; 0 0; 0 0; 0.5 0.5]; %cyclical amplitude in u_betaFS
%Cu_betaFS=[0.2 0.2]; %cyclical amplitude in u_betaFS: 0.25
%DDu_betaFSs=[0.001 1]; %temporal variance in u_betaFS
var_Lambda_c=[0 0; 0 0; 0.001 1; 0 0]; %temporal variance in u_betaFS
%du_betaFSdts=[rates(1)*ones(1,reps) rates(2)*ones(1,reps)]; %directional rate of change in u_betaFS
%cycu_betaFSs=[Cu_betaFS(1)*sin(2*pi*t/cycleFreq(1))*ones(1,reps)  Cu_betaFS(2)*sin(2*pi*t/cycleFreq(2))*ones(1,reps)]; %cyclical change in u_betaFS
%endTimes=abs(Du_betaFS./rates([2 1])); %constant environment time after change
endTimes=[10 10];

for scenario=1:size(R,1)
    Ds=var_F(scenario,:);
    rates=R(scenario,:);
    cycleFreq=f(scenario,:);
    initu_betaFS=initLambda_c(scenario);
    Du_betaFS=DLambda_c(scenario); %total directional u_betaFS change
    Cu_betaFS=A(scenario,:); %cyclical amplitude in u_betaFS
    DDu_betaFSs=var_Lambda_c(scenario,:); %temporal variance in u_betaFS

    harvestEvolPlotMixedCostBenefit_rateChange_2Dstochastic(initu_betaFS,rates,cycleFreq,endTimes,Ds,DDu_betaFSs,Du_betaFS,Cu_betaFS);
end