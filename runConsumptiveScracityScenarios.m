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