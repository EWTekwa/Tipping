# Tipping
Code for Tekwa &amp; Junquera 2022, A generalized adaptive harvesting model exhibits cusp bifurcation, noise, and rate-associated tipping pathways

Matlab code

1. Run "runConsumptiveScracityScenarios.m" to produce Figure 1 (consumptive scarcity as bifurcation parameter). Script calls the function "harvestEvolPlotMixedCostBenefit_solutions.m" to obtain numerical solutions for each scenario.
2. Run "runNonConsumptiveScarcityScenarios.m" to produce Figure 2 (non-consumptive scarcity as bifurcation parameter). Script calls the function "harvestEvolPlotUs_solutions.m" to obtain numerical solutions for each scenario.
3. Run "run2DCusp.m" to produce Figure 3 (cusp bifurcation). Load 2DCusp_solutions.mat to skip simulations and plot data.
4. Run "runTippingScenarios.m" to produce Figures 4 and 5 (different tipping phenomena). Script calls the function "harvestEvolPlotMixedCostBenefit_rateChange_2Dstochastic.m" to obtain stochastic solutions for each scenario.
