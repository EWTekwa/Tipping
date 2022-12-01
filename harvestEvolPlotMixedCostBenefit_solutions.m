function []=harvestEvolPlotMixedCostBenefit_solutions(initV,rates,cycleFreq,endTimes,Ds,DDVs,DV,CV,params)

% - find solutions to harvest rate evolution and plot corresponding stock
% equilibria for Pella-Tomlnson growth curve with all kinds of costs and
% benefits
% - assumes slow institutional adaptation and fast ecological dynamics
% - overlay with stochastic simulations of trajectories under different
% rates of change in reference benefit/cost ratio
% E.W. Tekwa, Apr 22, 2022

syms dS u(F) p E F Cfs V V2 Vs Vfs a r b xi R p w alph bet Vt(t) Seq Cf C3 L4 E alph t
syms S positive

%write resource stock dynamics and equilibrium under fast dynamics
%dS=S*(r-a*S-F);
dS=S*(r/p-((r^(1-p))*(a*S)^p)/p-F); %Pella-Tomlinson growth
%dS=S*(r/p-((r^(1-p))*(a*S)^p)/p-(E^alph)); %Pella-Tomlinson growth
%p=1; %rho, (-1 to Inf, with 1 being logistic growth and <1 moving Smsy closer to R=0) *note values other than 0.5 and 1 take longer to compute
p=params(1);
r=2; %intrinsic growth rate
a=1; %competition
alph=params(5);

Seq=eval(solve(eval(dS),S,'real',true))
Smsy=eval(subs(solve(eval(diff(dS,S)),'real',true),0));
Smax=eval(subs(eval(Seq),0));
%Smax=eval(solve(subs(eval(dS),F,0)));
%u=V*log(w*F*Seq)+V2*(F*Seq)-C*F*Seq-C2*F-C3*log(F*Seq)+Vs*Seq; %all types of cost and benefit

%bet=-1; %risk aversion controlling diminishing returns
bet=params(4);
% if bet==1
%     u=V*log(w*F*Seq)-Cfs*F*Seq+Vfs*F*Seq-Cf*F+Vs*Seq; %all types of cost and benefit
% else
%     u=V*(((w*F*Seq)^(1-bet))-1)/(1-bet)-Cfs*F*Seq+Vfs*F*Seq-Cf*F+Vs*Seq; %all types of cost and benefit
% end

if bet==1
    u=V*log(w*F*Seq)-Cfs*F*Seq-Cf*F^(1/alph)+Vs*Seq; %all types of cost and benefit
else
    u=V*(((w*F*Seq)^(1-bet))-1)/(1-bet)-Cfs*F*Seq-Cf*F^(1/alph)+Vs*Seq; %all types of cost and benefit
end

%Cobb-Douglas: let F be harvest rate, then:
%F*S=q*E^beta*S^alph = F0*S^alph = F0*S^(alph-1)*S
%F=F0*S^(alph-1) (if alph=1, F=F0)
%F now means total factor productivity x scaled labour = "harvest strategy"

w=1; %shape parameter for diminishing returns (larger=faster diminishing returns)
Cfs=1; %process marginal cost (marginal cost with catch)
Cf=params(3);
%Cf=0; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.01
C3=0; %multiplier for diminishing marginal cost with catch (economy of scales): (can be taken out)
%V=1; %reference price (multiplier for diminishing marginal benefit) - this
%is being varied in the plot
%V2=0; %constant marginal benefit with catch: 0.1
%Vfs=params(4);
%Vfs=0;
Vs=params(2);
%Vs=0; %linearly diminishing cost (increasing benefit) with stock size ("stock ecosystem service")
%D=1; %gaussian process noise variance in F: 2
Vsteps=200; %number of steps in V to plot

Fmax=eval(solve(Seq,F,'PrincipalValue',true)); %find maximum F after which extinction occurs
du=eval(simplify(diff(u,F))) %dF/dt is assumed proportional to du/dF in an adaptive process. Adaptive rate is set at 1.
ddu=simplify(diff(du,F));
F_symsols=solve(du,F)

%run trajectories with parameter change
reps=100; %replicates per rate of change in V
plotReps=3; %number of trajectories to plot
%rates=[0.1, 1]; %rates of change in V
%cycleFreq=[0.1 1]; %frequency of cycles in V
%initVs=0.9*ones(1,reps*2);
%DV=0; %total directional V change
%CV=[0.2 0.2]; %cyclical amplitude in V: 0.25
%DDV=0; %temporal variance in V
%dVdts=[rates(1)*ones(1,reps) rates(2)*ones(1,reps)]; %directional rate of change in V
cycVs=[CV(1)*(cos(2*pi*t/cycleFreq(1))-1), CV(2)*(cos(2*pi*t/cycleFreq(2))-1)]; %cyclical change in V
%endTimes=DV./rates([2 1]); %constant environment time after change
%endTimes=[50 50];
simTimeStep=0.01; %time step in Euler–Maruyama stochastic simulation

Vrange=linspace(0,2,Vsteps); %reference price

% set(0,'defaulttextinterpreter','tex');
% set(0, 'defaultAxesTickLabelInterpreter','tex');
% set(0, 'defaultLegendInterpreter','tex');
% set(0,'defaultaxeslinewidth',2)
% set(0,'DefaultAxesFontSize',16)
% scrsz = get(0,'ScreenSize');
% 
% %plot stock levels
% figs=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/3.5]);
% set(figs,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
%subplot(1,2,1)

hold on
Fsols=nan(16,length(Vrange));
Ssols=nan(16,length(Vrange));
Fstab=nan(16,length(Vrange));
for i=1:length(Vrange)
    V=Vrange(i);
    Fsol=eval(vpasolve(eval(du),F)); %compute equilibria
    Fsol((imag(Fsol)~=0))=NaN; %take out imaginary solutions
    Fsol(Fsol<=0)=NaN; %take out negative solutions
    Fsol(Fsol>=Fmax)=Fmax; %take out F>Fmax solutions (extinction)
    Fsols([1:length(Fsol)],i)=Fsol;
    Ssols(:,i)=eval(subs(eval(Seq),Fsols(:,i))); %compute resource biomass
    %Fstab(:,i)=eval(subs(eval(ddu),Fsols(:,i))); %compute stability (negative is stable)
    for j=1:sum(~isnan(Ssols(:,i)))
        Fstab(j,i)=eval(limit(eval(ddu),F,Fsols(j,i))); %compute stability (negative is stable)
    end
end
for curve=1:size(Ssols,1) %plot each curve (each row of equilibrium values)
    Ssols_stable=Ssols;
    Ssols_unstable=Ssols;
    Ssols_stable(curve,Fstab(curve,:)>=0)=NaN; %make unstable ranges NaN
    Ssols_unstable(curve,Fstab(curve,:)<0)=NaN; %make stable ranges NaN
    plot(Vrange,Ssols_stable(curve,:)'/Smsy,'-k','LineWidth',2); %plot stable parts of the curve
    plot(Vrange,Ssols_unstable(curve,:)'/Smsy,'--k','LineWidth',2); %plot unstable parts of the curve
end

text(1,1.8,{['\beta=' num2str(bet) ', \rho=' num2str(p) ', \alpha=' num2str(alph)];
    ['u_S=' num2str(Vs) ', u_E=' num2str(-Cf)];
    },'FontSize',16);

xlabel '\lambda_c'
ylabel 'biomass (S/S_{MSY})'
ylim([-0.01 Smax/Smsy])
ytickrecord=yticks;
