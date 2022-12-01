function []=harvestEvolPlotMixedCostBenefit_rateChange_2Dstochastic(initV,rates,cycleFreq,endTimes,Ds,DDVs,DV,CV)

% - find solutions to harvest rate evolution and plot corresponding stock
% equilibria for Pella-Tomlnson growth curve with all kinds of costs and
% benefits
% - assumes slow institutional adaptation and fast ecological dynamics
% - overlay with stochastic simulations of trajectories under different
% rates of change in reference benefit/cost ratio
% E.W. Tekwa, Apr 22, 2022

syms dS u(F) p E F Cfs V V2 Vs a r b xi R p w alph bet Vt(t) Seq Cf C3 L4 t
syms S positive

%write resource stock dynamics and equilibrium under fast dynamics
%dS=S*(r-a*S-F);
%dS=S*(r/p-((r^(1-p))*(a*S)^p)/p-F); %Pella-Tomlinson growth
dS=S*(r/p-((r^(1-p))*(a*S)^p)/p-F); %Pella-Tomlinson growth
p=1; %rho, (-1 to Inf, with 1 being logistic growth and <1 moving Smsy closer to R=0) *note values other than 0.5 and 1 take longer to compute
r=2; %intrinsic growth rate
a=1; %competition

Seq=eval(solve(eval(dS),S,'real',true))
Smsy=eval(subs(solve(eval(diff(dS,S)),'real',true),0));
Smax=eval(subs(eval(Seq),0));
%Smax=eval(solve(subs(eval(dS),F,0)));
%u=V*log(w*F*Seq)+V2*(F*Seq)-C*F*Seq-C2*F-C3*log(F*Seq)+Vs*Seq; %all types of cost and benefit

bet=1; %risk aversion controlling diminishing returns
if bet==1
    u=V*log(w*F*Seq)-Cfs*F*Seq-Cf*F+Vs*Seq; %all types of cost and benefit
else
    u=V*(((w*F*Seq)^(1-bet))-1)/(1-bet)-Cfs*F*Seq-Cf*F+Vs*Seq; %all types of cost and benefit
end

%Cobb-Douglas: let F be harvest rate, then:
%F*S=q*E^beta*S^alpha = F0*S^alpha = F0*S^(alpha-1)*S
%F=F0*S^(alpha-1) (if alpha=1, F=F0)
%F now means total factor productivity x scaled labour = "harvest strategy"

w=1; %shape parameter for diminishing returns (larger=faster diminishing returns)
Cfs=1; %process marginal cost (marginal cost with catch)
Cf=0; %effort-based marginal cost (marginal cost of each unit of effort x catchability): 0.01
C3=0; %multiplier for diminishing marginal cost with catch (economy of scales): (can be taken out)
%V=1; %reference price (multiplier for diminishing marginal benefit) - this
%is being varied in the plot
V2=0; %constant marginal benefit with catch: 0.1
Vs=0; %linearly diminishing cost (increasing benefit) with stock size ("stock ecosystem service")
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

set(0,'defaulttextinterpreter','tex');
set(0, 'defaultAxesTickLabelInterpreter','tex');
set(0, 'defaultLegendInterpreter','tex');
set(0,'defaultaxeslinewidth',2)
set(0,'DefaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');

%plot stock levels
figs=figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.2 scrsz(4)/3.5]);
set(figs,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
subplot(1,2,1)
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

text(1,1.8,{'economics:';['     u_S=' num2str(Vs,1) ', u_F=' num2str(Cf,1)   ', \beta=' num2str(bet,1)];
    'ecology:';['    \rho=' num2str(p,2)]
    },'FontSize',16);

xlabel 'price/(cost x MSY) [X=u_{\betaFS}/(u_{FS}MSY)]'
ylabel 'biomass (S/S_{MSY})'
ylim([-0.01 Smax/Smsy])
ytickrecord=yticks;

%simulate stochastic trajectories
nBiomass_int_slowChange=zeros(1,reps); %record biomass for slow change at the end of change
nBiomass_slowChange=zeros(1,reps); %record biomass for slow change after settling period
nBiomass_int_fastChange=zeros(1,reps); %record biomass for fast change at the end of change
nBiomass_fastChange=zeros(1,reps); %record biomass for fast change after settling period
for i=1:2 %low or high treatment
    dVdt=rates(i);
    %initV=initVs(i);
    cycV=cycVs(i);
    D=Ds(i);
    DDV=DDVs(i);
    Vpos=find(Vrange>=initV,1); %find position in Vrange closest to initV
    [maxS maxSpos]=max(Ssols_stable(:,Vpos));
    initF=Fsols(maxSpos,Vpos); %find initF at initV that leads to largest stable biomass
    %initF=initFs(i);
    t1=abs(DV/dVdt);
    if i==1
        times=[0:simTimeStep:t1+endTimes(1)];
    else
        times=[0:simTimeStep:t1+endTimes(2)];
    end
    for rep=1:reps
        [ts,Ft,Vt]=Euler_dFdt_slowInst(du,dVdt,cycV,initV,D,DDV,t1,times,initF,Fmax);
        t1pos=find(ts==t1);
        if i==1
            if rep<=plotReps
                plot(Vt,eval(subs(eval(Seq),Ft))/Smsy,'Color',[0 0 1 0.5]);
                startPt=scatter(initV,eval(subs(eval(Seq),initF))/Smsy,200,'^k','filled','LineWidth',2);
                %alpha(startPt,0.1);
                if DV~=0
                    intPt=scatter(Vt(t1pos),eval(subs(eval(Seq),Ft(t1pos)))/Smsy,100,'^b','LineWidth',2);
                    alpha(intPt,0.2);
                end
                endPt=scatter(Vt(end),eval(subs(eval(Seq),Ft(end)))/Smsy,400,'ob','LineWidth',2);
                alpha(endPt,0.2);
            end
            nBiomass_int_slowChange(rep)=eval(subs(eval(Seq),Ft(t1pos)))/Smsy;
            nBiomass_slowChange(rep)=eval(subs(eval(Seq),Ft(end)))/Smsy;
        else
            if rep<=plotReps
                plot(Vt,eval(subs(eval(Seq),Ft))/Smsy,'Color',[1 0 0 0.5]);
                startPt=scatter(initV,eval(subs(eval(Seq),initF))/Smsy,200,'^k','filled','LineWidth',2);
                %alpha(startPt,0.1);
                if DV~=0
                    intPt=scatter(Vt(t1pos),eval(subs(eval(Seq),Ft(t1pos)))/Smsy,100,'^r','LineWidth',2);
                    alpha(intPt,0.2);
                end
                endPt=scatter(Vt(end),eval(subs(eval(Seq),Ft(end)))/Smsy,400,'or','LineWidth',2);
                alpha(endPt,0.2);
            end
            nBiomass_int_fastChange(rep)=eval(subs(eval(Seq),Ft(t1pos)))/Smsy;
            nBiomass_fastChange(rep)=eval(subs(eval(Seq),Ft(end)))/Smsy;
        end
    end
end

% xlabel 'reference price/(cost x MSY) [V_{dFS}/(C_{FS}MSY)]'
% ylabel 'biomass (S/R_{MSY})'
% ylim([-0.01 Smax/Smsy])
% ytickrecord=yticks;
% title({['[C_{FS},C_F,C_{dFS},V_{FS},r,a,\rho]=[' num2str(C,1) ',' num2str(C2,1) ',' num2str(C3,1) ',' num2str(V2,1) ',' num2str(r,1) ',' num2str(a,1) ',' num2str(p,2) ']']; ...
%     ['[\sigma^2,dV_{dFS}/dt_{slow},dV_{dFS}/dt_{fast},dt]=['  num2str(D,1) ',' num2str(dVdts(1),1) ',' num2str(dVdts(end),1) ',' num2str(simTimeStep,1) ']']})
%title({['[C_{FS},C_F,C_{dFS},V_{FS},r,a,\rho]=[' num2str(C,1) ',' num2str(C2,1) ',' num2str(C3,1) ',' num2str(V2,1) ',' num2str(r,1) ',' num2str(a,1) ',' num2str(p,2) ']']; ...
%    ['[\sigma^2,dV_{dFS}/dt_{slow},dV_{dFS}/dt_{fast},dt]=['  num2str(D,1) ',' num2str(dVdts(1),1) ',' num2str(dVdts(end),1) ',' num2str(simTimeStep,1) ']']})
% text(1.3,1.8,{'ecology:';['    \rho=' num2str(p,2)];
%     'economics:';['    C_{FS}=' num2str(C,1) ' ,C_F=' num2str(C2,1)];['    V_{FS}=' num2str(V2,1) ' ,V_{S}=' num2str(Vs,1)];
%     },'FontSize',16);
if Ds(1)==Ds(2)
    text(1.3,1.32, {'process noise:';['    \sigma_F^2=' num2str(Ds(1))]
        },'FontSize',16);
else
    text(1.3,1.32, {'process noise:';['    \sigma_F^2_{ low}=' num2str(Ds(1))]},'Color','b','FontSize',16);
    text(1.3,1.09,['    \sigma_F^2_{ high}=' num2str(Ds(2))],'Color','r','FontSize',16);
end

trackY=0.65; %position of next text line
text(1.3,0.95,['X change:'],'Color','k','FontSize',16);
if DV==0
    text(1.3,0.8,['    dX/dt=0'],'Color','k','FontSize',16);
else
    text(1.3,0.8,['    dX/dt_{low}=' num2str(rates(1),1)],'Color','b','FontSize',16);
    text(1.3,0.65,['    dX/dt_{high}=' num2str(rates(2),1)],'Color','r','FontSize',16);
    trackY=trackY-0.15;
end
if sum(CV)==0
    text(1.3,trackY,['    cycle freq=0'],'Color','k','FontSize',16);
    trackY=trackY-0.15;
else
    text(1.3,trackY,['    cycle freq_{low}=' num2str(cycleFreq(1),2)],'Color','b','FontSize',16);
    text(1.3,trackY-0.15,['    cycle freq_{high}=' num2str(cycleFreq(2),2)],'Color','r','FontSize',16);
    trackY=trackY-0.3;
end
if DDVs(1)==DDVs(2)
    text(1.3,trackY,['    \sigma_X^2=0'],'Color','k','FontSize',16);
else
    text(1.3,trackY,['    \sigma_X^2_{ low}=' num2str(DDVs(1))],'Color','b','FontSize',16);
    text(1.3,trackY-0.15,['    \sigma_X^2_{ high}=' num2str(DDVs(2))],'Color','r','FontSize',16);
end

%title({[num2str(plotReps) ' trajectories per rate of change'];''},'fontweight','normal')
%text(minSpecies+5*(maxSpecies-minSpecies)/100,ylimits(2)*0.9,['Chao1 slope=' num2str(slope_Chao1,2) '\pm' num2str(slopeSD_Chao1,2) ', R^2*=' num2str(R2_Chao1,2) ', pairwise_{' num2str(maxRichDiff_all(1)) ',' num2str(maxRichDiff_all(2)) '}=' num2str(numCorrect_Chao1(1)/numTests,2) ',' num2str(numCorrect_Chao1(2)/numTests,2)],'Color','b')

percBelowSmsy_int_slow=100*sum(nBiomass_int_slowChange<1)/reps;
percBelowSmsy_int_fast=100*sum(nBiomass_int_fastChange<1)/reps;
percBelowSmsy_slow=100*sum(nBiomass_slowChange<1)/reps;
percBelowSmsy_fast=100*sum(nBiomass_fastChange<1)/reps;

subplot(1,2,2)
yyaxis left
hold on
histogram(nBiomass_slowChange,linspace(0,2,20),'FaceColor','b','EdgeAlpha',0,'FaceAlpha',0.2)
histogram(nBiomass_fastChange,linspace(0,2,20),'FaceColor','r','EdgeAlpha',0,'FaceAlpha',0.2)
xlabel 'biomass (S/S_{MSY})'
ylabel 'frequency'
xlim([0 2])
yyaxis right
hold on
[f_slow_int,xi_slow_int]=ksdensity(nBiomass_int_slowChange);
[f_fast_int,xi_fast_int]=ksdensity(nBiomass_int_fastChange);
[f_slow,xi_slow]=ksdensity(nBiomass_slowChange);
[f_fast,xi_fast]=ksdensity(nBiomass_fastChange);
plot(xi_slow,f_slow,'-b','LineWidth',2);
plot(xi_fast,f_fast,'-r','LineWidth',2);
curYlim=ylim;

trackY=1; %position of next text line
if DV~=0
    plot(xi_slow_int,f_slow_int,':b','LineWidth',2);
    plot(xi_fast_int,f_fast_int,':r','LineWidth',2);
    text(0.2,1*curYlim(2),'end of directional change (...):','FontSize',16);
    text(0.2,0.93*curYlim(2),['    ' num2str(percBelowSmsy_int_slow) '% below R_{MSY}'],'Color','b','FontSize',16);
    text(0.2,0.86*curYlim(2),['    ' num2str(percBelowSmsy_int_fast) '% below R_{MSY}'],'Color','r','FontSize',16);
    trackY=0.79;
end
text(0.2,trackY*curYlim(2),['after settling for ' num2str(endTimes(1)) ' time units (-):'],'FontSize',16);
text(0.2,(trackY-.07)*curYlim(2),['    ' num2str(percBelowSmsy_slow) '% below R_{MSY}'],'Color','b','FontSize',16);
text(0.2,(trackY-.14)*curYlim(2),['    ' num2str(percBelowSmsy_fast) '% below R_{MSY}'],'Color','r','FontSize',16);
ylabel 'probability density'
%title({[num2str(reps) ' replicates per rate of change'];''},'fontweight','normal')
