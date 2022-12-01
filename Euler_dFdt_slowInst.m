function [ts,Ft,Vt]=Euler_dFdt_slowInst(dF,dVdt,cycV,initV,D,DDV,t1,times,initF,Fmax)
%function [ts,Ft,Vt]=Euler_dFdt_slowInst(dF,L,L2,a,r,v,dVdt,initV,D,t1,times,initF)

%Euler–Maruyama stochastic simulation
%E.W. Tekwa Apr 21, 2022


Ft(1)=initF;
ts(1)=times(1);
Vt(1)=initV;
dt=times(2)-times(1); %timestep
for t_index=2:length(times)
    t=times(t_index);
    if t<t1
        V=max(initV+dVdt*t+sqrt(DDV)*sqrt(dt)*randn+eval(cycV),0);
    else
        V=max(initV+dVdt*t1+sqrt(DDV)*sqrt(dt)*randn+eval(cycV),0);
    end
%     if t<t1
%         V=max(initV+dVdt*t+dt*sqrt(DDV)*randn,0);
%     else
%         V=max(initV+dVdt*t1+dt*sqrt(DDV)*randn,0);
%     end
    F=Ft(t_index-1);
    %dFdt=v+((F*L)/a + (L*(F - r))/a + (V*a*(F/a + (F - r)/a))/(F*(F - r)))+sqrt(D)*randn;
    %dFdt=eval(dF)+sqrt(D)*randn;
    dFdt=eval(dF);
    Ft(t_index)=Ft(t_index-1)+dt*dFdt+sqrt(D)*sqrt(dt)*randn;
    ts(t_index)=times(t_index);
    Vt(t_index)=V;
    if Ft(t_index)>Fmax
        Ft(t_index)=Fmax;
        return
    elseif Ft(t_index)<0
        Ft(t_index)=0;
    end
end