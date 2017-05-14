clear all
clc;

init_popS=10000;
init_popR=0;
init_popI=5;

n = 1000;
t_end=10;
tspan = linspace(0,t_end,n+1);

Init_pop = [init_popS init_popI init_popR];

[t, y] = ode45(@SIR,tspan,Init_pop);

figure(1)

plot(t,y(:,1),'b','LineWidth',1.5); hold on;
plot(t,y(:,2),'--r','LineWidth',1.5); hold on;
plot(t,y(:,3),':g','LineWidth',1.5); 
xlabel('Time (in years)')
legend('susceptible','infected','recovered')
title('zombie model')

 function [dzom] = SIR(t,zom)
   

birthRate=0.6;
natdeathRate=0.01;
zomtransRate=0.05;
zominfecRate=0.08;
zomrecRate=0.1;

dzom = zeros(3,1); 

S = zom(1);
I = zom(2);
R  = zom(3);

dzom(1) = birthRate-(zomtransRate*S*I)-(zominfecRate*S)-(natdeathRate*S);
dzom(2) = (zomtransRate*S*I)-(zomrecRate*R)-(zominfecRate+natdeathRate)*I;
dzom(3) = (zominfecRate-natdeathRate)*R;



 end