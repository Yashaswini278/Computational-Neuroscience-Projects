clear 
clc
close all 

warning off

% PART1 - MORRIS LECAR MODEL
% declaring variables
global gcal;
global gk; 
global gl;
global vca;
global vk;
global vl;
global phi;
global v1;
global v2;  
global v3; 
global v4; 
global v5; 
global v6;
global c;
global Iext; 

gcal=4.4; 
gk=8.0; 
gl=2; 
vca=120; 
vk=-84; 
vl=-60;
phi=0.02; 
v1=-1.2; 
v2=18;  
v3=2; 
v4=30; 
v5=2; 
v6=30;
c=20;
Iext=0;

% 2. plotting nullclines (V and w*100 axes)
V_nullcline = @(V,w_scaled) (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w_scaled*(vk-V)/100 + gl*(vl-V)+Iext);
w_nullcline = @(V,w_scaled) phi*((0.5*(1+tanh((V-v3)/v4)))-w_scaled/100)/(1/cosh((V-v3)/v4));
fimplicit(@(V,w_scaled) w_nullcline(V,w_scaled), [-100 100 0 100], 'Linewidth',4); 
hold on
fimplicit(@(V,w_scaled) V_nullcline(V,w_scaled), [-100 100 0 100], 'Linewidth',4);
xlabel('V');
ylabel('w*100');
title('Quiver plot with nullclines');
legend('w nullcline', 'V nullcline');
hold on 

% 2. quiver plot
x = linspace(-100,100,100);
y = linspace(0,1,100);
[V_quiv,w_quiv] = meshgrid(x, y);
U = (1/c)*(gcal*(0.5*(1+tanh((V_quiv-v1)/v2))).*(vca-V_quiv) + gk*w_quiv .*(vk-V_quiv) + gl*(vl-V_quiv) + Iext);
V = phi*((0.5*(1+tanh((V_quiv-v3)/v4)))-w_quiv).*cosh((V_quiv-v3)/(2*v4));
quiver(x,y*100,U,V*100,5,'k');


% 2. finding intersection point of the two nullclines i.e, equilibrium point
syms V;
syms w;
V_nullcline =(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*(0.5*(1+tanh((V-v3)/v4)))*(vk-V) + gl*(vl-V)+Iext);
Veq = vpasolve(V_nullcline,V);
weq = 0.5*(1+tanh((Veq-v3)/v4));
plot(Veq, weq*100, 'o');
text(Veq, weq*100, '(-60.853, 0.0149)');
legend('w nullcline', 'V nullcline', 'quiver', 'equilibrium point')
hold off

% 3. computing jacobian to find eigenvalue of eqm point and deduce stability
fprintf('3. computing jacobian to find eigenvalue of eqm point and deduce stability\n');
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {Veq, weq});
eigenv = eig(jaceq);
stab = check_stability(eigenv);
disp('The jacobian is:')
disp(jaceq)
disp('The eigen values are:');
disp(eigenv);
fprintf('Therefore the stability = %s\n\n',stab); 

%4. ODE solver options
options = odeset('RelTol',1e-3,'AbsTol',1e-6, 'refine',5, 'MaxStep', 1);

%5. generating action potential and plotting phase plane trajectories for
%different phi values 

Iext = 0;
weq = double(weq);
Veq = double(Veq);
tSpan = [0, 500];
initial1 = [0,weq];
initial2 = [-20,weq];

phi = 0.01;
[t1_1, S1_1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial1, options);
[t1_2, S1_2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial2, options);

phi = 0.02;
[t2_1, S2_1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial1, options);
[t2_2, S2_2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial2, options);

phi = 0.04;
[t3_1, S3_1] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial1, options);
[t3_2, S3_2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial2, options);

V_nullcline = @(V,w) (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nullcline = @(V,w) phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));

figure
fimplicit(@(V,w) w_nullcline(V,w), [-100 100 0 1]); 
hold on
fimplicit(@(V,w) V_nullcline(V,w), [-100 100 0 1]);
hold on 
plot(S1_1(:,1),S1_1(:,2), 'LineWidth', 5)
hold on
plot(S1_2(:,1),S1_2(:,2),'LineWidth', 5)
hold on
plot(S2_1(:,1), S2_1(:,2),'LineWidth', 5)
hold on 
plot(S2_2(:,1), S2_2(:,2),'LineWidth', 5)
hold on 
plot(S3_1(:,1), S3_1(:,2),'LineWidth', 5)
hold on
plot(S3_2(:,1), S3_2(:,2),'LineWidth', 5)
hold off
xlabel('V');
ylabel('w');
title('Phase plane trajectories for different phi values and initial conditions');
legend('w Nullcline', 'V nullcline', 'Action Potential Phi = 0.01', 'No Action Potenial Phi = 0.01', 'Action Potential Phi = 0.02', 'No Action Potential Phi = 0.02','Action Potential Phi = 0.04', 'No Action Potential Phi = 0.04');

% 6. plotting Voltage versus t depending upon initial value of V
% and maximum voltage reached vs initial value of V
phi = 0.02; 
V_initial = linspace(-60,0,100);
V_max = zeros(100, 1);
figure
for i = 1:100
    tSpan = [0, 300];
    initial = [V_initial(i),weq];
    [t2, S2] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan, initial, options);
    V_max(i) = max(S2(:,1));
    plot(t2, S2(:,1))
    hold on
end
xlabel('t');
ylabel('V');
title('Voltage vs t depending upon initial value of V');

figure
plot(V_initial, V_max, 'LineWidth', 3);
xlabel('V intial');
ylabel('V max');
title('Maximum amplitude of V reached versus intital Value of V');

% 7. Investigating effect of Iext on nullclines, equilibrium point and
% trajectories with different initial conditions
fprintf('7. Investigating effect of Iext on nullclines, equilibrium point and trajectories with different initial conditions'); 

phi = 0.02;
Iext = 86.0;
V_nullcline1 = @(V,w) (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nullcline1 = @(V,w) phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));

% find new equlibrium point
syms V;
syms w;
V_nullcline =(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*(0.5*(1+tanh((V-v3)/v4)))*(vk-V) + gl*(vl-V)+Iext);
Veq_new = vpasolve(V_nullcline,V);
weq_new = 0.5*(1+tanh((Veq_new-v3)/v4));
Veq_new = double(Veq_new);
weq_new = double(weq_new);

% check stability of new equilirbium point
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {Veq_new, weq_new});
eigenv2 = eig(jaceq);
stab2 = check_stability(eigenv2);

fprintf('\nFor Iext = %d, stability = %s\n', Iext, stab2);

tSpan_new = [0,1000];
initial_1 = [Veq, weq];
initial_2 = [Veq_new, weq_new];
initial_3 = [-27.9, 0.17];
[t4, S4] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan_new, initial_1, options);
[t5, S5] = ode15s(@(t,S)morris_lecar_ddt(t,S),tSpan_new, initial_2, options);
[t6, S6] = ode15s(@(t,S)morris_lecar_ddt(t,S), tSpan_new, initial_3, options);

figure
fimplicit(@(V,w) w_nullcline1(V,w), [-100 100 0 1]); 
hold on
fimplicit(@(V,w) V_nullcline1(V,w), [-100 100 0 1]);
hold on 
p3=plot(S4(:,1), S4(:,2));
set(p3,'color','k');
hold on
p4=plot(S5(:,1), S5(:,2));
set(p4,'color','b');
hold on
p5=plot(S6(:,1), S6(:,2));
set(p5,'color','g');
plot(Veq_new, weq_new, 'o');
title('Trajectories for different starting points for Iext = 86uA/cm^2');
legend('w Nullcline', 'V nullcline', 'start at [Veq, weq] for Iext = 0', 'start at [Veq, weq] for Iext = 86', 'start at [-27.9,0.17]', 'equilibrium point');


% 8. Unstable periodic orbit plotted backwards in time
[t7, S7] = ode15s(@(t,S)morris_lecar_ddt(t,S), [0,-500], initial_3, options);
figure
fimplicit(@(V,w) w_nullcline1(V,w), [-100 100 0 1]); 
hold on
fimplicit(@(V,w) V_nullcline1(V,w), [-100 100 0 1]);
hold on 
p3=plot(S4(:,1), S4(:,2));
set(p3,'color','k');
hold on
p4=plot(S5(:,1), S5(:,2));
set(p4,'color','b');
hold on
p5=plot(S6(:,1), S6(:,2));
set(p5,'color','g');
p6 = plot(S7(:,1), S7(:,2));
set(p6, 'color', 'r');
plot(Veq_new, weq_new, 'o');
title('Trajectories for different starting points for Iext = 86uA/cm^2');
legend('w Nullcline', 'V nullcline', 'start at [Veq, weq] for Iext = 0', 'start at [Veq, weq] for Iext = 86', 'start at [-27.9,0.17]', 'UPO', 'equilibrium point');


% 9. Checking the equilibrium points for 
%1) Iext = 80
%2) Iext = 86
%3) Iext = 90
fprintf('\n9. Checking the equilibrium pointd for various Iext\n');
Iext = 80.0;
% find new equlibrium point
syms V;
syms w;
V_nullcline =(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*(0.5*(1+tanh((V-v3)/v4)))*(vk-V) + gl*(vl-V)+Iext);
Veq_80 = vpasolve(V_nullcline,V);
weq_80 = 0.5*(1+tanh((Veq_80-v3)/v4));
Veq_80 = double(Veq_80);
weq_80 = double(weq_80);
% find stability of equilibrium point 
% check stability of new equilirbium point
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {Veq_80, weq_80});
eigenv2 = eig(jaceq);
stab80 = check_stability(eigenv2);
%Results: 
fprintf('Iext = %d, Equilibrium point [V,n] = [%.2f,%.2f], Stability = %s\n',Iext,Veq_80,weq_80,stab80);

Iext = 86.0;
% find new equlibrium point
syms V;
syms w;
V_nullcline =(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*(0.5*(1+tanh((V-v3)/v4)))*(vk-V) + gl*(vl-V)+Iext);
Veq_86= vpasolve(V_nullcline,V);
weq_86 = 0.5*(1+tanh((Veq_86-v3)/v4));
Veq_86 = double(Veq_86);
weq_86 = double(weq_86);
% find stability of equilibrium point 
% check stability of new equilirbium point
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {Veq_86, weq_86});
eigenv2 = eig(jaceq);
stab86 = check_stability(eigenv2);
%Results: 
fprintf('Iext = %d, Equilibrium point [V,n] = [%.2f,%.2f], Stability = %s\n',Iext,Veq_86,weq_86,stab86);

Iext = 90.0;
% find new equlibrium point
syms V;
syms w;
V_nullcline =(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*(0.5*(1+tanh((V-v3)/v4)))*(vk-V) + gl*(vl-V)+Iext);
Veq_90= vpasolve(V_nullcline,V);
weq_90 = 0.5*(1+tanh((Veq_90-v3)/v4));
Veq_90 = double(Veq_90);
weq_90 = double(weq_90);
% find stability of equilibrium point 
% check stability of new equilirbium point
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {Veq_90, weq_90});
eigenv2 = eig(jaceq);
stab90 = check_stability(eigenv2);
%Results: 
fprintf('Iext = %d, Equilibrium point [V,n] = [%.2f,%.2f], Stability = %s\n',Iext,Veq_90,weq_90,stab90);

frequency = zeros(21);
curr = 80:1:100;
for i = 1:21
    Iext = curr(i);
    syms V w
    Vnc1_eqn = (1/c)*(Iext - gcal*(0.5*(1+tanh((V-v1)/v2)))*(V-vca) - gk*w*(V-vk) - gl*(V-vl)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    eq_pt_1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w]);
    V_eq1 = double(eq_pt_1.V);
    w_eq1 = double(eq_pt_1.w);
    initpoint = [V_eq1+0.01, w_eq1+0.001];
    tspan = [0 2000];
    [t1,S1]=ode15s(@(t,S)morris_lecar_ddt(t,S), tspan, initpoint, options);
    frequency(i) = rate_of_ap(t1, S1);
end
figure
plot(curr, frequency);
xlabel('Iext'); 
ylabel('frequency'); 
title('frequency vs current plot'); 

% 10. MLE for parameter set-2
fprintf('\n10. MLE for parameter set-2\n')
gcal=4; 
gk=8.0; 
gl=2; 
vca=120; 
vk=-84; 
vl=-60;
phi=0.0667; 
v1=-1.2; 
v2=18;  
v3=12; 
v4=17.4; 
v5=12; 
v6=17.4;
c=20;
Iext=30;

%plotting the nullclines
V_nullcline1 = @(V,w) (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nullcline1 = @(V,w) phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));
fimplicit(@(V,w) w_nullcline1(V,w), [-100 100 -0.5 1], 'Linewidth',2); 
hold on
fimplicit(@(V,w) V_nullcline1(V,w), [-100 100 -0.5 1], 'Linewidth',2);
xlabel('V');
ylabel('w');
title('MLE parameter set-2 nullclines for Iext=30');
hold on

% Find the 3 equilibrium points 
syms V;
syms w;
V_nc = (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nc = phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));
eqns = [V_nc == 0, w_nc == 0]; 
eqm1 = vpasolve(eqns, [V, w], [-40, 0]);
eqm2 = vpasolve(eqns, [V, w], [-19, 0.02]);
eqm3 = vpasolve(eqns, [V, w], [4, 0.2]);
V_eqm1 = double(eqm1.V); 
V_eqm2 = double(eqm2.V); 
V_eqm3 = double(eqm3.V); 
w_eqm1 = double(eqm1.w); 
w_eqm2 = double(eqm2.w); 
w_eqm3 = double(eqm3.w);
plot(V_eqm1, w_eqm1, 'o');
plot(V_eqm2, w_eqm2, 'o');
plot(V_eqm3, w_eqm3, 'o');
hold on 

% stability analysis
syms V w
V_eqm = [V_eqm1, V_eqm2, V_eqm3];
w_eqm = [w_eqm1, w_eqm2, w_eqm3];
for i = 1:3
    V1 = V_eqm(:,i);
    w1 = w_eqm(:,i);
    jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
    jaceq = subs(jac, {sym('V'), sym('w')}, {V1, w1});
    eigenv = eig(jaceq);
    stab = check_stability(eigenv);
    fprintf('Equilibrium point %d, [Veq, weq] = [%.2f, %.2f], Stability = %s\n',i,V1,w1,stab);
end


[t1,S1] = ode15s(@(t,S)morris_lecar_ddt(t,S), [0, 1000], [-19.25, 0.024], options);
[t2,S2] = ode15s(@(t,S)morris_lecar_ddt(t,S), [0, 1000], [-19.6, 0.028], options);
[t3,S3] = ode15s(@(t,S)morris_lecar_ddt(t,S), [0, -100], [-20, 0.03], options);
[t4,S4] = ode15s(@(t,S)morris_lecar_ddt(t,S), [0, -60], [-19.56, 0.024], options);
plot(S1(:,1),S1(:,2), 'k');
hold on 
plot(S2(:,1), S2(:,2), 'k');
hold on
plot(S3(:,1),S3(:,2), 'm');
hold on 
plot(S4(:,1), S4(:,2), 'm');
legend('w nullcline', 'V nullcline', 'eqm point 1 (stable)', 'eqm point 2 (saddle)', 'eqm point 3 (unstable)', 'stable manifold','stable manifold', 'unstable manifold', 'unstable manifold'); 
fprintf('\n');

% 11. To investigate change in Iext 
fprintf('\n 11. To investigate change in Iext\n');
curr = [30, 38, 39, 39.1, 39.4, 39.9, 40, 45, 48, 50]; 
for i = 1:10
Iext = curr(i);
V_nullcline1 = @(V,w) (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nullcline1 = @(V,w) phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));
figure
fimplicit(@(V,w) w_nullcline1(V,w), [-100 100 -0.5 1]); 
hold on
fimplicit(@(V,w) V_nullcline1(V,w), [-100 100 -0.5 1]);
xlabel('V');
ylabel('w');
title(strcat('MLE parameter set-2 nullclines for Iext = ', num2str(Iext)));
end 

Iext = 45;
V_nc = (1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
w_nc = phi*((0.5*(1+tanh((V-v3)/v4)))-w)/(1/cosh((V-v3)/v4));
eqns = [V_nc == 0, w_nc == 0]; 
eqm = vpasolve(eqns, [V, w]);
V_eqm45 = double(eqm.V);
w_eqm45 = double(eqm.w); 
jac = jacobian([(1/c)*(gcal*(0.5*(1+tanh((V-v1)/v2)))*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext), phi*((0.5*(1+tanh((V-v3)/v4)))-w)*cosh((V-v3)/v4)],[V,w]);
jaceq = subs(jac, {sym('V'), sym('w')}, {V_eqm45, w_eqm45});
eigenv = eig(jaceq);
stab45 = check_stability(eigenv);

fprintf('For Iext = 30 - 39, there is 3 equilibrium points, stability = stable spiral, saddle, unstable\n');
fprintf('For Iext = 40 - 50, there is 1 equilibrium point, stability = %s\n',stab45);

frequency = zeros(16);
curr = 30:1:45;
for i = 1:10
    Iext = curr(i);
    [t1,S1]=ode15s(@(t,S)morris_lecar_ddt(t,S), tspan, [V_eqm2+0.1, w_eqm2+0.01], options);
    frequency(i) = rate_of_ap(t1, S1);
end
for i = 11:16
    Iext = curr(i); 
    [t1,S1]=ode15s(@(t,S)morris_lecar_ddt(t,S), tspan, [V_eqm2+0.1, w_eqm2+0.01], options);
    frequency(i) = rate_of_ap(t1, S1);
end
figure
plot(curr, frequency);
xlabel('Iext'); 
ylabel('frequency'); 
title('frequency vs current plot'); 


% Part-2 Hodgkin Huxley model


%declaring variables
global Gk;
global Ek;
global Gna;
global Ena;
global Gl;
global El;
global C;
global Q;
global T;
global phi;
global Iext;
global Vrest;
global h_inf_rest;
global f;
global n_inf_rest;

Iext = 0; 
phi = 1;
Gk = 36;
Gna = 120; 
Gl = 0.3;
Ek = -72;
Ena = 55; 
Q = 3; 
T = 6.3; 
C = 1; 

options = odeset('RelTol',1e-3,'AbsTol',1e-6, 'refine',5, 'MaxStep', 1);

% 13. To find Eleak so that Vrest = -60mV
fprintf('13. To find Eleak so that Vrest = -60mV\n');
Vrest = -60;
an = -0.01*(Vrest+50)/(exp(-(Vrest+50+1e-12)/10)-1);
bn = 0.125*phi*exp(-(Vrest+60)/80);
am = -0.1*phi*(Vrest+35)/(exp(-(Vrest+35+1e-12)/10)-1);
bm = 4*phi*exp(-(Vrest+60)/18);
ah = 0.07*phi* exp(-(Vrest+60)/20);
bh = phi/(exp(-(Vrest+30+1e-12)/10)+1);
n_inf = an/(an+bn);
m_inf = am/(am+bm); 
h_inf = ah/(ah+bh);
h_inf_rest = h_inf;
n_inf_rest = n_inf; 
Vleak = Vrest - (1/C)*(Iext-Gk*n_inf^4*(Vrest-Ek)-Gna*m_inf^3*h_inf*(Vrest-Ena))/(Gl);
El = Vleak;
fprintf('Eleak = %.2f\n\n', El); 
% To show action potentials with Iext = 10
Iext = 10;
tspan = [0,100];
[t1, S1] = ode15s(@(t,S)hh_eqn(t,S), tspan, [-60, n_inf, m_inf, h_inf], options);
figure
plot(t1, S1(:,1));
xlabel('t');
ylabel('V'); 
title('Simulating Hodgkin Huxley model for Iext = 10muA/cm^2');

% 14. To check stability of model at rest with Iext = 0 and to find
% threshold voltage
fprintf('14. To check stability of model at rest with Iext = 0 and to find threshold voltage\n');

Iext = 0; 
stab = check_stability_I(Iext);
fprintf('Iext = %d, stability = %s\n',Iext,stab); 

impulse_curr = linspace(0,10,100); 
V_max = zeros(100,1);
for i = 1:100
    [t2, S2] = ode15s(@(t,S)hh_eqn(t,S), tspan, [-60+impulse_curr(i)/C, n_inf, m_inf, h_inf], options);
    V_max(i) = max(S2(:,1));
end

figure
plot(impulse_curr, V_max); 
xlabel('Impulse current'); 
ylabel('Maximum voltage reached'); 
title('Maximum membrane voltage vs brief current pulses'); 


incr = V_max(2)-V_max(1);
for j = 3:100
    if(V_max(j)-V_max(j-1)>5*incr)
        break;
    end
end
V_thresh = -60+impulse_curr(j)/C ;

fprintf('Impulse current threshold = %.2f Therefore, Voltage threshold = %.2f\n\n',impulse_curr(j), V_thresh); 

% 15. To check stability of equilibrium point for various Iext
fprintf('15. To check stability of equilibrium point for various Iext\n');
for i=8:12
    Iext = i; 
    stab = check_stability_I(i);
    fprintf('Iext = %d, stability = %s\n',Iext, stab);
end

% 16. Myotonia HH model
Iext = 10;
f = 0;
[t3, S3] = ode15s(@(t,S)hh_myotonia(t,S), tspan, [-52, n_inf, m_inf, h_inf], options);
f = 0.1;
[t4, S4] = ode15s(@(t,S)hh_myotonia(t,S), tspan, [-52, n_inf, m_inf, h_inf], options);
f = 0.17;
[t5, S5] = ode15s(@(t,S)hh_myotonia(t,S), tspan, [-52, n_inf, m_inf, h_inf], options);
f = 0.2;
[t6, S6] = ode15s(@(t,S)hh_myotonia(t,S), tspan, [-52, n_inf, m_inf, h_inf], options);
figure
plot(t3, S3(:,1)); 
hold on 
plot(t4, S4(:,1)); 
hold on
plot(t5, S5(:,1));
hold on
plot(t6, S6(:,1), 'k');
xlabel('t');
ylabel('V');
title('Hodgkin Huxley model for Myotonia when Iext = 10');
legend('f=0', 'f=0.1', 'f=0.17', 'f=0.2');

% 17. V-n reduced model of HH 
fprintf('17. V-n reduced model of HH\n'); 
Iext = 10;
f = 0;
tspan = [0,100];
[t7, S7] = ode15s(@(t,S)hh_nV(t,S), tspan, [-60, n_inf], options);

Iext = 0; 
impulse_curr = linspace(0,10,100); 
V_max_new = zeros(100);
for i = 1:100
    [t2, S2] = ode15s(@(t,S)hh_nV(t,S), tspan, [-60+impulse_curr(i)/C, n_inf], options);
    V_max_new(i) = max(S2(:,1));
end

figure
plot(impulse_curr, V_max_new);
xlabel('Impulse current');
ylabel('Max membrane voltage');
title('V max vs Impulse current');

incr = V_max_new(2)-V_max_new(1);
for j = 3:100
    if(V_max_new(j)-V_max_new(j-1)>5*incr)
        break;
    end
end

V_thresh_new = -60+impulse_curr(j)/C ;
fprintf('Impulse current threshold = %.2f μA/cm^2, Therefore, V threshold = %.2f mV\n\n',impulse_curr(j),V_thresh_new); 


% 18. n-V phase plane analysis depending on parameter f
Iext = 0;     
an = @(V) -0.01 * (V + 50)/(exp(-(V + 1e-12 + 50)/10)-1);
am = @(V) -0.1 * (V  + 35)/(exp(-(V + 1e-12 + 35)/10)-1);
bn = @(V) 0.125 * exp(-(V + 60)/80);
bm = @(V) 4 * exp(-(V + 60+1e-12)/18);
m_inf_eqn = @(V) am(V) / (am(V) + bm(V));

fs = linspace(0.02, 0.4, 10);
for i = 1:10
f = fs(i);
Vnc = @(V) ( (Iext - Gna*(1-f)*(m_inf_eqn(V)^3)*h_inf_rest*(V-Ena) - Gna*f*(m_inf_eqn(V)^3)*(V-Ena) - Gl*(V-El))/(Gk * (V-Ek)) )^(1/4);
nnc = @(V) an(V)/(an(V) + bn(V));

fplot(@(V) Vnc(V), [-100 100]);
hold on 
fplot(@(V) nnc(V), [-100 100]);
xlabel('V(in mV)');
ylabel('n');
title('Phase Plane Plot (HH-V-n reduced)');
hold on;
end
    

% 19. Hodgkin Huxley model demonstrates Anode break excitation
Iext = -3;
tspan = [0,20]; 
[t1, S1] = ode15s(@(t,S)hh_eqn(t,S), tspan, [-60,n_inf,m_inf,h_inf], options);
Iext = 0; 
tspan1 = [20,80] ;
[t2,S2] = ode15s(@(t,S)hh_eqn(t,S),tspan1, [S1(end,1), S1(end,2), S1(end,3), S1(end,4)], options); 
tnet = [t1;t2];
Vnet = [S1(:,1);S2(:,1)];
nnet = [S1(:,2);S2(:,2)];
mnet = [S1(:,3);S2(:,3)];
hnet = [S1(:,4);S2(:,4)];
figure
plot(tnet, Vnet); 
hold on
xline(20, '--k', {'t_{end}'});
xlabel('t'); 
ylabel('V'); 
title('Action potential in Anode break excitation'); 
figure
plot(tnet, nnet); 
hold on 
plot(tnet, mnet); 
hold on 
plot(tnet, hnet);
hold on 
xline(20, '--k', {'t_{end}'}); 
xlabel('t');  
title('n,m,h parameters as a function of time'); 
legend('n','m','h'); 



% 20. Anode break excitation in V-m reduced HH (in phase plane)

fprintf('20. Anode break excitation in V-m reduced HH (in phase plane)\n'); 

an = @(V) -0.01 * (V + 50)/(exp(-(V + 1e-12 + 50)/10)-1);
am = @(V) -0.1 * (V + eps + 35)/(exp(-(V + 1e-12 + 35)/10)-1);
ah = @(V) 0.07 * exp(-(V + 60)/20);
bn = @(V) 0.125 * exp(-(V + 60)/80);
bm = @(V) 4 * exp(-(V + 60)/18);
bh = @(V) 1/(exp(-(V + 30)/10) + 1);

minf = am(Vrest)/(am(Vrest)+bm(Vrest)); 
ninf = an(Vrest)/(an(Vrest)+bn(Vrest)); 
hinf = ah(Vrest)/(ah(Vrest)+bh(Vrest)); 

% Starting point = [-Vrest, minf]
Iext = -3;
tSpan1 = [0, 20];
[t3, S3] = ode15s(@(t,S)hh_mV(t,S,ninf,hinf), tSpan1, [-60,minf], options);
Vnc1 = @(V) (((Iext - Gk * (V -Ek) * (ninf^4)- Gl*(V-El))/(Gna*hinf*(V-Ena)))^(1/3));
mnc = @(V) am(V)/(am(V) + bm(V));
figure;
hold on;
fplot(@(V) Vnc1(V), [-80 100]);
fplot(@(V) mnc(V), [-80 100]);

% Starting point = At the end of anodal stimulus 
Iext = 0;
tSpan2 = [20, 100];
ninf1 = S1(end,2); 
hinf1 = S1(end,4);
    
[t4, S4] = ode15s(@(t,S)hh_mV(t,S,ninf1,hinf1), tSpan2, [S3(end,1), S3(end,2)], options);
Vnc2 = @(V) ( (Iext - Gk * (V -Ek) * (ninf1^4)- Gl*(V-El))/(Gna*hinf1*(V-Ena)))^(1/3);

fplot(@(V) Vnc2(V), [-80 100]);
                                                              
ylabel('m');
xlabel('Voltage(in mV)');
title('Phase Plane for anode break');
legend('V nullcline starting [-60, minf(Vr)]', 'm nullcline', 'V nullcline starting after anodal stimulus'); 

Iext = -3;
syms V m;
    
V_nc = (1/C) * (Iext - Gk * ninf^4 * (V - Ek) - Gna * m^3 * hinf* (V - Ena) - Gl * (V - El));
m_nc = am *(1-m) - bm*m;
eqns = [V_nc == 0, m_nc == 0]; 
eqm1 = vpasolve(eqns, [V, m], [-70, 0]);
eqm2 = vpasolve(eqns, [V, m], [-54, 0.08]);
eqm3 = vpasolve(eqns, [V, m], [60, 1]);
V_eqm1_1 = double(eqm1.V);
V_eqm1_2 = double(eqm2.V); 
V_eqm1_3 = double(eqm3.V); 
m_eqm1_1 = double(eqm1.m); 
m_eqm1_2 = double(eqm2.m); 
m_eqm1_3 = double(eqm3.m);

V_eqm = [V_eqm1_1, V_eqm1_2, V_eqm1_3];
m_eqm = [m_eqm1_1, m_eqm1_2, m_eqm1_3];

fprintf('For part1\n');

for i = 1:3
    V1 = V_eqm(:,i);
    m1 = m_eqm(:,i);
    jac = jacobian([V_nc, m_nc],[V,m]);
    jaceq = subs(jac, {sym('V'), sym('m')}, {V1, m1});
    eigenv = eig(jaceq);
    stab = check_stability(eigenv);
    fprintf('Equilibrium point %d = %.2f, %.2f\tStability = %s\n',i,V1,m1,stab');
end
        
syms V m;
Iext=0;
V_nc = (1/C) * (Iext - Gk * ninf1^4 * (V - Ek) - Gna * m^3 * hinf1* (V - Ena) - Gl * (V - El));
m_nc = am *(1-m) - bm*m;
eq_pt = solve([V_nc == 0, m_nc == 0], [V, m]);
V_eq2 = double(eq_pt.V);
m_eq2 = double(eq_pt.m);
jac = jacobian([V_nc, m_nc],[V,m]);
jaceq = subs(jac, {sym('V'), sym('m')}, {V_eq2, m_eq2});
eigenv = eig(jaceq);
stab = check_stability(eigenv); 
fprintf('\nFor part-2\n');
fprintf('Equilibrium point 1 = %.2f, %.2f\tStability = %s\n',V_eq2,m_eq2,stab');


function dS = hh_eqn(t,S)
global Gk;
global Ek;
global Gna;
global Ena;
global Gl;
global El;
global C;
global phi;
global Iext;

V = S(1); 
n = S(2); 
m = S(3); 
h = S(4); 

an = -0.01*(V+50)/(exp(-(V+50+1e-12)/10)-1);
bn = 0.125*phi*exp(-(V+60)/80);
am = -0.1*phi*(V+35)/(exp(-(V+35+1e-12)/10)-1);
bm = 4*phi*exp(-(V+60)/18);
ah = 0.07*phi* exp(-(V+60)/20);
bh = phi/(exp(-(V+30+1e-12)/10)+1);

dV_dt = (1/C)*(Iext - Gk*n^4*(V-Ek) - Gna*m^3*h*(V-Ena) - Gl*(V-El));
dn_dt = an*(1-n) - bn*n; 
dm_dt = am*(1-m) - bm*m;
dh_dt = ah*(1-h) - bh*h;

dS = [dV_dt; dn_dt; dm_dt; dh_dt];
end


function stab = check_stability_I(I)
global Gk;
global Ek;
global Gna;
global Ena;
global Gl;
global El;
global C;
global phi;

syms V
syms n
syms m 
syms h

an = -0.01*(V+50)/(exp(-(V+50+1e-12)/10)-1);
bn = 0.125*phi*exp(-(V+60)/80);
am = -0.1*phi*(V+35)/(exp(-(V+35+1e-12)/10)-1);
bm = 4*phi*exp(-(V+60)/18);
ah = 0.07*phi* exp(-(V+60)/20);
bh = phi/(exp(-(V+30+1e-12)/10)+1);

V_nullcline = (1/C)*(I - Gk*n^4*(V-Ek) - Gna*m^3*h*(V-Ena) - Gl*(V-El));
n_nullcline = an*(1-n) - bn*(n); 
m_nullcline = am*(1-m) - bm*(m);
h_nullcline = ah*(1-h) - bh*(h);
eq = vpasolve([V_nullcline == 0, n_nullcline==0, m_nullcline==0, h_nullcline==0], [V,n,m,h]);
Veq0 = double(eq.V);
neq0 = double(eq.n);
meq0 = double(eq.m); 
heq0 = double(eq.h); 
jac = jacobian([V_nullcline, n_nullcline, m_nullcline, h_nullcline],[V,n,m,h]);
jaceq = subs(jac, {sym('V'), sym('n'), sym('m'), sym('h')}, {Veq0, neq0, meq0, heq0});
eigenv = eig(jaceq);

neg = 0;
pos = 0;

for i = 1:4
    if(real(eigenv(i))<0)
        neg = neg + 1;
    elseif (real(eigenv(i))>0)
        pos = pos + 1; 
    end
end

if (neg == 4)
    stab = 'Equilibrium point is stable'; 
elseif (pos == 4) 
    stab = 'Equilibrium point is unstable'; 
else
    stab = 'Cant say'; 
end
end

function dS = hh_myotonia(t,S)
    global C;
    global Iext;
    global Gk;
    global Gna;
    global Gl;
    global Ek;
    global Ena;
    global El;
    global f;
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    
    an =  -0.01 * (V + 50)/(exp(-(V + 50+1e-12)/10)-1);
    am =  -0.1 * (V + 35)/(exp(-(V + 35+1e-12)/10)-1);
    ah = 0.07 * exp(-(V + 60)/20);
    bn = 0.125 * exp(-(V + 60)/80);
    bm = 4 * exp(-(V + 60)/18);
    bh = 1/(exp(-(V + 30 + 1e-12)/10) + 1);
    
    dV_dt = (1/C) * (Iext - Gk * n^4 * (V - Ek) - Gna * (1-f) * m^3 * h * (V-Ena) - Gna * f * m^3 * (V-Ena)  - Gl * (V - El));
    dn_dt = an * (1 - n) - bn * n;
    dm_dt = am * (1 - m) - bm * m;
    dh_dt = ah * (1 - h) - bh * h;
    
    dS = [dV_dt; dn_dt; dm_dt; dh_dt];
end

function dS = hh_nV(t,S)
    global C;
    global Iext;
    global Gk;
    global Gna;
    global Gl;
    global Ek;
    global Ena;
    global El;
    global f;
    global h_inf_rest;
    
    V = S(1);
    n = S(2);
    
    an =  -0.01 * (V + 50)/(exp(-(V + 50+1e-12)/10)-1);
    am =  -0.1 * (V + 35)/(exp(-(V + 35+1e-12)/10)-1);
    bn = 0.125 * exp(-(V + 60)/80);
    bm = 4 * exp(-(V + 60)/18);
    
    m_inf = am/(am+bm);
    
    dV_dt = (1/C) * (Iext - Gk * n^4 * (V - Ek) - Gna * (1-f) * m_inf^3 * h_inf_rest * (V-Ena) - Gna * f * m_inf^3 * (V-Ena)  - Gl * (V - El));
    dn_dt = an * (1 - n) - bn * n;
    
    dS = [dV_dt; dn_dt];
end

function dS = hh_mV(t,S, ninf, hinf)
    global C;
    global Iext;
    global Gk;
    global Gna;
    global Gl;
    global Ek;
    global Ena;
    global El;
    
    V = S(1);
    m = S(2);
    
    am =  -0.1 * (V + 35)/(exp(-(V + 35+1e-12)/10)-1);
    bm = 4 * exp(-(V + 60)/18);
        
    dV_dt = (1/C) * (Iext - Gk * ninf^4 * (V - Ek) - Gna * m^3 * hinf * (V-Ena) - Gl * (V - El));
    dm_dt = am * (1 - m) - bm * m;
    
    dS = [dV_dt; dm_dt];
end

function dS = morris_lecar_ddt(t,S)

global gcal;
global gk; 
global gl;
global vca;
global vk;
global vl;
global phi;
global v1;
global v2;  
global v3; 
global v4; 
global v5; 
global v6;
global c;
global Iext; 

V=S(1);
w=S(2);

m_inf = (0.5*(1+tanh((V-v1)/v2)));
w_inf = (0.5*(1+tanh((V-v3)/v4)));

ddt_V = (1/c)*(gcal*m_inf*(vca-V) + gk*w*(vk-V) + gl*(vl-V)+Iext);
ddt_w = phi*(w_inf-w)*cosh((V-v3)/(2*v4));

dS=[ddt_V; ddt_w];

end

function stab = check_stability(eigenv)
if(imag(eigenv(1)) == 0)
    if(eigenv(1) < 0 && eigenv(2) < 0)
        stab = 'Equilibrium point is stable';
    elseif (eigenv(1)>0 && eigenv(2) > 0)
        stab = 'Equilibrium point is unstable';
    elseif (eigenv(1)*eigenv(2) < 0 )
        stab = 'Equilibrium point is saddle point';
    end
else
    if(real(eigenv(1)) > 0)
        stab = 'Equilibrium point is unstable and spiral';
    else
        stab = 'Equilibrium point is stable and spiral';
    end
end
end

function freq = rate_of_ap(t, S)
S = S(:,1); 
n = size(S);
has_spikes = 0; 
for i = 1:n
    if(S(i)>0) 
        has_spikes = 1; 
        break; 
    end
end
if(has_spikes == 0)
    freq = 0; 
    return 
end
if(i>n) 
    freq = 0; 
    return
end
while(S(i)-S(i-1)>0)
    i = i+1;
    if(i>n) 
        freq = 0; 
        return
    end
end
max1t = t(i); 
while(S(i)>0)
    i = i+1; 
    if(i>n) 
        freq = 0; 
        return
    end
end
while(S(i)<0)
    i = i+1; 
    if(i>n) 
        freq = 0; 
        return
    end
end
while(S(i)-S(i-1)>0)
    i = i+1; 
    if(i>n) 
        freq = 0; 
        return
    end
end
max2t = t(i); 
freq = 1000/(max2t-max1t); 
end

