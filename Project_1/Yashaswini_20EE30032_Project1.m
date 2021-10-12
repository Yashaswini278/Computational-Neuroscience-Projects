
mu1 = 0.1; 
mu2 = 1; 
mu3 = 100;

f1 = @(t, y) VanDerPol(t, y, mu1);
f2 = @(t, y) VanDerPol(t, y, mu2);
f3 = @(t, y) VanDerPol(t, y, mu3);

t_span1 = [0, 50];
t_span2 = [0, 100];
t_span3 = [0, 500]; 

%initial condition
y0 = [1; 0]; 

[t, y] = ode15s(f1, t_span1, y0);
[t1, y1] = ode15s(f2, t_span2, y0);
[t2, y2] = ode15s(f3, t_span3, y0); 

%[t, y] = ode45(f1, t_span1, y0);
%[t1, y1] = ode45(f2, t_span2, y0);
%[t2, y2] = ode45(f3, t_span3, y0);

% time plots
figure
subplot(3,1,1);
plot(t, y(:,1));
title('mu = 0.1')
xlabel('t'); 
ylabel('y(t)');
subplot(3,1,2);
plot(t1, y1(:,1));
title('mu = 1')
xlabel('t'); 
ylabel('y(t)');
subplot(3,1,3);
plot(t2, y2(:,1));
title('mu = 100')
xlabel('t'); 
ylabel('y(t)');

% phase plane plots
figure
subplot(3,1,1);
plot(y(:,1),y(:,2));
title('mu = 0.1');
xlabel('y(t)'); 
ylabel("y'(t)");
subplot(3,1,2);
plot(y1(:,1), y1(:,2));
title('mu = 1');
xlabel('y(t)'); 
ylabel("y'(t)");
subplot(3,1,3);
plot(y2(:,1), y2(:,2));
title('mu = 100');
xlabel('y(t)'); 
ylabel("y'(t)");


function dy = VanDerPol(t, y, mu)
dy = [y(2); -y(1) + mu*(1-y(1)^2)*y(2)];
end


