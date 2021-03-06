clc;
clear;
close all 

warning off

load('data_cn_project_iii.mat'); 

%% Question 1
R_xx = zeros(101,1);
for tau = -50:50
    if (tau > 0)
        for i = 1:20000-tau
            R_xx(tau+51) = R_xx(tau+51) + (Stimulus(i)*Stimulus(i+tau));
        end
        R_xx(tau+51) = R_xx(tau+51)/(20000-tau);
    end
    if (tau<=0)
        for i = -1*tau+1:20000
             R_xx(tau+51) = R_xx(tau+51) + (Stimulus(i)*Stimulus(i+tau));
        end
        R_xx(tau+51) = R_xx(tau+51)/(20000+tau-1);
    end
end

t = linspace(-50, 50, 101); 

figure
plot(t, R_xx); 
xlabel('\tau');
ylabel('R(\tau)');
title('Autocorrelation function of Stimulus');

var_wn = max(R_xx); 


%% Question 2
psth = zeros(4, 20000);
nbins = 20*1000;

for i = 1:4
    for rep = 1:50
        f = figure('visible', 'off'); 
        psth(i,:) = psth(i,:) + histcounts(All_Spike_Times{i,rep}*1000, 0:20000)*1000/50;
    end
end

figure
subplot(2,2,1)
plot(psth(1,:))
xlabel('time (ms)');
ylabel('Rate (spikes/sec)');
title('PSTH for neuron-1');
subplot(2,2,2)
plot(psth(2,:))
xlabel('time (ms)');
ylabel('Rate (spikes/sec)');
title(['PSTH for neuron-2']);
subplot(2,2,3)
plot(psth(3,:))
xlabel('time (ms)');
ylabel('Rate (spikes/sec)');
title(['PSTH for neuron-3']);
subplot(2,2,4)
plot(psth(4,:))
xlabel('time (ms)');
ylabel('Rate (spikes/sec)');
title(['PSTH for neuron-4']);

smallpsth=zeros(4,100);    
for i = 1:100
    for neur_no=1:4
        smallpsth(neur_no,i)=mean(psth(neur_no,15000+i*50-49:15000+50*i));
    end
end

%% Question 3
time_bins = [10, 20, 50, 100, 200, 500];
Mean = zeros(50);
varince = zeros(50); 

for j = 1:6
    figure (2+j)
    bin = time_bins(j); 
    for i = 1:4
        for rep = 1:50
            num_spikes = histcounts(All_Spike_Times{i,rep}*1000, 0:bin:20000);
            Mean(rep+(i-1)*50) = mean(num_spikes);
            varince((i-1)*50+rep) = var(num_spikes); 
        end
        if(i == 1)
            scatter(Mean, varince, 'b');
        elseif(i == 2)
            scatter(Mean, varince, 'g');
        elseif(i == 3)
            scatter(Mean, varince, 'r');
        elseif(i == 4)
            scatter(Mean, varince, 'k');
        end
        hold on 
    end
    plot([min(varince), max(varince)],[min(varince), max(varince)], 'c');
    xlabel('Mean'); 
    ylabel('Variance');
    title(sprintf('Scatter plot of mean vs variance for time binsize = %d for all 4 neurons', bin));
end

%% Question 4

sta = zeros(4,100);

for neur = 1:4
    spikes = 0;
    for j = 1:50
        a = All_Spike_Times{neur,j}<=15;
        spikes = spikes + nnz(a);
        All_Spike_Times_train = All_Spike_Times{neur,j}(a);
        [m,n] = size(All_Spike_Times_train);
        for i = 1:n
            spike_time = max(1,ceil(1000*(All_Spike_Times_train(i))));
            stim = Stimulus(max(spike_time-99,1):spike_time);
            sta(neur,:) = sta(neur,:) + [zeros(1,100-length(stim)), stim]; 
        end
    end
    sta(neur,:) = sta(neur,:)/spikes;
end

sta(:,:) = sta/var_wn; 

figure (11)
subplot(2,2,1)
plot(0:99, sta(1,:));
ylim([-0.8,0.8])
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 1');
subplot(2,2,2)
plot(0:99, sta(2,:));
ylim([-0.8,0.8])
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 2');
subplot(2,2,3)
plot(0:99, sta(3,:));
ylim([-0.8,0.8])
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 3');
subplot(2,2,4)
plot(0:99, sta(4,:));
ylim([-0.8,0.8])
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 4');
sgtitle('Spike Triggered Average without correction'); 

R_xx_new = zeros(1,100); 

for tau = 0:99
    for i=1:20000-tau
        R_xx_new(tau+1)=R_xx_new(tau+1)+(Stimulus(i)*Stimulus(tau+i));
    end
    R_xx_new(tau+1)=R_xx_new(tau+1)/(20000-tau);
end
 
R_xx_new_mat = toeplitz(R_xx_new);


h = ((R_xx_new_mat\sta')');

figure (20)
subplot(2,2,1)
plot(0:99, h(1,:));
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 1');
subplot(2,2,2)
plot(0:99, h(2,:));
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 2');
subplot(2,2,3)
plot(0:99, h(3,:));
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 3');
subplot(2,2,4)
plot(0:99, h(4,:));
xlabel('\tau');
ylabel('R_{px}(\tau)');
title('Neuron 4');
sgtitle('Spike Triggered Average with correction')

%xlim([0,100]);

%% Question 5
pred(1, :) = conv(Stimulus(1:15000), h(1, :)); 
y(1,1:15000) = pred(1, 1:15000); 
pred(2, :) = conv(Stimulus(1:15000), h(2, :)); 
y(2,1:15000) = pred(2, 1:15000);
pred(3, :) = conv(Stimulus(1:15000), h(3, :)); 
y(3,1:15000) = pred(3, 1:15000);
pred(4, :) = conv(Stimulus(1:15000), h(4, :)); 
y(4,1:15000) = pred(4, 1:15000);
    
psth_train = psth(:, 1:15000);   
bin_size = 30;
for i = 1:ceil(15000/bin_size)
    e = i*bin_size;
    if e>15000
        e = 15000;
   end
        x1(i) = mean( y(1, (1+(i-1)*bin_size):e) );
        x2(i) = mean( y(2, (1+(i-1)*bin_size):e) );
        x3(i) = mean( y(3, (1+(i-1)*bin_size):e) );
        x4(i) = mean( y(4, (1+(i-1)*bin_size):e) );
    
        y1(i) = mean( psth_train(1, (1+(i-1)*bin_size):e) );
        y2(i) = mean( psth_train(2, (1+(i-1)*bin_size):e) );
        y3(i) = mean( psth_train(3, (1+(i-1)*bin_size):e) );
        y4(i) = mean( psth_train(4, (1+(i-1)*bin_size):e) );
 end

[param1]=sigm_fit(x1, y1);
[param2]=sigm_fit(x2, y2); 
[param3]=sigm_fit(x3, y3); 
[param4]=sigm_fit(x4, y4); 

fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4))); 

x_1 = linspace(-0.4, 0.5);
y_1 = fsigm(param1, x_1)

x_2 = linspace(-5, 5);
y_2 = fsigm(param2, x_2)

x_3 = linspace(-4, 4);
y_3 = fsigm(param3, x_3) 

x_4 = linspace(-0.3, 0.5);
y_4 = fsigm(param4, x_4) 

pred1 = conv(Stimulus(15001:20000), h(1, :)); 
pred1 = pred1(1:5000);
pred1 = fsigm(param1,pred1) + 0.0005*randn(size(pred1));

pred2 = conv(Stimulus(15001:20000), h(2, :)); 
pred2 = pred2(1:5000);
pred2 =fsigm(param2,pred2) + 0.0005*randn(size(pred2));

pred3 = conv(Stimulus(15001:20000), h(3, :)); 
pred3 = pred3(1:5000);
pred3 =fsigm(param3,pred3) + 0.0005*randn(size(pred3));

pred4 = conv(Stimulus(15001:20000), h(4, :)); 
pred4 = pred4(1:5000);
pred4 =fsigm(param4,pred4) + 0.0005*randn(size(pred4));

figure(12)
subplot(2,2,1)
scatter(x1, y1, '*')
hold on 
plot(x_1, y_1, 'linewidth', 2)
xlabel('y(t)')
ylabel('\lambda(t)');
title('Neuron 1');
subplot(2,2,2)
scatter(x2, y2, '*')
hold on 
plot(x_2, y_2, 'linewidth', 2)
xlabel('y(t)')
ylabel('\lambda(t)');
title('Neuron 2');
subplot(2,2,3)
scatter(x3, y3, '*')
hold on 
plot(x_3, y_3, 'linewidth', 2)
xlabel('y(t)')
ylabel('\lambda(t)');
title('Neuron 3');
subplot(2,2,4)
scatter(x4, y4, '*')
hold on 
plot(x_4, y_4, 'linewidth', 2)
xlabel('y(t)')
ylabel('\lambda(t)');
title('Neuron 4');
sgtitle('Estimating output non linearity');

%% Question 6
clc
 R1=corrcoef(psth(1, 15001:20000),pred1); 
 R_squared1=R1(2)^2 

 R2=corrcoef(psth(2, 15001:20000),pred2) ;
 R_squared2=R2(2)^2 
 
 R3=corrcoef(psth(3, 15001:20000),pred3) ;
 R_squared3=R3(2)^2 
 
 
 R4=corrcoef(psth(4, 15001:20000),pred4) ;
 R_squared4=R4(2)^2 
 
    A=zeros(100,1);
    B=zeros(100,1);
    initial_r_sq=R_squared2;
    old_r_sq=R_squared2;
    new_r_sq=old_r_sq;
    pos=zeros(1,1);
    count=0;
    h_new = h(:,:);
    while (old_r_sq - new_r_sq) < 0.001 && count < 100
        old_r_sq=new_r_sq;
        min_=1000000000;
        for i = 1:100 
            if abs(h_new(2,i))<min_ && abs(h_new(2,i))~=0
                min_=abs(h(2,i));
                pos=i;
            end
        end
        h_new(2,pos)=0;
        pos(1,1)=0;
        pred2 = conv(Stimulus(15001:20000), h_new(2, :)); pred2 = pred2(1:5000);
        pred2 =fsigm(param2,pred2) + 0.0005*randn(size(pred2));
        new_r=corrcoef(psth(2, 15001:20000),pred2); new_r_sq=new_r(2)^2;
        count=count+1;
        A(count,1)=count;
        B(count,1)=new_r_sq;
    end

 
 figure
 scatter(A,B, '*')
 xlabel('iterations')
 ylabel('r^2'); 
 title('Pruning parameters for neuron 2')
 
 max_correlation_coeff2 = max(B)
 
 
 f_t2 = fft(h_new(2,:));
 f_t1 = fft(h(1,:));

 vect = [-50:49];
 vect = vect';
 
 figure(98)
 subplot(2,2,1)
 stem(h(1,:));
 title('Neuron 1');
 subplot(2,2,2)
 stem(h_new(2,:));
 title('Neuron 2');
  subplot(2,2,3)
 stem(h(3,:));
 title('Neuron 3');
  subplot(2,2,4)
 stem(h(4,:));
 title('Neuron 4');
 sgtitle('Linear filters h(t)')
 
 figure(100)
 subplot(2,2,1)
 plot(vect, fft(h(1,:)));
 title('Neuron 1');
 subplot(2,2,2)
plot(vect, f_t2);
 title('Neuron 2');
  subplot(2,2,3)
plot(vect, fft(h(3,:)));
 title('Neuron 3');
  subplot(2,2,4)
plot(vect, fft(h(4,:)));
 title('Neuron 4');
 sgtitle('Linear filters h(t) in frequency domain')
 
 %% Question 7

 
%% Borowed function from internet to estimate sigmoid parameters  
function [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
% Optimization of parameters of the sigmoid function
%
% Syntax:
%       [param]=sigm_fit(x,y)       
%
%       that is the same that
%       [param]=sigm_fit(x,y,[],[],[])     % no fixed_params, automatic initial_params
%
%       [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
%       [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor
%       [param]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
% param = [min, max, x50, slope]
%
% if fixed_params=[NaN, NaN , NaN , NaN]        % or fixed_params=[]
% optimization of "min", "max", "x50" and "slope" (default)
%
% if fixed_params=[0, 1 , NaN , NaN]
% optimization of x50 and slope of a sigmoid of ranging from 0 to 1
%
%
% Additional information in the second output, STAT
% [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
%
% Example:
% %% generate data vectors (x and y)
% fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
% param=[0 1 5 1];  % "min", "max", "x50", "slope"
% x=0:0.1:10;
% y=fsigm(param,x) + 0.1*randn(size(x));
%
% %% standard parameter estimation
% [estimated_params]=sigm_fit(x,y)
%
% %% parameter estimation with forced 0.5 fixed min
% [estimated_params]=sigm_fit(x,y,[0.5 NaN NaN NaN])
%
% %% parameter estimation without plotting
% [estimated_params]=sigm_fit(x,y,[],[],0)
%
%
% Doubts, bugs: rpavao@gmail.com
% Downloaded from http://www.mathworks.com/matlabcentral/fileexchange/42641-sigmoid-logistic-curve-fit
% warning off
x=x(:);
y=y(:);
if nargin<=1 %fail
    fprintf('');
    help sigm_fit
    return
end
automatic_initial_params=[quantile(y,0.05) quantile(y,0.95) NaN 1];
if sum(y==quantile(y,0.5))==0
    temp=x(y==quantile(y(2:end),0.5));    
else
    temp=x(y==quantile(y,0.5));
end
automatic_initial_params(3)=temp(1);
if nargin==2 %simplest valid input
    fixed_params=NaN(1,4);
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==3
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==4
    plot_flag=1;    
end
if exist('fixed_params','var')
    if isempty(fixed_params)
        fixed_params=NaN(1,4);
    end
end
if exist('initial_params','var')
    if isempty(initial_params)
        initial_params=automatic_initial_params;
    end
end
if exist('plot_flag','var')
    if isempty(plot_flag)
        plot_flag=1;
    end
end
%p(1)=min; p(2)=max-min; p(3)=x50; p(4)=slope como em Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
%f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));
f_str='f = @(param,xval)';
free_param_count=0;
bool_vec=NaN(1,4);
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        f_str=[f_str ' param(' num2str(free_param_count) ')'];
        bool_vec(i)=1;
    else
        f_str=[f_str ' ' num2str(fixed_params(i))];
        bool_vec(i)=0;
    end
    if i==1; f_str=[f_str ' + (']; end
    if i==2;
        if isnan(fixed_params(1))            
            f_str=[f_str '-param(1) )./ (   1 + 10.^( (']; 
        else
            f_str=[f_str '-' num2str(fixed_params(1)) ')./ (1 + 10.^((']; 
        end
    end    
    if i==3; f_str=[f_str ' - xval ) *']; end
    if i==4; f_str=[f_str ' )   );']; end
end
eval(f_str)
[BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,initial_params(bool_vec==1));
stat.param=BETA';
% confidence interval of the parameters
stat.paramCI = nlparci(BETA,RESID,'Jacobian',J);
% confidence interval of the estimation
[stat.ypred,delta] = nlpredci(f,x,BETA,RESID,'Covar',COVB);
stat.ypredlowerCI = stat.ypred - delta;
stat.ypredupperCI = stat.ypred + delta;
% plot(x,y,'ko') % observed data
% hold on
% plot(x,ypred,'k','LineWidth',2)
% plot(x,[lower,upper],'r--','LineWidth',1.5)
free_param_count=0;
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        param(i)=BETA(free_param_count);
    else
        param(i)=fixed_params(i);
    end    
end
    
end





        
    




    
    
 
 