% Variable setting
set(groot, 'defaultLineLineWidth', 0.5, 'defaulttextinterpreter', 'latex', ...
    'defaultAxesFontSize', 18, 'defaultLegendInterpreter', 'latex');
rng(9);
antdata = readtable('UVA_antibody_data.xlsx','Sheet','Antibody');
demodata = readtable('UVA_antibody_data.xlsx','Sheet','Timedep_meas_ONLY');

%% Use demodata1 as first group you'd like to compare to demodata2 as the second
%demodata2 = balanceData(demodata,'Diabetes','T2','n');
%demodata2 = balanceData(demodata,'Diabetes','Non','n');
%demodata2 = balanceData(demodata,'All','~','y');
%demodata2 = balanceData(demodata,'Sex','M','y');
%demodata2 = balanceData(demodata,'Sex','F','y');
demodata1 = balanceData(demodata,'Age','Below','y');
%demodata2 = balanceData(demodata,'Age','Above','y');
%demodata2 = balanceData(demodata,'Smoker','Never','n');
demodata2 = balanceData(demodata,'Smoker','Smoker','n');


% Naive data points
naiData = antdata(strcmp(antdata.Visit_ID,'Base'),:);

% Gathering tables for infected then vaxinated and infected data points in 
% infData and vaxData respectively
vax_logical = logical(strcmp(demodata1.Vaxed, 'y') & ~isnan(demodata1.Days_Start));
vaxData = demodata1(vax_logical, :);
inf_logical = logical(~strcmp(demodata1.Visit_ID, 'Base') & strcmp(demodata1.Vaxed, 'n'));
infData = demodata1(inf_logical, :);

% Gathering tables for correspinding time
timeI = infData.Days_Start; % for inf
time = demodata.Days_Start;
timeV = vaxData.Days_LastVax;
timeV = str2double(timeV); % Necessary conversion
timeSS = vaxData.Days_Start; % Time since start, 
t_scale_fact = 100; % Scaling factor to assist with convergence of algorithm


% Log scaling and gathering spike data points
nai = log2(naiData.Spike_IgG + 2) - 1;
spike = log2(demodata1.Spike_IgG + 2) - 1;
vaxData = log2(vaxData.Spike_IgG + 2) - 1;
infData = log2(infData.Spike_IgG + 2) - 1;

% Claeaning nans
nai(isnan(nai)) = []; 
time(isnan(spike)) = [];
spike(isnan(spike)) = [];
timeI(isnan(infData)) = [];
infData(isnan(infData)) = [];

timeV(isnan(vaxData)) = [];
timeSS(isnan(vaxData)) = [];
vaxData(isnan(vaxData)) = [];
vaxData(isnan(timeV)) = [];
timeSS(isnan(timeV)) = [];
timeV(isnan(timeV)) = [];

% Add noise to help converge if necessary
% n = 0;
% vaxNoise = vaxData + n*randn(1);
% infNoise = infData + n*randn(1);
% timeIN = timeI + n*randn(1);
% timeVN = timeV + n*randn(1);
% timeSSN = timeSS +n*randn(1);
% 
% infData = [infData;infNoise];
% vaxData = [vaxData;vaxNoise];
% timeI = [timeI;timeIN];
% timeV = [timeV;timeVN];
% timeSS=[timeSS;timeSSN];


%% MLE for naive data %%%%%%%%

F = @(x) (-sum(log(gampdf(nai,x(1),x(2)))));
x0 = [1,1];
[naiPg1,naiNLL] = fminsearch(F,x0);

%% Creating the two event model%%%%%%%%%

% Shape param 
aFunc = @(theta1,theta2,t,k) ((theta1 .* (t)) ./ ( 1 + (theta2 .* ((t).^k))));

% First event. Using MLE to find parameters.
E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t/t_scale_fact,x(3)) + naiPg1(1),naiPg1(2));
%E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t,x(3)) + naiP(1),naiP(2)); %
%without time scaling

G = @(x) (-sum(log(E1(infData,timeI,x))));
g1 = 0.2;
g2 = 0.0002;
g3=2;
x0 = [g1*100,g2*(10^4),g3];
%x0 = [g1,g2,g3]; % without time scaling
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[infPg1,infNLL] = fmincon(G,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
k1g1 = infPg1(3);

%% Additional function for two event model, combines E1 with the second
% event. Uses parameters from first event. 

E2 = @(r,t,tr,x) gampdf(r,aFunc(infPg1(1),infPg1(2),(t-tr)/t_scale_fact,k1g1) + aFunc(x(1),x(2),tr/t_scale_fact,x(3)) + naiPg1(1),naiPg1(2));
%E2 = @(r,t,tr,x) gampdf(r,aFunc(infP(1),infP(2),(t-tr),k1) + aFunc(x(1),x(2),tr,x(3)) + naiP(1),naiP(2));

H = @(x) (-sum(log(E2(vaxData,timeSS,timeV,x))));
g1 = 3;
g2 = 0.003;
g3 = 2;
x0 =  [g1*100,g2*(10^4),g3];
%x0 =  [g1,g2,g3];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[twoEvPg1,twoEvNLL] = fmincon(H,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
%[twoEvP,twoEvNLL] = fmincon(H,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[]);

k2g1 = twoEvPg1(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Group 2
%% MLE for naive data %%%%%%%%

% Gathering tables for infected then vaxinated and infected data points in 
% infData and vaxData respectively
vax_logical = logical(strcmp(demodata2.Vaxed, 'y') & ~isnan(demodata2.Days_Start));
vaxData = demodata2(vax_logical, :);
inf_logical = logical(~strcmp(demodata2.Visit_ID, 'Base') & strcmp(demodata2.Vaxed, 'n'));
infData = demodata2(inf_logical, :);

% Gathering tables for correspinding time
timeI = infData.Days_Start; % for inf
time = demodata.Days_Start;
timeV = vaxData.Days_LastVax;
timeV = str2double(timeV); % Necessary conversion
timeSS = vaxData.Days_Start; % Time since start, 
t_scale_fact = 100; % Scaling factor to assist with convergence of algorithm


% Log scaling and gathering spike data points
nai = log2(naiData.Spike_IgG + 2) - 1;
spike = log2(demodata2.Spike_IgG + 2) - 1;
vaxData = log2(vaxData.Spike_IgG + 2) - 1;
infData = log2(infData.Spike_IgG + 2) - 1;

% Claeaning nans
nai(isnan(nai)) = []; 
time(isnan(spike)) = [];
spike(isnan(spike)) = [];
timeI(isnan(infData)) = [];
infData(isnan(infData)) = [];

timeV(isnan(vaxData)) = [];
timeSS(isnan(vaxData)) = [];
vaxData(isnan(vaxData)) = [];
vaxData(isnan(timeV)) = [];
timeSS(isnan(timeV)) = [];
timeV(isnan(timeV)) = [];

% Add noise to help converge if necessary
% n = 0;
% vaxNoise = vaxData + n*randn(1);
% infNoise = infData + n*randn(1);
% timeIN = timeI + n*randn(1);
% timeVN = timeV + n*randn(1);
% timeSSN = timeSS +n*randn(1);
% 
% infData = [infData;infNoise];
% vaxData = [vaxData;vaxNoise];
% timeI = [timeI;timeIN];
% timeV = [timeV;timeVN];
% timeSS=[timeSS;timeSSN];


%% MLE for naive data %%%%%%%%

F = @(x) (-sum(log(gampdf(nai,x(1),x(2)))));
x0 = [1,1];
[naiPg1,naiNLL] = fminsearch(F,x0);

%% Creating the two event model%%%%%%%%%

% Shape param 
aFunc = @(theta1,theta2,t,k) ((theta1 .* (t)) ./ ( 1 + (theta2 .* ((t).^k))));

% First event. Using MLE to find parameters.
E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t/t_scale_fact,x(3)) + naiPg1(1),naiPg1(2));
%E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t,x(3)) + naiP(1),naiP(2)); %
%without time scaling

G = @(x) (-sum(log(E1(infData,timeI,x))));
g1 = 0.2;
g2 = 0.0002;
g3=2;
x0 = [g1*100,g2*(10^4),g3];
%x0 = [g1,g2,g3]; % without time scaling
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[infPg2,infNLL] = fmincon(G,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
k1g2 = infPg2(3);

%% Additional function for two event model, combines E1 with the second
% event. Uses parameters from first event. 

E2 = @(r,t,tr,x) gampdf(r,aFunc(infPg2(1),infPg2(2),(t-tr)/t_scale_fact,k1g2) + aFunc(x(1),x(2),tr/t_scale_fact,x(3)) + naiPg2(1),naiPg2(2));
%E2 = @(r,t,tr,x) gampdf(r,aFunc(infP(1),infP(2),(t-tr),k1) + aFunc(x(1),x(2),tr,x(3)) + naiP(1),naiP(2));

H = @(x) (-sum(log(E2(vaxData,timeSS,timeV,x))));
g1 = 3;
g2 = 0.003;
g3 = 2;
x0 =  [g1*100,g2*(10^4),g3];
%x0 =  [g1,g2,g3];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[twoEvPg2,twoEvNLL] = fmincon(H,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
%[twoEvP,twoEvNLL] = fmincon(H,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[]);

k2g2 = twoEvPg2(3);

%% Population functions
P1 = @(r,t) time_dep_two_plot(r,t/t_scale_fact, t_cross/t_scale_fact,[infPg1(1),infPg1(2), twoEvPg1(1),twoEvPg1(2), naiPg1],k1g1,k2g1);
P2 = @(r,t) time_dep_two_plot(r,t/t_scale_fact, t_cross/t_scale_fact,[infPg2(1),infPg2(2), twoEvPg2(1),twoEvPg2(2), naiPg2],k1g2,k2g2);

%% Pmin function
Pmin = @(r,t) min(P1(r,t),P2(r,t));

%% Double integral
OVL = (integral2(Pmin,0,10,0,t_cross)+integral2(Pmin,0,10,t_cross,600)+ integral2(Pmin,0,10,600,800))/800;


%% Functions

function shape_fun = shape_fun(t, params) 
%gives the shape parameter for time-dependent distribution (inf or vax) 
% at time period t.

    shape = params(3);
    theta = params([1 2]);
    shape_fun = (theta(1).*t)./(1+(theta(2).*t.^params(4))) + shape;


end

function W = time_dep_two_plot(r, t, t_cross, p,k1,k2)
% two-event distribution


if t < t_cross % cross-over time (from 1st to second event)
    shape_1 = shape_fun(t, [p(1), p(2), p(5),k1]);
    shape = shape_1;
    %RL: we want this model to be for both 1 and 2 events. So, if we're
    %before the "cross-over time", we want to plot the 1-event model.
end
if t >= t_cross
    t_since = t - t_cross;
    shape_1 = shape_fun(t, [p(1), p(2), p(5),k1]);
    shape_2 = shape_fun(t_since, [p(3), p(4), p(5),k2]);
    shape = shape_1 + shape_2 - p(5);
    %RL: plot the second event (the second increase in antibody response)
end
scale = p(6);
W = gampdf(r, shape, scale);

end