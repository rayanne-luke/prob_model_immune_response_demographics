% Code for creating models of two event antibody responce over time.
% Code by Kaitlyn Sullivan
% 11/29/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FINALTWOEVCODE(demoLabel,demoCat,pairTag,hist,single,two)
    % demoCat - deomographic category
    %   Smoker, Sex, Age, Race, Health, Diabetes, All
    % demoLabel - subpopulation
    %   Smoker: Smoker or Never
    %   Age: Below or Above (55)
    %   Sex: M or F
    %   Diabetes: T2, Non
    %   Race: White Non Hispanic, Black Non Hispanic, Hispanic
    %   Health: Healthy, Unhealthy
    %       Unhealthy: BMI > 30, Smoker/ Former, Age > 60 
    %       Healthy: BMI<30, Never smoker, Age < 45
    % pairTag - wheter or not data is paired
    %   'y' or 'n'
    % hist - 'y' or 'n' for if user wants a histogram 
    % single - 'y' or 'n' for if user wants a single event model
    % two - 'y' or 'n' for if user wants a two event model
    
    
    % Variable setting
    set(groot, 'defaultLineLineWidth', 0.5, 'defaulttextinterpreter', 'latex', ...
        'defaultAxesFontSize', 18, 'defaultLegendInterpreter', 'latex');
    rng(9);
    antdata = readtable('UVA_antibody_data.xlsx','Sheet','Antibody');
    demodata = readtable('UVA_antibody_data.xlsx','Sheet','Timedep_meas_ONLY');
    
    % Data table of desired demogrpahic antibody data, ids, etc. Use 'n' or 'y'
    % for desired pairing.
    demodata = balanceData(demodata,demoCat,demoLabel,pairTag);
    
    % Naive data points
    naiData = antdata(strcmp(antdata.Visit_ID,'Base'),:);
    
    % Gathering tables for infected then vaxinated and infected data points in 
    % infData and vaxData respectively
    vax_logical = logical(strcmp(demodata.Vaxed, 'y') & ~isnan(demodata.Days_Start));
    vaxData = demodata(vax_logical, :);
    inf_logical = logical(~strcmp(demodata.Visit_ID, 'Base') & strcmp(demodata.Vaxed, 'n'));
    infData = demodata(inf_logical, :);
    
    % Gathering tables for correspinding time
    timeI = infData.Days_Start; % for inf
    time = demodata.Days_Start;
    timeV = vaxData.Days_LastVax;
    timeV = str2double(timeV); % Necessary conversion
    timeSS = vaxData.Days_Start; % Time since start, 
    t_scale_fact = 100; % Scaling factor to assist with convergence of algorithm
    
    % Log scaling and gathering spike data points
    nai = log2(naiData.Spike_IgG + 2) - 1;
    spike = log2(demodata.Spike_IgG + 2) - 1;
    vaxData = log2(vaxData.Spike_IgG + 2) - 1;
    infData = log2(infData.Spike_IgG + 2) - 1;
    
    % Removing nans to ensure usable data
    nai(isnan(nai)) = []; 
    time(isnan(spike)) = [];
    spike(isnan(spike)) = [];
    timeI(isnan(infData)) = [];
    infData(isnan(infData)) = [];
    
    % Accounting for NA in days last vax
    timeV(isnan(vaxData)) = [];
    timeSS(isnan(vaxData)) = [];
    vaxData(isnan(vaxData)) = [];
    vaxData(isnan(timeV)) = [];
    timeSS(isnan(timeV)) = [];
    timeV(isnan(timeV)) = [];
    
    
    % Use if necessary to locate best parameters. Can help when there is a
    % lack of data points and model is not converging. 
    % 
    % Adding noise to help converge
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
    
    % Using maximum likelihood estimation to find best fit parameters for data
    % points
    F = @(x) (-sum(log(gampdf(nai,x(1),x(2))))); 
    x0 = [1,1];
    [naiP,naiNLL] = fminsearch(F,x0);
    
    % % New figure for histogram of the naive data
    if strcmp('y',hist)
        figure('Name','Naive Antibody Histogram');
    
        histogram(nai,15,'Normalization','pdf'); 
        hold on
        y = 0:0.01:max(nai);
        dist = gampdf(y,naiP(1),naiP(2));
        plot(y,dist,LineWidth=2);
        xticks(0:1:7);
        ylabel('Scaled Histogram');
        xlabel('Naive Spike Measurement');
    end
    
    %% Creating the two event model%%%%%%%%%
    
    % Shape parameter of gamma dist
    aFunc = @(theta1,theta2,t,k) ((theta1 .* (t)) ./ ( 1 + (theta2 .* ((t).^k))));
    
    % First event. Using MLE to find parameters. t_scale_fact assists with
    % convergence
    E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t/t_scale_fact,x(3)) + naiP(1),naiP(2));
    % E1 = @(r,t,x) gampdf(r,aFunc(x(1),x(2),t,x(3)) + naiP(1),naiP(2)); is E1
    % without time_scale_factor to assist with covergence
    
    G = @(x) (-sum(log(E1(infData,timeI,x))));
    g1 = 0.1;
    g2 = 0.0005;
    g3= 2;
    x0 = [g1*100,g2*(10^4),g3]; % Initial guesses with time_scale_factor
    % x0 = [g1,g2,g3]; is initial guess without time_scale_factor
    
    % infP - array containing parameters of the single event infected antibody
    % model
    % infNLL - negative log likelihood of found parameters fitting the dataset
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [infP,infNLL] = fmincon(G,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
    k1 = infP(3);
    
    % Create matrix of results of E1 at each point with calculated parameters
    s = linspace(0,13,300);
    t = linspace(0,800,300);
    if strcmp('y',single)
        z1 = zeros(300,300);
        for i = 1:300
            for j = 1:300
                z1(i,j) = E1(s(i),t(j),infP); 
            end
        end
    end
    
    %% Additional function for two event model
    % Combines E1 with the second event. Uses parameters from first event. 
    % Repeats process from above.
    
    E2 = @(r,t,tr,x) gampdf(r,aFunc(infP(1),infP(2),(t-tr)/t_scale_fact,k1) ...
        + aFunc(x(1),x(2),(tr)/t_scale_fact,x(3)) + naiP(1),naiP(2));
    %E2 = @(r,t,tr,x) gampdf(r,aFunc(infP(1),infP(2),(t-tr),k1) +
    %aFunc(x(1),x(2),tr,x(3)) + naiP(1),naiP(2)); is E2 without assistance from
    % time scaling
    
    H = @(x) (-sum(log(E2(vaxData,timeSS,timeV,x))));
    g1 =0.1;
    g2 = 0.0005;
    g3 = 2;
    x0 =  [g1*100,g2*(10^4),g3]; % with time scaling
    % x0 =  [g1,g2,g3]; % without time scaling
    
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [twoEvP,twoEvNLL] = fmincon(H,x0,[],[],[],[],[sqrt(eps),sqrt(eps),1 + sqrt(eps)],[],[],options);
    k2 = twoEvP(3);
    
    % Choosing the crossover point (the point in time in which the infected 
    % data points visually become infected then vaccinated data points in our 
    % model) as 120. 120 is chosen to visually best respresent the respective 
    % peaks and decays of infection and infection then vaccination
    
    %t_cross = mean(timeV);
    t_cross = 120;
    
    % Creating matrix of resulting two event model data points using calculated 
    % parameters
    if strcmp('y',two)
        z = zeros(300);
        for i = 1:300
            for j = 1:300
                z(i,j) = time_dep_two_plot(s(i),t(j)/t_scale_fact, t_cross/t_scale_fact,[infP(1),infP(2), twoEvP(1),twoEvP(2), naiP],k1,k2);
            end
        end
    end
    
    %% Plotting
    if strcmp('y',single)
        figure('Name','Contour of Infected Model');
        
        % Creating the contour plot of single event
        contourf(t,s,z1,30,'EdgeColor','None','FaceAlpha',0.6);
        hold on
        [C1,h1] = contour(t,s,z1,[0.001,0.01,0.1,0.25,0.5],'EdgeColor','black', ...
            'FaceAlpha',0.5, 'ShowText','on');
        h1.LineWidth = 1;
        
        % Scatter over contour
        scatter(timeI,infData,45,'MarkerFaceColor',"r",'MarkerEdgeColor','#cc0808');
        
        % Labeling
        legend('','','Infected');
        title(strcat(demoLabel,'Spike Single Event Model'));
        txt = {strcat('{Infected params: Peak  = }',num2str(infP(1)),'{ Decay = }',num2str(infP(2)),'{ k = }',num2str(k1))};
        subtitle(txt);
        xlabel("Days in Relative Time");
        ylabel("Log Scaled Antibody Measurments");
    end
    
    % Creating two event model contour
    if strcmp('y',two)
        figure('Name', 'Contour of Infected + Vaccinated Model')
        contourf(t,s,z,30,'EdgeColor','None','FaceAlpha',0.6);
        hold on
        [C,h] = contour(t,s,z,[0.001,0.01,0.1,0.25,0.5],'EdgeColor','black', ...
            'FaceAlpha',0.5, 'ShowText','on');
        h.LineWidth = 1.3;
        
        % Scatter personal trajectories over contour
        id = demodata.Subject_ID;
        anti = log2(demodata.Spike_IgG + 2) - 1;
        time = demodata.Days_Start;
        vax = demodata.Vaxed;
        
        id(isnan(anti)) = [];
        time(isnan(anti)) = [];
        anti(isnan(anti)) = [];
        
        hold on
        i = 1;
        while i <= size(id,1)
            j = i;
            a = true;
            while a
                % if strcmp(vax(i),'y') % If loop to color based on vax status
                % 
                %     plot(time(i),anti(i),'o','MarkerFaceColor','none','MarkerEdgeColor','none');
                % else
                %     plot(time(i),anti(i),'o','MarkerFaceColor','none','MarkerEdgeColor','none');
                % end
                if i + 1 > size(id,1) % If loop, lines through same subject IDs
                    a = false;
                elseif id(i) == id(i + 1)
                    i = i+1;
                else
                    a = false;
                end
            end
            plot(time(j:i),anti(j:i),'Color','blue','LineWidth',1);
            i = i + 1;
        end
        
        % JUST Scatter over contour NO TRAJECTORIES
        scatter(timeI,infData,45,'MarkerFaceColor','r','MarkerEdgeColor','#cc0808');
        scatter(timeSS,vaxData,45,'MarkerFaceColor','g','MarkerEdgeColor','#18c418');
        
        %Line at ave time since vax
        xline(t_cross,'--c',sprintf('%.4f',t_cross));
        
        % Labeling
        L1 = plot(nan,nan,'o','MarkerFaceColor','g','MarkerEdgeColor','g');
        L2 = plot(nan, nan,'o','MarkerFaceColor','r','MarkerEdgeColor','r');
        %L3 = plot(nan, nan,'Color','b');
        %legend([L3,L1, L2], {'SubjectID','Infected then Vaccinated', 'Infected'})
        legend([L1, L2], {'Infected then Vaccinated', 'Infected'})
        title(strcat(demoLabel,'Two Event Model'));
        txt = {strcat('{Infected params: Peak  = }',num2str(infP(1)),'{ Decay = }',num2str(infP(2)),'{ k = }',num2str(infP(3))),strcat('{Infected then Vaccination Params: Peak  = }',num2str(twoEvP(1)),'{ Decay = }',num2str(twoEvP(2)),'{ k = }',num2str(twoEvP(3)))};
        subtitle(txt);
        xlabel("Days in Relative Time");
        ylabel("Log Scaled Antibody Measurments");
    end
end

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