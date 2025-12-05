% Graphics
set(groot, 'defaultLineLineWidth', 1, ...
           'defaultTextInterpreter', 'latex', ...
           'defaultAxesFontSize', 18, ...
           'defaultLegendInterpreter', 'latex');

% Preallocate storage for both models params
params = struct();

% Lambda and Time Scale
tScale    = 100;

load('OCSS_Vars.mat')

naive_spike = readtable("NaiveSpikeData.xlsx");
naive_spike = naive_spike(strcmp(naive_spike.Visit_ID, 'Base'), :);
naive_spike = naive_spike.Log_Scaled_Spike;
naive_spike_norm = (naive_spike - min(naive_spike))/(max(naive_spike) - min(naive_spike));

% Set Tables/Vars
Table1    = Normalized_Cytokines;
Table2    = NovaxData;
varMap.Table1.DependentVar   = "Acute_Cytokine_Signature";  
varMap.Table1.IndependentVar = "Days_Start";  
varMap.Table2.DependentVar   = "Normal_Spike_Score";  
varMap.Table2.IndependentVar = "Days_Start";  
lambdaMap.Table1 = 0.1; 
lambdaMap.Table2 = 0; 

% Create separate models for both tables
datasets = {"Table1","Table2"};
figure(1)
h1 = tiledlayout('flow');
for idx = 1:numel(datasets)
    nexttile, hold on;

    % Variable names
    DependentVar   = varMap.(datasets{idx}).DependentVar;   
    IndependentVar = varMap.(datasets{idx}).IndependentVar; 

    % Load the current table
    DataTable  = eval(datasets{idx});     

    % Split into naive and regular data tables
    if idx == 1
        TitleStr = "Cytokine";
        NaiveData   = DataTable(DataTable.Subject_ID >= 200, :);
        % MLE on naive for alpha/beta
        NaiveDependentVar       = NaiveData.(DependentVar);
        NaiveLogLikelihoodFcn   = @(param) -sum(log(gampdf(NaiveDependentVar, param(1), param(2))));
        OptimizedParameters     = fminsearch(NaiveLogLikelihoodFcn, [1,1]);
        a = OptimizedParameters(1);
        b = OptimizedParameters(2);
    else
        TitleStr = "Spike";
        NaiveData = naive_spike_norm;
        OptimizedParameters = fitdist(NaiveData, 'Gamma');
        a = OptimizedParameters.a;
        b = OptimizedParameters.b;
    end
    RegularData = DataTable(DataTable.Subject_ID < 200, :);


    % Store a and b
    if strcmp(datasets{idx}, 'Table1')
        params.neg_a_table1 = a;
        params.b_table1 = b;
    else
        params.neg_a_table2 = a;
        params.b_table2 = b;
    end

    % Shape model and penalty times
    ShapeFcn      = @(t,t1,t2,k) ((t1*(t/tScale))./(1 + t2*(t/tScale).^k)) + a;
    if idx == 1
        InitialGuess  = [1500, 500, 1.8];
    else
        InitialGuess = [100, 10, 2];
    end
    smallT        = 0:0.1:20;

    % Penalized negative log‐likelihood
    LogLikelihoodFcn = @(th) -sum(log( ...
                      gampdf( RegularData.(DependentVar), ...
                              ShapeFcn(RegularData.(IndependentVar), th(1), th(2), th(3)), ...
                              b ) ));
    regLambda = lambdaMap.(datasets{idx});
    RegularizedLLF   = @(th) LogLikelihoodFcn(th) + regLambda * max(ShapeFcn(smallT, th(1), th(2), th(3)));

    % Fit parameters [theta1, theta2, k]
    [OptimizedParameters2, fval, exit_flag] = fmincon(RegularizedLLF, ...
        InitialGuess, [], [], [], [], [sqrt(eps),sqrt(eps),1+sqrt(eps)], []);
    theta1               = OptimizedParameters2(1);
    theta2               = OptimizedParameters2(2);
    k                    = OptimizedParameters2(3);

    % Store for overlap calculation
    if strcmp(datasets{idx}, "Table1")
        params.fval_table1 = fval;
        params.exit_flag_table1 = exit_flag;
        params.theta1_table1 = theta1;
        params.theta2_table1 = theta2;
        params.k_table1      = k;
    else
        params.fval_table2 = fval;
        params.exit_flag_table2 = exit_flag;
        params.theta1_table2 = theta1;
        params.theta2_table2 = theta2;
        params.k_table2      = k;
    end

    % Contours
    [Rmat,Tmat] = meshgrid(0:0.1:9, 0:0.1:700);
    levels      = [0.01, 0.05, 0.1, 0.2, 0.3];
    pdfMat      = gampdf(Rmat, ShapeFcn(Tmat, theta1, theta2, k), b);
    % Add Text Labels
    [C,hf] = contourf(Tmat, Rmat, pdfMat, levels, 'LineColor','none');
    hold on;
    [C2,hc] = contour(Tmat, Rmat, pdfMat, levels, 'LineColor','k');
    clabel(C2, hc, 'LabelSpacing', 1500);
    hText = findobj(gca, 'Type', 'text');

    % Style them (bold, black, readable)
    set(hText, 'Color','k', 'FontWeight','bold', ...
        'FontSize',16, 'BackgroundColor','w', 'Margin',2);




    % Subject Trajectories
    groups  = findgroups(RegularData.Subject_ID);
    uniqGrp = unique(groups);
    for ii = 1:numel(uniqGrp)
        grpData = RegularData(groups==uniqGrp(ii), :);
        plot(grpData.(IndependentVar), grpData.(DependentVar), 'k');
    end

    % Plots
    if idx == 1
        h3 = scatter(RegularData.(IndependentVar), RegularData.(DependentVar), 'filled', 'mv');
    else
        h3 = scatter(RegularData.(IndependentVar), RegularData.(DependentVar), 'filled', 'v',  'MarkerFaceColor', "#ffa500");
    end

    if idx == 1
       h4 = scatter(NaiveData.(IndependentVar),   NaiveData.(DependentVar),   'filled','k');
    else
      h4 = scatter(zeros(length(NaiveData),1), NaiveData, 'filled', 'k');
    end
    axis([0 700 0 6]);
    legend([h3, h4], {"Post Infection","Control"},"Location","northeast", 'fontsize', 20);
    title(TitleStr);
    if idx == 1
        subtitle(sprintf('$\\theta_1=%3.0f$, $\\theta_2=%3.0f$, $k=%.2f$, $\\lambda=%.2f$', theta1, theta2, k, regLambda));
    else
        subtitle(sprintf('$\\theta_1=%3.1f$, $\\theta_2=%3.2f$, $k=%.2f$', theta1, theta2, k));
    end
end
xlabel(h1, "Days Since Infection", 'fontsize', 24, 'interpreter', 'latex');
ylabel(h1, "Signature", 'fontsize', 24, 'interpreter', 'latex');

% Overlap at times =[x,y]
times   = [10, 30, 100, 700];
x_vals  = 0:0.01:4;

figure(2)
h2 = tiledlayout('flow');
for i = 1:numel(times)
    t         = times(i);
    a_table1  = ((theta1*(t/tScale))./(1 + theta2*(t/tScale).^k)) + params.neg_a_table1;
    a_table2  = ((theta1*(t/tScale))./(1 + theta2*(t/tScale).^k)) + params.neg_a_table2;

    % Overlap Metric
    f_table1  = gampdf(x_vals, a_table1, params.b_table1);
    f_table2  = gampdf(x_vals, a_table2, params.b_table2);
    ov        = integral(@(r) min(gampdf(r,a_table1,params.b_table1), ...
        gampdf(r,a_table2,params.b_table2)), 0, Inf);

    % Plot
    nexttile, hold on;
    plot(x_vals, f_table1, 'm', 'LineWidth',2, 'DisplayName','Cytokine');
    plot(x_vals, f_table2, 'MarkerFaceColor', "#ffa500", 'LineWidth',2, 'DisplayName','Spike');
    if i == 4
        plot(x_vals, gampdf(x_vals, params.neg_a_table1, params.b_table1), ...
            'm--', 'LineWidth',2, 'DisplayName','Cytokine Na{\"i}ve');
        plot(x_vals, gampdf(x_vals, params.neg_a_table2, params.b_table2), '--',  ...
            'MarkerFaceColor', "#ffa500", 'LineWidth',2, 'DisplayName','Spike Na{\"i}ve');
    end
    yy = min(f_table1, f_table2);
    area(x_vals, yy, 'FaceColor',[0.7 0.7 0.7], 'DisplayName','Overlap region');
    title(sprintf('t = %d (Overlap = %.3f)', t, ov));
    if i == 4
        legend('Location','northeast', 'fontsize', 20);
    end
    ylim([0, 4])
end
xlabel(h2, 'Signature', 'fontsize', 24, 'interpreter', 'latex');
ylabel(h2, 'Probability Density', 'fontsize', 24, 'interpreter', 'latex');