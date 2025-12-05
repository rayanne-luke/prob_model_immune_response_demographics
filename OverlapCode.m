% Graphics
set(groot, 'defaultLineLineWidth', 1, ...
           'defaultTextInterpreter', 'latex', ...
           'defaultAxesFontSize', 18, ...
           'defaultLegendInterpreter', 'latex');

% Preallocate storage for both models params
params = struct();

% Time-Scale
tScale    = 100;

% Create seperate models for both tables
datasets = {"Normal","Abnormal"};

figure; 
    hh1 = tiledlayout('flow');

for idx = 1:numel(datasets)
    % Set Lambda Values For Penalty
    if idx == 1
        regLambda = 0.03;
    else
        regLambda = 0.06;
    end 

    % Variable names
    DependentVar = "Acute_Cytokine_Signature";
    IndependentVar = "Days_Start";

    % Load the current table
    DataTable  = eval(datasets{idx});     % "Normal" first, then "Abnormal"
    TitleStr = sprintf("%s Subjects", datasets{idx});

    % Split into naive and regular data tables
    NaiveData = Normalized_Cytokines(Normalized_Cytokines.Subject_ID >= 200, :);
    RegularData  = DataTable(DataTable.Subject_ID < 200, :);

    % MLE on naive for alpha/beta
    NaiveDependentVar   = NaiveData.(DependentVar);
    NaiveLogLikelihoodFcn = @(param) -sum(log(gampdf(NaiveDependentVar, param(1), param(2))));
    OptimizedParameters     = fminsearch(NaiveLogLikelihoodFcn, [1,1]);
    a = OptimizedParameters(1);
    b = OptimizedParameters(2);

    % Shape model and penalty times
    ShapeFcn    = @(t,t1,t2,k) ((t1*(t/tScale))./(1 + t2*(t/tScale).^k)) + a;
    InitialGuess = [1*tScale, 0.01*(tScale^2), 1];
    smallT      = 0:0.1:20;

    % Penalized negative log‐likelihood
    LogLikelihoodFcn = @(th) -sum(log( ...
                  gampdf( RegularData.(DependentVar), ...
                          ShapeFcn(RegularData.(IndependentVar), th(1), th(2), th(3)), ...
                          b ) ));
    RegularizedLLF  = @(th) LogLikelihoodFcn(th) + regLambda * max(ShapeFcn(smallT, th(1), th(2), th(3)));

    % Fit parameters [theta1, theta2, k]
    OptimizedParameters2 = fmincon(RegularizedLLF, InitialGuess, [], [], [], [], [0,0,1], []);
    theta1   = OptimizedParameters2(1);
    theta2   = OptimizedParameters2(2);
    k = OptimizedParameters2(3);

    % Store for overlap calculation
    if strcmp(datasets{idx}, "Normal")
        params.theta1_norm = theta1;
        params.theta2_norm = theta2;
        params.k_norm  = k;
    else
        params.theta1_abn = theta1;
        params.theta2_abn = theta2;
        params.k_abn  = k;
    end

    nexttile, hold on
    % Contours
    [Rmat,Tmat] = meshgrid(0:0.1:9, 0:0.1:700);
    levels      = [0.01, 0.05, 0.2, 0.4];
    pdfMat      = gampdf(Rmat, ShapeFcn(Tmat, theta1, theta2, k), b);
    % Add Text Labels
    [C,hf] = contourf(Tmat, Rmat, pdfMat, levels, 'LineColor','none');
    hold on;
    [C2,hc] = contour(Tmat, Rmat, pdfMat, levels, 'LineColor','k');
    clabel(C2, hc, 'LabelSpacing', 1500);
    hText = findobj(gca, 'Type', 'text');

    % Subject Trajectories
    groups  = findgroups(RegularData.Subject_ID);
    uniqGrp = unique(groups);
    for ii = 1:numel(uniqGrp)
        grpData = RegularData(groups==uniqGrp(ii), :);
        plot(grpData.(IndependentVar), grpData.(DependentVar), 'k');
    end
    
    % Plots
    if idx == 1
        h3 = scatter(Normal.(IndependentVar),   Normal.(DependentVar),   'filled','b', 'Marker','diamond');
    else
        h3 = scatter(Abnormal.(IndependentVar), Abnormal.(DependentVar), 'filled','r', 'marker', '^');
    end
    h4 = scatter(NaiveData.(IndependentVar),  NaiveData.(DependentVar),  'filled','k');

    ax = gca; ax.FontSize = 22; 
    axis([0 700 0 6])
    if idx == 1
        legend([h3, h4], {"Normal","Na{\""i}ve","Location","best"}, 'fontsize', 24);
    else
        legend([h3, h4], {"Abnormal","Na{\""i}ve","Location","best"}, 'fontsize', 24);
    end
    title(TitleStr);
    subtitle(sprintf('$\\theta_1=%3.0f$, $\\theta_2=%3.0f$, $k$=%.2f, $\\lambda=%.2f$', theta1, theta2, k, regLambda));
end

xlabel(hh1, "Days Since Infection", 'fontsize', 24, 'interpreter', 'latex');
ylabel(hh1, "Cytokine Signature", 'fontsize', 24, 'interpreter', 'latex');

% Overlap at times =[x,y]
times = [10, 30, 100, 700];
x_vals  = 0:0.01:4;

% Shape Function
shapeAt = @(theta1,theta2,k,t) ((theta1*(t/tScale))./(1 + theta2*(t/tScale).^k)) + a;

figure; 
hh = tiledlayout('flow'); 

for i = 1:numel(times)
    t     = times(i);
    a_Norm = shapeAt(params.theta1_norm, params.theta2_norm, params.k_norm, t);
    a_Abn  = shapeAt(params.theta1_abn,  params.theta2_abn,  params.k_abn,  t);

    % Overlap Metric
    f_naive = gampdf(x_vals, a, b);
    f_Norm = gampdf(x_vals, a_Norm, b);
    f_Abn = gampdf(x_vals, a_Abn,  b);
    ov = integral(@(r) min(gampdf(r,a_Norm,b), gampdf(r,a_Abn,b)), 0, Inf);
    ov_naive_abn = integral(@(r) min(gampdf(r,a,b), gampdf(r,a_Abn,b)), 0, Inf); % abnormal-naive
    ov_naive_norm = integral(@(r) min(gampdf(r,a,b), gampdf(r,a_Norm,b)), 0, Inf); % normal-naive

    % Plot
    nexttile, hold on
    if t ~= 700
        plot(x_vals, f_Norm, 'b', 'LineWidth',1.5);
        plot(x_vals, f_Abn, 'r', 'LineWidth',1.5);
        yy = min(f_Norm, f_Abn);
        area(x_vals, yy, 'FaceColor',[0.7 0.7 0.7], 'DisplayName','Overlap region');
        title(sprintf('t = %d (Overlap = %.3f)', t, ov));
        ylim([0, 2.2])
        ax = gca; ax.FontSize = 24;
    elseif t == 700
         plot(x_vals, f_Norm, 'b', 'LineWidth',1.5, 'DisplayName','Normal PDF');
        plot(x_vals, f_Abn, 'r', 'LineWidth',1.5, 'DisplayName','Abnormal PDF');
        plot(x_vals, f_naive, 'k--', 'LineWidth',1.5, 'DisplayName','Naive PDF');
        yy = min(f_Norm, f_Abn);
        area(x_vals, yy, 'FaceColor',[0.7 0.7 0.7], 'DisplayName','Overlap region');
        title(sprintf('t = %d (Overlap = %.3f)', t, ov));
        legend('Location','best', 'fontsize', 28);
        ylim([0, 2.2])
        ax = gca; ax.FontSize = 24;
    end
end
  
 xlabel(hh, 'Cytokine Signature', 'fontsize', 30, 'interpreter', 'latex');
 ylabel(hh, 'Probability Density', 'fontsize', 30, 'interpreter', 'latex');
