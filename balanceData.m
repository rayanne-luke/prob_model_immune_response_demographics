function G = balanceData(data,group,subgroup,pairedTag)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data - data table of antibody data, ids, etc
    % group - is population: Smoker, Sex, Age, Race, Health, Diabetes, All
    % subgroup is subpopulation: 
    %   Smoker: Smoker or Never
    %   Age: Below or Above (55)
    %   Sex: M or F
    %   Diabetes: T2, Non
    %   Race: White Non Hispanic, Black Non Hispanic, Hispanic
    %   Health: Healthy, Unhealthy
    %       Unhealthy: BMI > 30, Smoker/ Former, Age > 60 
    %       Healthy: BMI<30, Never smoker, Age < 45
    
    % pairedTag - indicated users desire to pair the data
    %   'y' or 'n'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng(9);

    % Clearing nans in spike antibody data
    data(isnan(data.Spike_IgG),:) = [];

    %% Balancing Sex
    if strcmp(group,'Sex')
        % G1 is WOMEN
        % G2 is MEN

        %% Infected
        infWL = logical(strcmp(data.Sex,'F') & strcmp(data.Vaxed,'n'));
        infW = data(infWL,:);
        infML = logical(strcmp(data.Sex,'M') & strcmp(data.Vaxed,'n'));
        infM = data(infML,:);

        L = logical(infW.Age < 55);
        infWbelow55 = infW(L,:);
        L = logical(infW.Age >= 55);
        infWabove55 = infW(L,:);
        L = logical(infM.Age < 55);
        infMbelow55 = infM(L,:);
        L = logical(infM.Age >= 55);
        infMabove55 = infM(L,:);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'F')
                [infWbelow55,infWabove55]=balance(infWbelow55,infWabove55);
            elseif strcmp(subgroup,'M')
                [infMbelow55,infMabove55]=balance(infMbelow55,infMabove55);
            end 
        end

        %% Doing the same for infected + vaccinated data

        %Removing nans from necessary vaccinated data
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        vaxWL = logical(strcmp(data.Sex,'F') & strcmp(data.Vaxed,'y'));
        vaxW = data(vaxWL,:);
        vaxML = logical(strcmp(data.Sex,'M') & strcmp(data.Vaxed,'y'));
        vaxM = data(vaxML,:);

        L = logical(vaxW.Age < 55);
        vaxWbelow55 = vaxW(L,:);
        L = logical(vaxW.Age >= 55);
        vaxWabove55 = vaxW(L,:);
        L = logical(vaxM.Age < 55);
        vaxMbelow55 = vaxM(L,:);
        L = logical(vaxM.Age >= 55);
        vaxDMabove55 = vaxM(L,:);

        % Balancing each group

        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'F')
                [vaxWbelow55,vaxWabove55]=balance(vaxWbelow55,vaxWabove55);
            elseif strcmp(subgroup,'M')
                [vaxMbelow55,vaxDMabove55]=balance(vaxMbelow55,vaxDMabove55);
            end 
        end
    
        if strcmp(subgroup,'F')
            G1 = combineSorted(infWbelow55, infWabove55);
            G2 = combineSorted(vaxWbelow55, vaxWabove55);
            G = combineSorted(G1, G2);
        else
            G1 = combineSorted(infMbelow55, infMabove55);
            G2 = combineSorted(vaxMbelow55, vaxDMabove55);
            G = combineSorted(G1, G2);
        end

    elseif strcmp(group,'Age')

        %% Infected
        infWL = logical(strcmp(data.Sex,'F') & strcmp(data.Vaxed,'n'));
        infW = data(infWL,:);
        infML = logical(strcmp(data.Sex,'M') & strcmp(data.Vaxed,'n'));
        infM = data(infML,:);

        L = logical(infW.Age <= 45);
        infWbelow55 = infW(L,:);
        L = logical(infW.Age >= 60);
        infWabove55 = infW(L,:);
        L = logical(infM.Age <= 45);
        infMbelow55 = infM(L,:);
        L = logical(infM.Age >= 60);
        infMabove55 = infM(L,:);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'Below')
                [infWbelow55,infMbelow55]=balance(infWbelow55,infMbelow55);
            elseif strcmp(subgroup,'Above')
                [infWabove55,infMabove55]=balance(infWabove55,infMabove55);
            end 
        end
        %% Doing the same for infected + vaccinated data

        %Removing nans from necessary vaccinated data
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        vaxWL = logical(strcmp(data.Sex,'F') & strcmp(data.Vaxed,'y'));
        vaxW = data(vaxWL,:);
        vaxML = logical(strcmp(data.Sex,'M') & strcmp(data.Vaxed,'y'));
        vaxM = data(vaxML,:);

        L = logical(vaxW.Age <= 45);
        vaxWbelow55 = vaxW(L,:);
        L = logical(vaxW.Age >= 60);
        vaxWabove55 = vaxW(L,:);
        L = logical(vaxM.Age <= 45);
        vaxMbelow55 = vaxM(L,:);
        L = logical(vaxM.Age >= 60);
        vaxDMabove55 = vaxM(L,:);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'Below')
                [vaxWbelow55,vaxMbelow55]=balance(vaxWbelow55,vaxMbelow55);
            elseif strcmp(subgroup,'Above')
                [vaxWabove55,vaxDMabove55]=balance(vaxWabove55,vaxDMabove55);
            end 
        end
    
        if strcmp(subgroup,'Below')
            G1 = combineSorted(infWbelow55, infMbelow55);
            G2 = combineSorted(vaxWbelow55, vaxMbelow55);
            G = combineSorted(G1, G2);
        elseif strcmp(subgroup,'Above')
            G1 = combineSorted(infWabove55, infMabove55);
            G2 = combineSorted(vaxWabove55, vaxDMabove55);
            G = combineSorted(G1, G2);
        end    
    elseif strcmp(group,'Smoker')
        data(strcmp(data.Smoker,'NA'),:)=[];
        %% Infected
        infSL = logical(~strcmp(data.Smoker,'Never') & strcmp(data.Vaxed,'n'));
        infS = data(infSL,:);
        infNL = logical(strcmp(data.Smoker,'Never') & strcmp(data.Vaxed,'n'));
        infN = data(infNL,:);

        L = logical(infS.Age < 55 & strcmp(infS.Sex,'F'));
        infSWbelow55 = infS(L,:);
        L = logical(infS.Age >= 55 & strcmp(infS.Sex,'F'));
        infSWabove55 = infS(L,:);
        L = logical(infS.Age < 55 & strcmp(infS.Sex,'M'));
        infSMbelow55 = infS(L,:);
        L = logical(infS.Age >= 55 & strcmp(infS.Sex,'M'));
        infSMabove55 = infS(L,:);
 
        L = logical(infN.Age < 55 & strcmp(infN.Sex,'F'));
        infNWbelow55 = infN(L,:);
        L = logical(infN.Age >= 55 & strcmp(infN.Sex,'F'));
        infNWabove55 = infN(L,:);
        L = logical(infN.Age < 55 & strcmp(infN.Sex,'M'));
        infNMbelow55 = infN(L,:);
        L = logical(infN.Age >= 55 & strcmp(infN.Sex,'M'));
        infNMabove55 = infN(L,:);

        listSA = {infSMbelow55, infSWbelow55, infSMabove55, infSWabove55};
        listNA = {infNMbelow55, infNWbelow55, infNMabove55, infNWabove55};
        listSM = {infSMbelow55, infSMabove55};
        listNM = {infNMbelow55, infNMabove55};
        smallestSA = findSmallest(listSA);
        smallestNA = findSmallest(listNA);
        smallestSM = findSmallest(listSM);
        smallestNM = findSmallest(listNM);

        % Balancing each group

        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'Never')
            % Non smoker
            [smallestNA,infNWbelow55]=balance(smallestNA,infNWbelow55);
            [smallestNA,infNMbelow55]=balance(smallestNA,infNMbelow55);
            [smallestNA,infNWabove55]=balance(smallestNA,infNWabove55);
            [smallestNA,infNMabove55]=balance(smallestNA,infNMabove55);
            elseif strcmp(subgroup,'Smoker')
            % Smoker
            [smallestSA,infSWbelow55]=balance(smallestSA,infSWbelow55);
            [smallestSA,infSMbelow55]=balance(smallestSA,infSMbelow55);
            [smallestSA,infSWabove55]=balance(smallestSA,infSWabove55);
            [smallestSA,infSMabove55]=balance(smallestSA,infSMabove55);
            elseif strcmp(subgroup,'Never Men')
            % Male Non Smokers
    
            [smallestNM,infNMbelow55]=balance(smallestNM,infNMbelow55);
    
            [smallestNM,infNMabove55]=balance(smallestNM,infNMabove55);
            elseif strcmp(subgroup,'Smoker Men')
            % Male Smoker
            [smallestSM,infSMbelow55]=balance(smallestSM,infSMbelow55);
            [smallestSM,infSMabove55]=balance(smallestSM,infSMabove55);
            end 
        end

        %% Vaccinated
        %Removing nans
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        vaxSL = logical(~strcmp(data.Smoker,'Never') & strcmp(data.Vaxed,'y'));
        vaxS = data(vaxSL,:);
        vaxNL = logical(strcmp(data.Smoker,'Never') & strcmp(data.Vaxed,'y'));
        vaxN = data(vaxNL,:);

        L = logical(vaxS.Age < 57 & strcmp(vaxS.Sex,'F'));
        vaxSWbelow55 = vaxS(L,:);
        L = logical(vaxS.Age >= 57 & strcmp(vaxS.Sex,'F'));
        vaxSWabove55 = vaxS(L,:);
        L = logical(vaxS.Age < 57 & strcmp(vaxS.Sex,'M'));
        vaxSMbelow55 = vaxS(L,:);
        L = logical(vaxS.Age >= 57 & strcmp(vaxS.Sex,'M'));
        vaxSMabove55 = vaxS(L,:);
 
        L = logical(vaxN.Age < 57 & strcmp(vaxN.Sex,'F'));
        vaxNWbelow55 = vaxN(L,:);
        L = logical(vaxN.Age >= 57 & strcmp(vaxN.Sex,'F'));
        vaxNWabove55 = vaxN(L,:);
        L = logical(vaxN.Age < 57& strcmp(vaxN.Sex,'M'));
        vaxNMbelow55 = vaxN(L,:);
        L = logical(vaxN.Age >= 57 & strcmp(vaxN.Sex,'M'));
        vaxNMabove55 = vaxN(L,:);

        listSA = {vaxSMbelow55, vaxSWbelow55, vaxSMabove55, vaxSWabove55};
        listNA = {vaxNMbelow55, vaxNWbelow55, vaxNMabove55, vaxNWabove55};
        listSM = {vaxSMbelow55, vaxSMabove55};
        listNM = {vaxNMbelow55, vaxNMabove55};
        smallestSA = findSmallest(listSA);
        smallestNA = findSmallest(listNA);
        smallestSM = findSmallest(listSM);
        smallestNM = findSmallest(listNM);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'Never')
            % Nonsmoker
            [smallestNA,vaxNWbelow55]=balance(smallestNA,vaxNWbelow55);
            [smallestNA,vaxNMbelow55]=balance(smallestNA,vaxNMbelow55);
            [smallestNA,vaxNWabove55]=balance(smallestNA,vaxNWabove55);
            [smallestNA,vaxNMabove55]=balance(smallestNA,vaxNMabove55);
            elseif strcmp(subgroup, 'Smoker')
            % Smoker
            [smallestSA,vaxSWbelow55]=balance(smallestSA,vaxSWbelow55);
            [smallestSA,vaxSMbelow55]=balance(smallestSA,vaxSMbelow55);
            [smallestSA,vaxSWabove55]=balance(smallestSA,vaxSWabove55);
            [smallestSA,vaxSMabove55]=balance(smallestSA,vaxSMabove55);
            elseif strcmp(subgroup, "Never Men")
            % Nonsmoker men
            [smallestNM,vaxNMbelow55]=balance(smallestNM,vaxNMbelow55);
            [smallestNM,vaxNMabove55]=balance(smallestNM,vaxNMabove55);
            elseif strcmp(subgroup,'Smoker Men')
            % Smoker men
            [smallestSM,vaxSMbelow55]=balance(smallestSM,vaxSMbelow55);
            [smallestSM,vaxSMabove55]=balance(smallestSM,vaxSMabove55);
            end 
        end
        
        if strcmp(subgroup, 'Smoker')
            G1 = combineSorted(infSMabove55,infSWabove55);
            G2 = combineSorted(infSMbelow55,infSWbelow55);
            G3 = combineSorted(vaxSWbelow55,vaxSMbelow55);
            G4 = combineSorted(vaxSWabove55,vaxSMabove55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup,'Never')
            G1 = combineSorted(infNWbelow55, infNMbelow55);
            G2 = combineSorted(infNMabove55,infNWabove55);
            G3 = combineSorted(vaxNWbelow55,vaxNMbelow55);
            G4 = combineSorted(vaxNWabove55,vaxNMabove55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup,'Smoker Men')
            G1 = combineSorted(infSMabove55,infSMbelow55);
            G3 = combineSorted(vaxSMbelow55,vaxSMabove55);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup,'Never Men')
            G1 = combineSorted(infNMbelow55,infNMabove55);
            G3 = combineSorted(vaxNMbelow55,vaxNMabove55);
            G = combineSorted(G1, G3);
        end

    elseif strcmp(group,'Diabetes')
        data(strcmp(data.Diabetes,''),:)=[];
        
        %% Infected
        infDL = logical(strcmp(data.Diabetes,'T2') & strcmp(data.Vaxed,'n'));
        infD = data(infDL,:);
        infNL = logical(~strcmp(data.Diabetes,'T2') & strcmp(data.Vaxed,'n'));
        infN = data(infNL,:);


        L = logical(infD.Age < 55 & strcmp(infD.Sex,'F'));
        infDWbelow55 = infD(L,:);
        L = logical(infD.Age >= 55 & strcmp(infD.Sex,'F'));
        infDWabove55 = infD(L,:);
        L = logical(infD.Age < 55 & strcmp(infD.Sex,'M'));
        infDMbelow55 = infD(L,:);
        L = logical(infD.Age >= 55 & strcmp(infD.Sex,'M'));
        infDMabove55 = infD(L,:);

        L = logical(infN.Age < 55 & strcmp(infN.Sex,'F'));
        infNWbelow55 = infN(L,:);
        L = logical(infN.Age >= 55 & strcmp(infN.Sex,'F'));
        infNWabove55 = infN(L,:);
        L = logical(infN.Age < 55 & strcmp(infN.Sex,'M'));
        infNMbelow55 = infN(L,:);
        L = logical(infN.Age >= 55 & strcmp(infN.Sex,'M'));
        infNMabove55 = infN(L,:);

        listD = {infDMbelow55, infDWbelow55,infDMabove55, infDWabove55};
        smallestD = findSmallest(listD);
        listN = {infNMbelow55, infNWbelow55,infNMabove55, infNWabove55};
        smallestN = findSmallest(listN);

        % Balancing each group

        if strcmp(pairedTag,'y')
            if strcmp(subgroup,'T2')
                % Healthy (< 25)
                [smallestD, infDWbelow55] = balance(smallestD, infDWbelow55);
                [smallestD, infDMbelow55] = balance(smallestD, infDMbelow55);
                [smallestD, infDWabove55] = balance(smallestD, infDWabove55);
                [smallestD, infDMabove55] = balance(smallestD, infDMabove55);
            else
                % Obese (>30)
                [smallestN, infNWbelow55] = balance(smallestN, infNWbelow55);
                [smallestN, infNMbelow55] = balance(smallestN, infNMbelow55);
                [smallestN, infNWabove55] = balance(smallestN, infNWabove55);
                [smallestN, infNMabove55] = balance(smallestN, infNMabove55);
            end
        end

        %% Vaccinated
        %Removing nans
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        % vaxHL = logical(data.BMI < 25 & strcmp(data.Vaxed,'y'));
        % vaxH = data(vaxHL,:);
        % vaxOvL = logical(data.BMI >= 25 & data.BMI < 30  & strcmp(data.Vaxed,'y'));
        % vaxOv = data(vaxOvL,:);
        % vaxObL = logical(data.BMI >= 30 & strcmp(data.Vaxed,'y'));
        % vaxOb = data(vaxObL,:);
        vaxDL = logical(strcmp(data.Diabetes,'T2') & strcmp(data.Vaxed,'y'));
        vaxD = data(vaxDL,:);
        vaxNL = logical(~strcmp(data.Diabetes,'T2')  & strcmp(data.Vaxed,'y'));
        vaxN = data(vaxNL,:);

        L = logical(vaxD.Age < 55 & strcmp(vaxD.Sex,'F'));
        vaxDWbelow55 = vaxD(L,:);
        L = logical(vaxD.Age >= 55 & strcmp(vaxD.Sex,'F'));
        vaxDWabove55 = vaxD(L,:);
        L = logical(vaxD.Age < 55 & strcmp(vaxD.Sex,'M'));
        vaxDMbelow55 = vaxD(L,:);
        L = logical(vaxD.Age >= 55 & strcmp(vaxD.Sex,'M'));
        vaxDMabove55 = vaxD(L,:);
 
        L = logical(vaxN.Age < 55 & strcmp(vaxN.Sex,'F'));
        vaxNWbelow55 = vaxN(L,:);
        L = logical(vaxN.Age >= 55 & strcmp(vaxN.Sex,'F'));
        vaxNWabove55 = vaxN(L,:);
        L = logical(vaxN.Age < 55 & strcmp(vaxN.Sex,'M'));
        vaxNMbelow55 = vaxN(L,:);
        L = logical(vaxN.Age >= 55 & strcmp(vaxN.Sex,'M'));
        vaxNMabove55 = vaxN(L,:);

        listD = {vaxDMbelow55, vaxDWbelow55,vaxDMabove55, vaxDWabove55};
        smallestD = findSmallest(listD);
        listN = {vaxNMbelow55, vaxNWbelow55,vaxNMabove55, vaxNWabove55};
        smallestN = findSmallest(listN);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'T2')
                % Healthy (< 25)
                [smallestD, vaxDWbelow55] = balance(smallestD, vaxDWbelow55);
                [smallestD, vaxDMbelow55] = balance(smallestD, vaxDMbelow55);
                [smallestD, vaxDWabove55] = balance(smallestD, vaxDWabove55);
                [smallestD, vaxDMabove55] = balance(smallestD, vaxDMabove55);
            else
                % Overweight (25 - 30)
                [smallestN, vaxNWbelow55] = balance(smallestN, vaxNWbelow55);
                [smallestN, vaxNMbelow55] = balance(smallestN, vaxNMbelow55);
                [smallestN, vaxNWabove55] = balance(smallestN, vaxNWabove55);
                [smallestN, vaxNMabove55] = balance(smallestN, vaxNMabove55);
            end
        end

        if strcmp(subgroup, 'T2')
            G1 = combineSorted(infDWbelow55, infDMbelow55);
            G2 = combineSorted(infDWabove55, infDMabove55);
            G3 = combineSorted(vaxDWbelow55,vaxDMbelow55);
            G4 = combineSorted(vaxDWabove55,vaxDMabove55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        else
            G1 = combineSorted(infNWabove55,infNMabove55);
            G2 = combineSorted(infNMbelow55,infNWbelow55);
            G3 = combineSorted(vaxNWabove55,vaxNMabove55);
            G4 = combineSorted(vaxNWbelow55,vaxNMbelow55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        end

    elseif strcmp(group,'Health')
    %% Infected
    infHML = logical(strcmp(data.Sex,'M') & strcmp(data.Smoker,'Never') & data.Age<=45 & strcmp(data.Vaxed,'n'));
    infHM = data(infHML,:);
    infHWL = logical(strcmp(data.Sex,'F') & strcmp(data.Smoker,'Never') & data.Age<=45 & strcmp(data.Vaxed,'n'));
    infHW = data(infHWL,:);    
    infUnML = logical(strcmp(data.Sex,'M') & ~strcmp(data.Smoker,'Never') & data.Age>=60 & strcmp(data.Vaxed,'n'));
    infUnM = data(infUnML,:);
    infUnWL = logical(strcmp(data.Sex,'F') & ~strcmp(data.Smoker,'Never') & data.Age>=60 & strcmp(data.Vaxed,'n'));
    infUnW = data(infUnWL,:);

    % Balancing each group
    if strcmp(pairedTag,'y')
        if strcmp(subgroup,'Healthy')
            % Healthy (< 25)
            [infHM,infHW] = balance(infHM,infHW);
        elseif strcmp(subgroup,'Unhealthy')
            [infUnM,infUnW] = balance(infUnM,infUnW);
        end
    end

    %% Vaccinated
    %Removing nans
    data(isnan(data.Days_Start),:) = [];
    data(isnan(str2double(data.Days_LastVax)),:) = [];

    % Removing people who were vaccinated before getting infected
    timeV = str2double(data.Days_LastVax);
    timeSS = data.Days_Start; %time since start, 
    vax_inf_logical = logical(timeSS - timeV < 0);
    data(vax_inf_logical,:) = [];

    vaxHML = logical(strcmp(data.Sex,'M') & strcmp(data.Smoker,'Never') & data.Age<=45 & strcmp(data.Vaxed,'y'));
    vaxHM = data(vaxHML,:);
    vaxHWL = logical(strcmp(data.Sex,'F') & strcmp(data.Smoker,'Never') & data.Age<=45 & strcmp(data.Vaxed,'y'));
    vaxHW = data(vaxHWL,:);    
    vaxUnML = logical(strcmp(data.Sex,'M') & ~strcmp(data.Smoker,'Never') & data.Age>=60 & strcmp(data.Vaxed,'y'));
    vaxUnM = data(vaxUnML,:);
    vaxUnWL = logical(strcmp(data.Sex,'F')  & ~strcmp(data.Smoker,'Never') & data.Age>=60 & strcmp(data.Vaxed,'y'));
    vaxUnW = data(vaxUnWL,:);

    % Balancing each group
    if strcmp(pairedTag,'y')
        if strcmp(subgroup,'Healthy')
            [vaxHM,vaxHW] = balance(vaxHM,vaxHW);
        elseif strcmp(subgroup,'Unhealthy')
            [vaxUnM,vaxUnW] = balance(vaxUnM,vaxUnW);
        end
    end

    if strcmp(subgroup, 'Healthy')
        G1 = combineSorted(infHW,infHM);
        G2 = combineSorted(vaxHW,vaxHM);
        G = combineSorted(G1,G2);
    elseif strcmp(subgroup, 'Unhealthy')
        G1 = combineSorted(infUnW,infUnM);
        G2 = combineSorted(vaxUnW,vaxUnM);
        G = combineSorted(G1,G2);
    end
    elseif strcmp(group,'Race')
        %% Infected

        % White non hispanic
        infWNHL = logical(strcmp(data.Race,'W') & strcmp(data.Ethinicity,'Non') & strcmp(data.Vaxed,'n'));
        infWNH = data(infWNHL,:);
        % Black non hispanic
        infBNHL = logical(strcmp(data.Race,'B') & strcmp(data.Ethinicity,'Non') & strcmp(data.Vaxed,'n'));
        infBNH = data(infBNHL,:);
        % All Hispanic
        infHIS = data(strcmp(data.Ethinicity,'His'),:);
        infHIS = infHIS(strcmp(infHIS.Vaxed,'n'),:);
        % % Not hispanic, white, or black
        logicalOTHER = logical(~strcmp(data.Race,'W') & ~strcmp(data.Race,'B'));
        infOTHER = data(logicalOTHER,:);
        infOTHER = infOTHER(strcmp(infOTHER.Ethinicity,'Non'),:);

        L = logical(infWNH.Age < 55 & strcmp(infWNH.Sex,'F'));
        infWNHWbelow55 = infWNH(L,:);
        L = logical(infWNH.Age >= 55 & strcmp(infWNH.Sex,'F'));
        infWNHWabove55 = infWNH(L,:);
        L = logical(infWNH.Age < 55 & strcmp(infWNH.Sex,'M'));
        infWNHMbelow55 = infWNH(L,:);
        L = logical(infWNH.Age >= 55 & strcmp(infWNH.Sex,'M'));
        infWNHMabove55 = infWNH(L,:);
 
        L = logical(infBNH.Age < 55 & strcmp(infBNH.Sex,'F'));
        infBNHWbelow55 = infBNH(L,:);
        L = logical(infBNH.Age >= 55 & strcmp(infBNH.Sex,'F'));
        infBNHWabove55 = infBNH(L,:);
        L = logical(infBNH.Age < 55 & strcmp(infBNH.Sex,'M'));
        infBNHMbelow55 = infBNH(L,:);
        L = logical(infBNH.Age >= 55 & strcmp(infBNH.Sex,'M'));
        infBNHMabove55 = infBNH(L,:);

        L = logical(infHIS.Age < 55 & strcmp(infHIS.Sex,'F'));
        infHISWbelow55 = infHIS(L,:);
        L = logical(infHIS.Age >= 55 & strcmp(infHIS.Sex,'F'));
        infHISWabove55 = infHIS(L,:);
        L = logical(infHIS.Age < 55 & strcmp(infHIS.Sex,'M'));
        infHISMbelow55 = infHIS(L,:);
        L = logical(infHIS.Age >= 55 & strcmp(infHIS.Sex,'M'));
        infHISMabove55 = infHIS(L,:);
        
        L = logical(infOTHER.Age < 55 & strcmp(infOTHER.Sex,'F'));
        infOTHERWbelow55 = infOTHER(L,:);
        L = logical(infOTHER.Age >= 55 & strcmp(infOTHER.Sex,'F'));
        infOTHERWabove55 = infOTHER(L,:);
        L = logical(infOTHER.Age < 55 & strcmp(infOTHER.Sex,'M'));
        infOTHERMbelow55 = infOTHER(L,:);
        L = logical(infOTHER.Age >= 55 & strcmp(infOTHER.Sex,'M'));
        infOTHERMabove55 = infOTHER(L,:);

        listWNH = {infWNHMbelow55, infWNHWbelow55,infWNHMabove55, infWNHWabove55};
        smallestWNH = findSmallest(listWNH);
        listBNH = {infBNHMbelow55, infBNHWbelow55,infBNHMabove55, infBNHWabove55};
        smallestBNH = findSmallest(listBNH);
        listBNHW = {infBNHWbelow55, infBNHWabove55};
        smallestBNHW = findSmallest(listBNHW);
        listHIS = {infHISMbelow55, infHISWbelow55,infHISMabove55, infHISWabove55};
        smallestHIS = findSmallest(listHIS);
        listHISW = {infHISWbelow55,infHISWabove55};
        smallestHISW = findSmallest(listHISW);
        listOTHER = {infOTHERMbelow55, infOTHERWbelow55,infOTHERMabove55, infOTHERWabove55};
        smallestOTHER = findSmallest(listOTHER);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'White Non Hispanic')
                % White Non Hispanic
                [smallestWNH, infWNHWbelow55] = balance(smallestWNH, infWNHWbelow55);
                [smallestWNH, infWNHMbelow55] = balance(smallestWNH, infWNHMbelow55);
                [smallestWNH, infWNHWabove55] = balance(smallestWNH, infWNHWabove55);
                [smallestWNH, infWNHMabove55] = balance(smallestWNH, infWNHMabove55);
            elseif strcmp(subgroup, 'Black Non Hispanic')
                % Black Non Hispanic
                [smallestBNH, infBNHWbelow55] = balance(smallestBNH, infBNHWbelow55);
                [smallestBNH, infBNHMbelow55] = balance(smallestBNH, infBNHMbelow55);
                [smallestBNH, infBNHWabove55] = balance(smallestBNH, infBNHWabove55);
                [smallestBNH, infBNHMabove55] = balance(smallestBNH, infBNHMabove55);
    
            elseif strcmp(subgroup, 'Black Non Hispanic Women')
                % Black Non Hispanic
                [smallestBNHW, infBNHWbelow55] = balance(smallestBNHW, infBNHWbelow55);
    
                [smallestBNHW, infBNHWabove55] = balance(smallestBNHW, infBNHWabove55);
        
            elseif strcmp(subgroup, 'Hispanic')
                % Hispanic
                [smallestHIS, infHISWbelow55] = balance(smallestHIS, infHISWbelow55);
                [smallestHIS, infHISMbelow55] = balance(smallestHIS, infHISMbelow55);
                [smallestHIS, infHISWabove55] = balance(smallestHIS, infHISWabove55);
                [smallestHIS, infHISMabove55] = balance(smallestHIS, infHISMabove55);
            elseif strcmp(subgroup, 'Hispanic Women')
                % Hispanic
                [smallestHIS, infHISWbelow55] = balance(smallestHIS, infHISWbelow55);
                [smallestHIS, infHISWabove55] = balance(smallestHIS, infHISWabove55);
            elseif strcmp(subgroup, 'Other')
                % Other
                [smallestOTHER, infOTHERWbelow55] = balance(smallestOTHER, infOTHERWbelow55);
                [smallestOTHER, infOTHERMbelow55] = balance(smallestOTHER, infOTHERMbelow55);
                [smallestOTHER, infOTHERWabove55] = balance(smallestOTHER, infOTHERWabove55);
                [smallestOTHER, infOTHERMabove55] = balance(smallestOTHER, infOTHERMabove55);
            end
        end

        %% Vaccinated
        %Removing nans
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        % Vaxinated data
        % White non hispanic
        vaxWNHL = logical(strcmp(data.Race,'W') & strcmp(data.Ethinicity,'Non') & strcmp(data.Vaxed,'y'));
        vaxWNH = data(vaxWNHL,:);
        % Black non hispanic
        vaxBNHL = logical(strcmp(data.Race,'B') & strcmp(data.Ethinicity,'Non') & strcmp(data.Vaxed,'y'));
        vaxBNH = data(vaxBNHL,:);
        % All Hispanic
        vaxHIS = data(strcmp(data.Ethinicity,'His'),:);
        vaxHIS = vaxHIS(strcmp(vaxHIS.Vaxed,'y'),:);
        % % Not hispanic, white, or black
        logicalOTHER = logical(~strcmp(data.Race,'W') & ~strcmp(data.Race,'B'));
        vaxOTHER = data(logicalOTHER,:);
        vaxOTHER = vaxOTHER(strcmp(vaxOTHER.Ethinicity,'Non'),:);

        L = logical(vaxWNH.Age < 55 & strcmp(vaxWNH.Sex,'F'));
        vaxWNHWbelow55 = vaxWNH(L,:);
        L = logical(vaxWNH.Age >= 55 & strcmp(vaxWNH.Sex,'F'));
        vaxWNHWabove55 = vaxWNH(L,:);
        L = logical(vaxWNH.Age < 55 & strcmp(vaxWNH.Sex,'M'));
        vaxWNHMbelow55 = vaxWNH(L,:);
        L = logical(vaxWNH.Age >= 55 & strcmp(vaxWNH.Sex,'M'));
        vaxWNHMabove55 = vaxWNH(L,:);
 
        L = logical(vaxBNH.Age < 55 & strcmp(vaxBNH.Sex,'F'));
        vaxBNHWbelow55 = vaxBNH(L,:);
        L = logical(vaxBNH.Age >= 55 & strcmp(vaxBNH.Sex,'F'));
        vaxBNHWabove55 = vaxBNH(L,:);
        L = logical(vaxBNH.Age < 55 & strcmp(vaxBNH.Sex,'M'));
        vaxBNHMbelow55 = vaxBNH(L,:);
        L = logical(vaxBNH.Age >= 55 & strcmp(vaxBNH.Sex,'M'));
        vaxBNHMabove55 = vaxBNH(L,:);

        L = logical(vaxHIS.Age < 55 & strcmp(vaxHIS.Sex,'F'));
        vaxHISWbelow55 = vaxHIS(L,:);
        L = logical(vaxHIS.Age >= 55 & strcmp(vaxHIS.Sex,'F'));
        vaxHISWabove55 = vaxHIS(L,:);
        L = logical(vaxHIS.Age < 55 & strcmp(vaxHIS.Sex,'M'));
        vaxHISMbelow55 = vaxHIS(L,:);
        L = logical(vaxHIS.Age >= 55 & strcmp(vaxHIS.Sex,'M'));
        vaxHISMabove55 = vaxHIS(L,:);
        
        L = logical(vaxOTHER.Age < 55 & strcmp(vaxOTHER.Sex,'F'));
        vaxOTHERWbelow55 = vaxOTHER(L,:);
        L = logical(vaxOTHER.Age >= 55 & strcmp(vaxOTHER.Sex,'F'));
        vaxOTHERWabove55 = vaxOTHER(L,:);
        L = logical(vaxOTHER.Age < 55 & strcmp(vaxOTHER.Sex,'M'));
        vaxOTHERMbelow55 = vaxOTHER(L,:);
        L = logical(vaxOTHER.Age >= 55 & strcmp(vaxOTHER.Sex,'M'));
        vaxOTHERMabove55 = vaxOTHER(L,:);

        listWNH = {vaxWNHMbelow55, vaxWNHWbelow55,vaxWNHMabove55, vaxWNHWabove55};
        smallestWNH = findSmallest(listWNH);
        listBNH = {vaxBNHMbelow55, vaxBNHWbelow55,vaxBNHMabove55, vaxBNHWabove55};
        smallestBNH = findSmallest(listBNH);
        listBNHW = {vaxBNHWbelow55, vaxBNHWabove55};
        smallestBNHW = findSmallest(listBNHW);
        listHIS = {vaxHISMbelow55, vaxHISWbelow55,vaxHISMabove55, vaxHISWabove55};
        smallestHIS = findSmallest(listHIS);
                listHISW = {vaxHISWbelow55,vaxHISWabove55};
        smallestHISW = findSmallest(listHISW);
        listOTHER = {vaxOTHERMbelow55, vaxOTHERWbelow55,vaxOTHERMabove55, vaxOTHERWabove55};
        smallestOTHER = findSmallest(listOTHER);

        % Balancing each group
        if strcmp(pairedTag,'y')
            if strcmp(subgroup, 'White Non Hispanic')
                % White Non Hispanic
                [smallestWNH, vaxWNHWbelow55] = balance(smallestWNH, vaxWNHWbelow55);
                [smallestWNH, vaxWNHMbelow55] = balance(smallestWNH, vaxWNHMbelow55);
                [smallestWNH, vaxWNHWabove55] = balance(smallestWNH, vaxWNHWabove55);
                [smallestWNH, vaxWNHMabove55] = balance(smallestWNH, vaxWNHMabove55);
            elseif strcmp(subgroup, 'Black Non Hispanic')
                % Black Non Hispanic
                [smallestBNH, vaxBNHWbelow55] = balance(smallestBNH, vaxBNHWbelow55);
                [smallestBNH, vaxBNHMbelow55] = balance(smallestBNH, vaxBNHMbelow55);
                [smallestBNH, vaxBNHWabove55] = balance(smallestBNH, vaxBNHWabove55);
                [smallestBNH, vaxBNHMabove55] = balance(smallestBNH, vaxBNHMabove55);
            elseif strcmp(subgroup, 'Black Non Hispanic Women')
                % Black Non Hispanic
                [smallestBNHW, vaxBNHWbelow55] = balance(smallestBNHW, vaxBNHWbelow55);
              
                [smallestBNHW, vaxBNHWabove55] = balance(smallestBNHW, vaxBNHWabove55);
        
            
            elseif strcmp(subgroup, 'Hispanic')
                % Hispanic
                [smallestHIS, vaxHISWbelow55] = balance(smallestHIS, vaxHISWbelow55);
                [smallestHIS, vaxHISMbelow55] = balance(smallestHIS, vaxHISMbelow55);
                [smallestHIS, vaxHISWabove55] = balance(smallestHIS, vaxHISWabove55);
                [smallestHIS, vaxHISMabove55] = balance(smallestHIS, vaxHISMabove55);
            elseif strcmp(subgroup, 'Hispanic Women')
                % Hispanic
                [smallestHIS, vaxHISWbelow55] = balance(smallestHIS, vaxHISWbelow55);
                [smallestHIS, vaxHISWabove55] = balance(smallestHIS, vaxHISWabove55);
            elseif strcmp(subgroup, 'Other')
                % Other
                [smallestOTHER, vaxOTHERWbelow55] = balance(smallestOTHER, vaxOTHERWbelow55);
                [smallestOTHER, vaxOTHERMbelow55] = balance(smallestOTHER, vaxOTHERMbelow55);
                [smallestOTHER, vaxOTHERWabove55] = balance(smallestOTHER, vaxOTHERWabove55);
                [smallestOTHER, vaxOTHERMabove55] = balance(smallestOTHER, vaxOTHERMabove55);
            end
        end

        if strcmp(subgroup, 'White Non Hispanic')
            G1 = combineSorted(infWNHWbelow55, infWNHMbelow55);
            G2 = combineSorted(infWNHWabove55, infWNHMabove55);
            G3 = combineSorted(vaxWNHWbelow55,vaxWNHMbelow55);
            G4 = combineSorted(vaxWNHWabove55,vaxWNHMabove55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup, 'Black Non Hispanic')
            G1 = combineSorted(infBNHWabove55,infBNHMabove55);
            G2 = combineSorted(infBNHMbelow55,infBNHWbelow55);
            G3 = combineSorted(vaxBNHWabove55,vaxBNHMabove55);
            G4 = combineSorted(vaxBNHWbelow55,vaxBNHMbelow55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup, 'Black Non Hispanic Women')
            G1 = combineSorted(infBNHWabove55,infBNHWbelow55);
            G3 = combineSorted(vaxBNHWabove55,vaxBNHWbelow55);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup, 'Hispanic')
            G1 = combineSorted(infHISWabove55,infHISMabove55);
            G2 = combineSorted(infHISMbelow55,infHISWbelow55);
            G3 = combineSorted(vaxHISWabove55,vaxHISMabove55);
            G4 = combineSorted(vaxHISWbelow55,vaxHISMbelow55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup, 'Hispanic Women')
            G1 = combineSorted(infHISWabove55,infHISWbelow55);
            G3 = combineSorted(vaxHISWabove55,vaxHISWbelow55);
            G = combineSorted(G1, G3);
        elseif strcmp(subgroup, 'Other')
            G1 = combineSorted(infOTHERWabove55,infOTHERMabove55);
            G2 = combineSorted(infOTHERMbelow55,infOTHERWbelow55);
            G3 = combineSorted(vaxOTHERWabove55,vaxOTHERMabove55);
            G4 = combineSorted(vaxOTHERWbelow55,vaxOTHERMbelow55);
            G1 = combineSorted(G1,G2);
            G3 = combineSorted(G3,G4);
            G = combineSorted(G1, G3);
        end
    elseif strcmp(group,'All')
        %% Infected
        infL = logical(strcmp(data.Vaxed,'n'));
        inf = data(infL,:);

        L = logical(inf.Age < 55 & strcmp(inf.Sex,'F'));
        infWbelow55 = inf(L,:);
        L = logical(inf.Age >= 55 & strcmp(inf.Sex,'F'));
        infWabove55 = inf(L,:);
        L = logical(inf.Age < 55 & strcmp(inf.Sex,'M'));
        infMbelow55 = inf(L,:);
        L = logical(inf.Age >= 55 & strcmp(inf.Sex,'M'));
        infMabove55 = inf(L,:);

        list = {infMbelow55, infWbelow55, infMabove55, infWabove55};
        smallest = findSmallest(list);


        % Balancing each group
        if strcmp(pairedTag,'y')
        % All
        [smallest,infWbelow55]=balance(smallest,infWbelow55);
        [smallest,infMbelow55]=balance(smallest,infMbelow55);
        [smallest,infWabove55]=balance(smallest,infWabove55);
        [smallest,infMabove55]=balance(smallest,infMabove55);
        end
        %% Vaccinated
        %Removing nans
        data(isnan(data.Days_Start),:) = [];
        data(isnan(str2double(data.Days_LastVax)),:) = [];

        % Removing people who were vaccinated before getting infected
        timeV = str2double(data.Days_LastVax);
        timeSS = data.Days_Start; %time since start, 
        vax_inf_logical = logical(timeSS - timeV < 0);
        data(vax_inf_logical,:) = [];

        vaxL = logical(strcmp(data.Vaxed,'y'));
        vax = data(vaxL,:);

        L = logical(vax.Age < 55 & strcmp(vax.Sex,'F'));
        vaxWbelow55 = vax(L,:);
        L = logical(vax.Age >= 55 & strcmp(vax.Sex,'F'));
        vaxWabove55 = vax(L,:);
        L = logical(vax.Age < 55 & strcmp(vax.Sex,'M'));
        vaxMbelow55 = vax(L,:);
        L = logical(vax.Age >= 55 & strcmp(vax.Sex,'M'));
        vaxDMabove55 = vax(L,:);

        list = {vaxMbelow55, vaxWbelow55, vaxDMabove55, vaxWabove55};
        smallest = findSmallest(list);

        % Balancing each group
        if strcmp(pairedTag,'y')
        % All
        [smallest,vaxWbelow55]=balance(smallest,vaxWbelow55);
        [smallest,vaxMbelow55]=balance(smallest,vaxMbelow55);
        [smallest,vaxWabove55]=balance(smallest,vaxWabove55);
        [smallest,vaxDMabove55]=balance(smallest,vaxDMabove55);
        end
     
        G1 = combineSorted(infMabove55,infWabove55);
        G2 = combineSorted(infMbelow55,infWbelow55);
        G3 = combineSorted(vaxWbelow55,vaxMbelow55);
        G4 = combineSorted(vaxWabove55,vaxDMabove55);
        G1 = combineSorted(G1,G2);
        G3 = combineSorted(G3,G4);
        G = combineSorted(G1, G3);
    end
    
end

function [G1,G2] = balance(g1,g2)
    % Take in two arrays, randomly removes elements from the larger array 
    % and returns two arrays of the same size
        
        % Selecting the larger array
        if size(g1,1) > size(g2,1)
            N = size(g1,1); % generates vector of zero values
            random_num = zeros(N,1); % initialize with zeros
            num_true_val =  size(g2,1); % wanted number of true values
            spots = randperm(N,num_true_val);
            random_num(spots) = 1; 
            L = logical(random_num);
    
            G1 = g1(L,:);
            G2 = g2;
        elseif size(g1,1) < size(g2,1)
            N = size(g2,1); % generates vector of zero values
            random_num = zeros(N,1); % initialize with zeros
            num_true_val =  size(g1,1); % wanted number of true values
            spots = randperm(N,num_true_val);
            random_num(spots) = 1; 
            L = logical(random_num);
    
            G2 = g2(L,:);
            G1 = g1;
        else
            G1 = g1;
            G2 = g2;
        end
end

function matrix = combineSorted(g1, g2)
    % % Takes in infected data and vax data and creates an ordered matrix
    % % based on ID
    size1 = size(g1,1);
    size2 = size(g2,1);
    matrix = [g1;g2];
    g1ID = g1.Subject_ID;
    g2ID = g2.Subject_ID;
    i = 1;
    j = 1;
    n = 1;
    a = true;

    if(~size1 == 0 && ~size2 == 0)
        while a
            if (g1ID(i)== g2ID(j))
                if g1.Days_Start(i) < g2.Days_Start(j)
                    matrix(n,:) = g1(i,:);
                    n = n + 1;
                    i = i + 1;
                else
                    matrix(n,:) = g2(j,:);
                    n = n + 1;
                    j = j + 1;
                end
            elseif (g1ID(i) < g2ID(j))
                matrix(n,:) = g1(i,:);
                n = n + 1;
                i = i + 1;
            elseif (g1ID(i) > g2ID(j))
                matrix(n,:) = g2(j,:);
                n = n + 1;
                j = j + 1;
            end
            if (j > size2)
                a = false;
                matrix = [matrix(1:n-1,:); g1(i:size1,:)];
            elseif(i > size1)
                a = false;
                matrix = [matrix(1:n-1,:); g2(j:size2,:)];
            end
        end
    end
end

function smallest = findSmallest(list)
    % Returns the index of the smallest array in the given list. 
    b = zeros(length(list),1);
    for i = 1:length(list)
        b(i) = size(list{i},1);
    end
    index=find(b==min(b));
    smallest = list{index(1)};
end
