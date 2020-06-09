 function [xval,fval,varargout]=FalconObjFun(estim,k)
% FalconObjFun serves as the objective function for the optimisation.
% Apply the non-linear optimiser 'fmincon' with the default algorithm (interior-point)
% Return the optimised parameters values and fitting cost calculated from the sum-of-squared error (SSE)
% [xval,fval, varargout]=FalconObjFun(estim,k)

% :: Input ::
% estim      complete model definition
% k          initial guess for parameter values (initial condition)
% 
% :: Output ::
% fval       all fitting costs over different optimization rounds
% xval       all optimized parameter values over different optmization rounds
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

MSE = []; AIC = []; Nparams = []; BIC = [];

[xval,fval] = fmincon(@nestedfun,k,estim.A,estim.b,estim.Aeq,estim.beq,estim.LB,estim.UB,[],estim.options);
% global MSE; global AIC; global Nparams; global BIC;


    function [ Diff ] = nestedfun(k)
        
    n = estim.NrStates;
    N = numel(estim.Output) - sum(sum(isnan(estim.Output)));    
    
    if isfield(estim, 'Lambda')
        l = estim.Lambda;
    else
        l = 0;
    end
    
    if isfield(estim, 'RegMatrix')
        if isfield(estim.RegMatrix, 'Groups')
            RegGroups = estim.RegMatrix.Groups;
        end
        if isfield(estim.RegMatrix, 'Smooth')
            RegSmooth = estim.RegMatrix.Smooth;
        end
        if isfield(estim.RegMatrix, 'Cluster')
            RegCluster = estim.RegMatrix.Cluster;
        end
    end
    
    if isfield(estim, 'Reg')
        if strcmp(estim.Reg, 'none')
            Var = 0;
        elseif strcmp(estim.Reg, 'L1')
            Var = sum(abs(k));
        elseif strcmp(estim.Reg, 'L1Groups')
            Var = sum(sum(abs(k(RegGroups) - repmat(mean(k(RegGroups), 2), 1, size(RegGroups, 2)))));
        elseif strcmp(estim.Reg, 'L1Smooth')
            Var = 0;
            for v = 1:size(RegSmooth, 1)
                Var = Var + sum(abs(k(RegSmooth(v, 2:end)) - k(RegSmooth(v, 1:end-1))));
            end
        elseif strcmp(estim.Reg, 'L2')
            Var = sum(k.^2);
        elseif strcmp(estim.Reg, 'L1/2')
            Var = sum(k.^0.5);
        elseif strcmp(estim.Reg, 'Lx')
            Var = 1/sum(k.^2);
        elseif strcmp(estim.Reg, 'Ldrug')
            Var(1) = sum(k.^0.5);
            Var(2) = 0;
            for v = 1:size(RegGroups, 1)
                km = mean(k(RegGroups(v, :)));
                Var(2) = Var(2) + sum(abs(k(RegGroups(v, :)) - km));
            end
        elseif strcmp(estim.Reg, 'LTriple')
            Var(1) = sum(k.^0.5);
            Var(2) = 0; [S,T] = size(RegCluster);
            kg = sort(k(RegCluster), 2, 'ascend');
            for c = 1:T-1
                for cc = c+1:T
                    Var(2) = Var(2) + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
                end
            end
            Var(2) = S/Var(2);
            Var(3) = 0;
            for v = 1:size(RegSmooth, 1)
                Var(3) = Var(3) + sum(abs(k(RegSmooth(v, 2:end)) - k(RegSmooth(v, 1:end-1))));
            end
        elseif strcmp(estim.Reg, 'LCluster')
            Var = 0; [S,T] = size(RegCluster);
            kg = sort(k(RegCluster), 2, 'ascend');
            for c = 1:T-1
                for cc = c+1:T
                    Var = Var + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
                end
            end
            Var = S/Var;
        elseif strcmp(estim.Reg, 'PruneCluster')
            VarA = 0; [S, T] = size(RegCluster);
            kg = sort(k(RegCluster), 2, 'ascend');
            for c = 1:T-1
                for cc = c+1:T
                    VarA = VarA + sum(abs((kg(:, cc) - kg(:, c)) - repmat((cc - c)/T, S, 1)));
                end
            end
            Var(2) = S/VarA;
            Var(1) = sum(k.^0.5);            
        end
    else
        Var = 0;
    end
    
    
    %initial and successive number of steps for evaluation
    %this still needs to be worked on
    if n <= 25, initial_t = 10; step_t = 10;
    elseif n > 25 && n <= 100, initial_t = 100; step_t = 10;
    elseif n > 100, initial_t = 300; step_t = 150;
    end
    
    if isfield(estim, 'MaxTime')
        MaxTime = estim.MaxTime;
    else
        MaxTime = 0;
    end    
    
    ma = estim.ma;
    mi = estim.mi;
    param_index = estim.param_index;
    kmap = k(estim.kInd)'; %extend k including boolean gates. Now k2 has the same length as param_index
    

    pd = param_index(~estim.FixBool, :); %remove fixed Boolean gates from param_index
    kA = kmap(find((pd(:, 3) > 0))); %remap k for activating links
    kI = kmap(find((pd(:, 4) > 0))); %remap k for inhibiting links
    for i = 1:length(estim.IdxInAct) %for each activating link
        ma(estim.IdxOutAct(i), estim.IdxInAct(i)) = kA(i); %map to ma
    end
    for i = 1:length(estim.IdxInInh) %for each inhibiting link
        mi(estim.IdxOutInh(i), estim.IdxInInh(i)) = kI(i); %map to mi
    end

    % Get Boolean gate indices in for-loop
    BoolMax = max(param_index(:,5)); %number of Boolean gates
    if BoolMax > 0 %if there is at least one Boolean gate
        Gate_fill_indices = [];
        Gate_value_fill = [];

        for counter2 = 1:BoolMax %for each Boolean gate: trace the inputs and output
            gate_indices = find(ismember(param_index(:, 5), counter2))';
            current_ma_value = unique(ma(unique(param_index(gate_indices, 1)),param_index(gate_indices, 2)));
            Current_fill_indices = [unique(param_index(gate_indices, 1)), ... %output
                param_index(gate_indices(1), 2), ... %input 1
                param_index(gate_indices(2), 2), ... %input 2
                unique(param_index(gate_indices, 6))]; % type of gate (OR = 1, AND = 2)
            Gate_value_fill = [Gate_value_fill; current_ma_value];
            Gate_fill_indices = [Gate_fill_indices; Current_fill_indices];
        end

        Temp_BoolVal = zeros(n,size(estim.Output_idx, 1));
    end

    % Load input and output
    Input_index = estim.Input_idx;
    Output_index = estim.Output_idx;
    Inputs = estim.Input;
    Measurements = estim.Output;
    Variances = (estim.SD).^2;

    % Evaluation
    TimeSoFar = tic;
    
    x = rand(n,size(Measurements, 1)); %initial random values for the nodes
    x(Input_index(1, :), :) = Inputs';%(counter_exp,:); %fixing input nodes
    xmeas = Measurements; %expected values for output nodes
    break_point_ss = 1;
    runs = initial_t; % Initialize number of time steps

    while break_point_ss
        if MaxTime > 0
            if toc(TimeSoFar) > MaxTime
                break
            end
        end
        pre_x = x;
        for counter = 1:runs
            if BoolMax > 0 %calculate values for Boolean gates
                for counter2 = 1:size(Gate_fill_indices, 1)
                    if Gate_fill_indices(counter2, 4) == 1 % OR gate

                        Temp_BoolVal(Gate_fill_indices(counter2, 1), :) = ma(Gate_fill_indices(counter2, 1), Gate_fill_indices(counter2, 2)).*...
                            (1-(1-(x(Gate_fill_indices(counter2, 2), :))) .* (1-x(Gate_fill_indices(counter2, 3), :)));

                    elseif Gate_fill_indices(counter2,4) == 2 % AND gate

                        Temp_BoolVal(Gate_fill_indices(counter2, 1), :) = ma(Gate_fill_indices(counter2, 1), Gate_fill_indices(counter2, 2)).*...
                            (x(Gate_fill_indices(counter2, 2), :) .* x(Gate_fill_indices(counter2, 3), :));
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%
            %%% core equation %%%
            x = (ma * x) .* (ones(size(x)) - mi * x);
            %%%%%%%%%%%%%%%%%%%%

            if BoolMax > 0 %if there are Boolean gates
                x(Gate_fill_indices(:, 1), :) = Temp_BoolVal(Gate_fill_indices(:, 1), :);
            end                    
        end

        if any(sum(abs(pre_x - x)) > estim.SSthresh) %if the network did not reach steady-state
            runs = runs + step_t;
        else %if we are at steady-state
            break_point_ss = 0;
        end

    end
    xsim = x(Output_index(1, :), :)';
    mask = isnan(xmeas);
    xsim(mask) = 0; xmeas(mask) = 0;
    
    %calculate the sum-of-squared errors
    Fun = 1;
    if isfield(estim, 'ObjFunction')
        
        if strcmp(estim.ObjFunction, 'weighted')
            Fun = 2;
        elseif strcmp(estim.ObjFunction, 'unweighted')
            Fun = 1;
        end
    end
    
    if Fun == 2
        MSE = (sum(nansum(((xsim-xmeas).^2)./Variances)))/N;
    else
        MSE = (sum(sum((xsim-xmeas).^2)))/N;
    end
    
    Diff = MSE + sum(l.*Var);    
    
    Nparams = (sum(k > 0.01));

    if isfield(estim, 'Reg')
        if strcmp(estim.Reg, 'L1Groups')
            Std_group = std(k(RegGroups), 0, 2);
            Collapsed = Std_group < 0.01;
            Nparams = sum(Collapsed) + size(RegGroups, 2) * sum(~Collapsed);
        elseif strcmp(estim.Reg, 'L1Smooth')
            Smoothed = (max(k(RegSmooth), [], 2) - min(k(RegSmooth), [], 2)) < 0.1;
            RP =(mean(k(RegSmooth), 2)) > 0.01;
            Nparams = sum(~Smoothed .* RP) + size(RegSmooth, 2) * sum(~Smoothed .* RP);
        elseif strcmp(estim.Reg, 'Ldrug')
            Std_group = std(k(RegGroups), 0, 2);
            Collapsed = Std_group < 0.01;
            RP=(mean(k(RegGroups), 2)) > 0.01;
            Nparams = sum(Collapsed .* RP) + size(RegGroups, 2) * sum(~Collapsed .* RP);            
        elseif strcmp(estim.Reg, 'LTriple')
                %%% TO DO: check with example
            
            ParamIdentity = (abs(bsxfun(@minus, k, k')) < 0.01);
            ParamMatters = zeros(length(k));
            for kline = 1:size(RegCluster,1)
                for kcol1 = 1:size(RegCluster,2) - 1
                    for kcol2 = kcol1 + 1:size(RegCluster, 2)
                        ParamMatters(RegCluster(kline, kcol1), RegCluster(kline, kcol2)) = 1;
                    end
                end
            end
            for kline = 1:size(RegSmooth, 1)
                for kcol1 = 1:size(RegSmooth, 2) - 1
                    for kcol2 = kcol1 + 1:size(RegSmooth, 2)
                        ParamMatters(RegSmooth(kline, kcol1), RegSmooth(kline, kcol2)) = 1;
                    end
                end
            end
            ParamMatters = ParamMatters .* (~eye(length(k)));
            ParamSignif = k' > 0.01;
            NCollapsed = ceil(sum((sum(ParamMatters.*ParamIdentity, 2) > 0) .* ParamSignif));
            Nparams = estim.NrParams - NCollapsed;
                %%% TO DO: check with example
            
        elseif strcmp(estim.Reg, 'LCluster')
            Par = (k(RegCluster));
            Par = sort(Par, 2);
            Dist = Par(:, 2:end) - Par(:, 1:end-1);
            Collapsed = Dist < 0.01;
            Nparams = size(RegCluster, 1) + sum(~Collapsed(:));
        elseif strcmp(estim.Reg, 'PruneCluster')
            Par = (k(RegCluster));
            Par = sort(Par, 2);
            Sig = Par > 0.01;
            Dist = Par(:, 2:end) - Par(:, 1:end-1);
            Collapsed = Dist < 0.01;
            Sig = Sig .* [ones(size(Collapsed, 1), 1), ~Collapsed];
            Nparams = sum(Sig(:));
        end
    end
    
    AIC = N .* log(MSE) + 2 .* Nparams;
    BIC = N .* log(MSE) + Nparams * log(N);

    fprintf('MSE= %d \t reg cost= %d \t total= %d \t BIC= %d \n', MSE, sum(l.*Var), Diff, BIC);

    end

    varargout{1} = AIC;
    varargout{2} = MSE;
    varargout{3} = BIC;
    varargout{4} = Nparams;
end
