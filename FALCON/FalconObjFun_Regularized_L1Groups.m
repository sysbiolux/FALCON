 function [xval,fval]=FalconObjFun_Regularized_L1Groups(estim,k)
% FalconObjFun serves as the objective function for the optimisation.
% Apply the non-linear optimiser 'fmincon' with the default algorithm (interior-point)
% Return the optimised parameters values and fitting cost calculated from the sum-of-squared error (SSE)
% [xval,fval]=FalconObjFun(estim,k)

% :: Input ::
% estim      complete model definition
% k          initial guess for parameter values (initial condition)
% 
% :: Output ::
% fval       all fitting costs over different optmization rounds
% xval       all optimized parameter values over different optmization rounds
%
% :: Contact ::
% Prof. Thomas Sauter, University of Luxembourg, thomas.sauter@uni.lu
% Sebastien De Landtsheer, University of Luxembourg, sebastien.delandtsheer@uni.lu

[xval,fval]=fmincon(@nestedfun,k,estim.A,estim.b,estim.Aeq,estim.beq,estim.LB,estim.UB,[],estim.options);

    function [ diff ] = nestedfun(k)
        
    n=estim.NrStates;
    l=estim.Lambda;
    Reg=estim.RegMatrix;
    
    Var=0;
    for v=1:size(Reg,1)
        km=mean(k(Reg(v,:)));
        Var=Var+sum(abs(k(Reg(v,:))-km));
    end
    
    %initial and successive number of steps for evaluation
    %this still needs to be worked on
    if n<=25, initial_t=10; step_t=10;
    elseif n>25 && n<=100, initial_t=100; step_t=10;
    elseif n>100, initial_t=300; step_t=150;
    end
    
    if isfield(estim,'MaxTime')
        MaxTime=estim.MaxTime;
    else
        MaxTime=0;
    end    
    
    ma=estim.ma;
    mi=estim.mi;
    param_index=estim.param_index;
    kmap=k(estim.kInd)'; %extend k including boolean gates. Now k2 has the same length as param_index
    
    

    pd=param_index(~estim.FixBool,:); %remove fixed Boolean gates from param_index
    kA=kmap(find((pd(:,3)>0))); %remap k for activating links
    kI=kmap(find((pd(:,4)>0))); %remap k for inhibiting links
    for i=1:length(estim.IdxInAct) %for each activating link
        ma(estim.IdxOutAct(i),estim.IdxInAct(i))=kA(i); %map to ma
    end
    for i=1:length(estim.IdxInInh) %for each inhibiting link
        mi(estim.IdxOutInh(i),estim.IdxInInh(i))=kI(i); %map to mi
    end

    % Get Boolean gate indices in for-loop
    BoolMax=max(param_index(:,5)); %number of Boolean gates
    if BoolMax>0 %if there is at least one Boolean gate
        Gate_fill_indices=[];
        Gate_value_fill=[];

        for counter2=1:BoolMax %for each Boolean gate: trace the inputs and output
            gate_indices=find(ismember(param_index(:,5),counter2))';
            current_ma_value=unique(ma(unique(param_index(gate_indices,1)),param_index(gate_indices,2)));
            Current_fill_indices=[unique(param_index(gate_indices,1)), ... %output
                param_index(gate_indices(1),2), ... %input 1
                param_index(gate_indices(2),2), ... %input 2
                unique(param_index(gate_indices,6))]; % type of gate (OR = 1, AND = 2)
            Gate_value_fill=[Gate_value_fill; current_ma_value];
            Gate_fill_indices=[Gate_fill_indices; Current_fill_indices];
        end

        Temp_BoolVal=zeros(n,size(estim.Output_idx,1));
    end

    % Load input and output
    Input_index=estim.Input_idx;
    Output_index=estim.Output_idx;
    Inputs=estim.Input;
    Measurements=estim.Output;

    % Evaluation
    diff=0; % Initialize fitting cost
    TimeSoFar=tic;
    
    x=rand(n,size(Measurements,1)); %initial random values for the nodes
    x(Input_index(1,:),:)=Inputs';%(counter_exp,:); %fixing input nodes
    xmeas=Measurements; %expected values for output nodes
    break_point_ss=1;
    runs=initial_t; % Initialize number of time steps

    while break_point_ss
        if MaxTime>0
            if toc(TimeSoFar)>MaxTime
                break
            end
        end
        pre_x=x;
        for counter=1:runs
            if BoolMax>0 %calculate values for Boolean gates
                for counter2=1:size(Gate_fill_indices,1)
                    if Gate_fill_indices(counter2,4)==1 % OR gate

                        Temp_BoolVal(Gate_fill_indices(counter2,1),:)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2)).*...
                            (1-(1-(x(Gate_fill_indices(counter2,2),:))).*(1-x(Gate_fill_indices(counter2,3),:)));

                    elseif Gate_fill_indices(counter2,4)==2 % AND gate

                        Temp_BoolVal(Gate_fill_indices(counter2,1),:)=ma(Gate_fill_indices(counter2,1),Gate_fill_indices(counter2,2)).*...
                            (x(Gate_fill_indices(counter2,2),:).*x(Gate_fill_indices(counter2,3),:));
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%
            %%% core equation %%%
            x=(ma*x).*(ones(size(x))-mi*x);
            %%%%%%%%%%%%%%%%%%%%

            if BoolMax>0 %if there are Boolean gates
                x(Gate_fill_indices(:,1),:)=Temp_BoolVal(Gate_fill_indices(:,1),:);
            end                    
        end

        if any(sum(abs(pre_x-x))>estim.SSthresh) %if the network did not reach steady-state
            runs=runs+step_t;
        else %if we are at steady-state
            break_point_ss=0;
        end

    end
    xsim=x(Output_index(1,:),:)';
    mask=isnan(xmeas);
    xsim(mask)=0; xmeas(mask)=0;

    %calculate the sum-of-squared errors
    mse=(sum(sum((xsim-xmeas).^2)))/numel(estim.Output);
    diff=mse+l*Var;
    disp(['MSE: ', num2str(mse), ' ; reg cost: ',num2str(l*Var), ' ; Total: ', num2str(diff)])
    
    end

end
