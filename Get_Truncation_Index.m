%% 
function I = Get_Truncation_Index(Distr, Percentile, Resolution)

% This is a search function for homing in on index (denoted by "I") 
% associated with the given percentile (denoted by "Percentile") of some 
% distribution (denoted by "Distr"). 

I       = 1;                % the approximated index 
Sign	= 1;                % index is initially increasing
Step    = Resolution;       % Each increase is initially equal 
loop	= true;             % 
while loop
    
    % Stop condition: the search function starts with relatively big
    % Step changes to I, and decreases the "Step" over time, until it's less 
    % than 1. 
    if Step < 1         
        loop = false;
    end
    Step = ceil(Step);
    
    % update index
    switch Sign
        case -1      % Increase I untill the area under Distr exceeds the desired percentile    
            while sum(Distr(1:I)) > Percentile     
                I = max([1, I + (Sign*Step)]);
            end
        case 1      % Decrease I untill the area under Distr is less than the desired percentile
            while sum(Distr(1:I)) < Percentile      
                I = min([length(Distr), I + (Sign*Step)]);
            end
    end
    
    % Update search step
    Step  	= Step / 10;    % decrease the step size after each iteration
    Sign    = Sign * -1;    % flip the sign: if I was initially increasing, it will decrease next round, and vice versa. 
    
end