% Fuction for the non-linear constraint: Y^2sigma^2k+2thetaYk≤1
function [c, ceq] = nonlincon(param)
    % Estract parameters
    k_IG = param(1);
    sigma = param(2);
    theta = param(3);
    Y = param(4);
    
    %Impose the constraint
    c = Y^2 * sigma^2 * k_IG + 2 * theta * Y * k_IG - 1;
    
    % No equality constraint
    ceq = [];
end
