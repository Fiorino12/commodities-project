function [dS1,dS2] = fastInverseGaussian(mu, lambda, Nsim)
    % Simulate samples from an Inverse Gaussian distribution efficiently
    % Input:
    % mu: mean parameter
    % lambda: shape parameter
    % Nsim: number of samples to generate
    %
    % Output:
    % dS1: inverse gaussian sampled using U
    % dS2: inverse gaussian sampled using 1 - U
    
    % Step 1: sample Z ~ N(0, 1)
    Z = randn(Nsim, 1);
    
    % Step 2: Compute Y = Z^2
    Y = Z.^2;
    
    % Step 3: Compute X1
    sqrt_term = sqrt(4 * mu * lambda * Y + mu^2 * Y.^2);
    X1 = mu + (mu^2 .* Y) ./ (2 * lambda) - (mu ./ (2 * lambda)) .* sqrt_term;
    
    % Step 4: Compute X using uniform distribiution
    U = rand(Nsim, 1); % Uniform [0, 1]
    dS1 = X1; % Pre-allocation
    dS1(U > (mu ./ (mu + X1))) = (mu^2) ./ X1(U > (mu ./ (mu + X1)));

    dS2=X1;
    dS2(1-U > (mu ./ (mu + X1))) = (mu^2) ./ X1(1-U > (mu ./ (mu + X1)));
end
