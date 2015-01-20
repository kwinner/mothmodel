% test EM CI

%parameters learned for population A1 in the Gross et al paper
mu     = 6.9;
sigma  = 3.4;
lambda = 1/0.645;
alpha  = .32;
N      = 216;

realdata = sfs([fileparts(which(mfilename)),'/../../sfs2004.csv']);
t = [realdata(~isnan([realdata.A1])).t] + 1;

params = struct('N', N, 'alpha', alpha, 't', t);
theta_0 = struct('mu', mu, 'sigma', sigma, 'lambda', lambda);

use_ARS = true;
initial_MCMC_iterations = 5000;
initial_MCMC_burnin = 4000;

num_iters = 100;

thetas(num_iters) = theta_0;
theta_hessians = cell(num_iters,1);
states = cell(num_iters,1);
runtimes = cell(num_iters,1);
coverage = [0,0,0,0];
parfor iter = 1:num_iters
	y = sampleState(struct('mu',mu,'sigma',sigma,'lambda',lambda), struct('N',N,'t',t,'alpha',alpha));
	
    [theta, theta_hessian, state, runtime] = stochupEM(y, params, theta_0, 'use_ARS',true, 'initial_MCMC_iterations',initial_MCMC_iterations, 'initial_MCMC_burnin',initial_MCMC_burnin);
    
    thetas(iter) = theta;
    theta_hessians{iter} = theta_hessian;
    states{iter} = state;
    runtimes{iter} = runtime;
 
end

save EM_CI_data thetas theta_hessians states runtimes

theta_cov = inv(theta_hessians{end});
CI_width = abs(2.*sqrt(diag(theta_cov)))';

theta_error = abs(theta_zonn - [mu, sigma, lambda, N]);
coverage = coverage + (theta_error <= CI_width);