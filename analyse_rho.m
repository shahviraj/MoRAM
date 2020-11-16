% This file is for experimenting with various new ideas : mainly for EURASI
% P response
clear all;
clc;
tic;
rng ('shuffle')
pr = struct;
%Fixed parameters
pr.n = 1000; %length of the input signal
pr.b = 1; %number of blocks if signal is block-sparse; otherwise keep 1
pr.tol1 = 1e-8; %error tolerance for measurements
pr.tol2 = 1e-8;
pr.max_iter = 15;
pr.flip_frac = 0.01;
pr.R = 4; %period of the modulo function
pr.rho = 3;%spread of the true measurements, y =A*z
pr.del = 1; %truncation factor for supp estimation
pr.spgl_opts = spgSetParms('verbosity',0);
%Tuned parameters
%pr.mspan1 = [100:100:2000];
%pr.mspan2 = [600:100:1000];
%pr.mspan=[pr.mspan1,pr.mspan2];
pr.mspan =400:400:400;
pr.num_trials = 10;
pr.s_span = 6:6:24; % sparsity
pr.amp = 1; %amplification factor 
pr.del_p = 0.15; % ps = del*m (sparsity pertaining to error in p)
all_rho = zeros(length(pr.mspan),length(pr.s_span),pr.num_trials);
for j = 1:length(pr.mspan)
    m = pr.mspan(j);
    
    for k = 1:length(pr.s_span)
        s = pr.s_span(k);
        
        for l = 1:pr.num_trials
            %Generate the ground truth signal
            [z,z_ind] =  generate_signal(pr.n,s,pr.b, pr.amp);
            [y_mod, y, A] = modulo_measure_signal(m,z,pr.R);
            
            all_rho(j,k,l) = max(abs(y));
            
        end
    end
end

disp("Rho statistics")
disp("mean rho")
disp(mean2(all_rho))