function [y_mod,y, A] = multi_period_modulo_measure_signal(m,z,R)
%edited 2/15/2017
n = length(z);
%% signal measurement
A = randn(m,n);
y = A*z; %measurements
y_mod = mod(y,R);

end