function [y_mod,y_p,A] = modulo_measure_signal(m,z,R)
%edited 2/15/2017
n = length(z);
%% signal measurement
A = randn(m,n);
y = A*z; %measurements
%y_mod = mod(y,R);
y_p = (-sign(y)+1)/2; %actual phase

%modified modulo function
y_mod = y + y_p*R;
end