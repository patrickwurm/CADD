function [stress_new, state_new, dsde] = Umatiso(delta_eps, props, stress_old, state_old)

global global_data;
dim = global_data.dimension;

sig_yield = props(1); %Initial Yield Stress
K = props(2); %Isotropic hardening modulus
lambda = props(3); %First Lame constant
mu = props(4); %Shear modulus

delta = math.kronecker([dim dim]);
d_ijkl = ttt(delta,delta);

c_tensor = lambda * d_ijkl + 2*mu*permute(d_ijkl,[1 3 2 4]);

stress_trial = stress_old + ttt(c_tensor,delta_eps,[3 4],[1 2]);

stress_new = stress_trial;
state_new = state_old;
dsde = c_tensor;



end
