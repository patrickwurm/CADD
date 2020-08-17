function [stress_new, state_new, dsde] = Umatel(delta_eps, props, stress_old, state_old)

global global_data;
dim = global_data.dimension;

sig_yield = props(1); %Initial Yield Stress
K = props(2); %Isotropic hardening modulus
lambda = props(3); %First Lame constant
mu = props(4); %Shear modulus

delta = math.kronecker([dim dim]);
d_ijkl = ttt(delta,delta);

c_tensor = lambda * d_ijkl + 2*mu*permute(d_ijkl,[1 3 2 4]);

% old stress- and delta_eps deviator, old state
s_old         = stress_old - 1/3*ttt(stress_old,delta,[1 2],[1 2])*delta;
delta_eps_dev = delta_eps - 1/3*ttt(delta_eps,delta,[1 2],[1 2])*delta;
alpha_old     = state_old{1};
eps_p_old     = state_old{2};

% trial state
stress_trial = stress_old + ttt(c_tensor,delta_eps,[3 4],[1 2]);
s_trial      = s_old + 2*mu*delta_eps_dev;
eps_p_trial  = eps_p_old;
alpha_trial  = alpha_old;
f_trial      = norm(s_trial)-sqrt(2/3)*(sig_yield+K*alpha_old);

if(f_trial<=1e-12)
    stress_new   = stress_trial;
    state_new{1} = alpha_trial;
    state_new{2} = eps_p_trial;
    dsde         = c_tensor;
else    
    delta_gamma = f_trial/(2*mu+2/3*K);
    n_new       = s_trial/norm(s_trial);
    stress_new  = stress_trial - 2*mu*delta_gamma*n_new;
    eps_p_new   = eps_p_old + delta_gamma*n_new;
    alpha_new   = alpha_old + sqrt(2/3)*delta_gamma;
    
    % Einschub: Hilfsgrößen für dsde
    n_dyad_n = ttt(n_new,n_new);
    d_ijkl   = ttt(delta,delta); % eig. nicht notwendig, da oben
    d_ikjl   = permute(d_ijkl,[1 3 2 4]); % siehe oben
    
    dsde = c_tensor - 2*mu/(1+K/(3*mu))*n_dyad_n-2*mu*delta_gamma/norm(s_trial)*(2*mu*(d_ikjl-1/3*d_ijkl)-2*mu*n_dyad_n);
    state_new{1} = alpha_new;
    state_new{2} = eps_p_new;
end

end