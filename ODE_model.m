%ODE SIS open model code

function [Classes] = ODE_model(para,ICs,mintime,maxtime)
%Extra functions
f  = @(M,k,z) M/((1+M*(1-z)/k)^(k+1));
phi = @(M,k,z) 1-((1+M*(1-z)/k)/(1+M*(2-z)/k))^(k+1);
%pi = @(x);



%Run ODE using ODE45
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, pop] = ode45(@diff_SIR_open_model, [mintime:1:maxtime], [ICs.M_c ICs.M_a ICs.l], opts, para);

%Convert output to structure
Classes = struct('M_c',pop(:,1),'M_a',pop(:,2),'l',pop(:,3),'t',t);


%Diff equations
function dPop = diff_SIR_open_model(t,pop,para)

%Variables
M_c=pop(1);
M_a=pop(2);
l=pop(3);

%Model

dM_c = para.beta_c*l - para.sigma*M_c;
dM_a = para.beta_a*l - para.sigma*M_a;
dl   = (para.R0*para.sigma*para.mu)/(para.beta_c*para.rho*para.n_c + para.beta_a*(1-para.rho)*para.n_a)...
    *(phi(M_c,para.k,para.z)*f(M_c,para.k,para.z)*para.n_c*para.rho + ...
    phi(M_a,para.k,para.z)*f(M_a,para.k,para.z)*para.n_a*(1-para.rho)) ...
    - para.mu*l;

dPop = [dM_c; dM_a ; dl];

end

end
