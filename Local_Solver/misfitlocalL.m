function [ f, g ] = misfitlocalL( m,Dobs, model, LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr,Uback)

% compute forward wavefield in local domain and its contribution to
% receiver locations
[D, U] = FLL(m,model,LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr,Uback);

% compute residual
Res     = Dobs(:) - D(:);

% objective value
f = 0.5*norm(Res(:)).^2;

% gradient
J  = oppDFLL(m, model, U, LD, Gplus1, Gdc, U0dc,Gplusr,Gdcr);
g  = J'*conj(Res(:)); % this conjugate comes from the formulation, It in general appears in line 99 in DFL.m
end

