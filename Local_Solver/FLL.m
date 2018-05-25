function [D, U] = FLL(m,model,LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr,Uback)

nz       = LD(1,2) - LD(1,1)+1;
nx       = LD(2,2) - LD(2,1)+1;
index    = spiralfunction(nz,nx);
[~,ind]  = sort(vec(index));
nb       = 2*(nx+nz-2); % number of nodes on the boundary
nn       = nx*nz; % number of nodes in the truncated domain
% comp. grid
model.nb = [1 1 0];
model.n  = [nz nx 1];
ot       = model.o-model.nb.*model.d;
dt       = model.d;
nt       = model.n+2*model.nb;
nfreq    = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
% model parameter: slowness [s/m] on computational grid.
nu = 1e-3*Px*sqrt(m);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);

spmd
    codistr   = codistributor1d(3,[],[model.nrec,model.nsurf,nfreq]);
    codistr1  = codistributor1d(3,[],[nz*nx,model.nsurf,nfreq]);
    freqloc   = getLocalPart(freq);
    wloc      = getLocalPart(w);
    Gplus1loc = getLocalPart(Gplus1);
    Gdcloc    = getLocalPart(Gdc);
    Gplusrloc = getLocalPart(Gplusr);
    Gdcrloc   = getLocalPart(Gdcr);
    U0dcloc   = getLocalPart(U0dc);
    Ubackloc  = getLocalPart(Uback);
    nfreqloc  = length(freqloc);
    Dloc      = zeros(model.nrec,model.nsurf,nfreqloc);
    Uloc      = zeros(nz*nx,model.nsurf,nfreqloc);
    for k = 1:nfreqloc
        Hk                       = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb,model.mode);
        Hk                       = Helm_woPML(Hk,nz,nx);
        Iden                     = eye(nb);
        % intialize the matrix in equation 10
        nt                       = nn+nb;
        PDF                      = zeros(nt);

        % first relation in List1
        PDF(1:nb,1:2*nb)         = [-Iden Iden];
        % second relation in List 1
        PDF(nb+1:2*nb,1:3*nb-8)  = [Iden (1/model.d(1)^2)*Gplus1loc(:,:,k) -(1/model.d(1)^2)*Gdcloc(:,:,k)];

        % Third relation in List 1
        PDF(2*nb+1:end,nb+1:end) = -Hk(nb+1:end,:);
        
        % right hand side of equation 10
        RHS                      = zeros(nn+nb,model.nsurf);
        RHS(1:nb,:)              = wloc(:,k)*U0dcloc(:,:,k);
        
        % compute the forward wavefield in the local domain
        Ulocal                   = PDF\RHS;
        Ulocal                   = Ulocal(nb+1:end,:);
        Uloc(ind,:,k)            = Ulocal;
        
        % propogate forward wavefield from boundary to receiver locations
        Usc_Res                  = -(1/model.d(1)^2)*Gplusrloc(:,:,k)*Ulocal(1:nb,:) + (1/model.d(1)^2)*Gdcrloc(:,:,k)*Ulocal(nb+1:2*nb-8,:); 
        Dloc(:,:,k)              = wloc(:,k)*Ubackloc(:,:,k)+Usc_Res;
    end
    D = codistributed.build(Dloc,codistr,'noCommunication');
    U = codistributed.build(Uloc,codistr1,'noCommunication');
end

end



function Hk = Helm_woPML(Hk,nz,nx) % extract 
% extract the Helm part without PML since current code need atlest 1
% boundary point to compute helmholtz coefficients
nz            = nz+2; % number of spatial points in local domain
nx            = nx+2; % number of depth point in local domain
index         = spiralfunction(nz,nx);
[~,ind]       = sort(vec(index));
nb            = 2*(nx+nz-2); % number of nodes on the boundary
% extract hemlholtz coeffcients in the interior of the truncated domain
% make sure helmholtz operator is defined for spiral format. rotate the
% rows of helmholtz operator first to make sure each row point to the
% spiral scheme, then, within each row permute the indices so that they
% follow the spiral scheme
Hk            = Hk(ind,:);
Hk            = Hk(:,ind);
Hk            = Hk(nb+1:end,nb+1:end);
end

