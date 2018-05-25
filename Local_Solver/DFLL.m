function [output] = DFLL( m, input, model, U, LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr,flag)


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

Px       = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
Pe       = opKron(opExtension(model.n(2),model.nb(2),0),opExtension(model.n(1),model.nb(1),0));
% model parameter: slowness [s/m] on computational grid.
nu       = 1e-3*Px*sqrt(m);
dnu      = opDiag(.5*1e-6./(1e-3*sqrt(m)));

% distribute frequencies according to standard distribution
freq     = distributed(model.freq);

if flag==1

else
    spmd
        freqloc   = getLocalPart(freq);
        Gplus1loc = getLocalPart(Gplus1);
        Gdcloc    = getLocalPart(Gdc);
        Uloc      = getLocalPart(U);
        U0dcloc   = getLocalPart(U0dc);
        Resloc    = getLocalPart(reshape(input,[model.nsurf,model.nrec,nfreq]));
        nfreqloc  = length(freqloc);
        gloc      = zeros(nz*nx,1);
        Adloc     = zeros(nz*nx,model.nsurf);
        for k = 1:nfreqloc
            [Hk,dHk]                 = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb,model.mode);
            Hk                       = Helm_woPML(Hk,nz,nx);
            dHk                      = reshape(diag(dHk),nt(1),nt(2));
            dHk                      = spdiags(vec(dHk(2:end-1,2:end-1)),0,nz*nx,nz*nx);
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

            % moving this conjugate to outside DFL, in misfitlocal.m line
            % 15.
            RHS(1:nb,:)              = U0dcloc(:,:,k)*Resloc(:,:,k);
            
            % compute the forward wavefield in the local domain
            Ulocal                   = PDF\RHS;
            Ulocal                   = Ulocal(nb+1:end,:);
            
            Adloc(ind,:)             = -conj(Ulocal);
            r                        = real(sum(conj(Uloc(:,:,k)).*(dHk'*((2*pi*freqloc(k))'*dnu'*Adloc)),2));
            gloc                     = gloc + r(:);
        end
        output = pSPOT.utils.global_sum(gloc);  `11111111111111111111111114567
    end
end
output = output{1};
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
