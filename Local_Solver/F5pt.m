function D = F5pt(m,Q,model)

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q*model.W,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
% model parameter: slowness [s/m] on computational grid.
nu = 1e-3*Px*sqrt(m);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);

spmd
    codistr  = codistributor1d(2,[],[nsrc*nrec,nfreq]);
    freqloc  = getLocalPart(freq);
    wloc     = getLocalPart(w);
    nfreqloc = length(freqloc);
    Dloc     = zeros(nrec*nsrc,nfreqloc);
    for k = 1:nfreqloc
        Hk        = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb,model.mode);
        Uk        = Hk\(wloc(k)*(Ps'*Q*model.W));
        Dloc(:,k) = vec(Pr*Uk);
    end
    D = codistributed.build(Dloc,codistr,'noCommunication');
end

end
