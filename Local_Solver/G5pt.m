function D = G5pt(m,Q,model)

% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nfreq  = length(model.freq);

% mapping from source/receiver/physical grid to comp. grid
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));
Pe = opKron(opExtension(model.n(2),model.nb(1,2),0),opExtension(model.n(1),model.nb(1,1),0));

% model parameter: slowness [s/m] on computational grid.
nu = 1e-3*Px*sqrt(m);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);

spmd
    codistr  = codistributor1d(3,[],[prod(model.n),nsrc,nfreq]);
    freqloc  = getLocalPart(freq);
    nfreqloc = length(freqloc);
    Dloc     = zeros(prod(model.n),nsrc,nfreqloc);
    for k = 1:nfreqloc
        Hk          = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb,model.mode);
        Uk          = Hk\(Ps'*Q);
        Dloc(:,:,k) = Pe'*Uk;
    end
    D = codistributed.build(Dloc,codistr,'noCommunication');
end

end
