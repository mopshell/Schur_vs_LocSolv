function [Gdc, Gplus1, U0dc, Gdcr, Gplusr, Uback] = Back_GreenFcnL( m0, model, LD )

% added spiral function which works with rectangular or any dimension, also
% added the low-rank approximation of Green's function

model.boundG  = 2;
% generate source side Green's function
Q             = -1*speye(length(model.xsrc)); %source weight (simultaneous sources)


DSrc          = G5pt(m0,Q,model);
% compute the contribution of background source at receiver locations
Uback         = reshape(DSrc,model.n(1),model.n(2),model.nsurf,length(model.freq));
Uback         = distributed(reshape(squeeze(gather(Uback(model.zind,model.rind,:,:))),model.nrec,model.nsurf,length(model.freq)));
% receiver side Greens function
DRec          = DSrc;

% on local-domain---Left
model.zsrc   = model.zt(LD(1,1):LD(1,2));
model.xsrc   = model.xt(LD(2,1));
model.ns     = length(model.zsrc);
Q            = -1*speye(model.ns); %source weight (simultaneous sources)
DL           = G5pt(m0,Q,model);

% on local-domain---bottom
model.zsrc   = model.zt(LD(1,2));
model.xsrc   = model.xt(LD(2,1)+1:LD(2,2)-1);
model.ns     = length(model.xsrc);
Q            = -1*speye(model.ns); %source weight (simultaneous sources)
DB           = G5pt(m0,Q,model);

% on local-domain---right
model.zsrc   = model.zt(LD(1,1):LD(1,2));
model.xsrc   = model.xt(LD(2,2));
model.ns     = length(model.zsrc);
Q            = -1*speye(model.ns); %source weight (simultaneous sources)
DR           = G5pt(m0,Q,model);

% on local-domain---top
model.zsrc   = model.zt(LD(1,1));
model.xsrc   = model.xt(LD(2,1)+1:LD(2,2)-1);
model.ns     = length(model.xsrc);
Q            = -1*speye(model.ns); %source weight (simultaneous sources)
DT           = G5pt(m0,Q,model);


%% form the Green's function matrix
% distribute frequencies according to standard distribution
freq          = distributed(model.freq);
nz            = LD(1,2) - LD(1,1)+1;
nx            = LD(2,2) - LD(2,1)+1;
index         = spiralfunction(nz,nx);
[~,ind]       = sort(vec(index));
nb            = 2*(nx+nz-2); % number of nodes on the boundary
nfreq         = length(model.freq);
spmd
    codistrGdc    = codistributor1d(3,[],[nb,nb-8,nfreq]);
    codistrGplus1 = codistributor1d(3,[],[nb,nb,nfreq]);
    codistrGdcr   = codistributor1d(3,[],[model.nsurf,nb-8,nfreq]);
    codistrGplusr = codistributor1d(3,[],[model.nsurf,nb,nfreq]);
    codistrU0dc   = codistributor1d(3,[],[nb,model.nsurf,nfreq]);
    DSrcloc       = getLocalPart(DSrc);
    DRecloc       = getLocalPart(DRec);
    DLloc         = getLocalPart(DL);
    DRloc         = getLocalPart(DR);
    DBloc         = getLocalPart(DB);
    DTloc         = getLocalPart(DT);
    nfreqloc      = length(getLocalPart(freq));
    Gdc           = zeros(nb,nb-8,nfreqloc);
    Gplus1        = zeros(nb,nb,nfreqloc);
    Gdcr          = zeros(model.nsurf,nb-8,nfreqloc);
    Gplusr        = zeros(model.nsurf,nb,nfreqloc);
    U0dc          = zeros(nb,model.nsurf,nfreqloc);
    for k = 1:nfreqloc
            [Gdc(:,:,k), Gplus1(:,:,k), U0dc(:,:,k), Gdcr(:,:,k), Gplusr(:,:,k)] = Greenfnc_mat(DSrcloc(:,:,k), DRecloc(:,:,k), DLloc(:,:,k), DBloc(:,:,k), DRloc(:,:,k), DTloc(:,:,k), nx, nz, nb, ind, model, LD);
    end
    Gdc    = codistributed.build(Gdc,codistrGdc,'noCommunication');
    Gplus1 = codistributed.build(Gplus1,codistrGplus1,'noCommunication');
    Gdcr   = codistributed.build(Gdcr,codistrGdcr,'noCommunication');
    Gplusr = codistributed.build(Gplusr,codistrGplusr,'noCommunication');
    U0dc   = codistributed.build(U0dc,codistrU0dc,'noCommunication');
end

end


function [Gdc, Gplus1, U0dc, Gdcr, Gplusr] = Greenfnc_mat(DSrc, DRec, DL, DB, DR, DT, nx, nz, nb, ind, model, LD)

testGdc    = @(test)[test(2,:)+test(nb,:);test(3:nz-2,:);test(nz-1,:)+test(nz+1,:);...
              test(nz+2:(nz+nx)-3,:);test((nz+nx)-2,:)+test(nz+nx,:);test(nz+nx+1:2*nz+nx-4,:);...
              test(2*nz+nx-3,:)+test(2*nz+nx-1,:);test(2*nz+nx:nb-1,:)];

testGplus1 = @(test)[zeros(1,size(test,2)); test(1:nz-2,:); zeros(1,size(test,2));test(nz-2,:);...
             test(nz-1:nz+nx-5,:);zeros(1,size(test,2));test(nz+nx-5,:);...
             test(nz+nx-4:2*nz+nx-8,:);zeros(1,size(test,2));test(2*nz+nx-8,:);test(2*nz+nx-7:end,:);test(1,:)];

% form the right hand side of equation 10 for forward modelling
DSrc      = reshape(DSrc,model.n(1),model.n(2),size(DSrc,2));
U0dc      = reshape(DSrc(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nz*nx,size(DSrc,3));
U0dc      = U0dc(ind,:);
U0dc      = U0dc(1:nb,:);


% form the Gdc and Gdc_plus_1 matrix
% left boundary
DL        = reshape(DL,model.n(1),model.n(2),size(DL,2));
DL        = reshape(DL(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nx*nz,size(DL,3));
DL        = DL(ind,:);
test      = DL(1:nb,:);
Gdc       = testGdc(test);
Gdc       = transp(Gdc);
test      = DL(nb+1:2*nb-8,:);
Gplus1    = testGplus1(test);
Gplus1    = transp(Gplus1); 

% bottom boundary
DB        = reshape(DB,model.n(1),model.n(2),size(DB,2));
DB        = reshape(DB(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nx*nz,size(DB,3));
DB        = DB(ind,:);
test      = DB(1:nb,:);
test      = testGdc(test);
Gdc       = [Gdc;transp(test)];
test      = DB(nb+1:2*nb-8,:);
test      = testGplus1(test);
Gplus1    = [Gplus1;transp(test)]; 

% right boundary
DR        = reshape(DR,model.n(1),model.n(2),size(DR,2));
DR        = reshape(DR(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nx*nz,size(DR,3));
DR        = DR(:,end:-1:1);
DR        = DR(ind,:);
test      = DR(1:nb,:);
test      = testGdc(test);
Gdc       = [Gdc;transp(test)];
test      = DR(nb+1:2*nb-8,:);
test      = testGplus1(test);
Gplus1    = [Gplus1;transp(test)];

% top boundary
DT        = reshape(DT,model.n(1),model.n(2),size(DT,2));
DT        = reshape(DT(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nx*nz,size(DT,3));
DT        = DT(:,end:-1:1);
DT        = DT(ind,:);
test      = DT(1:nb,:);
test      = testGdc(test);
Gdc       = [Gdc;transp(test)];
test      = DT(nb+1:2*nb-8,:);
test      = testGplus1(test);
Gplus1    = [Gplus1;transp(test)];

% receiver side Green's function
DRec      = reshape(DRec,model.n(1),model.n(2),size(DRec,2));
DRec      = reshape(DRec(LD(1,1):LD(1,2),LD(2,1):LD(2,2),:),nx*nz,size(DRec,3));
DRec      = DRec(ind,:);
test      = DRec(1:nb,:);
Gdcr      = testGdc(test);
Gdcr      = transp(Gdcr);
test      = DRec(nb+1:2*nb-8,:);
Gplusr    = testGplus1(test);
Gplusr    = transp(Gplusr);
end
