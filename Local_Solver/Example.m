% standard FWI example
clear all;clc
%if sign(parpool_size)==0
 %   torque.parpool(120,1,5);
%end
%addpath(genpath('minConf'));
load marm.mat
v                = v(1:150,1:400);
v                = [1500*ones(10,size(v,2));v];
vbase            = v;
ind              = find(v==1028);
vbase(ind)       = 1836;
v                = vbase(60:100,100:250);
pert             = 200;
v(14,91:130)     = v(14,91:130) - pert;
v(15,77:126)     = v(15,77:126) - pert;
v(16,65:114)     = v(16,65:114) - pert;
v(17,52:101)      = v(17,52:101) - pert;
v(18,40:89)      = v(18,40:89) - pert;
v(19,23:82)      = v(19,23:82) - pert;
vfin             = vbase;
vfin(60:100,100:250) = v;
v = vfin;
clear vfin;
%% define model paramater
n                = size(v);
model.n          = [n 1]; % dimension of model in physical domain
model.o          = [0 0 0]; % origin
model.d          = [10 10 1]; % grid spacing
model.xt         = 0:model.d(2):(model.n(2)-1)*model.d(2); % spatial physical grid
model.zt         = 0:model.d(1):(model.n(1)-1)*model.d(1); % depth grid
model.nb         = [60 60 0]; % PML boundary points
freq             = [(3:0.25:4)' (4:0.25:5)' (5:0.25:6)' (6:0.25:7)' (7:0.25:8)' (8:0.25:9)' (9:0.25:10)'...
                    (10:0.25:11)' (11:0.25:12)' (12:0.25:13)' (13:0.25:14)' (14:0.25:15)' (15:0.25:16)' (16:0.25:17)',...
                    (17:0.25:18)' (18:0.25:19)' (19:0.25:20)' (20:0.25:21)' (21:0.25:22)' (22:0.25:23)' (23:0.25:24)' (24:0.25:25)'];
model.f0         = 8; %peak freq of ricker wavelet
model.t0         = 0; %phase shift of wavelet in seconds
model.zsrc       = model.d(2); % source depth
model.xsrc       = model.xt(1:end); % source position
model.zrec       = model.d(2); % receiver depth
model.xrec       = model.xt(1:end); % reciever position
model.ns         = length(model.xsrc);
model.nr         = length(model.xrec);
model.mtrue      = 1e6./v(:).^2; % slowness-square model
model.mbase      = 1e6./vbase(:).^2; % slowness-square model
load update_baseline1_21.mat
nsrc             = length(model.xsrc);
% define local dimension
LD               = [60 100;100 250];
model.mode       = 2; % 5 point stencil
model.nsurf      = length(model.xsrc);
model.nrec       = length(model.xrec);
model.zind       = 2;
model.rind       = 1:length(model.xt);
mpt              = reshape(m0,model.n(1),model.n(2));
mpt              = mpt(LD(1,1):LD(1,2),LD(2,1):LD(2,2));
mpt              = vec(mpt);
model.vmin       = min(vec(v(LD(1,1):LD(1,2),LD(2,1):LD(2,2))));
model.vmax       = max(vec(v(LD(1,1):LD(1,2),LD(2,1):LD(2,2))));
model.mmin       = 1e6./model.vmax.^2;
model.mmax       = 1e6./model.vmin.^2;
model.collocated = 1;
mback            = mpt;
nz               = LD(1,2) - LD(1,1)+1;
nx               = LD(2,2) - LD(2,1)+1;
mback            = reshape(mback,nz,nx);
model.testGreen  = 1;
%% run inversion
for i = 8:size(freq,2)
    % minimize over m
    fprintf('\n **** SOLVING FOR m : freq batch %d ****** \n',i);
    model.freq  = freq(:,i);

    %initialize stuff
    Q         = -1*speye(model.ns); %source weight (simultaneous sources)

    % generate Fully seismic data for inversion
    D         = F5pt(model.mtrue,Q,model);
    fprintf('\n **** data generated for freq batch %d ****** \n',i);

    % generate background Green's function
    [Gdc, Gplus1, U0dc, Gdcr, Gplusr, Uback] = Back_GreenFcnL( m0, model, LD );

    % run inversion
    fh = @(x)misfitlocalL(x,D, model, LD,Gplus1,Gdc,U0dc,Gplusr,Gdcr,Uback);

    % minimize over model for fixed data
    options.maxIter  = 5;
    boptions.LB      = model.mmin*ones(size(mpt));
    boptions.UB      = model.mmax*ones(size(mpt));
    boptions.l       = model.mmin;
    boptions.u       = model.mmax;
    funProj          = @(x) boundProjectstan(x,boptions);  % incorporates water velocity constraints
    mpt              = minConf_TMP(fh,mpt,boptions.LB,boptions.UB,options);
    mpt              = reshape(mpt,nz,nx);
    mpt(1,:)         = mback(1,:);
    mpt(:,1)         = mback(:,1);
    mpt(end,:)       = mback(end,:);
    mpt(:,end)       = mback(:,end);
    mpt              = vec(mpt);

    figure(1);imagesc(reshape(mpt,nz,nx));colorbar;caxis([0.25 0.44]);drawnow;
    save(['Local_FullPDE_FWI_iter5_' num2str(i)] , 'mpt');
end
