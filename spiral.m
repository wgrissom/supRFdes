function [rf,mse,sar,m] = spiral(maps,algp)

addpath util/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genfigs = algp.genfigs;

[dimb1(1),dimb1(2),Nc] = size(maps.b1);

% normalize b1 by median so that regularization params can be
% meaningfully reported
b1Scale = 1/median(abs(maps.b1(repmat(maps.mask,[1 1 Nc]))));
maps.b1 = maps.b1*b1Scale;
% apply mask to b1 maps, if not already done
maps.b1 = maps.b1.*repmat(maps.mask,[1 1 Nc]);

%% Get the spiral trajectory
dt = algp.dt;
deltax = algp.spres; % spatial resolution of trajectory
FOVsubsamp = algp.spfov; % cm, xfov of spiral
forwardspiral = 0;
dorevspiralramp = 0;
gmax = algp.gmax;
dgdtmax = algp.dgdtmax;
traj = 'spiral';
get_traj;
k = k(NN(2)+1:end,:);
Nt = size(k,1);

% design grid dimension
dim = ceil(maps.deltax*dimb1/algp.spres*algp.overSamp);
        
switch algp.dPatternSelect
    case 'unif'
        d = ones(dim);
    case 'square'
        [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
        d = double(abs(xi) <= algp.squareWidth/2/(maps.deltax/algp.overSamp) & abs(yi) <= algp.squareWidth/2/(maps.deltax/algp.overSamp));
end

% interpolate the B1 maps to the same grid as the target pattern
[xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
[xb1,yb1] = ndgrid((-dimb1(1)/2:dimb1(1)/2-1)/dimb1(1),(-dimb1(2)/2:dimb1(2)/2-1)/dimb1(2));
b1Int = zeros([dim Nc]);
for ii = 1:Nc
    b1Int(:,:,ii) = interp2(yb1,xb1,maps.b1(:,:,ii),yi/dim(2),xi/dim(1),'linear',0);
end
b1 = b1Int;
maskInt = interp2(yb1,xb1,double(maps.mask),yi/dim(2),xi/dim(1),'linear',0);
mask = maskInt > 0.9; % get new mask for higher res b1

% remove the first channels phase from the rest
b1 = b1.*repmat(exp(-1i*angle(b1(:,:,1))),[1 1 Nc]);

% find the middle point of the mask
xmi = round(mean(xi(logical(mask(:)))));
ymi = round(mean(yi(logical(mask(:)))));

% center the target pattern in the mask
d = circshift(d,[xmi ymi]);

% smooth it a bit
d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));

% Set target phase equal to phase of B1 sum
d = d.*exp(1i*angle(sum(b1,3)));


%% Set up system matrix and penalty matrix

% columnize and mask the b1 maps
b1 = permute(b1,[3 1 2]);b1 = b1(:,:).';
b1 = b1(mask,:);

%% Build system matrix
A = Gmri_SENSE(k,mask,'fov',maps.deltax*dimb1,'sens',conj(b1))';

%% Solve
if ~isfield(algp,'SARlimits')
    
    disp 'Building SAR reg matrix'
    R = sqrt(algp.beta*sum(mask(:))/size(maps.Sv,3))*buildSARReg(sum(maps.Sv,3),Nt);
    %R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
    %    sqrt(roughbeta)*speye(Nt);
    %Rfull = kron(speye(Nc),R);
    
    disp 'Starting CG'
    ncgiters = floor(length(k)/8);       % # CG Iterations per RF update
    rf = qpwls_pcg(zeros(length(k)*Nc,1),A,1,d(mask),0,R,1,ncgiters,mask); % CG
    rf = reshape(rf(:,end),[length(k) Nc]); % qpwls_pcg returns all iterates
    
    % calculate excitation pattern
    m = zeros(dim);m(mask) = A*rf(:);
    
    mse = norm(mask.*(m-d))^2;
    sar = norm(R*rf(:))^2;
    cost = mse + sar;
    fprintf('Final MSE: %0.2f. SAR: %0.2f.\n',mse,sar);

else %% problem is SAR-constrained; use fmincon

    % 2016b syntax
    %options = optimoptions(@fmincon,'Algorithm','interior-point',...
    %    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
    % 2015a syntax  
    R = 0; % no regularization here; could add roughness or integrated power later
    HessFcn = @(x,lambda,v)HessMultFcn(x,lambda,v,A,R,maps.Sv);
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'GradConstr','on','GradObj','on','Hessian','user-supplied',...
        'SubproblemAlgorithm','cg','HessMult',HessFcn,'Display','iter-detailed');
    
    fun = @(x)quadobj(x,A,d(mask),R);
    nonlconstr = @(x)quadconstr(x,maps.Sv,algp.SARlimits*algp.SARTR/dt);
    x0 = zeros(2*Nt*Nc,1); % column vector
    [x,fval,eflag,output,lambda] = fmincon(fun,x0,...
        [],[],[],[],[],[],nonlconstr,options);
    
    rf = reshape(x,[length(k) Nc]); 
    
    % calculate excitation pattern
    m = zeros(dim);m(mask) = A*rf(:);
    
    mse = norm(mask.*(m-d))^2;cost = mse;
    maxSAR = 0;
    for ii = 1:size(maps.Sv,3)
        SARtmp = 0;
        for jj = 1:Nt
            SARtmp = SARtmp + real(conj(rf(jj,:))*(maps.Sv(:,:,ii)*rf(jj,:).'));
        end
        if SARtmp > maxSAR
            maxSAR = SARtmp;
        end
    end
    fprintf('Final MSE: %0.2f%. Max SAR: %f.\n',mse,maxSAR);
    
end


%% Display final results

if genfigs
    maxamp = max(abs([d(:);m(:)]));
    figure
    subplot(221)
    imagesc(abs(d),[0 maxamp]);axis image;colorbar
    title 'Desired pattern'
    subplot(222)
    imagesc(abs(m),[0 maxamp]);axis image;colorbar
    title(sprintf('Final pattern',Nc));
    subplot(224)
    imagesc(abs(m-d));axis image;colorbar
    title(sprintf('Cost = %0.2f%%',cost));
end

%% Helper functions
function [y,grady] = quadobj(x,A,d,R)

x = x(1:length(x)/2)+1i*x(length(x)/2+1:end);
Ax = A*x;
res = Ax - d;
Rx = R*x;
y = 1/2*real(res'*res) + 1/2*real(Rx'*Rx);
if nargout > 1
    grady = A'*res + R'*Rx;
    grady = [real(grady);imag(grady)];
end

function [y,yeq,grady,gradyeq] = quadconstr(x,Sv,SARlimits)

x = x(1:length(x)/2)+1i*x(length(x)/2+1:end); % restore complex
[Nc,~,Nvop] = size(Sv);
Nt = length(x)/Nc;
x = reshape(x,[Nt Nc]).'; % reshape to Nt x Nc; transpose to get Nc x Nt

% precompute Sv*x
Svx = zeros([Nc Nt Nvop]);
for ii = 1:Nvop % loop over constraints
    Svx(:,:,ii) = Sv(:,:,ii)*x;
end

y = zeros(1,Nvop);
for ii = 1:Nvop % loop over constraints
    for jj = 1:Nt % sum over time; duty cycle is absorbed in SARlimits
        y(ii) = y(ii) + x(:,jj)'*Svx(:,jj,ii);
    end
    y(ii) = real(y(ii)) - SARlimits(ii);
end
yeq = []; % no equality constraints
    
% should be able to reduce computation by calculating grady first and 
% using that to get y!
if nargout > 2
    % collapse time and coil dims; stack real and imaginary
    grady = permute(2*Svx,[3 2 1]); % Nvop x Nt x Nc
    grady = grady(:,:).'; grady = [real(grady);imag(grady)];
end
gradyeq = []; % no equality constraint derivatives

function W = HessMultFcn(x,lambda,v,A,R,Sv)

x = v(1:length(x)/2)+1i*v(length(x)/2+1:end); % restore complex
[Nc,~,Nvop] = size(Sv);
Nt = length(x)/Nc;

if max(abs(x)) > 0
    
    % objective second derivatives
    W = A'*(A*x) + R'*(R*x);
    
    % constraint second derivatives
    x = reshape(x,[Nt Nc]).'; % reshape to Nc x Nt
    grady = zeros([Nc Nt]);
    for ii = 1:Nvop % loop over constraints
        grady = grady + lambda.ineqnonlin(ii)*2*(Sv(:,:,ii)*x);
    end
    W = W + col(grady.');W = [real(W);imag(W)];

else    
    W = zeros(size(v));
end


function out = col(in)

out = in(:);



