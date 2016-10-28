function [rf,cost,m] = spiral(maps,algp)

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

%if strcmp(algp.dPatternSelect,'unif')
    d = d.*exp(1i*angle(sum(b1,3)));
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up system matrix and penalty matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build RF waveform roughness penalty matrix for full
% and compressed arrays

disp 'Building SAR reg matrix'
R = sqrt(algp.beta*sum(mask(:))/size(maps.Sv,3))*buildSARReg(maps.Sv,Nt);
%R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
%    sqrt(roughbeta)*speye(Nt);
%Rfull = kron(speye(Nc),R);

% columnize and mask the b1 maps
b1 = permute(b1,[3 1 2]);b1 = b1(:,:).';
b1 = b1(mask,:);

% Build system matrix
A = Gmri_SENSE(k,mask,'fov',maps.deltax*dimb1,'sens',conj(b1))';

disp 'Starting CG'
ncgiters = floor(length(k)/8);       % # CG Iterations per RF update
rf = qpwls_pcg(zeros(length(k)*Nc,1),A,1,d(mask),0,R,1,ncgiters,mask); % CG
rf = reshape(rf(:,end),[length(k) Nc]); % qpwls_pcg returns all iterates

% calculate excitation pattern
m = zeros(dim);m(mask) = A*rf(:);

mse = norm(mask.*(m-d))^2;
sar = norm(R*rf(:))^2;
cost = mse + sar;
fprintf('Final MSE: %0.2f%. SAR: %f.\n',mse,sar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display final results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
