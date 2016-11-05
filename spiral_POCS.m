function [rf,mse,sar,m] = spiral_POCS(maps,algp)

addpath util/

%% Problem parameters

genfigs = algp.genfigs;

[dimb1(1),dimb1(2),Nsl,Nc] = size(maps.b1);

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

%% design grid dimension
dim = ceil(maps.deltax*dimb1/algp.spres*algp.overSamp);

%% Target pattern definition
switch algp.dPatternSelect
    case 'unif'
        d = ones([dim Nsl]);
    case 'square'
        [xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
        d = double(abs(xi) <= algp.squareWidth/2/(maps.deltax/algp.overSamp) & abs(yi) <= algp.squareWidth/2/(maps.deltax/algp.overSamp));
        % smooth it a bit
        d = ift2((hamming(dim(1))*hamming(dim(2))').^2.*ft2(d));
        % replicate to all slices
        d = repmat(d,[1 1 Nsl]);
end

%% Process b1 maps and set up all slice-specific design variables
mask = false([dim Nsl]);
[xi,yi] = ndgrid(-dim(1)/2:dim(1)/2-1,-dim(2)/2:dim(2)/2-1);
[xb1,yb1] = ndgrid((-dimb1(1)/2:dimb1(1)/2-1)/dimb1(1),(-dimb1(2)/2:dimb1(2)/2-1)/dimb1(2));
A = {};
R = {};
for ii = 1:Nsl
    
    % normalize b1 by median so that regularization params can be
    % meaningfully reported
    
    tmp = sqz(maps.b1(:,:,ii,:));
    b1Scale = 1/median(abs(tmp(repmat(maps.mask(:,:,ii),[1 1 Nc]))));
    b1 = sqz(maps.b1(:,:,ii,:));
    b1 = b1*b1Scale;
    % apply mask to b1 maps, if not already done
    b1 = b1.*repmat(maps.mask(:,:,ii),[1 1 Nc]);
    
    % interpolate the B1 maps to the same grid as the target pattern
    b1Int = zeros([dim Nc]);
    for jj = 1:Nc
        b1Int(:,:,jj) = interp2(yb1,xb1,b1(:,:,jj),yi/dim(2),xi/dim(1),'linear',0);
    end
    maskInt = interp2(yb1,xb1,double(maps.mask(:,:,ii)),yi/dim(2),xi/dim(1),'linear',0);
    mask(:,:,ii) = maskInt > 0.9; % get new mask for higher res b1
    
    % remove the first channels phase from the rest
    b1Int = b1Int.*repmat(exp(-1i*angle(b1Int(:,:,1))),[1 1 Nc]);
    
    % find the middle point of the mask
    xmi = round(mean(xi(logical(col(mask(:,:,ii))))));
    ymi = round(mean(yi(logical(col(mask(:,:,ii))))));
    
    % center the target pattern in the slice's mask
    d(:,:,ii) = circshift(d(:,:,ii),[xmi ymi]);
    
    % Set target phase equal to phase of B1 sum
    d(:,:,ii) = d(:,:,ii).*exp(1i*angle(sum(b1Int,3)));
    
    %% Set up system matrix and penalty matrix
    
    % columnize and mask the b1 maps
    b1Int = permute(b1Int,[3 1 2]);b1Int = b1Int(:,:).';
    b1Int = b1Int(col(logical(mask(:,:,ii))),:);
    
    %% Build system matrix
    A{ii} = Gmri_SENSE(k,logical(mask(:,:,ii)),'fov',maps.deltax*dimb1,'sens',conj(b1Int))';
    
    %% Build SAR reg matrices
    R{ii} = sqrt(algp.beta*sum(sum(mask(:,:,ii)))/size(maps.Sv(:,:,:,ii),3))*buildSARReg(sum(maps.Sv(:,:,:,ii),3),Nt);
    %R = sqrt(roughbeta)*spdiags([-ones(Nt,1) ones(Nt,1)],[0 1],Nt,Nt) + ...
    %    sqrt(roughbeta)*speye(Nt);
    %Rfull = kron(speye(Nc),R);
end


%% outer loop: continue until shims stop changing
flag = 1; % to get loop started
itr = 1; % iteration counter
m = zeros([dim Nsl]);
rf = zeros([Nt Nc Nsl]);
while flag
    
    rfOld = rf;
    
    %% update pulses with a few CG iters
    disp 'Updating pulses'
    for ii = 1:Nsl
        
        dm = d(:,:,ii);
        xS = qpwls_pcg(col(rf(:,:,ii)),A{ii},1,dm(mask(:,:,ii)),0,R{ii},1,algp.ncgiters,mask(:,:,ii)); % CG
        rf(:,:,ii) = reshape(xS(:,end),[length(k) Nc]); % qpwls_pcg returns all iterates
        
        % calculate excitation pattern
        tmp = zeros(dim);tmp(mask(:,:,ii)) = A{ii}*col(rf(:,:,ii));
        m(:,:,ii) = tmp;
        
        mse(ii) = norm(mask(:,:,ii).*(m(:,:,ii)-d(:,:,ii)))^2;
        sar(ii) = norm(R{ii}*col(rf(:,:,ii)))^2;
        cost(ii) = mse(ii) + sar(ii);
        
    end
    
    %% project pulses onto space of predictable pulses
    disp 'Projecting pulses onto subspace of predictable pulses'
    if isfield(algp,'Fproj')
        
        rft = permute(rf,[3 1 2]);rft = rft(:,:); % collapse time + coil dims
        
        % project real and imag parts of each slice's shim for this
        % coil onto the space of predictable shims, assuming that
        % the prediction weights are given by fitting the
        % features to this training data across slices
        rft = algp.Fproj*rft;
        
        rft = reshape(rft,[size(rf,3) size(rf,1) size(rf,2)]);
        rf = permute(rft,[2 3 1]);
        
    end
    
    itr = itr + 1;
        
    % check stopping criterion
    if norm(rf(:)-rfOld(:))/norm(rfOld(:)) < algp.tol || itr > algp.maxIters
        flag = 0;
    end
        
    if rem(itr,100) == 0;
        fprintf('Iteration %d. Err all %f. SAR all %f.\n',itr,sum(mse),sum(sar));
    end
    
end

function out = col(in)

out = in(:);



