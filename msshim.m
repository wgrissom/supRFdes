% Script to do multi-slice RF shimming
% Copyright Zhipeng Cao and Will Grissom, Vanderbilt University, 2015.

Nc = 36;
subj = {'Duke','Ella'};
b1Path = '';
Nsubj = 5;

%b1case = 'v1_36chv3'; %'headb1';
mode = 'all'; % 'CP','all'
maps.b1 = [];
maps.mask = [];
maps.R = {};
for subjI = 10%:length(subj)*Nsubj
    % these cases model a 3-row x 12 element array, arranged
    % circumferentially in B1p_out, starting with the bottom row and
    % moving up. The CP mode is obtained by multiplying each B1 map
    % with exp(-1i*[0:35]/12*2*pi);
    %loadString = [b1Path subj{floor((subjI-1)/Nsubj)+1} '/' 'v' num2str(subjI-floor((subjI-1)/Nsubj)*Nsubj) '_' num2str(Nc) 'ch'];
    loadString = [b1Path subj{floor((subjI-1)/Nsubj)+1} 'S' num2str(subjI-floor((subjI-1)/Nsubj)*Nsubj) '_' num2str(Nc) 'ch'];
    fprintf('Loading %s\n',loadString);
    load(loadString);
    switch mode
        case 'CP'
            tmp = zeros(size(B1p_out)); tmp = tmp(:,:,:,1:2);
            for ii = 1:size(B1p_out,4) % loop over coils
                tmp(:,:,:,1) = tmp(:,:,:,1) + B1p_out(:,:,:,ii)*exp(-1i*(ii-1)/12*2*pi); % CP mode
                tmp(:,:,:,2) = tmp(:,:,:,2) + B1p_out(:,:,:,ii)*exp(-1i*2*(ii-1)/12*2*pi); % anti-CP mode
            end
            B1p_out = tmp;
        case 'all'
            % apply CP phases but do not sum
            for ii = 1:size(B1p_out,4) % loop over coils
                B1p_out(:,:,:,ii) = B1p_out(:,:,:,ii)*exp(-1i*(ii-1)/12*2*pi); % CP mode
            end
        case 'volume'
            tmp = permute(B1p_out,[4 1 2 3]);
            tmp = tmp(:,:,:); % collapse slice and y dims
            B1p_out = reshape(permute(tmp,[2 3 1]),[size(B1p_out,1) size(B1p_out,2)*size(B1p_out,3) 1 size(B1p_out,4)]);
            mask_out = mask_out(:,:);
        case 'opt'
            % load the best volumetric shim
        case 'splitCP'
            tmp = zeros(size(B1p_out)); tmp = tmp(:,:,:,1:2);
            for ii = 1:size(B1p_out,4)/2 % loop over coils
                tmp(:,:,:,1) = tmp(:,:,:,1) + B1p_out(:,:,:,2*ii-1)*exp(-1i*((2*ii-1)-1)/12*2*pi); % CP mode
                tmp(:,:,:,2) = tmp(:,:,:,2) + B1p_out(:,:,:,2*ii)*exp(-1i*(2*ii-1)/12*2*pi); % anti-CP mode
            end
            B1p_out = tmp;
            
    end
    maps.b1 = cat(3,maps.b1,B1p_out*10^6); % bring up to uT
    maps.mask = logical(cat(3,maps.mask,maskAll_out));
    loadStringVOP = [loadString '_10gVOP'];
    load(loadStringVOP);
    % get the average SAR matrix, and append it to the rest
    Sv(:,:,end+1) = Sv(:,:,1)*nSv(1);
    for ii = 2:size(Sv,3)-1
        Sv(:,:,end) = Sv(:,:,end) + Sv(:,:,ii)*nSv(ii);
    end
    Sv(:,:,end) = Sv(:,:,end)/sum(nSv);
    % scale it up 5x since the max local can be e.g. 20 W/kg but the max ave can only be 4 W/kg
    Sv(:,:,end) = 5*Sv(:,:,end);
    R = {}; R{1} = buildSARReg(sum(Sv,3),1);
    maps.R = cat(1,maps.R,R(ones(1,size(B1p_out,3)))); % stick these matrices on the end
end

load('../v1_36chv3_shims_axial','randAmpPhs');

orientation = 'axial'; % 'axial', 'sagittal', 'coronal'
switch orientation
    case 'axial'
        % prune off empty slices
        notEmpty = squeeze(sum(sum(maps.mask,1),2)) > 0;
        notEmptyInds = find(notEmpty);
        maps.mask = maps.mask(:,:,notEmpty);
        maps.b1 = maps.b1(:,:,notEmpty,:);
    case 'sagittal'
        maps.b1 = permute(maps.b1,[3 2 1 4]);
        maps.mask = permute(maps.mask,[3 2 1]);
        % prune off empty slices
        notEmpty = squeeze(sum(sum(maps.mask,1),2)) > 0;
        notEmptyInds = find(notEmpty);
        maps.mask = maps.mask(:,:,notEmpty);
        maps.b1 = maps.b1(:,:,notEmpty,:);
    case 'coronal'
        maps.b1 = permute(maps.b1,[1 3 2 4]);
        maps.mask = permute(maps.mask,[1 3 2]);
        % prune off empty slices
        notEmpty = squeeze(sum(sum(maps.mask,1),2)) > 0;
        notEmptyInds = find(notEmpty);
        maps.mask = maps.mask(:,:,notEmpty);
        maps.b1 = maps.b1(:,:,notEmpty,:);
end

[dimxy(1),dimxy(2),Nsl,Nc] = size(maps.b1); % number of physical coils (Ncoils)

% normalize b1 by median so that regularization params can be
% meaningfully reported
maps.b1 = maps.b1./median(abs(maps.b1(repmat(maps.mask,[1 1 1 Nc]))));

algp.beta = 10^-2;
%algp.tol = 1-0.999;
algp.tol = 0.99999;
algp.nRandStart = 100;
algp.randAmpPhs = randAmpPhs;

[rf,costAll,mAll,randAmpPhs] = msShim_randStart(maps,algp);

% calculate percent deviation from unit amplitude
for ii = 1:Nsl
    for jj = 1:algp.nRandStart + 2
        tmp = mAll(:,:,jj,ii);
        percErr(jj,ii) = 100*sqrt(mean(abs(abs(tmp(maps.mask(:,:,ii)))-1).^2));
    end
end

% for each rf solution, normalize to first coil
rfNorm = rf./repmat(rf(1,:,:),[Nc 1 1]);

% examine "valid" solutions more closely
threshPerc = 5; % percent error threshold

% examine "best" rf pulses for each slice,
% to see how similar they are over slices
[~,inds] = min(percErr,[],1);
for ii = 1:Nsl
    rfBest(:,ii) = rfNorm(:,inds(ii),ii);
    mBest(:,:,ii) = mAll(:,:,inds(ii),ii);
end

save(['subj' num2str(subjI) '_shims_' orientation]);
