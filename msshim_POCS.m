% Script to do multi-slice RF shimming
% Copyright Zhipeng Cao and Will Grissom, Vanderbilt University, 2015.

% TO IMPLEMENT:
% 1) Load slice/subject indices for train cases as in compare_kNN
% 2) Get features for train cases
% 3) Normalize features for train cases; add some feature transformations
% 4) Build Feature projector for train cases, pass to POCS ms-shim

if ~exist('mapsTrain','var')
    Nc = 36;
    subj = {'Duke'};%{'Duke','Ella'};
    b1Path = '';
    Nsubj = 1;
    
    %b1case = 'v1_36chv3'; %'headb1';
    mode = 'all'; % 'CP','all'
    mapsTrain.b1 = [];
    mapsTrain.mask = [];
    mapsTrain.R = {};
    for subjI = 1:length(subj)*Nsubj
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
        mapsTrain.b1 = cat(3,mapsTrain.b1,B1p_out*10^6); % bring up to uT
        mapsTrain.mask = logical(cat(3,mapsTrain.mask,maskAll_out)); % maskAll_out = water + fat, no bone; maskROI_out = water, no fat, no bone
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
        mapsTrain.R = cat(1,mapsTrain.R,R(ones(1,size(B1p_out,3)))); % stick these matrices on the end 
%         Rfoo = buildSARReg(sum(Sv,3),1);
%         CtCfoo2 = Rfoo'*Rfoo;
%         CtCfoo = R{1}'*R{1};
%         xfoo = randn(Nc,100)+1i*randn(Nc,100);
%         for ii = 1:100
%             foo1 = xfoo(:,ii)'*(CtCfoo*xfoo(:,ii));
%             foo2 = xfoo(:,ii)'*(Rfoo'*(Rfoo*xfoo(:,ii)));
%             foo3 = xfoo(:,ii)'*(CtCfoo2*xfoo(:,ii));
%             norm(foo1-foo2)
%             norm(foo1-foo3)
%         end
    end
    
    [dimxy(1),dimxy(2),Nsl,Nc] = size(mapsTrain.b1); % number of physical coils (Ncoils)

    orientation = 'axial'; % 'axial', 'sagittal', 'coronal'
    switch orientation
        case 'axial'
            % prune off empty slices
            notEmpty = squeeze(sum(sum(mapsTrain.mask,1),2)) > 0;
            notEmptyInds = find(notEmpty);
            mapsTrain.mask = mapsTrain.mask(:,:,notEmpty);
            mapsTrain.b1 = mapsTrain.b1(:,:,notEmpty,:);
        case 'sagittal'
            mapsTrain.b1 = permute(mapsTrain.b1,[3 2 1 4]);
            mapsTrain.mask = permute(mapsTrain.mask,[3 2 1]);
            % prune off empty slices
            notEmpty = squeeze(sum(sum(mapsTrain.mask,1),2)) > 0;
            notEmptyInds = find(notEmpty);
            mapsTrain.mask = mapsTrain.mask(:,:,notEmpty);
            mapsTrain.b1 = mapsTrain.b1(:,:,notEmpty,:);
        case 'coronal'
            mapsTrain.b1 = permute(mapsTrain.b1,[1 3 2 4]);
            mapsTrain.mask = permute(mapsTrain.mask,[1 3 2]);
            % prune off empty slices
            notEmpty = squeeze(sum(sum(mapsTrain.mask,1),2)) > 0;
            notEmptyInds = find(notEmpty);
            mapsTrain.mask = mapsTrain.mask(:,:,notEmpty);
            mapsTrain.b1 = mapsTrain.b1(:,:,notEmpty,:);
    end
    
end

load('../v1_36chv3_shims_axial','randAmpPhs');

% normalize b1 by median so that regularization params can be
% meaningfully reported
mapsTrain.b1 = mapsTrain.b1./median(abs(mapsTrain.b1(repmat(mapsTrain.mask,[1 1 1 Nc]))));

%algp.beta = 0;%10^0; % RF power regularization
algp.noCG = false;
algp.step = 0.1;
algp.beta = 10^-2;%10^0; % SAR regularization
%algp.tol = 1-0.999;
algp.tol = 1-0.999; % stopping parameter
algp.nRandStart = -1; % only do CP initialization
algp.ncgiters = 3; % number of cg iters for each RF update
algp.dofigs = 0; % whether to show pattern figure every 100 iters
algp.randAmpPhs = randAmpPhs; % random starting phases and amplitudes
%algp.Nccomp = 9; % number of compressed channels to threshold down to
algp.holeThresh = 0.25; % hole detection threshold
algp.maxIters = 5000; % max iters before exiting
if exist('Fproj','var') 
    algp.Fproj = Fproj; % projector onto space of predictable solutions
end

[rfPOCS,errAllPOCS,powAllPOCS,mAllPOCS,randAmpPhsPOCS] = msShim_randStart_POCS_vopReg(mapsTrain,algp);

% % calculate percent deviation from unit amplitude
% for ii = 1:Nsl
%     for jj = 1:algp.nRandStart + 2
%         tmp = mAll(:,:,jj,ii);
%         percErrPOCSshim(jj,ii) = 100*sqrt(mean(abs(abs(tmp(mapsTrain.mask(:,:,ii)))-1).^2));
%     end
% end
% 
% % for each rf solution, normalize to first coil
% rfNorm = rf./repmat(rf(1,:,:),[Nc 1 1]);
% 
% % examine "valid" solutions more closely
% threshPerc = 5; % percent error threshold
% 
% % examine "best" rf pulses for each slice,
% % to see how similar they are over slices
% [~,inds] = min(percErr,[],1);
% for ii = 1:Nsl;rfBest(:,ii) = rfNorm(:,inds(ii),ii);end
% 
% %save([b1case '_shims_' orientation '_mode_' mode]);
