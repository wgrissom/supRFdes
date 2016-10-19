function [goodSols,goodFeatures,percErrs,bestSols,bestFeatures,bestPercErrs,totalPerSlice,maps] = loadData(indices,nHood)

cases = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10'}; %'headb1';
orientations = {'axial'}; % list of orientations to fit to, for each b1case

% get the set of 'acceptable' shims; associate with their slice's features
%nHood = 3; % B1+ DFT neighborhood
Nc = 36; % # tx coils
threshPerc = 5;
corrThresh = 0.999;

maps = {};

disp 'Loading All Solutions and Features...'
goodSols = []; goodFeatures = [];percErrs = [];totalPerSlice = [];
bestSols = []; bestFeatures = []; bestPercErrs = [];
for b1i = 1:length(cases)
    for ori = 1:length(orientations)
        % load this b1 case
        printf('Loading case %s...',[cases{b1i} '_shims_' orientations{ori}]);
        x = load([cases{b1i} '_shims_' orientations{ori}]);
        maps{b1i} = x.maps;
        for ii = 1:x.Nsl
            
            if any(indices == (b1i-1)*x.Nsl + ii) % include this slice if its on the list
                
                % extract this slice's features
                % include central DFT coeffs for each coil, normalized to coil 1's DC
                % coefficient
                sliceFeatures = ii; % store slice position
                %sliceFeatures = x.notEmptyInds(ii)*([streq(orientations{ori},'axial');...
                %    streq(orientations{ori},'sagittal');...
                %    streq(orientations{ori},'coronal')]);
                % calculate the 'width' of the mask in x and y
                [xi,yi] = meshgrid(1:x.dimxy(1),1:x.dimxy(2));
                sliceFeatures = [sliceFeatures;std(xi(x.maps.mask(:,:,ii)));std(yi(x.maps.mask(:,:,ii)))];
                % calculate centroid of mask
                sliceFeatures = [sliceFeatures;mean(xi(x.maps.mask(:,:,ii)));mean(yi(x.maps.mask(:,:,ii)))];
                if nHood > 0
                    for jj = 1:Nc
                        tmp = x.maps.b1(:,:,ii,jj);
                        tmp = fftshift(fft2(fftshift(tmp)));
                        tmp = col(tmp(ceil(x.dimxy(1)/2)+1-(nHood-1)/2:ceil(x.dimxy(1)/2)+1+(nHood-1)/2,...
                            ceil(x.dimxy(2)/2)+1-(nHood-1)/2:ceil(x.dimxy(2)/2)+1+(nHood-1)/2));
                        tmp = [real(tmp);imag(tmp)];
                        %tmp = tmp./sum(col(x.maps.b1(:,:,ii,1))); % normalize to coil 1's DC
                        sliceFeatures = [sliceFeatures; tmp];
                    end
                end
                sliceFeatures = [sliceFeatures; sign(sliceFeatures(end-2*Nc:end))];
                %sliceFeatures = [sliceFeatures;sliceFeatures.^2;sliceFeatures.^3];
                
                goodSolst = []; % solutions that may be added
                goodFeaturest = []; % features that may be added
                percErrst = []; % keep errors of solutions that may be added, for pruning
                for jj = 1:x.algp.nRandStart + 2
                    if x.percErr(jj,ii) < threshPerc
                        goodSolst = [goodSolst x.rf(:,jj,ii)]; % assume coil 1 is just 1; so we are looking at relative amplitude and phase
                        goodFeaturest = [goodFeaturest sliceFeatures];
                        percErrst = [percErrst x.percErr(jj,ii)];
                    end
                end
                
                [~,bestInd] = min(x.percErr(:,ii));
                bestSols = [bestSols x.rf(:,bestInd,ii)];
                bestFeatures = [bestFeatures sliceFeatures];
                bestPercErrs = [bestPercErrs x.percErr(bestInd,ii)];
                
                % prune redundant solutions
                remaining = 1:size(goodSolst,2); % list of remaining solutions to group
                %keepInds = remaining; % for comparison
                keepInds = [];
                % normalize so that max(corrs) = 1
                goodSolsNorm = goodSolst./repmat(sqrt(sum(abs(goodSolst).^2,1)),[Nc 1]);
                while ~isempty(remaining)
                    % find lowest error among remaining
                    [~,minInd] = min(percErrst(remaining));
                    % add this guy to the list to keep
                    keepInds = [keepInds remaining(minInd)];
                    % find all solutions close to this guy and remove them from
                    % the list
                    corrs = real(goodSolsNorm(:,remaining(minInd))'*goodSolsNorm(:,remaining))./...
                        real(goodSolsNorm(:,remaining(minInd))'*goodSolsNorm(:,remaining(minInd)));
                    remaining = remaining(corrs < corrThresh);
                end
                goodSolst = goodSolst(:,keepInds);
                goodFeaturest = goodFeaturest(:,keepInds);
                percErrst = percErrst(keepInds);
                
                
                %             % prune again, after phase rotating
                %             remaining = 1:size(goodSolst,2); % list of remaining solutions to group
                %             %keepInds = remaining; % for comparison
                %             keepInds = [];
                %             % normalize so that max(corrs) = 1
                %             goodSolsNorm = goodSolst./repmat(sqrt(sum(abs(goodSolst).^2,1)),[Nc 1]);
                %             while ~isempty(remaining)
                %                 % find lowest error among remaining
                %                 [~,minInd] = min(percErrst(remaining));
                %                 % add this guy to the list to keep
                %                 keepInds = [keepInds remaining(minInd)];
                %                 % find all solutions close to this guy and remove them from
                %                 % the list
                %                 corrs = Inf+zeros(length(remaining),1);
                %                 for jj = 1:length(remaining)
                %                     phs = 0:2*pi/1000:2*pi;
                %                     for kk = 1:length(phs)
                %                     corrs(jj) = min(corrs(jj),exp(1i*phs(kk))*goodSolsNorm(
                %                 end
                %                 corrs = abs(goodSolsNorm(:,remaining(minInd))'*goodSolsNorm(:,remaining))./...
                %                     abs(goodSolsNorm(:,remaining(minInd))'*goodSolsNorm(:,remaining(minInd)));
                %                 remaining = remaining(corrs < corrThresh);
                %             end
                %             goodSolst = goodSolst(:,keepInds);
                %             goodFeaturest = goodFeaturest(:,keepInds);
                %             percErrst = percErrst(keepInds);
                
                
                % append the final set of features and solutions to the running
                % list
                goodSols = [goodSols goodSolst];
                goodFeatures = [goodFeatures goodFeaturest];
                percErrs = [percErrs percErrst];
                
                % Store how many acceptable solutions we got in each slice of this
                % dataset
                totalPerSlice = [totalPerSlice length(keepInds)];
                
            end
        end
        
    end
end
Nsl = sum(totalPerSlice > 0);
totalPerSlice = totalPerSlice(totalPerSlice > 0);
startinds = [1 cumsum(totalPerSlice(1:end-1))+1]; % starting index of each slice's solutions
endinds = cumsum(totalPerSlice); % ending index of each slice's solutions
nTraining = size(goodSols,2);
nFeatures = size(goodFeatures,1);
printf('Loaded a total of %d solutions and %d features\n',nTraining,nFeatures);


