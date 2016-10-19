% TODO: Implement POCS in a coordinate-wise fashion; test regularization
% and different numbers of features


nHood = 3; % fourier neighborhood over which to extract B1+ coeff features
K = 10; % K-fold cross validation

Nsl = 31; % # slices in each subject's maps
Nsubj = 10; % # subjects
Nc = 36;

doPOCS = true;
POCSreg = 0.1; % regularization on POCS pseudoinverse 

lambdaKern = 10^0;

kNeighbors = 5; % # neighbors to interpolate over for kNN

getNewInds = false;
phsRotate = false;


if getNewInds
    kfoldInds = crossvalind('Kfold',Nsl*Nsubj,K);
    save ../kfoldInds kfoldInds
else
    load ../kfoldInds
end

% init solution + error arrays
mDirect = zeros(46,55,Nsubj*Nsl);
percErrmDirect = zeros(length(kfoldInds),1);

mNN = zeros(46,55,Nsubj*Nsl);
percErrmNN = zeros(length(kfoldInds),1);
solErrNN = zeros(length(kfoldInds),1);

mkNN = zeros(46,55,Nsubj*Nsl);
percErrmkNN = zeros(length(kfoldInds),1);
solErrkNN = zeros(length(kfoldInds),1);

mwNN = zeros(46,55,Nsubj*Nsl);
percErrmwNN = zeros(length(kfoldInds),1);
solErrwNN = zeros(length(kfoldInds),1);

mPLS = zeros(46,55,Nsubj*Nsl);
percErrmPLS = zeros(length(kfoldInds),1);
solErrPLS = zeros(length(kfoldInds),1);

mPCR = zeros(46,55,Nsubj*Nsl);
percErrmPCR = zeros(length(kfoldInds),1);
solErrPCR = zeros(length(kfoldInds),1);

mPCRmp = zeros(46,55,Nsubj*Nsl);
percErrmPCRmp = zeros(length(kfoldInds),1);
solErrPCRmp = zeros(length(kfoldInds),1);

mPOCS = zeros(46,55,Nsubj*Nsl);
percErrmPOCS = zeros(length(kfoldInds),1);
solErrPOCS = zeros(length(kfoldInds),1);

for kk = 1:K % loop over the K folds

  fprintf('Processing fold %d...\n',kk);
  
    testInds = find(kfoldInds == kk);
    trainInds = find(kfoldInds ~= kk);
    nTest = length(testInds);
    nTrain = length(trainInds);
    
    % load the shims and features of the slices we don't train on -
    % REDFLAG: Assumes testInds is sorted!
    disp 'Loading test features'
    [~,~,~,testSols,testFeatures,testPercErrs,~,maps] = loadData(testInds,nHood);
    
    % load the training shims and features - REDFLAG: Assumes trainInds is
    % sorted!
    disp 'Loading training features'
    [~,~,~,trainSols,trainFeatures,trainPercErrs] = loadData(trainInds,nHood);
    
    if phsRotate % phase rotate the training shims so we can better interpolate over them
        phiTest = 0:2*pi/1000:2*pi-(2*pi/1000);
        for ii = 2:nTrain
            errTmp = Inf;
            for ll = 1:length(phiTest)
                if norm(trainSols(:,ii-1) - trainSols(:,ii)*exp(1i*phiTest(ll))) < errTmp
                    errTmp = norm(trainSols(:,ii-1) - trainSols(:,ii)*exp(1i*phiTest(ll)));
                    bestPhsInd = ll;
                end
            end
            trainSols(:,ii) = trainSols(:,ii)*exp(1i*phiTest(ll));
        end
    end
    
    
    % get the patterns of the direct design solutions
    for ii = 1:nTest
        
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*testSols(jj,ii);
        end
        mDirect(:,:,testInds(ii)) = tmp;
        percErrmDirect(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % for fun, trace out the cost function as we vary the maximum RF
        % coefficient
        if 0
        [~,maxCoeff] = max(abs([real(testSols(:,ii)); imag(testSols(:,ii))]));
        if maxCoeff > Nc
            maxCoeffIsReal = false;
            maxCoeff = maxCoeff - Nc;
        else
            maxCoeffIsReal = true;
        end
        % get the range of values to be tested
        maxCoeffVals = (-50:0.1:50)*(maxCoeffIsReal*real(testSols(maxCoeff,ii)) ...
            + ~maxCoeffIsReal*imag(testSols(maxCoeff,ii)));
        cost = zeros(size(maxCoeffVals));
        mNoMax = tmp - maps{subjInd}.b1(:,:,slInd,maxCoeff)*(maxCoeffIsReal*real(testSols(maxCoeff,ii)) ...
            + ~maxCoeffIsReal*1i*imag(testSols(maxCoeff,ii)));
        for jj = 1:length(maxCoeffVals)
            cost(jj) = norm(abs(mDirect(:,:,testInds(ii))) - abs(mNoMax + maps{subjInd}.b1(:,:,slInd,maxCoeff)*(maxCoeffIsReal ...
                    + ~maxCoeffIsReal*1i)*maxCoeffVals(jj)))^2;
        end
        % plot the error 
        figure;plot(maxCoeffVals,cost,'k');
        axis([-0.75 1.75 -10 150]);
        % show the patterns at two local minima
        figure;im([tmp; mNoMax + maps{subjInd}.b1(:,:,slInd,maxCoeff)*(maxCoeffIsReal ...
                    + ~maxCoeffIsReal*1i)*1.15]);
        end
        
    end % testInds loop
    
    % run the POCS-shim design
    % first get all the maps for the design into one array
    mapsTrain.b1 = [];
    mapsTrain.mask = [];
    mapsTrain.R = {};
    for ii = 1:nTrain
        
        subjInd = floor(trainInds(ii)/Nsl)+1;
        slInd = rem(trainInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        mapsTrain.b1 = cat(3,mapsTrain.b1,maps{subjInd}.b1(:,:,slInd,:));
        mapsTrain.mask = logical(cat(3,mapsTrain.mask,maps{subjInd}.mask(:,:,slInd)));
        
        mapsTrain.R = cat(1,mapsTrain.R,maps{subjInd}.R);
        
    end
    
    % then get all the features for the design
    % get normalized training features for nearest-neighbors methods
    featMean = mean(trainFeatures,2);
    featStd = std(trainFeatures,0,2);
    trainFeaturesNorm = (trainFeatures-repmat(featMean,[1 nTrain]))...
        ./repmat(featStd,[1 nTrain]);
    trainFeaturesNorm(isinf(trainFeaturesNorm) | isnan(trainFeaturesNorm)) = 0; % get divide by zero for some coordinates
    % get normalized test features
    %featMean = mean(testFeatures,2);
    %featStd = std(testFeatures,0,2);
    testFeaturesNorm = (testFeatures-repmat(featMean,[1 nTest]))...
        ./repmat(featStd,[1 nTest]);
    testFeaturesNorm(isinf(testFeaturesNorm) | isnan(testFeaturesNorm)) = 0; % get divide by zero for some coordinates
    
    % test the POCS-shim design
    % add a 1 to each slice's feature vector to account for offset
    trainFeaturesPOCSshim = [trainFeaturesNorm;ones(1,nTrain)];
    testFeaturesPOCSshim = [testFeaturesNorm;ones(1,nTest)];
    % calculate projector matrix
    F = trainFeaturesPOCSshim.';
    if doPOCS
      Fproj = F*((F'*F + POCSreg*eye(size(F,2)))\F'); % coordinate-wise
    end
%     F = sparse(Nc*2*nTrain,size(trainFeaturesPOCSshim,1)*Nc*2);
%     for ii = 1:nTrain
%         F((ii-1)*Nc*2+1:ii*Nc*2,:) = sparse(kron(speye(Nc*2),trainFeaturesPOCSshim(:,ii).'));
%     end
    
    %Fproj = F*inv(F'*F)*F';
    
    msshim_POCS;
    % get prediction matrix back and apply to get the test RF
    rfPOCS = sqz(rfPOCS);
    % get A matrix
    APOCS = [];
    Finv = (F'*F + POCSreg*eye(size(F,2)))\F';
    for ii = 1:Nc
        APOCS = [APOCS; (Finv*real(rfPOCS(ii,:)).').'; (Finv*imag(rfPOCS(ii,:)).').'];
    end
    %APOCS = reshape(F\col([real(rfPOCS);imag(rfPOCS)]),[size(trainFeaturesPOCSshim,1) 2*Nc]).';
    for ii = 1:nTest
        
        predPOCS = APOCS*testFeaturesPOCSshim(:,ii);
        predPOCS = predPOCS(1:2:2*Nc) + 1i*predPOCS(2:2:2*Nc);
        % predPOCS(1:Nc) + 1i*predPOCS(Nc+1:end);
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predPOCS(jj);
        end
        mPOCS(:,:,testInds(ii)) = tmp;
        percErrmPOCS(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));

        % get error between best solutions and predicted solutions
        solErrPOCS(testInds(ii)) = norm(testSols(:,ii) - predPOCS)/norm(testSols(:,ii));
        
    end
    
    
    % get the patterns of the Nearest-Neighbor solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        % get predicted weights using nearest-neighbors (NN)
        % do this by selecting good solution with min norm between
        % features; if multiple with same norm, take the one with
        % lowest error
        nnNorm = [];
        for jj = 1:nTrain
            nnNorm(jj) = norm(testFeaturesNorm(:,ii) - trainFeaturesNorm(:,jj));
        end
        minInd = find(nnNorm == min(nnNorm));
        predNN = trainSols(:,minInd);
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predNN(jj);
        end
        mNN(:,:,testInds(ii)) = tmp;
        percErrmNN(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrNN(testInds(ii)) = norm(testSols(:,ii) - predNN)/norm(testSols(:,ii));
        
    end % testInds loop
    
    
    % get the patterns of the k-Nearest-Neighbor solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        % get predicted weights using nearest-neighbors (NN)
        % do this by selecting good solution with min norm between
        % features; if multiple with same norm, take the one with
        % lowest error
        knnNorm = [];
        for jj = 1:nTrain
            knnNorm(jj) = norm(testFeaturesNorm(:,ii) - trainFeaturesNorm(:,jj));
        end
        [~,kInds] = sort(knnNorm); % norms in ascending order
        predkNN = mean(trainSols(:,kInds(1:kNeighbors)),2);
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predkNN(jj);
        end
        mkNN(:,:,testInds(ii)) = tmp;
        percErrmkNN(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrkNN(testInds(ii)) = norm(testSols(:,ii) - predkNN)/norm(testSols(:,ii));
        
    end % testInds loop
    
    
    % get the patterns of the weighted-Nearest-Neighbor solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        % get predicted weights using nearest-neighbors (NN)
        % do this by selecting good solution with min norm between
        % features; if multiple with same norm, take the one with
        % lowest error
        wnnKern = [];
        for jj = 1:nTrain
            wnnKern(jj) = 1/sqrt(2*pi)*exp(-1/2*norm(testFeaturesNorm(:,ii) - trainFeaturesNorm(:,jj))^2/lambdaKern^2);
        end
        predwNN = sum(repmat(wnnKern,[Nc 1]).*trainSols,2)./sum(wnnKern);
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predwNN(jj);
        end
        mwNN(:,:,testInds(ii)) = tmp;
        percErrmwNN(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrwNN(testInds(ii)) = norm(testSols(:,ii) - predwNN)/norm(testSols(:,ii));
        
    end % testInds loop
    
    % get real and imaginary values of training solutions
    trainSolsReIm = trainSols.';
    trainSolsReIm = [real(trainSolsReIm) imag(trainSolsReIm)];
    trainSolsMagPhs = [abs(trainSolsReIm) angle(trainSolsReIm)];
    
    % use partial-least squares regression
    betaPLS = [];
    for ii = 1:size(trainSolsReIm,2) % loop over shim values to predict
        [~,~,~,~,betaPLS(:,ii)] = plsregress(trainFeaturesNorm.',trainSolsReIm(:,ii),min(10,size(trainFeaturesNorm,1)));
    end
    % get the patterns of the partial least squares regression solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        predPLS = [];
        for jj = 1:size(trainSolsReIm,2) % loop over shim values to predict
            predPLS(jj) = [1 testFeaturesNorm(:,ii).']*betaPLS(:,jj);
        end
        predPLS = predPLS(1:length(predPLS)/2) + 1i*predPLS(length(predPLS)/2+1:end);
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predPLS(jj);
        end
        mPLS(:,:,testInds(ii)) = tmp;
        percErrmPLS(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrPLS(testInds(ii)) = norm(testSols(:,ii) - predPLS(:))/norm(testSols(:,ii));
        
    end % testInds loop
    
    
    % use principle components regression
    betaPCRt = [];
    betaPCR = [];
    [PCALoadings,PCAScores,~] = pca(trainFeaturesNorm.','Economy',false);
    for ii = 1:size(trainSolsReIm,2) % loop over shim values to predict
        betaPCRt(:,ii) = regress(trainSolsReIm(:,ii)-mean(trainSolsReIm(:,ii)),PCAScores(:,1:min(10,size(trainFeaturesNorm,1))));
        tmp = PCALoadings(:,1:min(10,size(trainFeaturesNorm,1)))*betaPCRt(:,ii);
        betaPCR(:,ii) = [mean(trainSolsReIm(:,ii)) - mean(trainFeaturesNorm.')*tmp; tmp];
    end
    % get the patterns of the partial least squares regression solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        predPCR = [];
        for jj = 1:size(trainSolsReIm,2) % loop over shim values to predict
            predPCR(jj) = [1 testFeaturesNorm(:,ii).']*betaPCR(:,jj);
        end
        predPCR = predPCR(1:length(predPCR)/2) + 1i*predPCR(length(predPCR)/2+1:end);
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predPCR(jj);
        end
        mPCR(:,:,testInds(ii)) = tmp;
        percErrmPCR(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrPCR(testInds(ii)) = norm(testSols(:,ii) - predPCR(:))/norm(testSols(:,ii));
        
    end % testInds loop
    
    
    % use principle components regression with magnitude and phase of the
    % shim values
    betaPCRt = [];
    betaPCR = [];
    %[PCALoadings,PCAScores,~] = pca(trainFeaturesNorm.','Economy',false);
    for ii = 1:size(trainSolsReIm,2) % loop over shim values to predict
        betaPCRt(:,ii) = regress(trainSolsMagPhs(:,ii)-mean(trainSolsMagPhs(:,ii)),PCAScores(:,1:min(10,size(trainFeaturesNorm,1))));
        tmp = PCALoadings(:,1:min(10,size(trainFeaturesNorm,1)))*betaPCRt(:,ii);
        betaPCR(:,ii) = [mean(trainSolsMagPhs(:,ii)) - mean(trainFeaturesNorm.')*tmp; tmp];
    end
    % get the patterns of the partial least squares regression solution
    for ii = 1:nTest
        
        % get subject and slice indices
        subjInd = floor(testInds(ii)/Nsl)+1;
        slInd = rem(testInds(ii),Nsl);
        if slInd == 0;slInd = Nsl;subjInd = subjInd - 1;end;
        
        predPCR = [];
        for jj = 1:size(trainSolsReIm,2) % loop over shim values to predict
            predPCR(jj) = [1 testFeaturesNorm(:,ii).']*betaPCR(:,jj);
        end
        %predPCR = predPCR(1:length(predPCR)/2) + 1i*predPCR(length(predPCR)/2+1:end);
        predPCR = predPCR(1:length(predPCR)/2).*exp(1i*predPCR(length(predPCR)/2+1:end));
        
        tmp = zeros(size(maps{subjInd}.b1,1),size(maps{subjInd}.b1,2));
        for jj = 1:Nc
            tmp = tmp + maps{subjInd}.b1(:,:,slInd,jj)*predPCR(jj);
        end
        mPCRmp(:,:,testInds(ii)) = tmp;
        percErrmPCRmp(testInds(ii)) = 100*std(abs(tmp(maps{subjInd}.mask(:,:,slInd))))/mean(abs(tmp(maps{subjInd}.mask(:,:,slInd))));
        
        % get error between best solutions and predicted solutions
        solErrPCRmp(testInds(ii)) = norm(testSols(:,ii) - predPCR(:))/norm(testSols(:,ii));
        
    end % testInds loop
    
    
end

% make a box plot comparing the patterns
X = [percErrmDirect(:) percErrmNN(:) percErrmkNN(:) percErrmwNN(:) percErrmPCR(:) percErrmPCRmp(:) percErrmPOCS(:)];
figure;subplot(122)
boxplot(X,{'Direct','NN','kNN','KNN','PCR','PCRmp','POCSshim'});
ylabel '|B1+| standard deviation (%)'
axis([0.5000    7.5000   0   100]);

% make a box plot comparing the predictions and solutions
Y = 100*[solErrNN(:) solErrkNN(:) solErrwNN(:) solErrPCR(:) solErrPCRmp(:)];
subplot(121)
boxplot(Y,{'NN','kNN','KNN','PCR','PCRmp'});
ylabel 'Shim Value NRMSE (%)'

% show the best, worst, and median patterns
%mDirect = mDirect(:,:,2:end);
bestDirect = abs(mDirect(:,:,find(percErrmDirect == min(percErrmDirect))));
[~,medInd] = min(abs(percErrmDirect - median(percErrmDirect)));
medianDirect = abs(mDirect(:,:,medInd));
worstDirect = abs(mDirect(:,:,find(percErrmDirect == max(percErrmDirect))));

%mNN = mNN(:,:,2:end);
bestNN = abs(mNN(:,:,find(percErrmNN == min(percErrmNN))));
[~,medInd] = min(abs(percErrmNN - median(percErrmNN)));
medianNN = abs(mNN(:,:,medInd));
worstNN = abs(mNN(:,:,find(percErrmNN == max(percErrmNN))));

%mkNN = mkNN(:,:,2:end);
bestkNN = abs(mkNN(:,:,find(percErrmkNN == min(percErrmkNN))));
[~,medInd] = min(abs(percErrmkNN - median(percErrmkNN)));
mediankNN = abs(mkNN(:,:,medInd));
worstkNN = abs(mkNN(:,:,find(percErrmkNN == max(percErrmkNN))));

%mwNN = mwNN(:,:,2:end);
bestwNN = abs(mwNN(:,:,find(percErrmwNN == min(percErrmwNN))));
[~,medInd] = min(abs(percErrmwNN - median(percErrmwNN)));
medianwNN = abs(mwNN(:,:,medInd));
worstwNN = abs(mwNN(:,:,find(percErrmwNN == max(percErrmwNN))));

%mPCR = mPCR(:,:,2:end);
bestPCR = abs(mPCR(:,:,find(percErrmPCR == min(percErrmPCR))));
[~,medInd] = min(abs(percErrmPCR - median(percErrmPCR)));
medianPCR = abs(mPCR(:,:,medInd));
worstPCR = abs(mPCR(:,:,find(percErrmPCR == max(percErrmPCR))));

%mPCRmp = mPCRmp(:,:,2:end);
bestPCRmp = abs(mPCRmp(:,:,find(percErrmPCRmp == min(percErrmPCRmp))));
[~,medInd] = min(abs(percErrmPCRmp - median(percErrmPCRmp)));
medianPCRmp = abs(mPCRmp(:,:,medInd));
worstPCRmp = abs(mPCRmp(:,:,find(percErrmPCRmp == max(percErrmPCRmp))));

%mPOCS = mPOCS(:,:,2:end);
bestPOCS = abs(mPOCS(:,:,find(percErrmPOCS == min(percErrmPOCS))));
[~,medInd] = min(abs(percErrmPOCS - median(percErrmPOCS)));
medianPOCS = abs(mPOCS(:,:,medInd));
worstPOCS = abs(mPOCS(:,:,find(percErrmPOCS == max(percErrmPOCS))));


bests = [bestDirect./median(bestDirect(bestDirect > 0));...
    bestNN./median(bestNN(bestNN > 0));...
    bestkNN./median(bestkNN(bestkNN > 0));...
    bestwNN./median(bestwNN(bestwNN > 0));...
    bestPCR./median(bestPCR(bestPCR > 0));...
    bestPCRmp./median(bestPCRmp(bestPCRmp > 0));...
    bestPOCS./median(bestPOCS(bestPOCS > 0))];
    
medians = [medianDirect./median(medianDirect(medianDirect > 0));...
    medianNN./median(medianNN(medianNN > 0));...
    mediankNN./median(mediankNN(mediankNN > 0));...
    medianwNN./median(medianwNN(medianwNN > 0));...
    medianPCR./median(medianPCR(medianPCR > 0));...
    medianPCRmp./median(medianPCRmp(medianPCRmp > 0));...
    medianPOCS./median(medianPOCS(medianPOCS > 0))];

worsts = [worstDirect./median(worstDirect(worstDirect > 0));...
    worstNN./median(worstNN(worstNN > 0));...
    worstkNN./median(worstkNN(worstkNN > 0));...
    worstwNN./median(worstwNN(worstwNN > 0));...
    worstPCR./median(worstPCR(worstPCR > 0));...
    worstPCRmp./median(worstPCRmp(worstPCRmp > 0));...
    worstPOCS./median(worstPOCS(worstPOCS > 0))];

figure;
imagesc([fliplr(bests) fliplr(medians) fliplr(worsts)].',[0 2])
colormap gray;
axis image;axis off
