function [rf,errAll,vopAll,mAll,randAmpPhs] = msShim_randStart_POCS_vopReg(maps,algp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = algp.beta; % SAR regularization parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build system matrix for each slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dimxy(1),dimxy(2),Nsl,Nc] = size(maps.b1);

for idx = 1:Nsl
    
    % reorg b1 maps into tall skinny matrix (row: space, col: coil)
    tmp = permute(squeeze(maps.b1(:,:,idx,:)),[3 1 2]);
    tmp = tmp(:,:).';
    % mask them
    tmp = tmp(logical(maps.mask(:,:,idx)),:);
    % remove coil 1 phase
    A{idx} = bsxfun(@times,exp(-1i*angle(tmp(:,1))),tmp);
    if algp.noCG
        ARinv{idx} = (A{idx}'*A{idx} + beta*(maps.R{idx}'*maps.R{idx}))\A{idx}';
    end
    
    % init target phase pattern
    if ~isfield(maps,'phsinit')
        dphs{idx} = zeros(sum(col(maps.mask(:,:,idx))),1);
    else
        tmp = maps.phsinit(:,:,idx);
        dphs{idx} = tmp(col(maps.mask(:,:,idx)));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do compressed design for increasing # of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a set of random amplitudes and phases to use in initializing each
% slice's target phase pattern
if ~isfield(algp,'randAmpPhs')
    randAmpPhs = rand(Nc,2,algp.nRandStart);
else
    randAmpPhs = algp.randAmpPhs;
end
nRandStart = algp.nRandStart + 2; % also add zero phase and direct sum phase
rf = zeros(Nc,nRandStart,Nsl);
errAll = zeros(nRandStart,Nsl);
vopAll = zeros(nRandStart,Nsl);
mAll = zeros(dimxy(1),dimxy(2),nRandStart,Nsl);
tol = algp.tol;

for raIdx = 1:nRandStart
    
    % get initial phase patterns for each slice
    for slIdx = 1:Nsl
        
        % precalculate regularized pseudoinverse
        %Ainv{slIdx} = (A{slIdx}'*A{slIdx} + R)\A{slIdx}';
        
        % get random initial starting phase
        switch raIdx
            case 2
                dphs{slIdx} = zeros(size(A{slIdx},1),1);
            case 1
                dphs{slIdx} = angle(sum(A{slIdx},2));
            otherwise
                dphs{slIdx} = angle(A{slIdx}*(randAmpPhs(:,1,raIdx-2).*exp(1i*2*pi*randAmpPhs(:,2,raIdx-2))));
        end
        
    end
    
    % outer loop: continue until shims stop changing
    flag = 1; % to get loop started
    itr = 1; % iteration counter
    holesExist = false;
    while flag
        
        rfOld = sqz(rf(:,raIdx,:));
        
        % update each slice's shims to minimize shim error
        if algp.noCG
            for slIdx = 1:Nsl
                rf(:,raIdx,slIdx) = (1-algp.step)*rf(:,raIdx,slIdx) + ...
                    algp.step*(ARinv{slIdx}*exp(1i*dphs{slIdx}));
            end
        else
            parfor slIdx = 1:Nsl
                % take a couple CG iterations
                xS = qpwls_pcg(rf(:,raIdx,slIdx),A{slIdx},1,...
                    exp(1i*dphs{slIdx}),0,sqrt(beta)*maps.R{slIdx},1,algp.ncgiters,ones(Nc,1));
                rf(:,raIdx,slIdx) = xS(:,end);
            end
        end
        
        % optional: project shims onto space of low-rank shims
        % threshold the SVD of the weights matrix
        if isfield(algp,'Nccomp')
          [u,s,v] = svd(sqz(rf(:,raIdx,:)),'econ');
        
          % Hard truncation
          s((algp.Nccomp+1):end,(algp.Nccomp+1):end) = 0;
          rf(:,raIdx,:) = u*(s*v');
        end
        %compWts = u(:,1:algp.Nccomp);
        
        % project shims onto space of predictable shims
        if isfield(algp,'Fproj');
            % all-together
%             tmp = sqz(rf(:,raIdx,:));
%             tmp = col([real(tmp);imag(tmp)]);
%             tmp = algp.Fproj*tmp;
%             tmp = reshape(tmp,[2*Nc Nsl]);
%             rf(:,raIdx,:) = tmp(1:Nc,:) + 1i*tmp(Nc+1:end,:);
            % coordinate-wise
            for ii = 1:Nc
                % collect real and imag parts of each slice's shim for this
                % coil
                rfr = real(sqz(rf(ii,raIdx,:)));
                rfi = imag(sqz(rf(ii,raIdx,:)));
                % project each of them
                rfr = algp.Fproj*rfr;
                rfi = algp.Fproj*rfi;
                % put them back into the pulse array
                rf(ii,raIdx,:) = rfr+1i*rfi;
            end
        end
        
        % update target phase patterns and calculate error
        rfBest = median(sqz(rf(:,raIdx,:)),2);
        for slIdx = 1:Nsl
            m = A{slIdx}*rf(:,raIdx,slIdx);
            if any(abs(m) < algp.holeThresh) && itr > 100
                m = A{slIdx}*rfBest;
                holesExist = true;
            end
            dphs{slIdx} = angle(m);
            % save error
            errAll(raIdx,slIdx) = 1/2*norm(m-exp(1i*dphs{slIdx}))^2;
            vopAll(raIdx,slIdx) = 1/2*beta*norm(maps.R{slIdx}*rf(:,raIdx,slIdx))^2;
            % save shimmed pattern
            mAll(:,:,raIdx,slIdx) = embed(m,logical(maps.mask(:,:,slIdx)));
        end
        
        itr = itr + 1;
        
        % check stopping criterion
        if (norm(sqz(rf(:,raIdx,:))-rfOld)/norm(rfOld) < tol && ~holesExist) || itr > algp.maxIters
            %if Nccomp == Nc
                flag = 0;
            %else
            %    Nccomp = Nccomp + 1;
            %end
        end
        
        if rem(itr,100) == 0;
            fprintf('Initialization %d. Iteration %d. Err all %f. SAR all %f.\n',raIdx,itr,sum(errAll(raIdx,:),2),sum(vopAll(raIdx,:),2));
            if algp.dofigs
                figure(1);clf;im(mAll(:,:,raIdx,:));
                drawnow
            end
        end
        
    end % while loop over all updates
    
    fprintf('Final stats for Initialization %d: Iterations %d. Err all %f. SAR all %f.\n',raIdx,itr,sum(errAll(raIdx,:),2),sum(vopAll(raIdx,:),2));
    
end % loop over random starts


