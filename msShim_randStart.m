function [rf,costAll,mAll,randAmpPhs] = msShim_randStart(maps,algp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = algp.beta; % 0    % RF regularization parameter
R = maps.R;

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

%R = algp.R;%beta*eye(Nc);

% get a set of random amplitudes and phases to use in initializing each
% slice's target phase pattern
if ~isfield(algp,'randAmpPhs')
    randAmpPhs = rand(Nc,2,algp.nRandStart);
else
    randAmpPhs = algp.randAmpPhs;
end
nRandStart = algp.nRandStart + 2; % also add zero phase and direct sum phase
rf = zeros(Nc,nRandStart,Nsl);
costAll = zeros(nRandStart,Nsl);
mAll = zeros(dimxy(1),dimxy(2),nRandStart,Nsl);
tol = algp.tol;
parfor slIdx = 1:Nsl
    
    % precalculate regularized pseudoinverse
    RtR = beta*(R{slIdx}'*R{slIdx});
    Ainv = (A{slIdx}'*A{slIdx} + RtR)\A{slIdx}';

    for raIdx = 1:nRandStart

        % get random initial starting phase
        switch raIdx
            case 1
                dphs = zeros(size(A{slIdx},1),1); 
            case 2
                dphs = angle(sum(A{slIdx},2));
            otherwise
                dphs = angle(A{slIdx}*(randAmpPhs(:,1,raIdx-2).*exp(1i*2*pi*randAmpPhs(:,2,raIdx-2))));    
        end

        flag = 1; % to get loop started
        cost = Inf; % initial cost
        itr = 1; % iteration counter
        while flag
            
            % calculate RF weights
            rft = Ainv*exp(1i*dphs);
            
            % calculate the excitation pattern
            m = A{slIdx}*rft;
            
            % update the target phase pattern
            dphs = angle(m);
            
            % calculate error
            itr = itr + 1;
            Err = 1/2*norm(m-exp(1i*dphs))^2;
            SAR = 1/2*(rft'*(RtR*rft));
            cost(itr) = Err + SAR;
            
            % check stopping criterion
            if cost(itr) > tol*cost(itr-1)
                flag = 0;
            end
            
            if rem(itr,1000) == 0;
                fprintf('Slice %d. Init %d. Iter %d. Err %f. SAR %f\n',slIdx,raIdx,itr,Err,SAR);
            end
            
        end    
        
        % save cost
        costAll(raIdx,slIdx) = cost(end);
        
        % save RF
        rf(:,raIdx,slIdx) = rft;

        % save shimmed pattern
        mAll(:,:,raIdx,slIdx) = embed(m,logical(maps.mask(:,:,slIdx)));

    end
    
    
end