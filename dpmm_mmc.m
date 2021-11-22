% Copyright (C) 2013 Gang Chen, gangchen@buffalo.edu
% distributable under GPL, see README.txt

function [params, tElapsed] = dpmm_mmc(data, num_its, alpha, params)
%function params = dpmm(data, num_its, params)
%standard dirichlet process mixture model, with gaussian observations
%"rao-blackwellised" from, which does not store explicit means or covs


addpath('/home/gangchen/Downloads/classifier/cluster/Tcodes');
ln = @(x)(log(x));

[T, dimension] = size(data);

%some stats
debug = false;
allmean = mean(data,1);
% allcov = cov(data);

C = 1.8; % constant weight % set 3 for glass, wdbc, iris, wine, jain, aggregation,

C = 4.25 % mnist

flag = false;
% flag = true;

tStart = tic;  % TIC,

% minus mean
data = data- repmat(allmean, [T 1]);
% normalize it
% data = data./repmat(sqrt(sum(data.^2,2)), [1 dimension]);

if (~exist('params','var'))
    params(1).alpha = T / 50; %1 / wishrnd(1,1);
    params(1).kappa = .1; %T / 1000; %a pseudo count on the mean
    params(1).nu = 6; %a pseudo-count on the covariance
    params(1).initmean = allmean;
    if flag
        params(1).w = randn(1, dimension+1); % data(randi(200,1),:) - allmean; % initialize weight for the first class
    else
        params(1).w = randn(1, dimension);
    end
    params(1).bias = 0; % initialize the bias
    params(1).m = [];
    params(1).wpca = [];
    % params(1).initcov = allcov / 10;
    params(1).num_classes = 0;
    params(1).counts = 0;
    params(1).sums = [];
    % params(1).cholSSE = [];
    params(1).classes = ones(T,1);
  %  params(1).SSE = [];
    params(1) = addNewClass_mmc(params(1));
    params(1) = unhidupdateweight_mmc(params(1), 1, data);
    if debug, if ~checkParams (params(1), data), disp('no check'); end, end
end
if(exist('alpha', 'var'))
    params(1).alpha = alpha;
end

varnum=size(data,2);
start_it = 1+size(params,2);
oridata= data; 

%hyperparmeter
h1=10;
h0=100;
k1=0.27;
%k1=3;
del=3;
a=3;
b=0.2;
w= 0.8;
gamma = 1;
Vnn=Rn_2(T,gamma); %이건 나중에 MFM할때
%load('Rn_2.mat');


for it = start_it:(start_it+num_its-1)
  c = (params(it-1).classes)';
  r = randi([1 varnum],1,101);
  
  gam = zeros(1,varnum);  
  if it == 2
      for i = r
      gam(i)=1;
      end
      params(it-1).gam=gam;
  else
      for gami = find(params(it-1).gam)
      gam(gami)=1;
      end
  end
  
  t = 20;
    for i = 1:t
      gamcandi=gam;%gamma candidata
      ran=rand;
            if ran >0.5 || prod(gam)~=0 || sum(gam)==0 %1만있거나 0만있거나
                ii=ceil(varnum*rand);
                gamcandi(ii)=abs(gam(ii)-1);
            else
              
                    nonzero=find(gam);
                    in=ceil(length(nonzero)*rand);
                    ind1 = nonzero(in);

                    zero=find(~gam);
                    in2=ceil(length(zero)*rand);
                    ind2=zero(in2);

                    gamcandi(ind1)=0;
                    gamcandi(ind2)=1;
                    
            end 
        while isempty(find(gamcandi,1))==1%==여야함
            if ran >0.5
                ii=ceil(varnum*rand);
                gamcandi(ii)=abs(gam(ii)-1);
            else
                if prod(gam)==0 && sum(gam)>0
                    nonzero=find(gam);
                    in=ceil(length(nonzero)*rand);
                    ind1 = nonzero(in);

                    zero=find(~gam);
                    in2=ceil(length(zero)*rand);
                    ind2=zero(in2);

                    gamcandi(ind1)=0;
                    gamcandi(ind2)=1;
                end     
            end
        end  
        
        cc = unique(c);%클러스터 종류들
        k=length(cc);%클러스터 갯수

        ww= params(it-1).w;
        fXrc_new=GammagivenXc(oridata',c,T,k,gamcandi,ww,C);
        fXrc_old=GammagivenXc(oridata',c,T,k,gam,ww,C);
        
        Gammaprior_new=Gammaprior(varnum,w,gamcandi);
        Gammaprior_old=Gammaprior(varnum,w,gam);
        logprob=(fXrc_new)-(fXrc_old);
        probb=exp(logprob)*exp(Gammaprior_new-Gammaprior_old);
        accprob = min([1,probb]);
        if accprob==1
        gam=gamcandi;
        else
        br=binornd(1,accprob);
            if br == 1
                 gam=gamcandi;
            end
        end    
        
    end
  
    % find(gam)
%% 위에까지 variable selection      
    params(it) = params(it-1);
   params(it).gam=gam;    
%      clear data
%     data = oridata.*gam;
%     data( :, ~any(gam,1) ) = [];
%     
%     ppp= params(it).initmean;
%     params(it).initmean=zeros(1,dimension);
%     ii=1;
%     for i = find(params(it-1).gam)
%         try
%         params(it).initmean(i) = ppp(ii);
%         catch
%         end
%         ii=ii+1;
%     end      
%     params(it).initmean = params(it).initmean .*gam;
%     params(it).initmean( :, find(gam==0) ) = [];
%     
%     numcl=params(it).num_classes;
%     ppp= params(it).w;
%     params(it).w=zeros(numcl,dimension);
%     for j = 1:numcl
%         ii=1;
%         for i = find(params(it-1).gam)
%             params(it).w(j,i) = ppp(j,ii);
%             ii=ii+1;
%         end
%         params(it).w = params(it).w .*gam;        
%     end
%     params(it).w( :, find(gam==0) ) = [];
%     
%     
%     
%     numcl=params(it).num_classes;
%     ppp= params(it).sums;
%     params(it).sums=zeros(numcl,dimension);
%     for j = 1:numcl
%         ii=1;
%         for i = find(params(it-1).gam)
%             params(it).sums(j,i) = ppp(j,ii);
%             ii=ii+1;
%         end
%         params(it).sums = params(it).sums .*gam;
%     end
%     params(it).sums( :, find(gam==0) ) = [];   
 
    disp (strcat (sprintf('iterations = %i [%i]: ',it,params(it).num_classes), ...
        sprintf(' %i',params(it).counts)));
    fprintf('alpha = %f ', params(it).alpha);
    %% Split merge
    if it>2

    ccandi=c;
    claunch=c;
    csplitmerge=c;
    %index뽑기
    iii=randi([1,T]);
    jjj=randi([1,T]);
    if iii==jjj
        while iii==jjj
        iii=randi([1,T]);
        jjj=randi([1,T]);
        end
    end
    ci=c(iii);
    cj=c(jjj);
    kset0=[iii,jjj];
    kset1=find(c==ci);
    kset2=find(c==cj);
    ksett=union(kset1,kset2);
    kset=setdiff(ksett,kset0);    
    
    
        if isempty(kset)==1 %==여야함
            if c(iii)==c(jjj) %split procedure
                newc=max(unique(c))+1;
                ccandi(iii)=newc;    
                qL=splitprob(oridata',c,T,k,gam,iii,jjj,alpha,gamma,ccandi,ww,C);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MFM수정%%%%%%%%%%%%%%%%%%%%%%%
                accprob2 = min([1,qL]);
                if accprob2==1
                c=ccandi;
                else
                br=binornd(1,accprob2);
                    if br == 1
                    c=ccandi;
                    end
                end    
            else %merge procedure
            ccandi(iii)=ccandi(jjj);
                qL=splitprob_2(oridata',c,T,k,gam,iii,jjj,alpha,gamma,ccandi,ww,C);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MFM수정%%%%%%%%%%%%%%%%%%%%%%%
                accprob2 = min([1,1/qL]);
                if accprob2==1
                c=ccandi;
                else
                br=binornd(1,accprob2);
                    if br == 1
                    c=ccandi;
                    end
                end   
            end


        else
        %lauch state(i)
            sww=size(ww);
            if c(iii)==c(jjj)
                iiiwidx=sww(1);
                jjjwidx=c(jjj);
                newc=max(unique(c))+1;
                claunch(iii)=newc;    
            else
                iiiwidx=c(iii);
                jjjwidx=c(jjj);
            end
         %lauch state gibbs sampling scan   
         %lauch state(ii)


         %lauch state(iii)
            for i = 1:5 % t intermediate restricted Gibbs sampling        
                for ii=kset
                    prci = exp(-0.5*sum(ww(iiiwidx, :).^2)+C*ww(iiiwidx, :)*data(t,:)');
                    prcj = exp(-0.5*sum(ww(jjjwidx, :).^2)+C*ww(jjjwidx, :)*data(t,:)');
% 
%                     prci=FH(oridata',claunch,T,k,gam,iii,kset,kset1,ksett,ii,claunch(iii),ww); %kset은 i랑 l이과 같은라벨 같고있는 observation의 집합(i,l 제외)이고(claunch 말고 ci) ii는 그중하나, claunch 로 들어가는게 맞음
%                     
%                     prcj=FH(oridata',claunch,T,k,gam,jjj,kset,kset2,ksett,ii,claunch(jjj),ww);

                    R=binornd(1,prci/(prci+prcj));
                    if R ==1
                        claunch(ii)=claunch(iii);
                    else
                        claunch(ii)=claunch(jjj);
                    end

                end
            end

       %step b            
            
            if c(iii)==c(jjj) % Split
                csplitmerge(iii)=claunch(iii);    
                csplitmerge(jjj)=claunch(jjj);  
                
                f=1;
                fqbmat=[];
                for ii=kset % final gibbs sampling scan 
                    prci= exp(-0.5*sum(ww(iiiwidx, :).^2)+C*ww(iiiwidx, :)*data(t,:)');
                    prcj= exp(-0.5*sum(ww(jjjwidx, :).^2)+C*ww(jjjwidx, :)*data(t,:)');
                    R=binornd(1,prci/(prci+prcj));
                    if R ==1
                        csplitmerge(ii)=csplitmerge(iii);
                        fqbmat(f)=(prci/(prci+prcj));
                    else
                        csplitmerge(ii)=csplitmerge(jjj);
                        fqbmat(f)=(prcj/(prci+prcj));
                    end
                    f=f+1;
                end
                try
                prodfqb=sum(log(1./(fqbmat)));
                catch
                end
                qL=splitprob2(oridata',c,gam,iii,jjj,alpha,prodfqb,csplitmerge,gamma,ww,C);
                accprob3 = min([1,qL]);
                if accprob3==1
                c=csplitmerge;
                else
                br=binornd(1,accprob3);
                    if br == 1
                    c=csplitmerge;        
                    end
                end
            else % Merge
                csplitmerge(iii)=c(jjj);    
                csplitmerge(jjj)=c(jjj);  
                for ii=kset % final gibbs sampling scan 
                csplitmerge(ii)=c(jjj);
                end
                
                f=1;
                fqbmat2=[];
                
                for ii=kset % final gibbs sampling scan   나중에 merge된 c launch로부터 이전의 split된 c로 갈 확률 구한다 X 거꾸로인듯...
                    prci = exp(-0.5*sum(ww(iiiwidx, :).^2)+C*ww(iiiwidx, :)*data(t,:)');
                    prcj = exp(-0.5*sum(ww(jjjwidx, :).^2)+C*ww(jjjwidx, :)*data(t,:)');
                    if claunch(ii)==cj
                        fqbmat2(f)=(prcj/(prci+prcj));
                        f=f+1;
                    else
                        fqbmat2(f)=(prci/(prci+prcj));;
                        f=f+1;  
                    end
                end
                prodfqb=sum(log(fqbmat2));
                qL=splitprob3(oridata',c,gam,iii,jjj,alpha,prodfqb,csplitmerge,gamma,ww,C);
                accprob3 = min([1,qL]);
                if accprob3==1
                c=csplitmerge;
                else
                br=binornd(1,accprob3);
                    if br == 1
                    c=csplitmerge;        
                    end
                end
                
            end
            
            
            
        end  
    end
    
    
    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GIBBS SAMPLING %%%%%%%%%%%%%%%%%%%%%%%
    t_order = randperm(T);
    for sctr = 1:T        
        %t = sctr+1;%for debugging
        t = t_order(sctr);

        old_class = params(it).classes(t);
        % Gang: consider late on how to delete the current point's contribution to this class
        %% params(it) = hideObservations(params(it),old_class,data(t,:));
        params(it) = hidupdateweight(params(it),old_class,data(t,:));
        params(it) = handleRemovedClasses_mmc(params(it));
        if debug, if ~checkParams(params(it),data,t), disp('no check at hide'); end, end
        
        %these are the probabilities that we will sample from
        %note we add one to include the chance of adding a new class
        log_p_obs = -inf * ones(params(it).num_classes+1,1);

        p_prior = [];
                
        %% params(it) = addNewClass(params(it));  %it will be removed if nothing is put in it
        params(it) = addNewClass_mmc(params(it), flag);  %it will be removed if nothing is put in it
        if debug, if ~checkParams(params(it),data,t), disp('no check at add class'); end, end
        
        kappabar = params(it).counts + params(it).kappa;
        nubar = params(it).counts + params(it).nu;
        factor = (kappabar + 1) ./ (kappabar .* (nubar - dimension - 1));
        %p_prior = params(it).counts + params(it).alpha * (params(it).counts == 0);        
        p_prior = (params(it).counts+gamma) + (gamma*Vnn(params(it-1).num_classes)) * (params(it).counts == 0);        

        
        for i = 1:params(it).num_classes
% %            if (params(it).counts(i) == 0), p_prior(i) = params(it).alpha; 
% %            else p_prior(i) = params(it).counts(i); end
%             try
%                 %integrating over the parameters of a
%                 %normal-inverse-Wishart yields student-t.  
%                 %this can be approximated by a "moment-matched" Gaussian, 
%                 %see sudderth thesis p 47
%                 %kappabar = params(it).counts(i) + params(it).kappa;
%                 %nubar = params(it).counts(i) + params(it).nu;
%                 %factor = (kappabar + 1) / (kappabar * (nubar - dimension - 1));
%                 log_p_obs(i) = normpdfln(data(t,:)', ...
%                     params(it).sums(i,:)' / kappabar(i),...
%                     sqrt(factor(i))*params(it).cholSSE(:,:,i));
%             catch
%                 disp('mvnpdf throws error');
%             end
            log_p_obs(i) = 0.5*sum(params(it).w(i, :).^2);
            score = 0;
            try
               if flag
                   score = params(it).w(i, :)*[1 data(t,:)]';
               else
                   score =  params(it).w(i, :)*data(t,:)';
               end
               %if score < 1
               %    log_p_obs(i) = log_p_obs(i) + C*(1-score);
               %end
               log_p_obs(i) =  -log_p_obs(i) + C*score;% + sum(data(t,:).^2);
               % log_p_obs(i) =  score/log_p_obs(i);
               
               % log_p_obs(i) =ln( 1/(1+exp(-log_p_obs(i))));
            catch
               disp('mvnpdf throws error');
            end

            
        end
        % log_p_obs = exp(log_p_obs)./ repmat(sum(exp(log_p_obs),1),[params(it).num_classes 1]);


        %lightspeed sample normalizes automatically
        classprobs = p_prior'.*exp(log_p_obs-max(log_p_obs));% classprobs = p_prior'.*log_p_obs;
        
        try
            new_class = sample(classprobs);
            if (params(it).counts(new_class) == 0)
%                disp('adding a guy');
            end
            params(it).classes(t) = new_class;
        catch
            disp('could not sample');
        end
        if debug, if ~checkParams(params(it),data,t), disp('no check at sample'); end, end
        % Gang: consider on line updating for the current model w
        %% params(it) = unhideObservations(params(it),new_class,data(t,:));
        if new_class == 2
            stop = 1;
        end
        params(it) = unhidupdateweight_mmc(params(it),new_class,data(t,:),flag, t_order(1:sctr));
        
        
%         if new_class >1 % more than one class, then deploy SVM or LDA
%             %%add more updating based on the current label
%             %cdata = data(t_order(1:sctr), :);
%             %clabel  = params(it).classes(t_order(1:sctr));
%             %%updating the model based on the current training data
%             params(it) = onlearnmodel(params(it), data, t_order(1:sctr));
%         end
        if debug, if ~checkParams(params(it),data), disp('no check at hide'); end, end
        
    end
    
    %%%%%%%%%%%%%%%%%%%% PARAMETER UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alpha is the "pseudo-count" for new classes.  it is estimated using ARS
    k = params(it).num_classes;
    n = T; 
    
    %can show that derivative is guaranteed to be positive / negative at
    %these points
    deriv_up = 2 / (n - k + 3/2);
    deriv_down = k * n / (n - k + 1);
    
    %this is the version with a conjugate inverse gamma prior on alpha, as
    %in Rasmussen 2000
    params(it).alpha = ars(@logalphapdf, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);

    %this is the version with a totally non-informative prior
    %params(it).alpha = ars(@logalphapdfNI, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);

end

    % params(it) = relabel(params(it), data);
    
    tElapsed = toc(tStart);  % TOC, pair 2  
end

%checks a set of parameters to see if they are self-consistent.
%for debugging
function [total c_basic c_count c_sum] = checkParams(params,data,exclude)
    if exist('exclude','var')
        c_basic = min(params.classes([1:exclude-1 exclude+1:end]) > 0);
    else
        c_basic = min(params.classes > 0);
    end
    c_count = 1;
    c_sum = 1;
    for i = 1:params.num_classes
        statedata = data(find(params.classes == i),:);
        err_amount = params.sums(i,:) - sum(statedata) - params.kappa * params.initmean;
        statecount = size(statedata,1);
        if exist('exclude','var')
            if i == params.classes(exclude) 
                err_amount = err_amount - data(exclude,:); 
                statecount = statecount - 1;
            end
        end
        if (statecount ~= params.counts(i)), c_count = 0; end
        if (sum(err_amount) > .01), c_sum = 0; end
    end        
    total = c_basic * c_count * c_sum;
end