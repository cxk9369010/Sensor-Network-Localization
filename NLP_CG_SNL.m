function [R,Q]=NLP_CG_SNL(P0,D0,DD,outputdim,randSeed,pars);
%%
% NLP_CG_SNL: Solving the NonLinear Programming model by using CG method
%
% min \sum_(j,i)\in N_x (||R_(j,:)-R_(i,:)||^2-d_ji^2)^2+ 
%       \sum_(j,k)\in N_a (||R_(j,:)-Qa_k||^2-d_jk^2)^2
%
% where d_ji and d_jk are distances, Q is a dim-by-outputdim(dim+l) matrix 
% such that QQ^T=I_dim.
%
% Input:
% P0 -- dim-by-m matrix of coordinates of m anchors 
% D0 = [sensor - anchor (squared) distance];
% DD = [sensor - sensor (squared) distance];
% 
% dim:          The embedding dimensions (e.g. dim = 2 or 3)
% outputdim:    dim+l; (solving in dim+l dimensional space);
% randSeed:     the number of random problems you wish to solve
%
% pars: parameters and other information
%
% pars.tol_gk  = m -- the number of known anchors, =0 if none.
% pars.tol_fk  = PP = [PA, PS] -- Positions of Ancors/Sensors

% Output
%
% R    -- final localizations of the sensors in dim+l dimensional space;
% Q    -- QQ^T=I_dim
% QR^T -- projection to dim dimensional space;

%%
% Based on the paper:
% A New Unconstraint Nonlinear Programming Method for Solving 
% Sensor Network Localization
% by Xiaokai Chang, Sanyang Liu and Xunyang Wang (2016)
%
% Send your comments and suggestions to    %%%%%%
%          15293183303@163.com                %%%%%%
%
%%%%% Warning: Accuracy may not be guaranteed!!!!! %%%%%%%%

%%%%%%%%% This version:   October 7, 2016     %%%%%%%
if isfield(pars, 'tol_gk')
    tolerance_gk = pars.tol_gk;    % Stopping criterion. 
else
    tolerance_gk = 1e-06;
end
if isfield(pars, 'tol_fk')
    tolerance_fk = pars.tol_fk;    % Stopping criterion. 
else
    tolerance_fk = 1e-10;
end
if isfield(pars, 'showyes')
    showyes = pars.showyes;
else
    showyes = 0; % no plot
end

NumOfEdge = full(sum(sum(D0~=0))+sum(sum(DD~=0)));

dim = size(P0,1);
[stemp,mtemp] = size(D0);
%%
%                        Generating N_x and N_a (connection matrix)
[d1,d2,N_x,N_a] = initseries(D0,DD,stemp,mtemp);
%%
%                        Generating Initial Point
%[R,Q]=v0(stemp,outputdim,P0,D0',dim,randSeed);%迭代初始值
[R,Q] = goodv0(stemp,outputdim,P0,D0',DD,dim,randSeed);%迭代初始值
% figure(2)
% plotpositions(P0,PP,Q*R',[],1); 
% Xopt0=Q*R';
% errtrue = sum((Xopt0-PP).*(Xopt0-PP));
% RMSD = sqrt(sum(errtrue))/sqrt(stemp)

% Wolf Conditions and CG direction Parameter
c1 = 0.0001;     %  Wolf Conditions Parameter
ata = 0.01;      %  CG direction Parameter


maxiter   =2000;   % The maximum number of iterations.    



fprintf('@@@@  Starting Conjugate Gradient Method with Exact Linesearch \n');
 if showyes == 1  
    fprintf('iter     norm_gk       ave-fval  \n');
    disp(['-----------------------------------------------------------------' ])
 end
 for iter = 1:maxiter
%% (1). Compute the CG direction.

     if iter==1;
         %[D00]=gradient(outputdim,stemp,V0,A1,A2,d1,d2);
         [g0]=gradient1(outputdim,R,Q,D0,DD,P0);
         dir=-1*g0;
     else
         g0=g1;
         dir=-g0+belta*dir;
         fval0 = fval;
     end
%% (2). Compute step along the search direction
     %[f]=fun(V0,D00,D0,DD,PPtemp)
     [a,b,c,d,fval]=coefficient_SPARSE(Q,R,dir,P0,d1,d2,N_x,N_a); %%  Computing a,b,c,d,e
     if  d>0
         error('The direction is not decreasing.');  
     end
     [alpha]=linesearch(a,b,c,d,c1);    %%  Computing the step length (alpha)
     
%% (3). Move to the new point and Compute the Gradient.

      R = R+alpha*dir;
      [g1] = gradient1(outputdim,R,Q,D0,DD,P0);
      norm_gk = norm(g1,inf);
  
     if showyes == 1
        if iter==1|| mod(iter,50) == 0;
            fprintf('%3d   %0.3e     %0.5e   \n',iter,norm_gk,fval/NumOfEdge);
        end
     end
%% (4). Check the stopping criterion.
       if  norm_gk<=tolerance_gk
          %fprintf('%3d   %0.3e     %0.5e   \n',iter,norm_gk,fval/NumOfEdge);
           break
      %else if  norm_gk<=100*tolerance_gk & abs(fval0-fval)/(abs(fval0)+1)<tolerance_fk;
      elseif  iter>=2  
              if  abs(fval0-fval)/(abs(fval0)+1)<tolerance_fk;
      %else if  norm_gk<=100*tolerance_gk & alnn<tolerance
      %else if norm_gk<=50*tolerance_gk & alnn<1e-15*fval
            %     fprintf('%3d   %0.3e     %0.5e   \n',iter,norm_gk,fval/NumOfEdge);
              break
              end
       end
%% (5). Prepare for computing the CG direction
      if iter>=1
           y = g1-g0;
         dTy = mytrace(dir,y);
        ata1 = -1/(norm(dir)*min(ata,norm(g1)));
    %belta=1/dTy*mytrace(y-(2*mytrace(y,y)/dTy)*dir,g1);
       belta = 1/dTy*mytrace(y,g1);
       belta = max(ata1,belta);
    %belta=1/dTy*mytrace(g1,g1);
    %belta=1/mytrace(g0,g0)*mytrace(g1,g1);
      end
 end    

 fprintf('@@@@  Finished Conjugate Gradient Method with Exact Linesearch \n');
 fprintf('iter = %3d \n',iter);
 fprintf('norm_gk = %0.3e \n',norm_gk);
 fprintf('fval/NumOfEdge = %0.5e \n',fval/NumOfEdge);
 
 
 

   
 
   
   
       
  