clear  all;
results=[];
spse=0;
randSeed = 1; 
dim = 2;
l = 8;
outputdim = dim+l; %  d+l
noOfSensors = 900;
noOfAnchors = 100; 
rd = 0.1;                 %  distance radius 

nf = 0.1;                  %  nosiy factors
anchorType=2;              % 1--grid points 2---placed randomly 
SecondLocalization = 1;    % for noiseless problems
fprintf('####  noOfSensors = %3d \n',noOfSensors);
fprintf('####  noOfAnchors = %3d \n',noOfAnchors);
fprintf('####  nf = %1.1e \n',nf);
fprintf('####  radiorange = %1.1e \n',rd);
fprintf('####  dim = %3d \n',dim);
 
%   anchorType = 1
%       all anchors are placed at grid points in the interior of
%       [0,1]^{sDim}
% 	anchorType = 2
%       all anchors are placed randomly in [0,1]^{sDim}
%   anchorType = 3
%       sDim + 1 anchors on the origin and the coordinate axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: Generate problems as in SFSDP proposed in 
% Exploiting sparsity in SDP relaxation for sensor network localization
% by S. Kim, M. Kojima and H. Waki, (2009)
[xMatrix0,distanceMatrix0] = generateProblem(dim,nf,... 
    rd,noOfSensors,anchorType,noOfAnchors,randSeed);
 if spse==1
           pars.sparseSW = 10;
           pars.localMethodSW = 1;
          % [xMatrix0,distanceMatrix0] = SPARSEDATA(dim,noOfSensors,noOfAnchors,xMatrix0,distanceMatrix0,pars);
          degree=20;
           [distanceMatrix0] = gensparse(distanceMatrix0,degree);
 end


P0=xMatrix0(1:dim,noOfSensors+1:end);   %%   anchors
PP=xMatrix0(1:dim,1:noOfSensors);       %%   sensors

Dall=distanceMatrix0;                   %%   Dall=[DD D0];
DD=distanceMatrix0(1:noOfSensors,1:noOfSensors);   %%   sensor-sensor (squared) distance
D0=distanceMatrix0(1:noOfSensors,noOfSensors+1:end);   %%    sensor-anchor (squared) distance

NumOfEdge=full(sum(sum(D0~=0))+sum(sum(DD~=0)));   
fprintf('####  the number of distance (NumOfEdge) = %3d \n',NumOfEdge);

if dim==2
    pars.tol_gk = 1e-6;  % Stopping criterion. 
    if nf==0 
       pars.tol_fk = 1e-9;
    else
       pars.tol_fk = 1e-8;
    end
elseif dim==3
    pars.tol_gk = 1e-6;  % Stopping criterion. 
    if nf==0 
       pars.tol_fk = 1e-8;
    else
       pars.tol_fk = 5e-7;
    end
end
pars.showyes = 0;
pars.refinement = 1;
startingTime = tic; % cputime; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2: Solving the unconstraint NLP problems 
%%         By Conjugate Gradient Method with Exact Linesearch 
[R,Q]=NLP_CG_SNL(P0,D0,DD,outputdim,randSeed,pars);
plotpositions(P0,PP,Q*R',[],1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2':     Second localizatiin for the problem with exact distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if nf==0 && SecondLocalization==1; 
     pars.tol_gk = 1e-6;  % Stopping criterion. 
     if dim==2
       pars.tol_fk = 1e-10;
     elseif dim==3
         pars.tol_fk = 5e-10;
     end
       pars.plotyes = 0;
  [R,Q]=secondLocalization(Q,R,P0,D0,DD,outputdim,randSeed,dim,pars,PP);    
 end
Xopt0= Q*R';
%V0=[Q;R];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output information

fprintf('####  elapsed time for CG methods to solve NLP model= %8.2fs\n',toc(startingTime))  
BoxScale=1;
OPTIONS.alpha       = 1; %% regularization parameter 
OPTIONS.plotyes     = 1; 
OPTIONS.PP          = PP;   
OPTIONS.nf          = nf; 

plotyes  = 1;  
printyes = 1;
refinemaxit =3000;
   
if (plotyes)
     
     errtrue = sum((Xopt0-PP).*(Xopt0-PP));
     err=max(sqrt(errtrue));
	 RMSD1 = sqrt(sum(errtrue))/sqrt(noOfSensors);
     fprintf(' RMSD = %4.2e',full(RMSD1));
      fprintf(' \n');
     fprintf(' err = %4.2e',full(err)); 
      figure(101)
      plotpositions(P0,PP,Xopt0,[],BoxScale); 
      xlabel([' RMSD ',sprintf('%4.2e',full(RMSD1))]);
      title('Before refinement'); 
 end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3: Refinement step (code by Kim-Chuan)
% Refinement step (Step 4 in the code) is taken from Kim-Chaun Toh SFSDP
% solver
 
  if (refinemaxit)
      fprintf('\n-------------------------------------------------------')
      fprintf('\n        refine positions by steepest descent');
      fprintf('\n-------------------------------------------------------')
      tstart = clock; 
      plotidx = [];
      [XSD,Info] = refinepositions(Xopt0,P0,Dall,refinemaxit);
      %[XSD,Info] = refinepositions(Xopt0,P0,distanceMatrix1,refinemaxit);
      Xopt = XSD;
      RefineInfo = Info; 
      ttime = etime(clock,tstart);
      len = length(Info.objective);
      if (printyes)
         fprintf('\n####  objstart = %5.4e, objend = %5.4e',...
         Info.objective(1),Info.objective(len));
         fprintf('\n####  number of iterations = %2.1d, time = %4.1fs\n',len,ttime);
      end     
      if (plotyes)
         
         errtrue = sum((Xopt-PP).*(Xopt-PP));   
         err=max(sqrt(errtrue));
	      RMSD2 = sqrt(sum(errtrue))/sqrt(noOfSensors); 
          fprintf('####  err-refine = %4.2e',full(err));
          fprintf(' \n');
          fprintf('####  RMSD = %4.2e',full(RMSD2));
          fprintf(' \n');
          figure(103)
          plotpositions(P0,PP,Xopt,[],BoxScale); 
          xlabel(['RMSD %R = ',sprintf('%4.2e',full(RMSD2))]);
          title('After refinement'); 
      end
  end
 % plot(tvartemp,errr,'b*')

  result=[NumOfEdge eval(vpa(RMSD1,2)) eval(vpa(RMSD2,2)) toc(startingTime)];      
 
        
