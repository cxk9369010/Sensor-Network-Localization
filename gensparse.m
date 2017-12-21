function [distanceMatrix1] = gensparse(distanceMatrix0,degree);

% distanceMatrix0 :  a noOfSensors x (noOfSensors + noOfAnchors) upper
%                   triangular matrix to represent distances from sensors
%                   to sensors and anchors, where the (p,q)th element
%                   d_{pq} denotes the distance from the pth sensor to the qth
%                   sensor or anchor. When noisyFac = 0,
%                       d_{pq} = norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   When \sigma = noisyFac > 0,
%                       d_{pq} = max(1+epsilon,0.1)*(norm(xMatrix0(:,p) - xMatrix0(:,q))
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) <= radiorange,
%                       d_{pq} = 0 
%                           if norm(xMatrix0(:,p) - xMatrix0(:,q)) > radiorange,
%                   Here epsilon is a random number chosen from N(0,\sigma). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stemp,ntemp]=size(distanceMatrix0);
distanceMatrix1=zeros(stemp,ntemp);

D0=distanceMatrix0(1:stemp,stemp+1:ntemp);
DD=distanceMatrix0(1:stemp,1:stemp);
%DD=DD+DD';

   for i=1:stemp
       indx=find(DD(i,:));
       indx=indx(randperm(length(indx)));
       indx=indx(1:min(3+degree,length(indx)));
       distanceMatrix1(i,indx)=DD(i,indx);
   end
   
  for i=1:stemp
       indx=find(D0(i,:));
       indx=indx(randperm(length(indx)));
       indx=indx(1:min(3,length(indx)));
       distanceMatrix1(i,indx+stemp)=D0(i,indx);
   end
noOfEdges = nnz(distanceMatrix1); 
noOfEdges2 = nnz(distanceMatrix1(:,1:stemp)); 
degreeVector = sum(spones([distanceMatrix1(:,1:stemp)+distanceMatrix1(:,1:stemp)',...
    distanceMatrix1(:,stemp+1:ntemp)]),2); 
minDeg = full(min(degreeVector')); 
maxDeg = full(max(degreeVector')); 
averageDeg = full(sum(degreeVector')/stemp); 

fprintf('####  the number of dist. eq. between two sensors  = %d\n',noOfEdges2);
fprintf('####  the number of dist. eq. between a sensor & an anchor = %d\n',noOfEdges-noOfEdges2);
fprintf('####  the min., max. and ave. degrees over sensor nodes = %d, %d, %6.2f\n',minDeg,maxDeg,averageDeg);
  
end