function [R,Q]=secondLocalization(Q,R,P0,D0,DD,outputdim,randSeed,dim,pars,PP);  
 Xopt0= Q*R';
 Y_opt=R*R';
 noOfSensors=size(DD,2);
for j=1:noOfSensors
    tvartemp(j)=max(0,Y_opt(j,j)-R(j,:)*Q'*Q*R(j,:)');
end;

% err=(Xopt0-PP).*(Xopt0-PP);
% errr=sqrt(err(1,:)+err(2,:));
%    figure(555)
%    
%   plot(tvartemp,'gs')
%   hold on
%   plot(errr,'b*')
%   kkk=1:noOfSensors;
%   plot([kkk;kkk],[tvartemp(1,:);errr(1,:)],'-r')
%   
%   corrcoef(tvartemp,errr)
  if dim==2
     index_tr=find(tvartemp>1e-04);
  else
     index_tr=find(tvartemp>1e-04);
  end
% figure(111);
% if dim==2
%     plot(PP(1,index_tr),PP(2,index_tr),'r*');
% else
%     plot3(PP(1,index_tr),PP(2,index_tr),PP(3,index_tr),'r*');
% end

DDD=DD+DD';
tr=sparse(size(DDD));
[index1,index2]=find(DDD~=0);
for i=1:length(index1)  
    tr(index1(i),index2(i))=norm(R(index1(i),:)-R(index2(i),:))^2-norm((R(index1(i),:)-R(index2(i),:))*Q')^2;      
end

[index1,index2]=find(tr>1e-02);

index= union(union(index1,index2),index_tr);
index_setdiff=setdiff(1:noOfSensors,index);
  figure(222);
 if dim==2
     plot(PP(1,index),PP(2,index),'r*');
 else
     plot3(PP(1,index),PP(2,index),PP(3,index),'r*');
 end  
%PP_2round=PP(:,index);
DD_2=DD(index,index);
P0_2=[P0 Xopt0(:,index_setdiff)];
D0_2=[D0(index,:) DDD(index,index_setdiff)];
% Dall_2round=[D0_2round DD_2round];
% noOfAnchors_2round=noOfAnchors+length(index_setdiff);
% dij=[zeros(noOfAnchors_2round) D0_2round';
%     zeros(size(D0_2round)) DD_2round];
NumOfEdge=full(sum(sum(D0_2~=0))+sum(sum(DD_2~=0)));   

% noOfSensors_2round=length(index);
%pars.tol_gk =1e-6;  % Stopping criterion. 
%pars.tol_fk=1e-11;   % Stopping criterion.
fprintf('@@@@  Starting Second Localization! \n');
fprintf('####  the number of distance (NumOfEdge) = %3d \n',NumOfEdge);
pars.showyes = 0;

[R_2,Q] = NLP_CG_SNL(P0_2,D0_2,DD_2,outputdim,randSeed,pars);
R(index,:) = R_2;

