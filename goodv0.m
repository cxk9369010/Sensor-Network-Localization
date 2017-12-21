function [R,Q]=goodv0(npts,outputdim,P0,D0,DD,dim,randSeed);

%%
%                        Generating Initial Point
%%
rand('seed',randSeed);

Rtemp=zeros(dim,npts);
denote=zeros(1,npts);


xx=0.5*min(P0');
yy=0.5*max(P0');

DD=DD+DD';
set0=intersect(find(sum(DD)==0),find(sum(D0)==0));
if ~isempty (set0);
    fprintf('the graph is not connected!\n');
    %error('the graph is not connected!');
    Rtemp(:,set0)=(xx'+yy')*ones(1,length(set0));
    denote(set0)=ones(1,length(set0));
end



   for i=1:npts; 
       INDEX0=find(D0(:,i)~=0);
       if ~isempty(INDEX0);
        [a,b]=min(D0(INDEX0,i));
        Rtemp(:,i)=P0(:,INDEX0(b)); 
        denote(i)=INDEX0(b);
          %V0(i+sDim,:)=zeros(1,outputdim);  
          %V0(i+sDim,:)=rand(1,outputdim); 
          %V0(i+sDim,:)=0.5*ones(1,outputdim);
           
       end
   end
   INDEX1=find(denote==0);
   
while  ~isempty(INDEX1);
    INDEX11=INDEX1;
   for j=1:length(INDEX1)
        INDEX2=find(DD(:,INDEX1(j))~=0);
        INDEX0=find(denote(INDEX2)~=0);
        index=INDEX2(INDEX0);
%        KKKK
       if ~isempty(INDEX0);
        [a,b]=min(DD(index,INDEX1(j)));
        Rtemp(:,INDEX1(j))=P0(:,denote(index(b))); 
          %V0(i+sDim,:)=zeros(1,outputdim);  
          %V0(i+sDim,:)=rand(1,outputdim); 
          %V0(i+sDim,:)=0.5*ones(1,outputdim);
          denote(INDEX1(j))=denote(index(b));
       end
   end
  INDEX1=find(denote==0);
  if length(INDEX1)==length(INDEX11)
      fprintf('some sensor point is not connected, directly or indirectly, to an anchor point.\n');  
      Rtemp(:,INDEX1)=(xx'+yy')*ones(1,length(INDEX1));
      denote(INDEX1)=ones(1,length(INDEX1));
      break
  end
      
end
       %V0(K2(i,1)+sDim,:)=P0(K2(i,2),:); 
       %V0(K2(i,1)+2,:)=[0 0 0];
   %V0(3:4,1:2)=[-0.2,0.1;0.1 0.4];
   Qtemp=orth(rand(outputdim,outputdim));
  Q=Qtemp(1:dim,1:outputdim);
   %v=eye(sDim);
    %V0(1:sDim,1:outputdim)=v;
   %V0(1:sDim,1:sDim)=eye(sDim);
    R=(Q\Rtemp)';
   %V0(sDim+1:end,1:outputdim)=R';
   
   


end