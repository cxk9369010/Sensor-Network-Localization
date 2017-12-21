function [d_ji,djk,N_x,N_a]=initseries(D0,DD,npts,apts)
%%
%                        Generating N_x and N_a (connection matrix)
%%

[K1,H1]=find(DD>0);
   T=cat(2,K1,H1);
  DD=DD+DD'; 
   N_x=[];
   for i=1:npts
       indx=H1==i;
       k=T(indx,:);
       N_x=[N_x;k];
   end

   d_ji=DD((N_x(:,1)-1)*npts+N_x(:,2)).^2;

   
[K2,H2]=find(D0>0);

K=cat(2,K2,H2); 
   N_a=[];
   for i=1:apts
       indx= H2==i;
       k=K(indx,:);
       N_a=[N_a;k];
   end

   djk=D0((N_a(:,2)-1)*npts+N_a(:,1)).^2;


end