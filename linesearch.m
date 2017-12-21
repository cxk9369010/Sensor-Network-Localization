function [alpha]=linesearch(a,b,c,d,c1)
 
phi=[a b c d-c1*d 0];

rphi=roots(phi);
k=length(rphi);
alphahi=eps;
for i=1:k;
    t=isreal(rphi(i));
    if t==1&&rphi(i)>0&&alphahi<rphi(i);
         alphahi=rphi(i); 
    end
end

dphi=[4*a 3*b 2*c d];
rdphi=roots(dphi);

k=length(rdphi);
alpha=alphahi;
for i=1:k;
    t=isreal(rdphi(i));
    if t==1&&rdphi(i)>0&&rdphi(i)<=alphahi;
       alpha=rdphi(i); 
    end
end
%if d>=0
%      fprintf('wolfe_exact search: Error: Need a descent direction  \n');
%    alpha=0; 
  
% end
end







