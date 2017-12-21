function [a,b,c,d,fval]=coefficient_SPARSE(Q,R,D,P0,d1,d2,K1,K2)

%%  Computing a,b,c,d,e
P0=P0';

R_R=R(K1(:,1),:)-R(K1(:,2),:);
D_D=D(K1(:,1),:)-D(K1(:,2),:);

p0=sum(R_R.^2,2)-d1;
p1=2*sum(R_R.*D_D,2);
p2=sum(D_D.^2,2);


p14=sum(p2.*p2);
p13=sum(p2.*p1);
p12=sum(p1.*p1)+2*sum(p2.*p0);
p11=sum(p1.*p0);
p10=sum(p0.*p0);

R_alphaQ=R(K2(:,1),:)-P0(K2(:,2),:)*Q;


p0=sum(R_alphaQ.^2,2)-d2;
p1=2*sum(R_alphaQ.*D(K2(:,1),:),2);
p2=sum(D(K2(:,1),:).^2,2);


p24=sum(p2.*p2);
p23=sum(p2.*p1);
p22=sum(p1.*p1)+2*sum(p2.*p0);
p21=sum(p1.*p0);
p20=sum(p0.*p0);


    

a=0.5*(p14+p24);
b=p13+p23;
c=0.5*(p12+p22);
d=p11+p21;

fval=0.5*(p10+p20);
%while d>0
%    error('it is not desent direction');
   
%end


end