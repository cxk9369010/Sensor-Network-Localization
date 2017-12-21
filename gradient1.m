function[G]=gradient1(outputdim,R,Q,D0,DD,P0);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%         Compute the Gradient.
 P0=P0';
 D=DD+DD';
 D0=D0';
%D=DD;
 [apts,npts]=size(D0);
 %dim = size(V0,1);
  G0 = zeros(npts,outputdim);
  G1 = zeros(npts,outputdim);
  %F=0;
  if npts>1
  for j = 1:npts
      %xj = R(j,:); 
      %Dj = D(:,j);
      [idx,dummy,d] = find(D(:,j));
       tmp =ones(length(idx),1)*R(j,:) - R(idx,:);
      %tmp =ones(length(idx),1)*xj - PPtemp(idx+2,:)*V0(1:2,:);
      %nrm =(sum(tmp.*tmp,2))' ;
      %find(K1)
      alpha =(sum(tmp.*tmp,2))-d.^2;
      %F=F+alpha'*alpha;
      G0(j,:) = alpha'*tmp;
  end
  end

  
  for j = 1:npts
     %xj = R(j,:); 
      %Dj = D0(:,j);
      [idx,dummy,d] = find(D0(:,j));
  
      tmp =ones(length(idx),1)*R(j,:) - P0(idx,:)*Q;
      %nrm =(sum(tmp.*tmp,2))' ;
      %find(K1)
      alpha =(sum(tmp.*tmp,2))-d.^2;
      %F=F+alpha'*alpha;
    % alpha = (nrm-d.^2);
      G1(j,:) = alpha'*tmp;
  end
  
   G=2*(G0+G1);
   %F=0.5*F;
            
end
