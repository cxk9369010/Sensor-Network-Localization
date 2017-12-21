function val = mytrace(A,B);

if (nargin == 1)
   val = sum(diag(A));
elseif (nargin == 2)
   val = sum(sum(A.*B));
end
return