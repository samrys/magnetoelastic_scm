  function [D,x] = cheb(N)
  if N==0, D=0; x=1; return, end
  x = cos(mp('pi')*(0:N)/N)'; 
  c = [mp('2'); mp(ones(N-1,1)); mp('2')].*mp((-1).^(0:N)');
  X = mp(repmat(x,1,N+1));
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(mp(eye(N+1))));      % off-diagonal entries
  D  = D - diag(sum(D'));                 % diagonal entries