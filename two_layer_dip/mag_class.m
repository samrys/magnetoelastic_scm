function [class_vec] = mag_class(Xv,e,kv,N,par)
%Classifies the type of eigenvector found

%Inputs:
%    Xv = Matrix of eigenvectors
%    e = vector of eigenvalues

%Outputs:
%   class_vec = vector with wave classification
%       1: A1 (top)
%       2: A2 (top)
%       3: A3 (top)
%	4: A1 (bottom)
%	5: A2 (bottom)
%	6: A3 (bottom)
%	7: magnetic

class_vec = zeros(size(e));

ms = par.ms;
beta = par.g.*par.mu0.*par.lam.*ms;
rho = par.rho;
cs = par.cs;
g = par.g;
Astar = sqrt(ms*beta/(rho*g*cs^2));
C = ms./Astar;

parfor jj = 1:length(e)
    w = e(jj);
    X = Xv(:,jj);
    if real(w)<10000 && isinf(w)==0
        class_vec(jj) = 0;
    else        
    v1 = C*norm(X(1:N));        %A1 top material
    v2 = C*norm(X(N+1:2*N));    %A2 top material
    v3 = C*norm(X(2*N+1:3*N));  %A3 top material
    vm = norm(sqrt(real(X(3*N+1:4*N)).^2+imag(X(3*N+1:4*N)).^2+real(X(4*N+1:5*N)).^2+imag(X(4*N+1:5*N)).^2));%Magnetic elements
    v4 = C*norm(X(5*N+1:6*N));
    v5 = C*norm(X(6*N+1:7*N));
    v6 = C*norm(X(7*N+1:8*N));
    vt = [v1 v2 v3 v4 v5 v6 vm]; 

    %Energy concentration classification
    [~,z] = max(vt);
    class_vec(jj) = z;
    end
end
end

