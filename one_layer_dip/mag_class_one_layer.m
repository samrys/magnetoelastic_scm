function [class_vec] = mag_class_one_layer(Xv,e,kv,N,par)
%Classifies the type of eigenvector found

%Inputs:
%    Xv = Matrix of eigenvectors
%    e = vector of eigenvalues

%Outputs:
%   class_vec = vector with wave classification
%       1: longitudinal
%       2: magnetic
%       3: transverse


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
    if real(w)<10000 && isinf(w)==0
        class_vec(jj) = 0;
    else        
        X = Xv(:,jj)
    v1 = C*norm(X(1:N));        %A1 top material
    v2 = C*norm(X(N+1:2*N));    %A2 top material
    v3 = C*norm(X(2*N+1:3*N));  %A3 top material
    vm = norm(sqrt(real(X(3*N+1:4*N)).^2+imag(X(3*N+1:4*N)).^2+real(X(4*N+1:5*N)).^2+imag(X(4*N+1:5*N)).^2));%Magnetic elements
    vt = [v1 v2 v3 vm]; 

    %Energy concentration classification
    [~,z] = max(vt);
     if z==1
        class_vec(jj) = 1;
     elseif z==4
        class_vec(jj) = 2;
     else
        class_vec(jj) = 3;
     end
    end
end

