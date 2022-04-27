function [class_vec] = mag_class_one_layer_simp(Xv,e,kv,N,par)
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



parfor jj = 1:length(e)
    w = e(jj);
    if real(w)<10000 && isinf(w)==0
        class_vec(jj) = 0;
    else        
    vmg = par.g^2.*par.mu0^2.*kv.*par.lam.*par.ms.*(par.ms+2*par.lam*par.ms*kv^2-2*par.h0)/w;
    v3 = (abs(Xv(1:N,jj))*sqrt(par.cs*par.rho)*w).^2;  %A3 top material
    vm = (abs(Xv(N+1:3*N,jj))*sqrt(vmg*par.h0*.5/par.ms)).^2; %Magnetic elements
    vt = [sum(v3) sum(vm)]; 

    %Energy concentration classification
    [~,z] = max(vt);
     if z==100
        class_vec(jj) = 1;
     elseif z==2
        class_vec(jj) = 2;
     else
        class_vec(jj) = 3;
     end
    end
end

