function [vout,fnew] = mode_structure(f,k,par,N)
%Function that generates a vector with the energy of a mode given by 
% the wavenumber k and the frequency f
%
%Inputs:
%   f   = approximate frequency (Hz; not angular frequency)
%   k   = wavenumber
%   par = structure with physical parameters
%   N   = number of Chebyshev modes per layer
%
%Outputs:
%   vout = vector with energies of the given mode

%Generate Chebyshev matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Db = 2*D./par.ELL;
Dt2 = Dt*Dt;
Db2 = Db*Db;

%Convert to angular frequency
w_approx = 2*pi*f;

%Create field matrix and polynomial eigenvalue problem
[A,B,C]=cheb_mat_dip(k,Dt,Db,Dt2,Db2,N,par);
[Xv,eigval] = polyeig(A,B,C);

%Find mode closest to given frequency
[~,ind]=min(abs(real(eigval)-w_approx));
w = real(eigval(ind));

%err = abs((w-w_approx)/w_approx)

X = Xv(:,ind);

ms = par.ms;
beta = par.g.*par.mu0.*par.lam.*ms;
rho = par.rho;
cs = par.cs;
g = par.g;
Astar = sqrt(ms*beta/(rho*g*cs^2));
C = ms./Astar;

X(1:3*N) = C*X(1:3*N);        %Rescale elastic elements
X(5*N+1:8*N) = C*X(5*N+1:8*N);

vout(2) = norm(X(1:N));        %A1 top material
vout(3) = norm(X(N+1:2*N));    %A2 top material
vout(4) = norm(X(2*N+1:3*N));  %A3 top material
vout(1) = norm(sqrt(real(X(3*N+1:4*N)).^2+imag(X(3*N+1:4*N)).^2+real(X(4*N+1:5*N)).^2+imag(X(4*N+1:5*N)).^2)); %Magnetic elements
vout(5) = norm(X(5*N+1:6*N));
vout(6) = norm(X(6*N+1:7*N));
vout(7) = norm(X(7*N+1:8*N));

vout = vout/norm(vout);    

fnew = w/(2*pi);
end

