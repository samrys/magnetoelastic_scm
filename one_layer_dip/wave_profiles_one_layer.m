function [s,fnew] = wave_profiles_one_layer(f,k,par,N,xtop)
%Function that generates a vector with the structure of a mode given by 
% the wavenumber k and the frequency f
%
%Inputs:
%   f   = approximate frequency (Hz; not angular frequency)
%   k   = wavenumber
%   par = structure with physical parameters
%   N   = number of Chebyshev modes per layer
%   xtop = vector of points to project the function onto
%
%Outputs:
%   s = structure with wave profiles
%   fnew = the actual frequency at which the eigenvalue was found (Hz)

%Generate Chebyshev matrices
[D,~] = cheb(N-1);
Dt = 2*D./par.ell;
Dt2 = Dt*Dt;

%Convert to angular frequency
w_approx = 2*pi*f;

%Create field matrix and polynomial eigenvalue problem
[A,B,C]=cheb_mat_dip_one_layer(k,Dt,Dt2,N,par);
[Xv,eigval] = polyeig(A,B,C);

%Find mode closest to given frequency
[~,ind]=min(abs(real(eigval)-w_approx));
w = real(eigval(ind));

X = Xv(:,ind);

ms = par.ms;
beta = par.g.*par.mu0.*par.lam.*ms;
rho = par.rho;
cs = par.cs;
g = par.g;
Astar = sqrt(ms*beta/(rho*g*cs^2));
C = ms./Astar;

X(1:3*N) = C*X(1:3*N);        %Rescale elastic elements

s.a1 = cheby_int(X(1:N)',xtop);
s.a2 = cheby_int(X(N+1:2*N)',xtop);
s.a3 = cheby_int(X(2*N+1:3*N)',xtop);
s.m1 = cheby_int(X(3*N+1:4*N)',xtop);
s.m2 = cheby_int(X(4*N+1:5*N)',xtop);

s.m = sqrt(real(s.m1).^2+real(s.m2).^2+imag(s.m1).^2+imag(s.m2).^2);

%Sometimes the wave profiles have the interesting part in the complex part
%of the eigenvector
if norm(imag(s.a1))>norm(real(s.a1))
    s.a1 = imag(s.a1);
else
    s.a1 = real(s.a1);
end
if norm(imag(s.a2))>norm(real(s.a2))
    s.a2 = imag(s.a2);
else
    s.a2 = real(s.a2);
end
if norm(imag(s.a3))>norm(real(s.a3))
    s.a3 = imag(s.a3);
else
    s.a3 = real(s.a3);
end

%Normalize the output vectors so the largest element = 1
mm = max(abs([s.a2 s.a3 s.m]));
% s.a1 = real(s.a1)/mm;
s.a2 = sign(s.a2(end)).*s.a2/mm;
s.a3 = sign(s.a3(end)).*s.a3/mm;
s.m = s.m/mm;


fnew = w/(2*pi);
end

