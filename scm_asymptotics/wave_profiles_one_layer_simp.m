function [s,fnew] = wave_profiles_one_layer_simp(f,k,par,N,xtop)
%Function that generates a vector with the structure of a mode given by 
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
Dt2 = Dt*Dt;


%Convert to angular frequency
w_approx = 2*pi*f;

%Create field matrix and polynomial eigenvalue problem
[A,B,C]=cheb_mat_one_layer_simp(k,Dt,Dt2,N,par);
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

X(1:N) = C*X(1:N);        %Rescale elastic elements


s.a3 = cheby_int_test(X(1:N)',xtop);
s.m1 = cheby_int_test(X(N+1:2*N)',xtop);
s.m2 = cheby_int_test(X(2*N+1:3*N)',xtop);

s.m = sqrt(real(s.m1).^2+real(s.m2).^2+imag(s.m1).^2+imag(s.m2).^2);

if norm(imag(s.m1))>norm(real(s.m1))
    s.m1 = imag(s.m1);
else
    s.m1 = real(s.m1);
end
if norm(imag(s.m2))>norm(real(s.m2))
    s.m2 = imag(s.m2);
else
    s.m2 = real(s.m2);
end
if norm(imag(s.a3))>norm(real(s.a3))
    s.a3 = imag(s.a3);
else
    s.a3 = real(s.a3);
end

mm = max(abs([s.a3 s.m2]));
s.a3 = sign(s.a3(end)).*s.a3/mm;
%s.m1 = sign(s.m1(end)).*s.m1/mm;
s.m2 = sign(s.m2(end)).*s.m2/mm;


fnew = w/(2*pi);
end

