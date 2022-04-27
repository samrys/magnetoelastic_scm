function [s,fnew] = wave_profiles(f,k,par,N,xtop,xbot)
%Function that generates a vector with the structure of a mode given by 
% the wavenumber k and the frequency f
%
%Inputs:
%   f   = approximate frequency (Hz; not angular frequency)
%   k   = wavenumber
%   par = structure with physical parameters
%   N   = number of Chebyshev modes per layer
%   xtop = spatial variable for top layer
%   xbot = spatial variable for bottom layer

%Outputs:
%   s = structure with relevant wave profiles
%   fnew = actual frequency values determined via eigenvalue problem (Hz)

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


s.a1 = cheby_int(X(1:N)',xtop);
s.a2 = cheby_int(X(N+1:2*N)',xtop);
s.a3 = cheby_int(X(2*N+1:3*N)',xtop);
s.A1 = cheby_int(X(5*N+1:6*N)',xbot);
s.A2 = cheby_int(X(6*N+1:7*N)',xbot);
s.A3 = cheby_int(X(7*N+1:8*N)',xbot);
s.m1 = cheby_int(X(3*N+1:4*N)',xtop);
s.m2 = cheby_int(X(4*N+1:5*N)',xtop);

s.m = sqrt(real(s.m1).^2+real(s.m2).^2+imag(s.m1).^2+imag(s.m2).^2);


%Occasionally, the action is in the complex coefficient
if norm(imag([s.a1 s.A1]))>norm(real([s.a1 s.A1]))
    s.a1 = imag(s.a1);
    s.A1 = imag(s.A1);
else
    s.a1 = real(s.a1);
    s.A1 = real(s.A1);
end
if norm(imag([s.a2 s.A2]))>norm(real([s.a2 s.A2]))
    s.a2 = imag(s.a2);
    s.A2 = imag(s.A2);
else
    s.a2 = real(s.a2);
    s.A2 = real(s.A2);
end
if norm(imag([s.a3 s.A3]))>norm(real([s.a3 s.A3]))
    s.a3 = imag(s.a3);
    s.A3 = imag(s.A3);
else
    s.a3 = real(s.a3);
    s.A3 = real(s.A3);
end

%Rescale by max value
s.mm = max(abs([s.a2 s.a1 s.a3 s.A1 s.A2 s.A3 s.m]));
% s.a1 = real(s.a1)/mm;
s.a2 = sign(s.a2(end)).*s.a2/s.mm;
s.a3 = sign(s.a3(end)).*s.a3/s.mm;
s.A2 = sign(s.a2(1).*s.A2(end)).*s.A2/s.mm;
s.A3 = sign(s.a3(1).*s.A3(end)).*s.A3/s.mm;
s.m = s.m/s.mm;



fnew = w/(2*pi);
end

