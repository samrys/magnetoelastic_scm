function [A,B,C] =cheb_mat_one_layer_simp(k,Dt,Dt2,N,par)

cl2 = par.cl.^2;
cs2 = par.cs.^2;

b2 = par.b2;
ms = par.ms;
co1 = b2/(par.rho.*ms);
co2 = par.g.*b2;

z = zeros(N);

a33 = cl2*Dt2-k^2*cs2*eye(N);
a34 = co1.*(1i.*k.*eye(N));
a45 = par.g*par.h0.*eye(N)+par.beta.*(k^2.*eye(N)-Dt2);


a53 = -1i*co2.*k.*eye(N);
a54 = -a45;


A = [ a33 a34 z; ...
     z z a45; a53 a54 z ];


dv = Dt(1,:);
id = eye(N);
ev = 1i*k*id(1,:);
ew = 1i*k*id(N,:);
dw = Dt(N,:);
zv = z(1,:);
mv = id(1,:);
mw = id(N,:);

%Zero stress on top
%A(1,:) = [c44*ev b2/ms.*mv zv];
A(1,:) = [dv zv zv];

%Neumann conditions for magnetism
A(N,:) = [zv dv zv];
A(N+1,:) = [zv dw zv];
A(2*N,:) = [zv zv dv];
A(2*N+1,:) = [zv zv dw];

%zero stress at the bottom
%A(3*N,:) = [c44*ew b2/ms.*mw zv];
A(3*N,:) = [dw zv zv];

m = id;
m(1,1) = 0;
m(N,N) = 0;

%Magnetism equations have omega
B = 1i.*blkdiag(z,m,m);

%Elasticity equations include omega^2
C = blkdiag(m,z,z);


