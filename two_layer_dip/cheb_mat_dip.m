function [A,B,C] =cheb_mat_dip(k,Dt,Db,Dt2,Db2,N,par)

cl2 = par.cl.^2;
cs2 = par.cs.^2;
CL2 = par.CL.^2;
CS2 = par.CS.^2;
c11 = par.c11;
c12 = par.c12;
c44 = par.c44;
C11 = par.C11;
C12 = par.C12;
C44 = par.C44;
b2 = par.b2;
ms = par.ms;
co1 = b2/(par.rho.*ms);
co2 = par.g.*b2;

z = zeros(N);

a11 = cs2*Dt2-k^2*cl2*eye(N);
a13 = 1i*k*(cl2-cs2)*Dt;
a14 = co1.*Dt;
a22 = cs2.*(Dt2-k^2*eye(N));
a25 = a14;
a31 = a13;
a33 = cl2*Dt2-k^2*cs2*eye(N);
a34 = co1.*(1i.*k.*eye(N));
a42 = co2.*Dt;
a45 = par.g.*(par.h0-par.mu0*par.ms).*eye(N)+...
    par.g.*par.lam*par.mu0*ms.*(k^2.*eye(N)-Dt2);
a51 = -a42;
a53 = -1i*co2.*k.*eye(N);
a54 = -a45;




A11 = CS2*Db2-k^2*CL2*eye(N);
A13 = 1i*k*(CL2-CS2)*Db;
A22 = CS2.*(Db2-k^2*eye(N));
A31 = A13;
A33 = CL2*Db2-k^2*CS2*eye(N);

Top = [a11 z a13 a14 z ; z a22 z z a25 ; a31 z a33 a34 z ; ...
     z a42 z z a45 ; a51 z a53 a54 z  ];
 
Bot = [A11 z A13; z A22 z; A31 z A33];

%Dipole field
Dip =  Dt2 -k^2*eye(N);

A = blkdiag(Top,Bot,Dip);

%Also for dipole field
A(8*N+1:9*N,3*N+1:4*N) = -eye(N)*1i*k;
A(4*N+1:5*N,8*N+1:9*N) = -1i*par.g*par.mu0*par.ms*k*eye(N);

dv = Dt(1,:);
DV = Db(1,:);
id = eye(N);
ev = 1i*k*id(1,:);
ew = 1i*k*id(N,:);
dw = Dt(N,:);
DW = Db(N,:);
z1 = z(1,:);
z3 = [z1 z1 z1];
z5 = [z1 z1 z1 z1 z1];
mv = id(1,:);
mw = id(N,:);

%Zero stress on top
A(1,:) = [c44*dv z1 c44*ev b2/ms.*mv z1 z3 z1];
A(N,:) = [c12*ev z1 c11*dv z1 z1 z3 z1];
A(N+1,:) = [z1 c44.*dv z1 z1 b2/ms.*mv z3 z1];

%Continuity of displacements at the boundary
A(2*N,:) = [z1 mw z1 z1 z1 z1 -mv z1 z1];
A(2*N+1,:) = [mw z1 z1 z1 z1 -mv z1 z1 z1];
A(3*N,:) = [z1 z1 mw z1 z1 z1 z1 -mv z1];

%Neumann conditions for magnetism
A(3*N+1,:) = [z3 dv z1 z3 z1];
A(4*N,:) = [z3 dw z1 z3 z1];
A(4*N+1,:) = [z3 z1 dv z3 z1];
A(5*N,:) = [z3 z1 dw z3 z1];

%Continuity of stress at the boundary
A(5*N+1,:) = [c44*dw z1 c44*ew b2/ms.*mw z1 -C44*DV z1 -C44*ev z1];
A(6*N,:) = [c12*ew z1 c11*dw z1 z1 -C12*ev z1 -C11*DV z1];
A(6*N+1,:) = [z1 c44.*dw z1 z1 b2/ms.*mw z1 -C44.*DV z1 z1];

%Zero stress at the bottom
A(7*N,:) = [z5 DW z1  ew z1];
A(7*N+1,:) = [z5 C12*ew z1  C11*DW z1];
A(8*N,:) = [z5 z1 DW z1 z1]; 

%Neumann conditions for dipole field
A(8*N+1,:) = [z5 z3 dv+k*mv];
A(9*N,:) = [z5 z3 dw-k*mw];

m = id;
m(1,1) = 0;
m(N,N) = 0;

%Magnetism equations have omega
B = -1i.*blkdiag(z,z,z,m,m,z,z,z,z);

%Elasticity equations include omega^2
C = blkdiag(m,m,m,z,z,m,m,m,z);

