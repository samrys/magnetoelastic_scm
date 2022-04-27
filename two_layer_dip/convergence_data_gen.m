%Generates data for convergence figure

%Initialize arbitrary precision
addpath('C:\Users\samry\Documents\Multiprecision Computing Toolbox\')
p = 34; %Choose precision
mp.Digits(p); %Quadruple precision p = 34

k_init = .2;

%Choose initial frequency
%f_init = 137.2174;
% f_init = 167.7488;
% f_init = 300.5974;
% f_init = 249.8971;
% f_init = 199.5158;
f_init = 341.3463;
%f_init = 185.1;

k = k_init*10^9;
w = 2*pi*10^9*f_init;

save_data = ['conv_data_v2_',num2str(k/10^6,'%0.f'),'_f_',...
    num2str(w/(2*pi*10^9),'%0.f')];

%Elastic parameters for the top material (YIG)
par.cs = mp('3800'); %transverse wave speed m/sec
par.cl = mp('7200'); %Longitudinal wave Speed m/sec
par.rho = mp('5170'); %Density kg/m^3
par.ell = mp('30*10^-9');
par.c11 = par.rho*par.cl^2;
par.c12 = par.rho*(par.cl^2-2*par.cs^2);
par.c44 = par.rho*par.cs^2;

%Magnetic parameters for the top material
par.h0 = mp('.25');     %applied field
par.ms = mp('140000'); %saturation magnetization
par.lam = mp('3')*mp('10^-16'); %exchange length squared
par.mu0 = 4*mp('pi')*mp('10^-7'); %Vacuum permeability
par.g = mp('1.760859644')*mp('10^11'); %gyromagnetic ratio of an electron
par.b2 = mp('7000000'); %Magnetoelastic coupling constant


%Physical parameters for the bottom material (GGG)
par.CL = mp('6400');
par.CS = mp('3500');
par.RH = mp('7080');
par.C11 = par.RH*par.CL^2;
par.C12 = par.RH*(par.CL^2-2*par.CS^2);
par.C44 = par.RH*par.CS^2;
par.ELL = mp('50*10^-9');

Ntest = 5:20;
freq_vec = mp(zeros(size(Ntest)));
freq_conv_err = mp(zeros(size(Ntest)));
count = zeros(size(Ntest));
count(1) = nan;

freq_vec(1) = w/(2*mp('pi'));
freq_conv_err(1) = NaN;

ell = par.ell;
ELL = par.ELL;

parfor hh = 2:length(Ntest)
 
       Nh = Ntest(hh);
     
        [Dh,~] = cheb(Nh-1);
        Dth = 2*Dh./ell;
        Dbh = 2*Dh./ELL;
        Dt2h = Dth*Dth;
        Db2h = Dbh*Dbh;
        
 
        %Create field matrix and polynomial eigenvalue problem
        [Ah,Bh,Ch]=cheb_mat_dip(k,Dth,Dbh,Dt2h,Db2h,Nh,par);
        wh = polyeig(Ah,Bh,Ch);
        [~,hind] = min(abs(real(wh)-w));
        w2 = wh<Inf;
        w3 = (real(wh)/(2*pi*10^9))>0;
        count(hh) = sum(w2.*w3);
        wout = wh(hind);
        freq_vec(hh) = real(wout)/(2*mp('pi'));
        
end   
% 
for hh = 2:length(Ntest)
        freq_conv_err(hh) = abs(freq_vec(hh)-freq_vec(hh-1))/freq_vec(hh-1);
end
     data_dir = ['C:\Users\samry\Documents\MATLAB\Elastic\magnetoelastics_paper\data_v2\',save_data,'.mat'];
 %   save(data_dir,'Ntest','k','w','k_init','f_init','freq_vec','freq_conv_err','count');
