function [odeData,st,en,Peter,Uf_idx,Ropt] = Peter_aerodynamic_module(structural,panelc,sim,segc,N,aero,R_OC,h_sign,aero_ON_OFF,idx_seg_input)
% clear panel.PSIa panel.PHIh panel.PHIh_dl panel.PSIa_dl
% keyboard
IN_min = 8;
IN_max = 11;
narginchk(IN_min,IN_max);
laser_loc_span_global = sum(segc.l)*sim.lv;
if nargin <= IN_min
    
    cum_length     = cumsum(segc.l);
    idx_seg        = find(cum_length >= laser_loc_span_global,1); %non dimensional index location
    idx_start      = 1; %75
    idx_end        = numel(panelc.rhoA); %309
    Ltot           = sum(segc.l); %m, length of this segment
    laser_loc_span = laser_loc_span_global;
else
    nsON = nnz(aero_ON_OFF);
    if nsON > numel(segc.l)
        error('segments to apply aerodynamics must be <= total structural segments')
    end
%     laser_loc_span_global = sum(segc.l)*sim.lv;
% keyboard
    idx_seg         = idx_seg_input;
    cum_length      = cumsum(segc.l);
    idx_start       = sum(segc.Npanel(1:idx_seg - 1)) + 1; %75
    idx_end         = idx_start - 1 + segc.Npanel(idx_seg); %309
    Ltot            = segc.l(idx_seg); %m, length of this segment
    laser_loc_span  = laser_loc_span_global - cum_length(idx_seg - 1); %laser_loc_span local within the segment
end
% keyboard

%% INPUTS
% N        = 6; %total aerodynamic panels
Peter    = Peter_states(N); %
Peter.N  = N;

rho      = aero.rho;
T        = sim.T;

X0_disp  = sim.X0_disp;
X0_velo  = sim.X0_velo;

%% NON INPUTS
if aero.Np > idx_end - idx_start + 1
    Np      = idx_end - idx_start + 1;%500
    aero.Np = Np;
else
    Np = aero.Np;%500
end

idx           = round(linspace(idx_start,idx_end,Np)); %vector, index values of all panels included for applying aerodynamic loads

panel.a       = panelc.a(idx);
panel.b       = panelc.b(idx);


% keyboard
V_panel = round(laser_loc_span*numel(idx)/Ltot);
laser_c = panelc.chord(V_panel)*sim.v;
% keyboard
nLam      = N*Np;
aero.N_Np = nLam;


aero.idx       = idx;
sim.modes      = sim.nbending + sim.ntorsion;


PHI = structural.bending.PHI(idx,:);
PSI = structural.torsion.PSI(idx,:);


panel.PHIh    = PHI;
panel.PSIa    = PSI;


panel.PHIh_dl = structural.bending.PHI(idx,:).*panelc.dr(idx,:);
panel.PSIa_dl = structural.torsion.PSI(idx,:).*panelc.dr(idx,:);

panel.dPHI_dx   = structural.bending.dPHI_dx;
panel.d2PHI_dx2 = structural.bending.d2PHI_dx2;

panel.dPSI_dx   = structural.torsion.dPSI_dx;
panel.d2PSI_dx2 = structural.torsion.d2PSI_dx2;

panel.Et        = panelc.Et;
panel.dr        = panelc.dr;

Nb     = 100; %disxretization in chorwise direction

ylower = (1 - panelc.e(1))*panelc.b(1);
yupper = (1 + panelc.e(1))*panelc.b(1);

y      = (linspace(ylower, yupper, Nb))';

panel.y   = y;
dy        = (y(2) - y(1))*ones(numel(y),1);

int_dy    = sum(dy);    %integral of dy
int_y_dy  = (dy'*y);    %integral of ydy
int_y2_dy = (dy'*y.^2); %integral of y^2dy
int_y3_dy = (dy'*y.^3); %integral of y^3dy
int_y4_dy = (dy'*y.^4); %integral of y^4dy

% keyboard
panel.int.int_dy    = int_dy;
panel.int.int_y_dy  = int_y_dy;
panel.int.int_y2_dy = int_y2_dy;
panel.int.int_y3_dy = int_y3_dy;
panel.int.int_y4_dy = int_y4_dy;

panel.L  = sum(panelc.dr);
panel.t0 = (0.5.*panel.Et./panel.L);
panel.dy = dy;

bmax     =  panel.b(1).*(1 + panel.a(1));
bmin     = -panel.b(1).*(1 - panel.a(1));

panel.x        = fliplr(linspace(bmin,bmax,Nb));
panel.dx       = gradient(panel.x);
panel.Nb       = Nb;
panel.T0_c     = panelc.T0_c;


NLh_geom = 0;
NLa_geom = 0;

Peter.A_inv = inv(Peter.A);        %NxN
Peter.A_c   = Peter.A_inv*Peter.c; %Nx1

% keyboard

%% Time simulation
M_heave     = structural.matxRR.M_bend;
K_heave     = structural.matxRR.K_bend_EI + structural.matxRR.K_bend_T(T);

M_torsional = structural.matxRR.M_torsion;
K_torsional = structural.matxRR.K_tor_Gg + structural.matxRR.K_tor_T(T);

matrix.K = blkdiag(K_heave,K_torsional);
matrix.M = blkdiag(M_heave,M_torsional);

M_struct = [          M_heave        ,structural.Icoup.Tor_Bend';
            structural.Icoup.Tor_Bend,        M_torsional       ];
        
K_struct = blkdiag(K_heave, K_torsional);

Cp.L             = structural.CpL;
Cp.R             = structural.CpR;
theta_coupling.L = structural.theta_coupL ;
theta_coupling.R = structural.theta_coupR ;
% keyboard
%% Structural damping accounding for added air mass too

alphaM = 0.002;
betaK  = 0.005;


if strcmp(sim.include_damping,'N') == 1
    C_h(1:sim.nbending,1:sim.nbending) = zeros;
    C_a(1:sim.ntorsion,1:sim.ntorsion) = zeros;
else
    [C_param_h, C_param_a] = PTFE_AR18_Damping_Data(0.5);
%     M_aero_h = M_aero(1:sim.nbending,1:sim.nbending);
    [C_h,alphah,betah] = Damping_struct_exp(M_heave,K_heave,0*1,C_param_h);
%     M_aero_a = M_aero(sim.nbending+1:sim.nbending+sim.ntorsion,sim.nbending+1:sim.nbending+sim.ntorsion);
    [C_a,alphaa,betaa] = Damping_struct_exp(M_torsional,K_torsional,0*1,C_param_a);
end
% C_struct = 1.*double(blkdiag(C_h,C_a));
C_struct = alphaM.*M_struct + betaK.*K_struct;
% C_struct = structural.C_struct;%(blkdiag(C_h,C_a));
% keyboard

%% Initial conditions
st.rh   = 1;
en.rh   = st.rh - 1 + sim.nbending;
st.ra   = en.rh + 1;
en.ra   = st.ra - 1 + sim.ntorsion;

st.rh_d = en.ra + 1; %st.Vh   = sim.modes + 1;
en.rh_d = st.rh_d - 1 + sim.nbending;
st.ra_d = en.rh_d + 1;
en.ra_d = st.ra_d - 1 + sim.ntorsion;

st.lamb = en.ra_d + 1; %st.lamb = en.Va + 1;
en.lamb = st.lamb - 1 + nLam; %en.lamb = en.Va + nLam;

st.qL   = en.lamb + 1;
st.qR   = st.qL + 1;
% keyboard

flag   = 0;
flag_c = 0;
Xr     = [];

% keyboard
[Xr0,h0_chk,a0_chk] = Initial_condition_X0_V2(N,Ltot,Np,PHI(V_panel,:),PSI(V_panel,:),X0_disp,X0_velo,sim.nbending,sim.ntorsion,1); %sim,N,Np,PHI_h_x,PSI_a_x,x, odeS,loc_x,flag
% (sim,N,L,Np,PHI,PSI,X0_disp,X0_velo,nh,na,flag)
% keyboard
Force_RHS = zeros(sim.nbending + sim.ntorsion,1);
% h0_chk
% a0_chk

% keyboard

Xr0_ode = [Xr0,0,0]; %Extra zero to account for charge q [Xr0,0,0]

lambda0_alloc  = zeros(1,Np);
Ueff_alloc     = zeros(1,Np);
tic
disp(['Wind range = ',num2str(sim.U_speed(1)),' - ',num2str(sim.U_speed(end)),' m/s',])
Uf = 0;
Uf_out = [];
tafter = 0;
for i = 1:numel(sim.U_speed)
    Ui    = sim.U_speed(i);
    
    M     = double(M_struct );
    K     = double(K_struct );
    C     = double(0*C_struct);
    
    matrices.M = M;
    matrices.C = C;
    matrices.K = K;
    
%     keyboard
    tim = sim.tspan_L;
    t   = tim;
%     t1  = t(1:round(numel(t)/2));
%     keyboard                                                 f,t,panel,U,matrices,aero,sim,st,en,Force_RHS,Peter,theta_coupling,Cp,R_OC,lambda0_alloc,Ueff_alloc,h_sign,X0
    [Xr,F_track]        = ode4_pchatte2_PetersLCO(@Peter_3D_ode4_V6,sim.tspan_L,panel,Ui,matrices,aero,sim,st,en,Force_RHS,Peter,theta_coupling,Cp,R_OC,lambda0_alloc,Ueff_alloc,h_sign,Xr0_ode);
    [Div, h_V1,~,th_V1] = check_amplitude_growth(sim.tspan_L,Xr,[sim.nbending, sim.ntorsion, 0],st,en, PHI,PSI,panel.a,panel.b,V_panel,laser_c,Ui,flag,0,h_sign) ; %(tim,Xr,PHI_x,PSI_x,sim,geom,Uf,st,en,sim.U_speed(i),flag,x)
    figure
    subplot(2,1,1)
    plot(sim.tspan_L,h_V1.*1000,'.b','LineWidth',1);
    title(['U_\infty=',num2str(Ui),'[m/s]'])
    grid on
    
    subplot(2,1,2)
    plot(sim.tspan_L,rad2deg(th_V1),'.r','LineWidth',1);
    grid on
%     keyboard
    
    [~,~,fp1]   = get_fft(t,h_V1);
    if isempty(fp1) == 1
        Ropti = 0;
    else
        Ropti = 1/(2*pi*fp1*Cp.L);  Ropti(isnan(Ropti)) = 0;
    end
     Ropt(i)   = Ropti;
    [Xr,F_track2,alpha_eff,theta,Ueff_mag,lambda0] = ode4_pchatte2_PetersLCO(@Peter_3D_ode4_V6,sim.tspan_NL,panel,Ui,matrices,aero,sim,st,en,Force_RHS,Peter,theta_coupling,Cp,Ropti,lambda0_alloc,Ueff_alloc,h_sign,Xr0_ode);%Ropti
%     keyboard Xout,Force_track,a_eff,theta,Ueff,lamda0
    [Div, h_V,tnew,th_V] = check_amplitude_growth(sim.tspan_NL,Xr,[sim.nbending, sim.ntorsion, 0],st,en, PHI,PSI,panel.a,panel.b,V_panel,laser_c,Ui,flag,tafter,h_sign) ; %(tim,Xr,PHI_x,PSI_x,sim,geom,Uf,st,en,sim.U_speed(i),flag,x)
    
    subplot(2,1,1)
    hold on
    plot(sim.tspan_NL,h_V.*1000,'-b');
    ylabel('h_P[mm]')
    
    subplot(2,1,2)
    hold on
    plot(sim.tspan_NL,rad2deg(th_V),'-r');
    ylabel('\theta[^\circ]')
    xlabel('time[s]')
%     keyboard
    
    [fx,P,fp2] = get_fft(tnew,h_V);
    fprintf('OC(%3.2f KOhm) freq. = %3.4f [Hz] and after Ropt(%3.2f KOhm) = %3.4f [Hz]\n',R_OC/1e3,fp1,Ropti/1e3,fp2);
%     keyboard
    
    if  Div == 1 || flag == 1 || flag_c > 0
        flag    = 1; 
        flag_c  = flag_c + 1; 
        Uf      = sim.U_speed(i);
        Uf_out  = [Uf_out,Uf];
%         Xr0     = double(Xr(end,:));
%         Xr0_ode = Xr0;
    end
%     plot(t,Xr(:,1))
%     keyboard
    odeData(i).time_U    = tnew(1:end-1);
    odeData(i).Xr_U      = Xr;
    odeData(i).alpha_eff = alpha_eff;%(:,1:end-1)
    odeData(i).Ueff_mag  = Ueff_mag; %(:,1:end-1)
    odeData(i).lambda0   = lambda0;
    odeData(i).theta     = theta;%(:,1:end-1)
    odeData(i).ffx       = fx;
    odeData(i).ffy       = P;
    cprintf('string','%d. ODE solution Completed at U = %3.2f m/s\n',i,sim.U_speed(i));

end

Uf_idx = find(sim.U_speed == Uf_out(1));

toc



end