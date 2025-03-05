
clc;
clear;
%%
%moduli
K_gr=56*10^9; % bulk modulus of rock matrix
mu_gr=33*10^9; % shear modulus of rock matrix
K_f=3.8100*10^8; % bulk modulus of fluid

%aspect ratio 
r_p=0.9999;%aspect ratio of stiff pore
r_c=0.0001;%aspect ratio of crack

for i = 1:4
for j = 1:3
for ii=1:1:1000
if i==1||i==3
%porosity
if i==1
phi=0.005;
else
phi=0.2;
end
epsilon=0.1;%crack density
phi_c=epsilon*4*pi*r_c/3; % crack content
phi_p=phi-phi_c; % stiff pore content
% fluid pressure
% A_w: stiff pore pressure; B_w: crack pressure.
A_w = 0.001*10^(j-1);
B_w = (ii-1)*0.0003;
elseif i==2||i==4
%porosity
if i==2
phi=0.005;
else
phi=0.2;
end
epsilon=0.1;
phi_c=epsilon*4*pi*r_c/3; 
phi_p=phi-phi_c; 
% fluid pressure
A_w = (ii-1)*0.0001;
B_w = 0.1*j;
end  

% 6*6 identity matrix I
I=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1];

%stiffness tensor L0 of the matrix
L0_11=K_gr+4/3*mu_gr;
L0_12=K_gr-2/3*mu_gr;
L0_44=mu_gr;
L0=[L0_11, L0_12, L0_12, 0, 0, 0;
    L0_12, L0_11, L0_12, 0, 0, 0;
    L0_12, L0_12, L0_11, 0, 0, 0;
    0, 0, 0, L0_44, 0, 0;
    0, 0, 0, 0, L0_44, 0;
    0, 0, 0, 0, 0, L0_44];

% Sp_ani and Sc_ani are Eshelby's tensors for pores and cracks respectively based on:David E C, Zimmerman R W 2011 
% Sp_ani is the function of r_p and Poisson's ratio of matrix.
Sp_ani = 0.1*[5.1198    0.2399    0.2397         0         0         0
          0.2399    5.1198    0.2397         0         0         0
          0.2399    0.2399    5.1203         0         0         0
          0         0         0    4.8801         0         0
          0         0         0         0    4.8801         0
          0         0         0         0         0    4.8799];
% Sc_ani is the function of r_c and Poisson's ratio of matrix.
Sc_ani = 0.1*[0.0014    0.0001   -0.0003         0         0         0
          0.0001    0.0014   -0.0003         0         0         0
          3.3984    3.3984    9.9995         0         0         0
          0         0         0    9.9982         0         0
          0         0         0         0    9.9982         0
          0         0         0         0         0    0.0013];
% dry rock modulus tensor
Ld_anp = phi_p*(I-Sp_ani)^(-1);
Ld_anc = phi_c*(I-Sc_ani)^(-1);
Ld_anp_oa = oa(Ld_anp); % orientation average rule, Dvorak (2012)
Ld_anc_oa = oa(Ld_anc);
Ld=L0-L0*(Ld_anp_oa+Ld_anc_oa)*((1-phi)*I+Ld_anp_oa+Ld_anc_oa)^(-1); 

% fluid modulus tensor
a=K_f;
b=K_f*10^(-7);
c=a-2*b;
Lv=[a c c 0 0 0;
    c a c 0 0 0;
    c c a 0 0 0;
    0 0 0 b 0 0;
    0 0 0 0 b 0;
    0 0 0 0 0 b];

% modulus tensor of saturated rocks
L1_anp = (I-Sp_ani)^(-1)*Sp_ani/L0;
L1_anc = (I-Sc_ani)^(-1)*Sc_ani/L0;
L1_anp_oa = oa(L1_anp);
L1_anc_oa = oa(L1_anc);  

L1 = phi_p*L1_anp_oa*A_w+phi_c*L1_anc_oa*B_w;
L2 = phi_p*A_w+phi_c*B_w;
L3 = I-L0/Lv;

L_lef=Ld*L1+L2*I;
L_rig=L3*L2;
L_sat1 = (L_lef/L_rig-I)\(L_lef/L_rig*L0-Ld);
L_sat2 = (L_rig/L_lef-I)\(L_rig/L_lef*Ld-L0);
L_sat3 = (L_rig-L_lef)\(L_rig*Ld-L_lef*L0);
L_sat = L_sat3;

% moduli for saturated rocks
mu_n(ii)=L_sat(4,4)/10^9;%shear modulus
M_n(ii)=L_sat(1,1)/10^9;%P-wave modulus

Stiff_P(ii) = A_w;
Crack_P(ii) = B_w;

end

%plot
if i==1||i==3
figure (1)
if i==1
subplot(221)
else
subplot(223)
end
if j==1
plot(Crack_P,M_n,'k--','LineWidth', 1.5); hold on
elseif j==2
plot(Crack_P,M_n,'k-.','LineWidth', 1.5); hold on
else
plot(Crack_P,M_n,'k','LineWidth', 1.5); hold on
end
xlabel('σ_c^r'); 
ylabel('M (Gpa)');
legend('σ_p^r=0.001', 'σ_p^r=0.01', 'σ_p^r=0.1','Box', 'off', 'FontSize', 10, 'Location', 'east');

elseif i==2||i==4
figure (1)
if i==2
subplot(222)
else
subplot(224)
end
if j==1
plot(Stiff_P,M_n,'k--','LineWidth', 1.5); hold on
elseif j==2
plot(Stiff_P,M_n,'k-.','LineWidth', 1.5); hold on
else
plot(Stiff_P,M_n,'k','LineWidth', 1.5); hold on
end
xlabel('σ_p^r'); 
ylabel('M (Gpa)'); 
legend('σ_c^r=0.1', 'σ_c^r=0.2', 'σ_c^r=0.3','Box', 'off', 'FontSize', 10, 'Location', 'northeast');
end
if i==1||i==2
text(0.05, 0.95, 'Por=0.0005', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'norm', 'Color', 'black');
else
text(0.05, 0.95, 'Por=0.2', 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'norm', 'Color', 'black');
end

end
end
