% Program pre výpočet a zobrazenie TF pomocou stavového modelu el. obvodu. I. Zeman, 25.4.202
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Riešenie sústavy 3 lin.alg.rovníc v symb.MATLABe
% pre daný elektrický obvod:
%
% R1.x1 + L1.dx1 + R1.C.dx2 = Uin
% x2 - L1.dx1 + L2.dx3 = 0
% -R2.x3 + R2.C.dx2 - L2.dx3 = 0
%
% Neznámymi sú premenné dx1, dx2, dx3 predstavujúce derivácie stavových veličín
% dx1=dx1/dt, dx2=dx2/dt, dx3=dx3/d
%
% Znázornenie PrCH a LFCH viacslučkového obvodu 
% Zmeny veľkosti popisov osi, farby a hrúbky čiary
clear, clc, clf, format compact

% parametre obvodu
syms x1 x2 x3 dx1 dx2 dx3;
syms R1 R2 L1 L2 C Uin
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
Tstep=5e-3; wmin=1e1; wmax=1e5;              % parametre pre Step a Bode
Tchirp=0.04; wminchirp=1e1; wmaxchirp=1e3;   % parametre pre Chirp
Tsim=0.05                                    % doba simulácie v Simulinku (experimentálne určená)
color='r';                                   % farba grafu b,r,y,m,c,

% Zápis rovníc obvodu
eq1 = x1*R1+L1*dx1+R1*C*dx2==Uin
eq2 = x2-L1*dx1+L2*dx3==0
eq3 = -R2*x3+R2*C*dx2-L2*dx3==0

disp('Riešenie')
[dx1,dx2,dx3] = solve(eq1,eq2,eq3,dx1,dx2,dx3)
pretty(dx1), pretty(dx2), pretty(dx3) % úprava výpisu zlomkov

%% Stavový model (prepísaný z výsledkov riešenia alg. rovníc a doplnený výstupnou rovnicou
disp('Stavový model v symbolickom tvare:')

A=[(-R1*R2)/(L1*R1+L1*R2)    R2/(L1*R1+L1*R2)     (-R1*R2)/(L1*R1+L1*R2)
 -R1/(C*R1+C*R2)            -1/(C*R1+C*R2)        R2/(C*R1+C*R2)  
 (-R1*R2)/(L2*R1+L2*R2)     -R2/(L2*R1+L2*R2)   (-R1*R2)/(L2*R1+L2*R2)]
b=[R2/(L1*R1+L1*R2); 1/(C*R1+C*R2); R2/(L2*R1+L2*R2)]
cT=[(-R1*R2)/(R1+R2) (-R2)/(R1+R2) (-R1*R2)/(R1+R2)]
d=[-R2/(R1+R2)]

%% Náhrada symb.premenných hodnotami
R1=R1x; R2=R2x; C=Cx; L1=L1x; L2=L2x;   Uin=Uinx; %
disp('Stavový po dosadení hodnôt parametrov:')
A=[(-R1*R2)/(L1*R1+L1*R2)    R2/(L1*R1+L1*R2)     (-R1*R2)/(L1*R1+L1*R2)
 -R1/(C*R1+C*R2)            -1/(C*R1+C*R2)        R2/(C*R1+C*R2)  
 (-R1*R2)/(L2*R1+L2*R2)     -R2/(L2*R1+L2*R2)   (-R1*R2)/(L2*R1+L2*R2)]
b=[R2/(L1*R1+L1*R2); 1/(C*R1+C*R2); R2/(L2*R1+L2*R2)]
cT=[(-R1*R2)/(R1+R2) (-R2)/(R1+R2) (-R1*R2)/(R1+R2)]
d=[R2/(R1+R2)]

%% Výstupy
disp('Výpis stavového modelu:')
 printsys(A,b,cT,d)
disp('Výpis prenosovej funkcie:')
 [num,den]=ss2tf(A,b,cT,d);
 F=tf(num/den(end),den/den(end))
disp('Vlastné hodnoty matice A:')
 format long
 eig(A)
disp('Póly prenosovej funkcie:')
 roots(den)
 format short

%vykreľovanie priebehov 
subplot(1,2,1),step(A,b,cT,d),grid on
     title('Prechodová charakteristika','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',14)
subplot(1,2,2),bode(A,b,cT,d),grid on
     title('Frekvenčná charakteristika','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',2) % inštrukcia pre zmenu hrúbky čiary

sim('Stavovy_model_S')