% Program pre výpočet TF pomocou symb.MATLABu (MSP). I. Zeman, 27.3.2022
% Simulácia a znázornenie PrCH a LFCH viacslučkového obvodu HPF LCL
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Zmeny veľkosti popisov osi, farby a hrúbky čiary
% Výpočet núl a pólov pre LCL filter
% Simulácia a znázornenie PrCH a LFCH z blok. modelu v Simulinku

clear, clc, clf, format compact
syms s R1 R2 L1 L2 C Uin Uoutsim %deklarácia symbolických premenných

% Zadanie vstupných hodnôt
disp('Hornopriepustný LCL filter so záťažou - riešenie obvodu MSP')
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
Tstep=10e-3; wmin=1; wmax=1e5;                      % parametre pre Step a Bode
Tchirp=0.04; wminchirp=1e1; wmaxchirp=1e3;          % parametre pre Chirp
Tsim=0.02;                                          % doba simulácie v Simulinku (experimentálne určená)
color='r';                                          % farba grafu b,r,y,m,c,

% Záspis systému a výpočet TF v symbolickom tvare
Z=[R1+s*L1            -s*L1              0          %matica impedancií obvodu podľa MSP:
   -s*L1         1/(s*C)+s*L2+s*L1      -s*L2  
       0              -s*L2              R2+s*L2   ];
u=[Uin; 0; 0];                                        %vektor napätí obvodu
ZI3=[Z(:,1) Z(:,2) u];                              %submatica pre I3
I3=det(ZI3)/det(Z);                                 %výpočet slučkového prúdu I3 Cramerovým pravidlom 
Uout=R2*I3;                                          %výstupné napätie UR2 (OZ)
F=Uout/Uin                                             %nájdenie TF F(s) = Y/U v symbolickom tvare

% Spracovanie údajov TF v symbolickom tvare pre prechod do num. MATLABu
[cit,men]=numden(F);                     % oddelenie polynómov čitateľa a menovateľa
cit=subs(cit,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men=subs(men,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b=sym2poly(cit);                         % b - koeficienty polynómu čitateľa b(s)
a=sym2poly(men);                         % a - koeficienty polynómu menovateľa a(s)
b=double(b);                             % Prechod do numerickeho MATLABu
a=double(a); 
F=tf(b,a)                                % Výsledná TF v numerickom MATLABe
F=tf(b/a(end),a/a(end))                  % TF upravená pre a0=1 (normovanie TF)

figure(1)   % Vykresľovanie a popis priebehov PrCh a LFCh
subplot(121); step(F,Tstep,color), grid on,
   title('Prechodová charakteristika','FontSize',16)
     xlabel('\rightarrow T','FontSize',14)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',14)
subplot(122); bode(F,{wmin,wmax},color), grid on
     title('Frekvenčná charakteristika','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',14),ylabel('\rightarrow\phi','FontSize',14)
     set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary
         
%výpočet núl a pólov pre TF
roots([2.94e-10])                               %výpočet koreňov čitateľa
roots([3.234e-10 7.07e-07 0.00294 1])           %výpočet koreňov menovateľa