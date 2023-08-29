% Program pre výpočet TF pomocou symb.MATLABu (MUN). I. Zeman, 27.3.2022
% Znázornenie PrCH a LFCH viacslučkového obvodu HPF LCL
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Zmeny veľkosti popisov osi, farby a hrúbky čiary

clear, clc, clf, format compact
syms s R1 R2 L1 L2 C Uin %deklarácia symbolických premenných

% Zadanie vstupných hodnôt
disp('Hornopriepustný LCL filter so záťažou - riešenie obvodu MSP')
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
Tstep=10e-3; wmin=1; wmax=1e5;               % parametre pre Step a Bode
color='r';                                   % farba grafu b,r,y,m,c,

% Záspis systému a výpočet TF v symbolickom tvare
Z=[(1/R1)+(s*C)+1/(s*L1)              -s*C                    %matica impedancií obvodu podľa MSP:
   -s*C                     (1/R2)+(s*C)+1/(s*L2)];
u=[Uin/R1; 0];                                 %vektor napätí obvodu
ZUA=[u Z(:,2)];                              %submatica pre UA
UA=det(ZUA)/det(Z);                          %výpočet uzlového napätia UA Cramerovým pravidlom 
ZUB=[Z(:,1) u ];                             %submatica pre UB
UB=det(ZUB)/det(Z);                          %výpočet uzlového napätia UB Cramerovým pravidlom 
I4=UB/R2;                                    %výpočet prúdu I4 (1.KZ)
Uout=R2*I4;                                   %výstupné napätie R2 (OZ)
F=Uout/Uin                                      %nájdenie TF F(s) = Y/U v symbolickom tvare

% Spracovanie údajov TF v symbolickom tvare pre prechod do num. MATLABu
[cit,men]=numden(F);                     % oddelenie polynómov čitateľa a menovateľa
cit=subs(cit,{R1,R2,L1,L2,C,Uin},{R1x,R2x,L1x,L2x,Cx,Uinx});        % dosadenie hodnôt do polynómu čitateľa
men=subs(men,{R1,R2,L1,L2,C,Uin},{R1x,R2x,L1x,L2x,Cx,Uinx});        % dosadenie do polynómu menovateľa
b=sym2poly(cit);                         % b - koeficienty polynómu čitateľa b(s)
a=sym2poly(men);                         % a - koeficienty polynómu menovateľa a(s)
b=double(b);                             % Prechod do numerickeho MATLABu
a=double(a);
F=tf(b,a)                                % Výsledná TF v numerickom MATLABe
F=tf(b/a(end),a/a(end))                  % TF upravená pre a0=1 (normovanie TF)

%figure(1)   % Vykresľovanie a popis priebehov PrCh a LFCh
subplot(121); step(F,Tstep,color), grid on,
    title('Prechodová charakteristika','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',16)
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';

subplot(122); bode(F,{wmin,wmax},color), grid on
     title('Frekvenčná charakteristika','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',2) % inštrukcia pre zmenu hrúbky čiary 
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';