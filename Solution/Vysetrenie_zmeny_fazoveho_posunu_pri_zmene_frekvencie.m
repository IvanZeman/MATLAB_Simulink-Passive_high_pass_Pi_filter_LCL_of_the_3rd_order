% Program pre vyšetrenie zmeny fázového posunu pri zmene frekvencie. I. Zeman, 24.4.2022
% Znázornenie PrCH a znázornenie LFCH viacslučkového obvodu pri 3 rôznych frekvenciách
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Zmeny veľkosti popisov osi, farby a hrúbky čiary

clear, clc, clf, format compact
syms s R1 R2 L1 L2 C Uin Uoutsim %deklarácia symbolických premenných

% Zadanie vstupných hodnôt
disp('Hornopriepustný LCL filter so záťažou - riešenie obvodu MSP')
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
Tstep=10e-3; wmin=1; wmax=1e5;                     % parametre pre Step a Bode

%Nájdená TF v symbolickom tvare
F=(C*L1*L2*R2*s^3)/(R1*R2 + L1*R2*s + L2*R1*s + L1*L2*s^2 + C*L1*L2*R1*s^3 + C*L1*L2*R2*s^3 + C*L1*R1*R2*s^2 + C*L2*R1*R2*s^2)

% Spracovanie údajov TF v symbolickom tvare pre prechod do num. MATLABu
[cit,men]=numden(F);                     % oddelenie polynómov čitateľa a menovateľa
cit=subs(cit,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men=subs(men,{R1,R2,L1,L2,C,Uinx,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b=sym2poly(cit);                         % b - koeficienty polynómu čitateľa b(s)
a=sym2poly(men);                         % a - koeficienty polynómu menovateľa a(s)
b=double(b);                             % Prechod do numerickeho MATLABu
a=double(a); 
F=tf(b,a)                                % Výsledná TF v numerickom MATLABe
F=tf(b/a(end),a/a(end))                  % TF upravená pre a0=1 (normovanie TF)

figure(1) % Vyšetrovanie fázvého posunu pri rôznych frekvenciách
color='r';
n=3; % počet zobrazených periód ~ napätia
subplot(1,2,1), bode(F), grid % Kreslenie LFCh
title('Frekvenčná charakteristika','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',14),ylabel('\rightarrow\phi','FontSize',14)
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';

     %frekvencia w1
subplot(3,2,2)
    set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary
    w1=1e3; T1=2*pi/w1; Tkon1=n*T1; % pre zvolenú w1: doba periody T1,doba kon.
    Tkon1
    [u1,t1]=gensig('sin',T1,Tkon1,1e-5);% doba periódy, trvanie signálu,vzorkovanie
    lsim (F,u1,t1),grid 
    w1str=['w1 = ',num2str(w1),' rad/s']; w1text=join(w1str); title (w1text)
    ax = gca        %úprava popisu osí - farba, veľkosť, bold 
    ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';

     %frekvencia w2
subplot(3,2,4)
    set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiar
    w2=5e3; T2=2*pi/w2; Tkon2=n*T2; % pre zvolenú w2: doba periody T2,
    Tkon2
    [u2,t2]=gensig('sin',T2,Tkon2,1e-5); % doba periódy, trvanie signálu,vzorkovanie
    lsim (F,u2,t2), grid
    w2str=['w2 = ',num2str(w2),' rad/s']; w2text=join(w2str); title (w2text)
    ax = gca        %úprava popisu osí - farba, veľkosť, bold 
    ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';

     %frekvencia w3
subplot(3,2,6)
    w3=15e3 ; T3=2*pi/w3; Tkon3=n*T3; % pre zvolenú w3: doba periody T3,
    Tkon3
    [u3,t3]=gensig('sin',T3,Tkon3,1e-6); % doba periódy, trvanie signálu,vzorkovanie
    lsim (F,u3,t3), grid,
    w3str=['w3 = ',num2str(w3),' rad/s']; w3text=join(w3str); title (w3text)
    set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiar
    ax = gca        %úprava popisu osí - farba, veľkosť, bold 
    ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';