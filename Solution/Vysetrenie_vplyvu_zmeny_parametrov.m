% Program pre vyšetrenie vplyvu zmeny parametrov jednotlivých prvkov obvodu. I. Zeman, 22.4.2022
% Simulácia a znázornenie PrCH a LFCH pri zmene parametrov obvodu
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Zmeny veľkosti popisov osi, farby a hrúbky čiary

clear, clc, clf, format compact
syms s R1 R2 L1 L2 C Uin Uoutsim %deklarácia symbolických premenných

% Zadanie vstupných hodnôt
disp('Hornopriepustný LCL filter so záťažou - riešenie obvodu MSP')
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
R11=0.1*R1x; R12=10*R1x; 
L11=0.1*L1x; L12=100*L1x; 
L21=0.1*L2x; L22=100*L1x;
C11=0.1*Cx; C12=10*Cx; 

Tstep=10e-3; wmin=1; wmax=1e5;                     % parametre pre Step a Bode
Tchirp=0.04; wminchirp=1e1; wmaxchirp=1e3;          % parametre pre Chirp
Tsim=0.02;                                          % doba simulácie v Simulinku (experimentálne určená)
color='r';                                          % farba grafu b,r,y,m,c,

%Nájdená TF v symbolickom tvare
F=(C*L1*L2*R2*s^3)/(R1*R2 + L1*R2*s + L2*R1*s + L1*L2*s^2 + C*L1*L2*R1*s^3 + C*L1*L2*R2*s^3 + C*L1*R1*R2*s^2 + C*L2*R1*R2*s^2)

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



%Zmena parametra R1 za R11 
%Úprava TF v symbolickom tvare
FR11 =(C*L1*L2*R2*s^3)/(R11*R2 + L1*R2*s + L2*R11*s + L1*L2*s^2 + C*L1*L2*R11*s^3 + C*L1*L2*R2*s^3 + C*L1*R11*R2*s^2 + C*L2*R11*R2*s^2)      
[cit2,men2]=numden(FR11);                     % oddelenie polynómov čitateľa a menovateľa
cit2=subs(cit2,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men2=subs(men2,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b2=sym2poly(cit2);                         % b - koeficienty polynómu čitateľa b(s)
a2=sym2poly(men2);                         % a - koeficienty polynómu menovateľa a(s)
b2=double(b2);                             % Prechod do numerickeho MATLABu
a2=double(a2);  
FR11=tf(b2,a2)                             % Výsledná TF v numerickom MATLABe
FR11=tf(b2/a2(end),a2/a2(end))             % Výsledná TF v numerickom MATLABe

%Zmena parametra R1 za R12 
%Úprava TF v symbolickom tvare
FR12=(C*L1*L2*R2*s^3)/(R12*R2 + L1*R2*s + L2*R12*s + L1*L2*s^2 + C*L1*L2*R12*s^3 + C*L1*L2*R2*s^3 + C*L1*R12*R2*s^2 + C*L2*R12*R2*s^2)      
[cit3,men3]=numden(FR12);                     % oddelenie polynómov čitateľa a menovateľa
cit3=subs(cit3,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men3=subs(men3,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b3=sym2poly(cit3);                         % b - koeficienty polynómu čitateľa b(s)
a3=sym2poly(men3);                         % a - koeficienty polynómu menovateľa a(s)
b3=double(b3);                             % Prechod do numerickeho MATLABu
a3=double(a3);  
FR12=tf(b3,a3)                             % Výsledná TF v numerickom MATLABe
FR12=tf(b3/a3(end),a3/a3(end))             % Výsledná TF v numerickom MATLABe

%PrCh a LFCh pre vplyv zmeny parametra R1
figure(1)   % Vykresľovanie a popis priebehov PrCh a LFCh pre vplyv zmeny parametra R1
subplot(1,2,1), step(FR11,Tstep,'r',F,'g',FR12,'b'),grid on,
   title('PrCh pri zmene R1','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',16)
     legend('R1/10','R','10*R1')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';
subplot(1,2,2), bode(FR11,{wmin,wmax},'r',F,'g',FR12,'b'),grid on,
     title('LFCh pri zmene R1','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary 
     legend('R1/10','R1','10*R1')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';

%Zmena parametra L1 za L11 
%Úprava TF v symbolickom tvare
FL11 =(C*L11*L2*R2*s^3)/(R1*R2 + L11*R2*s + L2*R1*s + L11*L2*s^2 + C*L11*L2*R1*s^3 + C*L11*L2*R2*s^3 + C*L11*R1*R2*s^2 + C*L2*R1*R2*s^2)      
[cit3,men3]=numden(FL11);                     % oddelenie polynómov čitateľa a menovateľa
cit3=subs(cit3,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men3=subs(men3,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b3=sym2poly(cit3);                         % b - koeficienty polynómu čitateľa b(s)
a3=sym2poly(men3);                         % a - koeficienty polynómu menovateľa a(s)
b3=double(b3);                             % Prechod do numerickeho MATLABu
a3=double(a3);  
FL11=tf(b3,a3)                             % Výsledná TF v numerickom MATLABe
FL11=tf(b3/a3(end),a3/a3(end))             % Výsledná TF v numerickom MATLABe

%Zmena parametra L1 za L12 
%Úprava TF v symbolickom tvare
FL12=(C*L12*L2*R2*s^3)/(R1*R2 + L12*R2*s + L2*R1*s + L12*L2*s^2 + C*L12*L2*R1*s^3 + C*L12*L2*R2*s^3 + C*L12*R1*R2*s^2 + C*L2*R1*R2*s^2)      
[cit4,men4]=numden(FL12);                     % oddelenie polynómov čitateľa a menovateľa
cit4=subs(cit4,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men4=subs(men4,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b4=sym2poly(cit4);                         % b - koeficienty polynómu čitateľa b(s)
a4=sym2poly(men4);                         % a - koeficienty polynómu menovateľa a(s)
b4=double(b4);                             % Prechod do numerickeho MATLABu
a4=double(a4);  
FL12=tf(b4,a4)                             % Výsledná TF v numerickom MATLABe
FL12=tf(b4/a4(end),a4/a4(end))             % Výsledná TF v numerickom MATLABe

%PrCh a LFCh pre vplyv zmeny parametra L1
figure(2)   % Vykresľovanie a popis priebehov PrCh a LFCh pre vplyv zmeny parametra L1
subplot(1,2,1), step(FL11,Tstep,'r',F,'g',FL12,'b'),grid on,
   title('PrCh pri zmene L1','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',16)
     legend('L1/10','L1','100*L1')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';
subplot(1,2,2), bode(FL11,{wmin,wmax},'r',F,'g',FL12,'b'),grid on,
     title('LFCh pri zmene L1','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary 
     legend('L1/10','L1','100*L1')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';



%Zmena parametra C za C11
%Úprava TF v symbolickom tvare
FC11=(C11*L1*L2*R2*s^3)/(R1*R2 + L1*R2*s + L2*R1*s + L1*L2*s^2 + C11*L1*L2*R1*s^3 + C11*L1*L2*R2*s^3 + C11*L1*R1*R2*s^2 + C11*L2*R1*R2*s^2)      
[cit5,men5]=numden(FC11);                     % oddelenie polynómov čitateľa a menovateľa
cit5=subs(cit5,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men5=subs(men5,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b5=sym2poly(cit5);                         % b - koeficienty polynómu čitateľa b(s)
a5=sym2poly(men5);                         % a - koeficienty polynómu menovateľa a(s)
b5=double(b5);                             % Prechod do numerickeho MATLABu
a5=double(a5);  
FC11=tf(b5,a5)                             % Výsledná TF v numerickom MATLABe
FC11=tf(b5/a5(end),a5/a5(end))             % Výsledná TF v numerickom MATLABe

%Zmena parametra C za C12 
%Úprava TF v symbolickom tvare
FC12=(C12*L1*L2*R2*s^3)/(R1*R2 + L1*R2*s + L2*R1*s + L1*L2*s^2 + C12*L1*L2*R1*s^3 + C12*L1*L2*R2*s^3 + C12*L1*R1*R2*s^2 + C12*L2*R1*R2*s^2)      
[cit6,men6]=numden(FC12);                     % oddelenie polynómov čitateľa a menovateľa
cit6=subs(cit6,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men6=subs(men6,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b6=sym2poly(cit6);                         % b - koeficienty polynómu čitateľa b(s)
a6=sym2poly(men6);                         % a - koeficienty polynómu menovateľa a(s)
b6=double(b6);                             % Prechod do numerickeho MATLABu
a6=double(a6);  
FC12=tf(b6,a6)                             % Výsledná TF v numerickom MATLABe
FC12=tf(b6/a6(end),a6/a6(end))             % Výsledná TF v numerickom MATLABe

%PrCh a LFCh pre vplyv zmeny parametra C
figure(3)   % Vykresľovanie a popis priebehov PrCh a LFCh pre vplyv zmeny parametra C
subplot(1,2,1), step(FC11,Tstep,'r',F,'g',FC12,'b'),grid on,
   title('PrCh pri zmene C','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',16)
     legend('C/10','C','10*C')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';
subplot(1,2,2), bode(FC11,{wmin,wmax},'r',F,'g',FC12,'b'),grid on,
     title('LFCh pri zmene C','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary (2 body)
     legend('C/10','C','10*C')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';


%Zmena parametra L2 za L21 
%Úprava TF v symbolickom tvare
FL21 =(C*L1*L21*R2*s^3)/(R1*R2 + L1*R2*s + L21*R1*s + L1*L21*s^2 + C*L1*L21*R1*s^3 + C*L1*L21*R2*s^3 + C*L1*R1*R2*s^2 + C*L21*R1*R2*s^2)      
[cit7,men7]=numden(FL21);                     % oddelenie polynómov čitateľa a menovateľa
cit7=subs(cit7,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men7=subs(men7,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b7=sym2poly(cit7);                         % b - koeficienty polynómu čitateľa b(s)
a7=sym2poly(men7);                         % a - koeficienty polynómu menovateľa a(s)
b7=double(b7);                             % Prechod do numerickeho MATLABu
a7=double(a7);  
FL21=tf(b7,a7)                             % Výsledná TF v numerickom MATLABe
FL21=tf(b7/a7(end),a7/a7(end))             % Výsledná TF v numerickom MATLABe

%Zmena parametra L2 za L22 
%Úprava TF v symbolickom tvare
FL22=(C*L1*L22*R2*s^3)/(R1*R2 + L1*R2*s + L22*R1*s + L1*L22*s^2 + C*L1*L22*R1*s^3 + C*L1*L22*R2*s^3 + C*L1*R1*R2*s^2 + C*L22*R1*R2*s^2)      
[cit8,men8]=numden(FL22);                     % oddelenie polynómov čitateľa a menovateľa
cit8=subs(cit8,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie hodnôt do polynómu čitateľa
men8=subs(men8,{R1,R2,L1,L2,C,Uin,Uoutsim},{R1x,R2x,L1x,L2x,Cx,Uinx,Uoutsim});        % dosadenie do polynómu menovateľa
b8=sym2poly(cit8);                         % b - koeficienty polynómu čitateľa b(s)
a8=sym2poly(men8);                         % a - koeficienty polynómu menovateľa a(s)
b8=double(b8);                             % Prechod do numerickeho MATLABu
a8=double(a8);  
FL22=tf(b8,a8)                             % Výsledná TF v numerickom MATLABe
FL22=tf(b8/a8(end),a8/a8(end))             % Výsledná TF v numerickom MATLABe

%PrCh a LFCh pre vplyv zmeny parametra L2
figure(4)   % Vykresľovanie a popis priebehov PrCh a LFCh pre vplyv zmeny parametra L2
subplot(1,2,1), step(FL21,Tstep,'r',F,'g',FL22,'b'),grid on,
   title('PrCh pri zmene L2','FontSize',16)
     xlabel('\rightarrow T','FontSize',16)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',16)
     legend('L2/10','L2','10*L2')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';
subplot(1,2,2), bode(FL21,{wmin,wmax},'r',F,'g',FL22,'b'),grid on,
     title('LFCh pri zmene L2','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',16),ylabel('\rightarrow\phi','FontSize',16)
     set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiary (2 body)
     legend('L2/10','L2','10*L2')
     ax = gca        %úprava popisu osí - farba, veľkosť, bold 
     ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';