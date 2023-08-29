% Program pre výpočet a zobrazenie TF na základe blokovej schémy prenosového článku. I. Zeman, 10.5.2022
% Znázornenie PrCH a LFCH viacslučkového obvodu HPF LCL
% Príklad: trojslučkový obvod - hornopriepustný LCL filter
% Zmeny veľkosti popisov osi, farby a hrúbky čiary
% Simulácia a znázornenie PrCH a LFCH z blok. modelu v Simulinku

clear, clc, clf, format compact
syms s R1 R2 L1 L2 C Uin %deklarácia symbolických premenných
R1x=20; R2x=50; L1x=14e-3; L2x=7e-3; Cx=15e-6; Uinx=10;
Tstep=10e-3; wmin=1; wmax=1e5;                      % parametre pre Step a Bode

[A,B,C,D]=linmod('LCL_filter_blok_schem_S');
[num,den]=ss2tf(A,B,C,D)
F=tf(num/den(end),den/den(end))

figure(1)   % Vykresľovanie a popis priebehov PrCh a LFCh
color='r'; 
subplot(121); step(F,Tstep,color), grid on,
   title('Prechodová charakteristika','FontSize',16)
     xlabel('\rightarrow T','FontSize',14)
     ylabel('\rightarrow U_{out}/U_{in}','FontSize',14)
        set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiar
        ax = gca        %úprava popisu osí - farba, veľkosť, bold 
        ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';
subplot(122); bode(F,{wmin,wmax},color), grid on
     title('Frekvenčná charakteristika','fontsize',16)
     xlabel('\rightarrow \omega','FontSize',14),ylabel('\rightarrow\phi','FontSize',14) 
        set(findall(gcf,'type','line'),'linewidth',1.5) % inštrukcia pre zmenu hrúbky čiar
        ax = gca        %úprava popisu osí - farba, veľkosť, bold 
        ax.YColor = 'k'; ax.XColor = 'k'; ax.FontSize = 12; ax.FontWeight = 'bold';