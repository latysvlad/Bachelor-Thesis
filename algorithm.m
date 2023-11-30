
close all;
WWW1=1;
while WWW1==1
    WWW2=1;
while WWW2==1
clear all;
close all;
WWW1=1;
WWW2=1;
%dbstop if error;
scrsz = get(0,'ScreenSize');

%ЗАДАЧА ИСХОДНЫХ ДАННЫХ
xn=400;                         %начальная длина волны сигнала
xk=700;                         %конечная длина волны сигнала
dx=1;                         %дискретный шаг
da=(xk-xn)/8;
ja=(xk-xn)/dx+1;                %количество итераций(шагов)
ImRe=2;                         %Выбор мат/реальных фильтров(1,2)       
ca=3;                           %(количество СФ)+1(исходный сигнал)
LLa=3;                          %число спектров
TriVos=2;                       %Трёх-/Восьми-/Восьми-(реальный)/Девятиканал(1,2,3,4)
auto=1;                         %Автоматическая простановка данных(!=0)
curve=2;                        %Номер кривой(1,2,3,4)
GreTih=3; % Гревиль(1)/Тихонов(кастомный)(2)/Тихонов(библиотечный)(3)/Вейвлет(базисные)(4)/Годунов(5)
kray=xk-xn+1;                  %Распределение графиков спектров по длинам волнам
prozent=0.005;  %    0,5%           %относительный процент погрешности 
%lambda=[1 2 3 4 6];
lambda=0.31;                       %параметр регуляризации
Tau=0.41;
aprox=1;                          %0-выключение моей аппроксимации

                                       %БЛОК СЛУЧАЙНОЙ СОСТАВЛЯЮЩЕЙ
                                       
LS=(rand(1,ja)-0.5)/100; 
r1=rand-0.5;
r2=r1+0.02;

r3=rand-0.5;
r4=r3+0.1;
for j=1:ja 
    LS(j)=rand-0.5;
    if LS(j)>r1 && LS(j)<r2
        rnd=rand;
       if rnd>0.5
        LS(j)=0-LS(j);
       end 
        LS(j)=(LS(j)+rand-0.5)/10;
        
    elseif LS(j)>r3 && LS(j)<r4 
        rnd=rand;
       if rnd>0.5
        LS(j)=0-LS(j);
       end 
        LS(j)=(LS(j)+rand-0.5)/15;
        
    else LS(j)=LS(j)/40;        
        
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    LS=0;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %БЛОК ФОРМИРОВАНИЯ ФУНКЦИЙ ВХОДНОЙ ЯРКОСТИ(L), ЧУВСТВИТЕЛЬНОСТИ(S) И ПРОПУСКАНИЯ СФ(T)

%%                               ЗАДАЧА ВХОДНОГО СИГНАЛА

x0 = xn:dx:xk;                      %Разбитый диапазон длин волн
%L = sin(4*2*pi*(x0-xn)/(xk-xn))*0.3+0.5;
%L = sin(2*pi*(L*400-0.5-xn)/(xk-xn))*0.3+0.5;
%L=cos(2*pi*(L-0.5)*10)+1;
%L = x0;  %Пробник
Lo = sin(2*pi*(x0-xn)/(xk-xn))*0.3+0.5+LS;   %Исходная функция яркости
%Реальные кривые
%Красный
RealRed=double(rgb2gray(imread('RealSpectraPink.png'))); RealR=pic2val(RealRed, 5, 60, 1, kray+50);LR=RealR(1:dx:kray)+LS;
%Зеленый
RealGreen=double(rgb2gray(imread('RealSpectraGreen.png')));RealG=pic2val(RealGreen, 5, 60, 1, kray+50);LG=RealG(1:dx:kray)+LS;
%Желтый
RealYellow=double(rgb2gray(imread('RealSpectraYellow.png')));RealY=pic2val(RealYellow, 5, 60, 1, kray+50); LY=RealY(1:dx:kray)+LS;

%Выбор кривой сигнала яркости 
if auto~=0
    inL=curve;
else
    inL=input('Number of curve(1,2,3,4):');
end

L=[LR;LG;LY];
%L=[Lo;LG;LY];
   %нахождение максимума для нормировки погрешности  
 maximum(1:LLa) = 0;
 for LL=1:LLa
 for j=1:ja 
     if L(LL,j) > maximum(LL)
         maximum(LL) = L(LL,j);
     end
 end
 end
%%                              ЗАДАЧА ФУНКЦИЙ ПРОПУСКАНИЯ СФ

t0=ones(1,ja);                      %Без светофильтра          
t1=x0./1000;                        %СФ1(линейка+)
t2=(x0-300)./500+2;                 %СФ2(линейка-)
t3=exp(x0*1.1/1000)-1.45;           %СФ3(экспонента+)
t4=3.2-exp(x0*1.1/1000);            %СФ4(экспонента-)
t5=sqrt(x0/1000);                   %СФ5(корень+)
t6=cos(x0/1000)./sqrt(x0/1000);     %СФ6(несусветная лабудень)
t7=1./(x0-390)+1;                   %СФ7(гипербола)
t8=(x0.*x0-400)/160000;             %СФ8(парабола)
t9=cos(x0/1000);                    %СФ9(косинус)

%РЕАЛЬНЫЕ ФИЛЬТРЫ
%Без СФ
tr0=ones(1,ja);
%ЖЗС-5 tr1
JZS5=double(rgb2gray(imread('JZS-5.png')));
jzs5=pic2val(JZS5,3,60,5.43,kray);
jzs5=jzs5(1:dx:kray);
tr1=jzs5;
%ЖЗС-18 tr2
JZS18=double(rgb2gray(imread('JZS-18.png')));jzs18=pic2val(JZS18,6,60,5.43,kray);jzs18=jzs18(1:dx:xk-xn+1);tr2=jzs18;
%CЗC-16 tr3
SZS16=double(rgb2gray(imread('SZS-16.png')));szs16=pic2val(SZS16,6,80,5.43,kray);szs16=szs16(1:dx:xk-xn+1);tr3=szs16;
%CC-1 tr4
SS1=double(rgb2gray(imread('SS-1.png')));ss1=pic2val(SS1,6,20,1,kray);ss1=ss1(1:dx:xk-xn+1);tr4=ss1;
%S1 tr5
S1=double(rgb2gray(imread('S1.png')));s1=pic2val(S1,6,20,1,kray);s1=s1(1:dx:xk-xn+1);tr5=s1;
%S2 tr6
S2=double(rgb2gray(imread('S2.png')));s2=pic2val(S2,6,20,1,kray);s2=s2(1:dx:xk-xn+1);tr6=s2;
%S3 tr7
S3=double(rgb2gray(imread('S3.png')));s3=pic2val(S3,6,20,1,kray);s3=s3(1:dx:xk-xn+1);tr7=s3;
%S4 tr8
S4=double(rgb2gray(imread('S4.png')));s4=pic2val(S4,6,20,1,kray);s4=s4(1:dx:xk-xn+1);tr8=s4;


                        %ВЫВОД ГРАФИКОВ ФУНКЦИЙ ПРОПУСКАНИЯ СФ
                        
T1=[t0;t1;t3;t5;t6;t7;t8;t9];                           %Массив функций пропускания условных СФ
T2=[tr0;tr2;tr4;tr1;tr3;tr5;tr6;tr7;tr8];                                  %Массив функций пропускания реальных СФ
%T2=[tr0;tr5;tr6;tr7;tr8];
if auto~=0
    inT=ImRe;
else
    inT=input('Choose T Image/Real(1/2):');
end

if inT==1
    T=T1;
    if ca>8 || ca==0
        ca=8;
    end
else
    T=T2;
    if ca>9 || ca==0
        ca=9;
    end 
end

f1=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on
%for c=1:ca    
    
    %plot(x0,T(2:ca,:)*100,'k','linewidth', 1);
    plot(x0,T(2,:)*100,'Color',[0.6,0.8,0],'linewidth', 2);
    plot(x0,T(3,:)*100,'Color',[0,0,1],'linewidth', 2);
    %plot(x0,T(4,:)*100,'k-.','linewidth', 1.5);
    %plot(x0,T(5,:)*100,'k:','linewidth', 2);
    set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman'); 
    %ylabel('\bf $$\tau,{\enskip}\%$$','interpreter','latex','fontsize',22,'rotation',90');
    %ylabel('\bf $$\mathrm{Transmittance {\enskip}\tau,\%}$$','interpreter','latex','fontsize',22,'rotation',90');
    %xlabel('\bf $$\mathrm{Wavelength {\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
    ylabel('Пропукание \tau, %','fontsize',28,'rotation',90');
    xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
    axis([xn xk 0 100]);
    %legend('1','2','3','4');
    legend('ЖЗС-18','СС-1');
%end
hold off
%grid on;

%plot(fcoord, P, 'k -', 'linewidth', 1.2);
%set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
%%grid on; 
 %axis([Fn1 St1 0 1.1]);
 %ylabel('\bf $${h}\tilde,{\enskip}arb. units$$','interpreter','latex','fontsize',15,'rotation',90');
 %xlabel('\bf $$\nu,{\enskip}mm^{-1}.$$','interpreter','latex','fontsize',15,'rotation',0');

%%                              ЗАДАЧА СПЕКТРАЛЬНОЙ ЧУВСТВИТЕЛЬНОСТИ

x=size(10);                            
a=xn+da/2;                                  %задание основной длины волны первого ф.Байера
SS=(rand(1,ja)-0.5)/100;                %Погрешность
if TriVos==1
    ka=3;                        %количество субпикселей    
        Red=double(rgb2gray(imread('R.png')));red=pic2val(Red,6,60,5.43,kray+50);red=red(1:dx:xk-xn+1);
        for j=1:ja
        x(1,j,1)=xn+(j-1)*dx;
        x(2,j,1)=red(j);
        end
        Green=double(rgb2gray(imread('G.png')));green=pic2val(Green,6,60,5.43,kray+50);green=green(1:dx:xk-xn+1);
        for j=1:ja
        x(1,j,2)=xn+(j-1)*dx;
        x(2,j,2)=green(j);
        end
        Blue=double(rgb2gray(imread('B.png')));blue=pic2val(Blue,6,60,5.43,kray+50);blue=blue(1:dx:xk-xn+1);
        for j=1:ja
        x(1,j,3)=xn+(j-1)*dx;
        x(2,j,3)=blue(j);
        end
        for k=1:ka 
            for j=1:ja
                tau(j,k)=1-x(2,j,k);
                Smatrm(j,k)=x(2,j,k);
            end 
        end
        tau=tau';
f21=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
plot(x0,tau(1,:),'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tau(2,:),'k-o','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tau(3,:),'k-*','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman'); 
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
axis([xn xk 0 1.2]);
legend('1', '2','3');
%grid on;
%%{
f22=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
plot(x0,tr1.*blue,'k-.','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tr1.*red,'k--','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tr2.*red,'k--','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr2.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr3.*blue,'k-.','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr3.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*blue,'k-.','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*red,'k--','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman'); 
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
axis([xn xk 0 1]);
%grid on;
%legend('1', '2','3');
%}
%{
f22=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
plot(x0,tr1.*blue,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tr1.*red,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,tr2.*red,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr2.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr3.*blue,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr3.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*blue,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*red,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,tr4.*green,'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman'); 
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
axis([xn xk 0 1]);
%}

elseif TriVos==2
    ka=8;
    %%{
for k=1:ka                              %порядковый номер фильтра в порядке возрастания ?max
    r1=rand-0.5;
    r2=r1+0.02;

    r3=rand-0.5;
    r4=r3+0.1;
    for j=1:ja                          %итерационный счётчик ? и S(?)
        x(1,j,k)=xn+(j-1)*dx;           %дискретные значения длин волн ?
        
        SS(j)=rand-0.5;
    if SS(j)>r1 && SS(j)<r2
        rnd=rand;
       if rnd>0.5
        SS(j)=0-SS(j);
       end 
        SS(j)=(SS(j)+rand-0.5)/10;
        
    elseif SS(j)>r3 && SS(j)<r4 
        rnd=rand;
       if rnd>0.5
        SS(j)=0-SS(j);
       end 
        SS(j)=(SS(j)+rand-0.5)/15;
        
    else SS(j)=SS(j)/40;        
        
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    SS(j)=0;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %значения ф-ии распр.спектр-й чувст-ти S(?)
        %x(2,j,k)=1/sqrt(1+(x(1,j,k)-a).*(x(1,j,k)-a)/80);                 %Лоренцева с корнем
        x(2,j,k)=1/(1+(x(1,j,k)-a).*(x(1,j,k)-a)/120)+SS(j);                      %Лоренцева
        Smatrm(j,k)=x(2,j,k);
        %x(2,j,k)=exp((-1)*(x(1,j,k)-a).*(x(1,j,k)-a)/400);                 %Гауссова
        %x(2,j,k)=(exp((-1)*0.8*(x(1,j,k)-a).*(x(1,j,k)-a)/400)+0.05)/1.05; %Гауссова расширенная поднятая
    end
    a=a+da;                             %переход на следующий фильтр
end
%}
elseif TriVos==3
    ka=8;
for j=1:ja                          %итерационный счётчик ? и S(?)
    for k=1:ka 
        x(1,j,k)=xn+(j-1)*dx;           %дискретные значения длин волн ?
    end
        x(2,j,1)=1/(1+(x(1,j,1)-424).*(x(1,j,1)-424)/150)*0.62;                      %Лоренцева
        Smatrm(j,1)=x(2,j,1);
        x(2,j,2)=1/(1+(x(1,j,2)-457).*(x(1,j,2)-457)/150)*0.61;
        Smatrm(j,2)=x(2,j,2);
        x(2,j,3)=1/(1+(x(1,j,3)-497).*(x(1,j,3)-497)/150)*0.58;
        Smatrm(j,3)=x(2,j,3);
        x(2,j,4)=1/(1+(x(1,j,4)-534).*(x(1,j,4)-534)/150)*0.58;
        Smatrm(j,4)=x(2,j,4);
        x(2,j,5)=1/(1+(x(1,j,5)-561).*(x(1,j,5)-561)/150)*0.58;
        Smatrm(j,5)=x(2,j,5);
        x(2,j,6)=1/(1+(x(1,j,6)-599).*(x(1,j,6)-599)/150)*0.56;
        Smatrm(j,6)=x(2,j,6);
        x(2,j,7)=1/(1+(x(1,j,7)-640).*(x(1,j,7)-640)/150)*0.56;
        Smatrm(j,7)=x(2,j,7);
        x(2,j,8)=1/(1+(x(1,j,8)-680).*(x(1,j,8)-680)/150)*0.56;
        Smatrm(j,8)=x(2,j,8);
        
end
else
    ka=9;
for k=1:ka                              %порядковый номер фильтра в порядке возрастания ?max
    r1=rand-0.5;
    r2=r1+0.02;

    r3=rand-0.5;
    r4=r3+0.1;
    for j=1:ja                          %итерационный счётчик ? и S(?)
        x(1,j,k)=xn+(j-1)*dx;           %дискретные значения длин волн ?
        
        SS(j)=rand-0.5;
    if SS(j)>r1 && SS(j)<r2
        rnd=rand;
       if rnd>0.5
        SS(j)=0-SS(j);
       end 
        SS(j)=(SS(j)+rand-0.5)/10;
        
    elseif SS(j)>r3 && SS(j)<r4 
        rnd=rand;
       if rnd>0.5
        SS(j)=0-SS(j);
       end 
        SS(j)=(SS(j)+rand-0.5)/15;
        
    else SS(j)=SS(j)/40;        
        
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    SS(j)=0;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k<=ka-1
    %значения ф-ии распр.спектр-й чувст-ти S(?)
        %x(2,j,k)=1/sqrt(1+(x(1,j,k)-a).*(x(1,j,k)-a)/80);                 %Лоренцева с корнем
        x(2,j,k)=1/(1+(x(1,j,k)-a).*(x(1,j,k)-a)/120)+SS(j);                      %Лоренцева
        Smatrm(j,k)=x(2,j,k);
        %x(2,j,k)=exp((-1)*(x(1,j,k)-a).*(x(1,j,k)-a)/400);                 %Гауссова
        %x(2,j,k)=(exp((-1)*0.8*(x(1,j,k)-a).*(x(1,j,k)-a)/400)+0.05)/1.05; %Гауссова расширенная поднятая
    elseif k==ka
       x(2,j,k)=exp((-x(1,j,k)+330)*0.01)+0.1+SS(j);                      %Лоренцева
       Smatrm(j,k)=x(2,j,k);
    end
    end
    a=a+da;                             %переход на следующий фильтр
end
end
                    %ВЫВОД ГРАФИКОВ ЯРКОСТИ И СПЕКТРАЛЬНЫХ ЧУВСТВИТЕЛЬНОСТЕЙ
                    
f2=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
plot(x0,L(1,:),'r','MarkerIndices',1:5:length(L(1,:)),'linewidth', 2); hold on;
plot(x0,L(2,:),'g','MarkerIndices',1:5:length(L(1,:)),'linewidth', 2); hold on;
plot(x0,L(3,:),'Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(1,:)),'linewidth', 2);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman'); 
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
ylabel('Спектральная плотности яркости L, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
axis([xn xk 0 1.2]);
%grid on;
%legend('1', '2','3');

f3=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
for k=1:ka
   if k==1
       plot(x(1,:,k),x(2,:,k),'Color',[0,0,0.6],'linewidth', 2);
   elseif k==2
       plot(x(1,:,k),x(2,:,k),'Color',[0.3,0.3,0.8],'linewidth', 2);
   elseif k==3
       plot(x(1,:,k),x(2,:,k),'Color',[0,0.4,0.3],'linewidth', 2);
   elseif k==4
       plot(x(1,:,k),x(2,:,k),'Color',[0,0.7,0],'linewidth', 2);
   elseif k==5
       plot(x(1,:,k),x(2,:,k),'Color',[0.3,0.8,0],'linewidth', 2);
   elseif k==6
       plot(x(1,:,k),x(2,:,k),'Color',[0.5,0.7,0],'linewidth', 2);
   elseif k==7
       plot(x(1,:,k),x(2,:,k),'Color',[0.5,0,0],'linewidth', 2);
   elseif k==8
       plot(x(1,:,k),x(2,:,k),'Color',[0.7,0.1,0],'linewidth', 2);
   elseif k==9
       plot(x(1,:,k),x(2,:,k),'Color',[0,0,0],'linewidth', 2);
   end
       
   xlim([xn xk]); 
   set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
   %ylabel('\bf $$\mathrm{Sensitivity{\enskip}S(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
   %xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
ylabel('Спектральная чувствительность S, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
end
%grid on;


%%                    БЛОК ВЫЧИСЛЕНИЯ ЗНАЧЕНИЙ СПЕКТРАЛЬНОЙ ЧУВСТВИТЕЛЬНОСТИ(МОДИФИЦИРОВАННОЙ)
                    
S0=zeros(ka,ja,ca);
for c=1:ca
    for k=1:ka
            for j=1:ja
                S0(k,j,c)=x(2,j,k).*T(c,j);
            end
    end   
end
S=size(10);
c=1;
k=1;
USa=ca*ka;
for s=1:USa
    for j=1:ja
        S(s,j)=S0(k,j,c);
    end
    k=k+1;
    if k>ka
        k=k-ka;
        c=c+1;
    end
end
                            %ВЫВОД ГРАФИКОВ МОДИФИЦИРОВАННЫХ СПЕКТР.ЧУВСТ-ТЕЙ
                            
f4=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
for c=1:ca          
    for k=1:ka    
        plot(x(1,:,k),S(k+ka*(c-1),:),'linewidth', 1.2);

        set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
        ylabel('\bf $$S(\lambda)\cdot\tau$$','interpreter','latex','fontsize',22,'rotation',90');
        xlabel('\bf $$\lambda,{\enskip}nm$$','interpreter','latex','fontsize',22,'rotation',0');


    end
end
hold off
%grid on;
%%                           БЛОК ВЫЧИСЛЕНИЯ ЗНАЧЕНИЙ СИГНАЛОВ(U) НА ВЫХОДЕ
                            
U0=zeros(LLa,ca,ka);
for LL=1:LLa
for c=1:ca
    for k=1:ka
        U0(LL,c,k)=0;
        for j=1:ja-1
               
            o=x(1,j,k);
            o0=dx.*(L(LL,j+1).*x(2,j+1,k).*T(c,j+1)+L(LL,j).*x(2,j,k).*T(c,j))./2;
            o1=(L(LL,j+1).*x(2,j+1,k).*T(c,j+1)+L(LL,j).*x(2,j,k).*T(c,j))./2;
            o2=L(LL,j);
            o3=L(LL,j+1);
            o4=x(2,j,k);
            o5=x(2,j+1,k);
            
            U0(LL,c,k)=U0(LL,c,k)+dx.*(L(LL,j+1).*x(2,j+1,k).*T(c,j+1)+L(LL,j).*x(2,j,k).*T(c,j))./2;%%%%%%
        end
              
    end
end
end
maximumU(1:LLa,1:ca) = 0;
 for LL=1:LLa
 for c=1:ca 
     for k=1:ka 
     if U0(LL,c,k) > maximumU(LL,c)
         maximumU(LL,c) = U0(LL,c,k);
     end
     end
 end
 end
                             %ВВОД ШУМА В СИГНАЛ
Pogreh=zeros(LLa,ca,ka);                           
for LL=1:LLa  
for c=1:ca 
%pog(LL)=srednee(LL)*prozent;
pog(LL,c)=maximumU(LL,c)*prozent;
%Pogreh(LL,c,:) = randi([-1 1], 1, ka)*pog(LL,c);
Pogreh(LL,c) = randi([-100 100], 1, 1)*pog(LL,c)/100;

end
end  
U0=U0+Pogreh;

U=zeros(LLa,USa);
c=1;
k=1;
srednee=zeros(1,LLa);
for LL=1:LLa
    for u=1:USa
        U(LL,u)=U0(LL,c,k);
        k=k+1;
        if k>ka
            k=k-ka;
            c=c+1;
        end
    srednee(LL)=srednee(LL)+U(LL,u);   
    end
    c=1;
    k=1;
end
srednee=srednee/USa;
U=U';                        
                             %ВЫВОД ГРАФИКОВ НАПРЯЖЕНИЙ
m=size(ka);
for mi=1:ka                             
    m(mi)=mi; 
end
f5=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
xlim([1 ka])
for LL=1:LLa
for c=1:ca
    
    plot(m,U(ka*(c-1)+1:ka*c,LL));
    set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
    ylabel('\bf $$U,{\enskip}V$$','interpreter','latex','fontsize',22,'rotation',90');
    %xlabel('№ of sensitive element','interpreter','latex','fontsize',22,'rotation',0');
    %xlabel('Номер чувствительного элемента','interpreter','latex','fontsize',22,'rotation',0');

end
end
hold off;
%grid on;

%%                              %ВЫВОД ФИНАЛЬНОГО ОТВЕТА

%S*Lx=U; где Lx-столбец неизвестных.
if auto~=0
    inM=GreTih;
else
    inM=input('Grevil/Tihonov(1/2):');
end

if inM==1               %Гревиль
    Lx=pinv(S)*U;
    U0=S*Lx;
    Lx=Lx'/dx;
    g = gausswin(20); 
    g = g/sum(g);
    for LL=1:LLa
    Lx(LL,:) = conv(Lx(LL,:), g, 'same');
    end
elseif inM==2                   %Тихонов кастомный
    Lx=pinv(S+lambda*eye(USa,ja))*U;
    Lx=Lx'/dx;
    %O=pow2(det(S*
elseif inM==3                   %Тихонов библиотечный
    [W,z,V] = csvd(S);
    for LL=1:LLa
        [Lx(LL,:),rho,eta] = tikhonov(W,z,V,U(:,LL),lambda);
    end
    Lx=Lx/dx;
    
elseif inM==4                           %Базисные функции !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x=-1:(1/((kray-1)/2)):1;
    BasicF1=ones(1,kray);
    BasicF2=erf(3.91*x);
    BasicF3=erf(3.86*(x+0.8));
    BasicF4=erf(4.14*(x-0.21));
    BasicF5=erf(5.39*(x+0.72));
    BasicF6=erf(2.24*(x+0.25));
    BasicF7=erf(4.19*(x-0.27));
    BasicF8=erf(6.98*(x+0.25));
    BasicF=[BasicF1', BasicF2', BasicF3', BasicF4', BasicF5', BasicF6', BasicF7', BasicF8'];
    
        Snew=S*BasicF;
    
    Lx=BasicF*pinv(Snew)*U;
    Lx=Lx'/dx;
else
    D=diag(-80*ones(kray,1))+diag(40*ones(kray-1,1),1)+diag(40*ones(kray-1,1),-1);
    
    %%{
    Lx(1,:)=(inv(((1-Tau).^2)*S'*S+(Tau^2)*D'*D))*((((1-Tau)^2)*S'* U(:,1)));

    Lx(2,:)=(inv(((1-Tau).^2)*S'*S+(Tau^2)*D'*D))*((((1-Tau)^2)*S'* U(:,2)));

    Lx(3,:)=(inv(((1-Tau).^2)*S'*S+(Tau^2)*D'*D))*((((1-Tau)^2)*S'* U(:,3)));
    %}
    %Lx=(inv(((1-Tau).^2)*S'*S+(Tau^2)*D'*D))*((((1-Tau)^2)*S'* U));
    %Lx=Lx';

end
    

                            %СРАВНЕНИЕ СИГНАЛА НА ВХОДЕ И НА ВЫХОДЕ
f6=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
%{
plot(x0,L(1,:),'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); 
plot(x0,Lx(1,:),'k:','linewidth', 1.5);
plot(x0,L(2,:),'k-o','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,Lx(2,:),'k:o','MarkerIndices',1:5:length(Lx(2,:)),'linewidth', 1.5);
plot(x0,L(3,:),'k-*','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); 
plot(x0,Lx(3,:),'k:*','MarkerIndices',1:5:length(Lx(3,:)),'linewidth', 1.5);
%}
plot(x0,L(1,:),'r','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,Lx(1,:),'r:','linewidth', 1.5);
plot(x0,L(2,:),'g','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,Lx(2,:),'g:','linewidth', 1.5);
plot(x0,L(3,:),'Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,Lx(3,:),':','Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(1,:)),'linewidth', 1.5);

set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
ylabel('Спектральная плотности яркости L, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
%legend('1', '2', '3', '4', '5', '6');
legend('Кр.исх.', 'Кр.рез.','Зел.исх.', 'Зел.рез.','Жел.исх.', 'Жел.рез.');
axis([xn xk 0 1.2]);
%grid on;


%%                          %АППРОКСИМАЦИЯ(ОПЦИОНАЛЬНО)
if aprox~=0
f9=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
%ylim([0.45 0.85]);
%xlim([xn xn+(xk-xn)/2]);
xlim([xn xk]);
ylim([0 1.2]);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
%ylabel('\bf $$L(\lambda)$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\lambda,{\enskip}nm$$','interpreter','latex','fontsize',22,'rotation',0');
ylabel('Спектральная плотности яркости L, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
%set(gca, 'FontName', 'Times New Roman');
%set(gca, 'FontSize', 20);
%%grid on;
%legend('Исходный сигнал', 'Результат', 'Аппроксимация');

f7=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
%grid on;
ylim([0 1.5]);

f8=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
hold on;
%grid on;
ylim([-1.5 0]);

for LL=1:LLa 
    LxL=Lx(LL,:);
    %{
[zpks,zlocs]=findpeaks(LxL);
    zpks1=zpks;
    zlocs1=zlocs;
    LxL=-Lx(LL,:);
[zpks,zlocs]=findpeaks(LxL);
    zpks2=zpks;
    zlocs2=zlocs;
    %}
[zpks1,zlocs1]=findpeaks(LxL);
figure(f7);
findpeaks(LxL);
hold on;
%grid on;
[zpks2,zlocs2]=findpeaks(-LxL);
figure(f8);
findpeaks(-LxL);
hold on;
%grid on;

ln1 = length(zlocs1);
ln2 = length(zlocs2);
ln = ln1 + ln2;
if ln~=15
    WWW2=0;
    break
end
for i=1:ln
    ost=mod(i,2);
    zel(LL)=fix(i/2);
    if ost~=0
        zlocs(i+1)=zlocs1(zel(LL)+1);    
    else
        zlocs(i+1)=zlocs2(zel(LL)); 
    end

end
if LL>1
if zel(LL)~=zel(LL-1)
    WWW2=0;
    break
end
end

zlocs(1)=1;
zlocs(ln+2)=ja;
%      Первый вариант
%{
for i=1:ka*2
    zspline(i)=round((zlocs(i)+zlocs(i+1))/2);
end
%}
%      Второй вариант
for i=1:ln+1
    sr=(LxL(zlocs(i))+LxL(zlocs(i+1)))/2;
    Lxsr(i)=LxL(zlocs(i));
    for i0=zlocs(i):zlocs(i+1)
        if abs(LxL(i0)-sr)<abs(Lxsr(i)-sr)
            Lxsr(i)=LxL(i0);
            zspline(i)=i0;
        end
    end
end
for i=1:ln+1
    if zspline(i)<=0
       WWW2=0;
       break
    end
end

if WWW2==0
    break
end
Lxz=LxL(zspline);
zspline=zspline+xn-1;
Lxx(LL,:)=spline(zspline,Lxz,x0);
figure(f9);
if LL==1
plot(x0,L(LL,:),'r','linewidth', 1); 
hold on;
plot(x0,Lxx(LL,:),'r:','linewidth', 1.5);

ylim([0 1.2]);
elseif LL==2
plot(x0,L(LL,:),'g','MarkerIndices',1:5:length(L(2,:)),'linewidth', 1);
plot(x0,Lxx(LL,:),'g:','MarkerIndices',1:5:length(Lxx(LL,:)),'linewidth', 1.5);
hold on;
ylim([0 1.2]);
elseif LL==3
plot(x0,L(LL,:),'Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(3,:)),'linewidth', 1); 
plot(x0,Lxx(LL,:),':','Color',[0.9,0.7,0],'MarkerIndices',1:5:length(Lxx(LL,:)),'linewidth', 1.5);
hold on;
ylim([0 1.2]);
end
legend('1', '2', '3', '4', '5', '6');
xlim([xn xk]);
ylim([0 1.2]);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
ylabel('Спектральная плотности яркости L, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
end
end
if WWW2==0
    break
end
legend('Кр.исх.', 'Кр.рез.','Зел.исх.', 'Зел.рез.','Жел.исх.', 'Жел.рез.');
%grid on;
%%                                     Другой вариант аппроксимации                               
f10=figure('Position',[scrsz(3)*2/3 scrsz(4)/3 scrsz(3)/3 scrsz(4)/3],'Color','w');
box on
%grid on;
%{
for LL=1:LLa 
sp=9;               %степень полинома
p=polyfit(x0,Lx(LL,:),sp);
f(LL,:)=polyval(p,x0);
plot(x0,L-LS,x0,Lx,':',x0,f(LL,:),'linewidth', 2);
hold on
end
%}
sp=9;

p1=polyfit(x0,Lx(1,:),sp);
f(1,:)=polyval(p1,x0);
plot(x0,L(1,:),'r','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,f(1,:),'r:','linewidth', 1.5);
p2=polyfit(x0,Lx(2,:),sp);
f(2,:)=polyval(p2,x0);
plot(x0,L(2,:),'g','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); hold on;
plot(x0,f(2,:),'g:','linewidth', 1.5);
p3=polyfit(x0,Lx(3,:),sp);
f(3,:)=polyval(p3,x0);
plot(x0,L(3,:),'Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,f(3,:),':','Color',[0.9,0.7,0],'MarkerIndices',1:5:length(L(1,:)),'linewidth', 1.5);
%{
p=polyfit(x0,Lx(1,:),sp);
f(1,:)=polyval(p,x0);
plot(x0,L(1,:),'k-','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); 
hold on
plot(x0,f(1,:),'k:','linewidth', 1.5);
p=polyfit(x0,Lx(2,:),sp);
f(2,:)=polyval(p,x0);
plot(x0,L(2,:),'k-o','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1);
plot(x0,f(2,:),'k:o','MarkerIndices',1:5:length(Lx(2,:)),'linewidth', 1.5);
p=polyfit(x0,Lx(3,:),sp);
f(3,:)=polyval(p,x0);
plot(x0,L(3,:),'k-*','MarkerIndices',1:5:length(L(1,:)),'linewidth', 1); 
plot(x0,f(3,:),'k:*','MarkerIndices',1:5:length(Lx(3,:)),'linewidth', 1.5);
legend('1', '2', '3', '4', '5', '6');
%}
%ylim([0.45 0.85]);
%xlim([xn xn+(xk-xn)/2]);
xlim([xn xk]);
ylim([0 1.2]);
set(0,'DefaultAxesFontSize',22,'DefaultAxesFontName','Times New Roman');
%ylabel('\bf $$\mathrm{Brightness{\enskip}L(\lambda),r.u.}$$','interpreter','latex','fontsize',22,'rotation',90');
%xlabel('\bf $$\mathrm{Wavelength{\enskip}\lambda,\left[nm\right]}$$','interpreter','latex','fontsize',22,'rotation',0');
%set(gca, 'FontName', 'Times New Roman');
%set(gca, 'FontSize', 20);
ylabel('Спектральная плотности яркости L, о.е.','fontsize',28,'rotation',90');
xlabel('Длина волны \lambda, [нм]','fontsize',28,'rotation',0');
%legend('Исходный сигнал', 'Результат', 'Аппроксимация');

U0=S*L';
 
U0=S*Lx';
%grid on;
%%                    %Расчёт процента отклонения начальный
for LL=1:LLa
Otkl=0;
maxotkl=0;
for j=ceil(da/2):ja-ceil(da/2)
    Otkl=Otkl+abs(L(LL,j)-Lx(LL,j))/L(LL,j);
    if abs(L(LL,j)-Lx(LL,j))/L(LL,j)>maxotkl
        maxotkl=abs(L(LL,j)-Lx(LL,j))/L(LL,j);
    end
end
Otkl=Otkl/ja;
fprintf( 'Otkl=%d\n',Otkl);
fprintf( 'maxotkl=%d\n',maxotkl);
end

if aprox~=0
                    %Расчёт процента отклонения: аппроксимация 1
for LL=1:LLa                    
Otkl1=0;
maxotkl1=0;
for j=ceil(da/2):ja-ceil(da/2)
    Otkl1=Otkl1+abs(L(LL,j)-Lxx(LL,j))/Lxx(LL,j);
    if abs(L(LL,j)-Lxx(LL,j))/Lxx(LL,j)>maxotkl1
        maxotkl1=abs(L(LL,j)-Lxx(LL,j))/Lxx(LL,j);
    end
end
Otkl1=Otkl1/ja;
fprintf( 'Otkl1=%d\n',Otkl1);
fprintf( 'maxotkl1=%d\n',maxotkl1);
end
end

                    %Расчёт процента отклонения: аппроксимация 2
for LL=1:LLa 
Otkl2=0;
maxotkl2=0;
for j=ceil(da/2):ja-ceil(da/2)
    Otkl2=Otkl2+abs(L(LL,j)-f(LL,j))/L(LL,j);
    if abs(L(LL,j)-f(LL,j))/L(LL,j)>maxotkl2
        maxotkl2=abs(L(LL,j)-f(LL,j))/L(LL,j);
    end
end
Otkl2=Otkl2/ja;
fprintf( 'Otkl2=%d\n',Otkl2);
fprintf( 'maxotkl2=%d\n',maxotkl2);
end
WWW2=0;
WWW1=0;
end                         %конец внутреннего вайла(проверка на ошибку)

end                         %конец внешнего вайла