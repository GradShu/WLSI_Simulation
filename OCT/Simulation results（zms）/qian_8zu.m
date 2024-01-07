clc;clear;close all;
load('D:\Desktop\zms3\Wavelength.mat','Wavelength'); 
z=0;
M=1:2048;
N=2048;
lam1=Wavelength';
P3=polyfit(M,lam1,3);
% P3=vpa(poly2sym(P3),3)%显示三阶多项式
C3=P3(1,1);C2=P3(1,2);C1=P3(1,3);C0=P3(1,4);
I=xlsread('D:\Desktop\zms3\FANGZHAN\pinghua.xlsx');  %100个点，FFT平滑滤波后的数据
I=I-min(I);I=I';
z0=[100 200 300 400 500 600 700 800];
r=[1.0102	1.0104	1.0101	1.0102	1.0102	1.0103	1.0102	1.0102 ];                            

for i=1:8
    %%%%%%%%%%%%%%%%%%产生八组位置不同的干涉信号
    dataz=z-z0(i);
   %  dataz=100;
    S1=I.*cos((1./lam1)*dataz*4*pi)+3*I;
    
    %%%%%%%%%%%%%%%%%%计算校准前的半高宽并求均值

    %pp=[1.2*C3 C2 C1 C0];
    pp=[-1.808e-12  1.442e-09   8.143e-05   C0 ];
    % pp=[C3 C2 C1 C0];
    lam2=polyval(pp,M);
    k2=1./lam2;
    k2=linspace(k2(end),k2(1),2048);
%     lam2=linspace(lam2(1),lam2(end),2048); 
    S2=interp1(1./lam2,S1,k2,'spline');
    figure(3)
    plot(k2,S2);%S2~K2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%半高宽
    FH2=fftshift(fft(S2));
    Am2=abs(FH2);Ph2=angle(FH2);
    figure(4)
    zx=1:N;
    zx=zx-N/2-1; 
    zx=zx./(2*N*(k2(2)-k2(1))); %caiyangjainge
    plot(zx,Am2);
  
    p=100;FH2(1:1040)=0;
    P2=find(abs(FH2)==max(abs(FH2)));
    FH2(1:P2-0.5*p)=0;FH2(P2+0.5*p:N)=0;
    Am2=abs(FH2);
    
    
    FN=ifft(fftshift(FH2));                                                      
    ReF=FN;
    linear=polyfit(k2,unwrap(angle(ReF)),1);                          
    a1(1,i)=linear(1,1);                       
    zs(1,i)=a1(1,i)/4/pi;
    za(i)=zx(Am2==max(Am2));%za
    zx=zx.*(zs(1,i)./za(i));%zs
    zx=zx*(z0(i)./za(i));%真实位置dataz
%     figure(5)
%     plot(zx,Am2);
    figure(6)
    plot(zx,Am2,'LineWidth',0.8);
    hold on
   
    %Z_error(i,:)=Am2(P2-49:P2+50);
    xlabel('z (μm)','Fontname','Times New Roman','FontSize',18);
    ylabel('Amplitude (a.u.)','Fontname','Times New Roman','FontSize',18);
    xlim([0 1000]);ylim([0.00 2.5e5]);
    set(gca,'XTick',[0:200:1000],'Fontname','Times New Roman','FontSize',18);
    set(gca,'YTick',[0.00:0.5e5:2.5e5],'Fontname','Times New Roman','FontSize',18);   
    Am2=Am2(P2-20:P2+20);
     %Xdata2 = 1:size(Am2,2); % X轴平移前的坐标区间
     Xdata2 = 2.7778*(1:size(Am2,2));
    %Xdata2 =(1:size(Am2,2));
    Ydata2 = Am2; % Y轴为数据点的值
    [fitresult, gof] = createFit4(Xdata2, Ydata2);
    sigma2(i) = fitresult.c1; % 标准差
    FWHM2(i)=2*sqrt(2*log(2))*sigma2(i);
    zaa2(i)= Ydata2(Ydata2==max(Ydata2));
end
FWHM_mean2=mean(FWHM2);
FWHM_std2=std(FWHM2);
