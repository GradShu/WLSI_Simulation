clc;clear;close all;
load('D:\Desktop\zms3\Wavelength.mat','Wavelength'); 
load('D:\Desktop\zms3\FANGZHAN\8zunon.mat');
non=non(1,:);
z=0;
M=1:2048;
N=2048;
lam1=Wavelength';
P3=polyfit(M,lam1,3);
% P3=vpa(poly2sym(P3),3)%显示三阶多项式
C3=P3(1,1);C2=P3(1,2);C1=P3(1,3);C0=P3(1,4);
I=xlsread('D:\Desktop\zms3\FANGZHAN\pinghua.xlsx');            %100个点，FFT平滑滤波后的数据
I=I-min(I);I=I';

z0=[100 200 300 400 500 600 700 800];

for ii=1:1
    dataz=z-z0(ii);
    S1=I.*cos((1./lam1)*dataz*4*pi)+3*I;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%失准的干涉信号%%%%%%%%%%%%%%%%%%%%%%%%%%
   % % pp=[C3 C2 1*C1 C0];
     pp=[-1.808e-12  1.442e-09   8.143e-05   C0 ];

    lam2=polyval(pp,M);
    k2=1./lam2;
    k2=linspace(k2(end),k2(1),2048);
    S2=interp1(1./lam2,S1,k2,'spline');
    %%%%%%%%%%%%%%%%%%%%%%%%%半高宽
    FH2=fftshift(fft(S2));
    Am2=abs(FH2);Ph2=angle(FH2);
    figure(4)
    plot(Am2);
    p=40;FH2(1:1040)=0;
    P2=find(abs(FH2)==max(abs(FH2)));
    FH2(1:P2-0.5*p)=0;FH2(P2+0.5*p:N)=0;
    Am2=abs(FH2);
    zx=1:N;
    zx=zx-N/2-1; 
    zx2=zx./(2*N*(k2(2)-k2(1))); %caiyangjainge
    FN=ifft(fftshift(FH2));                                                      
    ReF=FN;
    linear=polyfit(k2,unwrap(angle(ReF)),1);                          
    a1(1,ii)=linear(1,1);                       
    zs2(1,ii)=a1(1,ii)/4/pi;
    za2(ii)=zx2(Am2==max(Am2));%za
    figure(77)
      %zx2=zx2.*(zs2(1,ii)./za2(ii));%zs
    zx2=zx2*(z0(ii)./za2(ii)); %%真实位置dataz
    plot(zx2,Am2,'LineWidth',0.8);
    hold on
    Am2=Am2(P2-20:P2+20);
    Xdata2 =1*( 1:size(Am2,2)); % X轴平移前的坐标区间
    Ydata2 = Am2; % Y轴为数据点的值
    [fitresult, gof] = createFit4(Xdata2, Ydata2);
    sigma2(ii) = fitresult.c1; % 标准差
    FWHM2(ii)=2*sqrt(2*log(2))*sigma2(ii);
    zaa2(ii)= Ydata2(Ydata2==max(Ydata2));
    
     %%%%%%%%%%%%%%%%%%非线性相位 %%%%%%%%%%%%%%%%%%
    F=fftshift(fft(S2));
    Am=abs(F); Ph=angle(F);
    F(1:1040)=0;
    P=find(abs(F)==max(abs(F)));
    F(1:P-0.5*p)=0;F(P+0.5*p:N)=0;
    FN=ifft(fftshift(F));
    %%%%%%%%%%%%%%%%%%%%%%%%校准
    c=7.89e-4;
    lam3=lam2-c*non;
    k3=1./lam3;
    k3=linspace(k3(end),k3(1),2048);
    S3=interp1(1./lam3,S1,k3,'spline');%S1~k3,1/λ3附近插值成k3
    FH3=fftshift(fft(S3));
    Am3=abs(FH3);Ph3=angle(FH3);
    zx=1:N;
    zx=zx-N/2-1; 
    zx=zx./(2*N*(k3(2)-k3(1))); %caiyangjainge
   
    FH3(1:1040)=0;
    P3=find(abs(FH3)==max(abs(FH3))); 
    FH3(1:P3-0.5*p)=0;FH3(P3+0.5*p:N)=0;
    
    Am3=abs(FH3);
    FN=ifft(fftshift(FH3));                                                      
    ReF=FN;
    linear=polyfit(k3,unwrap(angle(ReF)),1);                          
    a1(1,ii)=linear(1,1);                       
    zs(1,ii)=a1(1,ii)/4/pi;
    za(ii)=zx(Am3==max(Am3));

%    %  jiaozhun_fw(ii,:)=Am3(P3-49:P3+50);
    figure(11)
      zx=zx.*(zs(1,ii)./za(ii));%zs
      % % zx=zx*(z0(ii)./za(ii));%真实值
    plot(zx,Am3,'LineWidth',0.8);
    hold on
    
    xlabel('z (μm)','Fontname','Times New Roman','FontSize',18);
    ylabel('Amplitude (a.u.)','Fontname','Times New Roman','FontSize',18);
    xlim([0 1000]);ylim([0.00 2.5e5]);
    set(gca,'XTick',[0:200:1000],'Fontname','Times New Roman','FontSize',18);
    set(gca,'YTick',[0.00:0.5e5:2.5e5],'Fontname','Times New Roman','FontSize',18); 
    

    Am3=Am3(P3-10:P3+10);
 
    %Xdata3 = 1:size(Am3,2); % X轴平移前的坐标区间
    Xdata3 = 1*(1:size(Am3,2)) ;
    Ydata3 = Am3; % Y轴为数据点的值
    [fitresult, gof] = createFit4(Xdata3, Ydata3);  
    sigma3(ii) = fitresult.c1; % 标准差
    FWHM3(ii)=2*sqrt(2*log(2))*sigma3(ii);
    zaa3(ii)= Ydata3(Ydata3==max(Ydata3));
end
FWHM_mean2=mean(FWHM2);
FWHM_std2=std(FWHM2);
FWHM_mean3=mean(FWHM3);
FWHM_std3=std(FWHM3);