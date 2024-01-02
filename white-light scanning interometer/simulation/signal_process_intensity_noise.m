clear;close all
%The simulation use for signal processing in wavenumber scanning domains
% of white-light interferometers. In the simulation, The IFT(Inverse Fourie
% -r Transform) of I can be get.
%Creat guass function:
%      iSigma =  a1*exp(-((sigma-b1)/c1)^2)
%Coefficients (with 95% confidence bounds):
%        a1 =       568.8
%        b1 =        1.45  
%        c1 =      0.3129
%Creat a phase functoion:
%      pSigma = -4*pi*z0*sigma;

%Related parameters
%intensityMark,phaseMark,ifftIntenMark,ifftPhaseMark:image identification
intensityMark = 221;
phaseMark = 223;
ifftIntenMark = 222;
ifftPhaseMark = 224;
%npoint:sampling point 
nPoint = 2^10;
z0 = 10;
sigma = 0:0.0083:(nPoint-1)*0.0083;
fSigma = fitfSigma(sigma);
pSigma1 = 313.5*sigma-905.1;
pSigma2 = (19-4*3.14*10)*sigma-153.6;
pSigma3 = (93.61-4*3.14*10)*sigma-768;
% pSigma1=-4*pi*((16+z0)-6*sigma).*sigma;
% iSigma = fSigma.*exp(1i*pSigma);
%calculate ifft of iSigma
% [intensity,phase,zDataInter]=ifftOfIsigma(iSigma,sigma(1),sigma(end),nPoint);
% zeroPhase = findZeroPhase(zDataInter,phase,0);

% figure(1);
% subplot(2,2,1);
% setImage(221,sigma,fSigma,'Wavenumber \sigma(\mum^{-1})','Intensity(a.u.)');
% subplot(2,2,3);
% setImage(223,sigma,pSigma,'Wavenumber \sigma(\mum^{-1})','Phase(rad)');
hold on
plot(sigma,pSigma1);
plot(sigma,pSigma2);
plot(sigma,pSigma3);

% subplot(2,2,2);
% setImage(222,zDataInter,intensity,'Position z(\mum)','Intensity(a.u.)');
% subplot(2,2,4);
% setImage(224,zDataInter,phase,'Position z(\mum)','Phase(rad)',zeroPhase);%

function [intensity,phase,zDataInter]=ifftOfIsigma(input,sigmaStart,sigmaEnd,nPoint)
%This function use for calculating ifft of iSigma.
%The output is intensity and phase after ifft.
deltaSigma = (sigmaEnd-sigmaStart)/(nPoint-1);
deltaZ = 1/(2*nPoint*deltaSigma);
zData = (1:nPoint)*deltaZ;
intensity = abs(ifft(input));
phase = unwrap(angle(ifft(input)));

%interpolation
zDataInter = linspace(zData(1),zData(end),nPoint*10);
 intensity = interpn(zData,intensity,zDataInter,'linear');
phase = wrapToPi(interpn(zData,phase,zDataInter,'linear'));
end

function zeroPhase = findZeroPhase(zData,phaseData,zeroPoint)
diff = 0.5;
%limit the range of zData,about (zp-lambda/8,zp+lambda/8)
LimZData = zData(zData>(35.3723-(1/1.45)/8)&zData<(35.3723+(1/1.45)/8));
phaseData = phaseData(zData>(35.3723-(1/1.45)/8)&zData<(35.3723+(1/1.45)/8));
% phaseData = phaseData()
positivePhase = phaseData(phaseData>0);
negativePhase = phaseData(phaseData<0);
while true 
    if length(find(abs(positivePhase-zeroPoint)<diff))>1
       lastDiff = diff;
       diff = diff/2;
    elseif (isempty(find(abs(positivePhase-zeroPoint)<diff, 1)))
        diff = (diff + lastDiff)/2;
    else
        zeroPhase1= positivePhase(abs(positivePhase-zeroPoint)<diff);
        break
    end
end
while true 
    if length(find(abs(negativePhase-zeroPoint)<diff))>1
       lastDiff = diff;
       diff = diff/2;
    elseif (isempty(find(abs(negativePhase-zeroPoint)<diff, 1)))
        diff = (diff + lastDiff)/2;
    else
        zeroPhase2= negativePhase(abs(negativePhase-zeroPoint)<diff);
        break
    end
end
zeroData1 = LimZData(phaseData==zeroPhase1);
zeroData2 = LimZData(phaseData==zeroPhase2);
zeroPhase = ((0-zeroPhase2)*(zeroData1-zeroData2))/(zeroPhase1-zeroPhase2)...
            +zeroData2;
end

function setImage(imageMark,xData,yData,xlabel,ylabel,zeroPhase)
%This function setImage use for setting image property, such as axes, line,
%colour and so on. 
%ampMaxZ:The abscissa where the maximum amplitude value is located.
global ampMaxZ
yDataLine = plot(xData,yData);
%Find the max of Intensity.
[yDataMax,yDataMaxindex]= max(yData);
%Set the location of the tick marks along the axis
if imageMark == 221
    hold on
    ampMaxZ = xData(yDataMaxindex);
    yMaxLine = plot([1.5438,1.5438],...
               [0,yData(187)]);
    set(gca,'Xlim',   [1.1 2.1]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(yMaxLine, 'color',       'k',...
                  'LineStyle',  '--',...
                  'LineWidth',   2);
elseif imageMark == 223
    set(gca,'Xlim',   [1.1 2.1]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
elseif imageMark == 222
    hold on
    ampMaxZ = xData(yDataMaxindex);
    yMaxLine = plot([xData(yDataMaxindex),xData(yDataMaxindex)],...
               [0,yDataMax]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(yMaxLine, 'color',       'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
    set(gca,'Xlim',   [xData(yDataMaxindex)-2 xData(yDataMaxindex)+2]);
else
    hold on
    zeroLine = plot([ampMaxZ-5,ampMaxZ+5],[0,0]);
    zeroPhaseLine = plot([zeroPhase,zeroPhase],[-3.14,0]);
    set(gca,'Xlim', [ampMaxZ-2 ampMaxZ+2]);
    set(gca,'Ylim', [-3.14 3.14]);
    set(yDataLine,'color',      'blue',...
                  'LineStyle',  '-',...
                  'LineWidth',   2);
    set(zeroLine, 'color',       'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
    set(zeroPhaseLine, 'color',  'k',...
                  'LineStyle',   '--',...
                  'LineWidth',    2);
end

set(gca,'Xcolor',     [0 0 0],...
        'Ycolor',     [0 0 0],...
        'Color' ,     [1 1 1],...
        'FontName',   'Times New Roman',...
        'FontAngle',  'normal',...
        'FontSize',    13);
%Set label names for X coordinates and Y coordinates
set(get(gca,'XLabel'), 'String',xlabel,...
                       'FontName','Times New Roman',...
                       'FontSize',14);
set(get(gca,'YLabel'), 'String',ylabel,...
                       'FontName','Times New Roman',...
                       'FontSize',14);
end

function fSigma = fitfSigma(x)
%Coefficients (with 95% confidence bounds):
       a1 =       88.19  ;%(51.68, 124.7)
       b1 =       1.378  ;%(1.347, 1.409)
       c1 =      0.1243  ;%(0.0879, 0.1606)
       a2 =      0.9322  ;%(-1.209, 3.073)
       b2 =       1.542  ;%(1.518, 1.565)
       c2 =     0.01341  ;%(-0.02512, 0.05195)
       a3 =       69.18  ;%(-6.223, 144.6)
       b3 =        1.27  ;%(1.228, 1.312)
       c3 =     0.07202  ;%(0.03205, 0.112)
       a4 =       528.7  ;%(493.5, 563.8)
       b4 =       1.488  ;%(1.483, 1.494)
       c4 =      0.2199  ;%(0.2153, 0.2245)
       a5 =       62.39  ;%(-14.96, 139.7)
       b5 =       1.207  ;%(1.188, 1.225)
       c5 =     0.05735  ;%(0.04772, 0.06698)
       a6 =         181  ;%(163.3, 198.6)
       b6 =        1.87  ;%(1.867, 1.873)
       c6 =       0.158  ;%(0.153, 0.1629)
       a7 =      -60.41  ;%(-80.21, -40.6)
       b7 =        1.92  ;%(1.789, 2.051)
       c7 =      0.4102  ;%(0.3363, 0.4841)
       a8 =       57.07  ;%(54.59, 59.55)
       b8 =       1.735  ;%(1.728, 1.742)
       c8 =       1.039  ;%(1.023, 1.054)
     fSigma = a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2) +... 
              a3*exp(-((x-b3)/c3).^2) + a4*exp(-((x-b4)/c4).^2) + ...
              a5*exp(-((x-b5)/c5).^2) + a6*exp(-((x-b6)/c6).^2) + ...
              a7*exp(-((x-b7)/c7).^2) + a8*exp(-((x-b8)/c8).^2);
end
