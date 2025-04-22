function [radX SpikeTimeX] = computeRadFiring(t_ECG,t_SPIKES,ind_SPIKES)
% t_ECG=t_events_plus,
% t_SPIKES,
% ind_SPIKES=1:length(t_SPIKES);


t_CORTOpeak=t_SPIKES(ind_SPIKES)';%ECG_RESULTS.t_events;  %BELT_Rpeak;

t_LUNGOpeak=t_ECG;%t_maxima(1:end);
N_ipol=50; %interpolation

% 
% for i=1:length(t_CORTOpeak)%floor(length(t_CORTOpeak)/2)
%     
%     if t_CORTOpeak(1)<t_LUNGOpeak(1)
%         t_CORTOpeak(1)=[];
%     end
%     
%     if t_CORTOpeak(end)>t_LUNGOpeak(end)
%         t_CORTOpeak(end)=[];
%     end
%     
% end


%RR interval 
rrint=[];
tint=[];
for i=1:length(t_CORTOpeak)-1 
    rrint(i)=t_CORTOpeak(i+1)-t_CORTOpeak(i); 
    tint(i)=t_CORTOpeak(i+1);
end

%Using Maximas
%Origin = Maxima
%Using interpolation - 50
%Low-pass filtering @ 1Hz
 

Tcnt=1;
Tmat=[]; rrmat=[];
Tmat(:,:)=0; rrmat(:,:)=0; tmct=1; 

for cct1=1:length(t_LUNGOpeak)-1%:2:2:2:2
    for cct2=Tcnt:length(tint)
        if tint(cct2)<t_LUNGOpeak(cct1) 
         continue;
        elseif tint(cct2)>=t_LUNGOpeak(cct1)... 
         & tint(cct2)<t_LUNGOpeak(cct1+1)
            Tmat(tmct,cct1)=tint(cct2); 
            rrmat(tmct,cct1)=rrint(cct2); 
            tmct=tmct+1;
        else
        Tcnt=cct2; 
        tmct=1; 
        break;

        end
    end
end
 
tmat_c = {[]};
rrmat_c ={[]};
for i =1:length(t_LUNGOpeak)-1
    tmat_c{i} = tint(tint>=t_LUNGOpeak(i) & tint<t_LUNGOpeak(i+1));
    rrmat_c{i} = rrint(tint>=t_LUNGOpeak(i) & tint<t_LUNGOpeak(i+1));
end



for i=1:length(t_LUNGOpeak)-1
    try
        DTresp=(t_LUNGOpeak(i+1)-t_LUNGOpeak(i));
    for j=1:length(Tmat(:,i))
        radX(j,i)=(Tmat(j,i)-t_LUNGOpeak(i))/DTresp*2*pi;
        if  radX(j,i)<=0
            radX(j,i)=0;
        end
    end
    catch
        radX(j,i)=0;
    end
end




for i=1:length(t_LUNGOpeak)-1
    try
    
        DTresp=(t_LUNGOpeak(i+1)-t_LUNGOpeak(i));
    for j=1:length(Tmat(:,i))
        SpikeTimeX(j,i)=(Tmat(j,i)-t_LUNGOpeak(i));
        if  SpikeTimeX(j,i)<=0
            SpikeTimeX(j,i)=0;
        end
    end
    catch
        radX(j,i)=0;
    end
end

radX_c = {[]};
SpikeTimeX_c = {[]};
DTresp_v = diff(t_LUNGOpeak);
for i = 1 : length(tmat_c)
    radX_c{i} = (tmat_c{i}-t_LUNGOpeak(i))/DTresp_v(i)*2*pi;
    SpikeTimeX_c{i} = tmat_c{i}-t_LUNGOpeak(i);
end


% for i=1:length(radX(:,1))
% alfaXrad(i)= mean(nonzeros(radX(i,:)));
% end

% for i =1 : length(tmat_c)
%     check_rad(i) = isequal(radX(radX(:,i) ~= 0,i),radX_c{i}');
%     check_stx(i) = isequal(SpikeTimeX(SpikeTimeX(:,i) ~= 0,i),SpikeTimeX_c{i}');
% end

[rro0,rcol0]=size(rrmat); 
rripolmat=[]; tmipolmat=[]; 
rripolmat(:,:)=0; tmipolmat(:,:)=0; 

for rct1=1:rcol0
rripolmattemp0=[]; 
tmipolmattemp0=[]; 
rpmct=1;

for rct2=2:rro0
    if rct2==1 | rrmat(rct2,rct1)>0 
        rripolmattemp0(rpmct)=rrmat(rct2,rct1); 
        tmipolmattemp0(rpmct)=Tmat(rct2,rct1); 
        rpmct=rpmct+1;    
    else
    
    break;

    end
end
 
 
tempA0=[];tempB0=[];
    if length(tmipolmattemp0)<2 
        continue;
    else
    tempA0=tmipolmattemp0(1):((tmipolmattemp0(end)-tmipolmattemp0(1))...
    /(N_ipol-1)):tmipolmattemp0(end); 
    tempB0=spline(tmipolmattemp0,rripolmattemp0,tempA0);
    end

    for tpct=1:length(tempA0) 
        tmipolmat(tpct,rct1)=tempA0(tpct); 
        rripolmat(tpct,rct1)=tempB0(tpct);
    end
end


%CANCELLA %CANCELLA %CANCELLA
% for i=1:length(rripolmat(1,:))
%     SuperaCheck=find(abs(rripolmat(:,i))>0.1)>0;
%     if length(SuperaCheck)>0
%         rripolmat(:,i)=0.001;
%     end
% end

% rripolmat(find(abs(rripolmat)>0.2))=0;%CANCELLA


[n,m]=size(rripolmat); 
x0=[]; x0(n,1)=0;
for rct1=1:n
x0ct=0;
    for rct2=1:m
        if rripolmat(rct1,rct2)==0 
            continue;
        else
        x0(rct1,1)=x0(rct1,1)+rripolmat(rct1,rct2); 
        x0ct=x0ct+1;
        end
    end
    
    if x0ct==0
    continue; 
    else
    x0(rct1,1)=x0(rct1,1)/x0ct;
    end
    
end

%phaseDom=0:(2*pi)/(iPol-1):2*pi; % Phase Domain 0-2pi 
phaseDom=0:(2*pi)/(n-1):2*pi;

% RRI Deviation 
rrmatMean=[];
 


rridevmat=[]; rridevmat(:,:)=0; 
for rrict1=1:m
rrmatTemp=[]; rmct=1; 
    for rrict2=1:n
        if rripolmat(rrict2,rrict1)==0 
            continue;
        else
        rrmatTemp(rmct)=rripolmat(rrict2,rrict1); 
        rmct=rmct+1;
        end
    end
rrmatMean(rrict1)=mean(rrmatTemp);

    for rrict3=1:length(rrmatTemp) 
        rridevmat(rrict3,rrict1)=rrmatTemp(rrict3)-rrmatMean(rrict1);
    end
    
end


x0dev(n,1)=0; 
for rct1=1:n
x0ct=0;
    for rct2=1:m
        if rridevmat(rct1,rct2)==0 
            continue;
        else
        x0dev(rct1,1)=x0dev(rct1,1)+(rridevmat(rct1,rct2)); 
        x0ct=x0ct+1;
        end
    end
    
    if x0ct==0
    continue; 
    else
    x0dev(rct1,1)=x0dev(rct1,1)/x0ct;
    end
    
end



RSAampInit=max(x0dev)-min(x0dev); % RSA amplitude before removing outliers

phaseDelayInit=0;
for xdct=1:length(x0dev)
    
    if x0dev(xdct)==min(x0dev)
        phaseDelayInit=phaseDom(xdct);
     break;
    end
    
end
 


phaseDelay2Init=0;
for xdct2=1:length(x0dev)
    if x0dev(xdct2)==max(x0dev)
        phaseDelay2Init=phaseDom(xdct2);
        break;
    end
end

debugFIG1=0;
if debugFIG1==1;
fig1=figure;

subplot(221)
plot(phaseDom,rripolmat); xlabel('Heart cycle phase [rad]');... 
    ylabel('Spike Intervals (s)');
subplot(222)
plot(phaseDom,rridevmat*1000); xlabel('Heart cycle phase [rad]');... 
    ylabel('Spike Intervals deviation (ms)');
subplot(223)
plot(phaseDom,x0); xlabel('Heart cycle phase [rad]');... 
    ylabel('Avg Spike Intervals (s)');
subplot(224)
plot(phaseDom,x0dev*1000); xlabel('Heart cycle phase [rad]');... 
    ylabel('Avg Spike Intervals deviation (ms)');
end

%Hilbert Transform RSA 
htAvRSA=hilbert(x0); 
magRSA=abs(htAvRSA); % magnitude 
phsRSA=angle(htAvRSA); % Phase
phi0=angle(max(htAvRSA)); % phase of maximum value



debugAFTER=0;
if debugAFTER==1

% Variance 
[x0R,x0C]=size(x0); 
Vj(1,m)=0;
for i1=1:m
    for i2=1:n
        if rripolmat(i2,i1)==0 
            continue;
        elseif i1>x0R
        Vj(1,i1)=0;

        else
        Vj(1,i1)=Vj(1,i1)+((x0(i1,1)-rripolmat(i2,i1)).^2);

        end 
    end
end
 
try 
%Phase Respiration
htRes=hilbert(RespFilter);
phiJ=[]; Rmct=1; %maxRPmat=[]; phase of jth respiration cycle max value 
Tcnt=1;
for cct3=1:length(t_LUNGOpeak)-1
mmHTmat=[]; Hmct=1;
    for cct4=Tcnt:length(time)
        if time(cct4)<t_LUNGOpeak(cct3) 
            continue;
        elseif time(cct4)>=t_LUNGOpeak(cct3)... 
                & time(cct4)<t_LUNGOpeak(cct3+1)
        mmHTmat(Hmct)=htRes(cct4); 
        Hmct=Hmct+1;
        else
        Tcnt=cct4; 
        break;
        end
    end

phiJ(Rmct)=angle(max(mmHTmat)); 
Rmct=Rmct+1;

end
 
 

%Phase difference 
deltaPHI=[];
for i=1:m
deltaPHI(i)=abs(phi0-phiJ(i));
end

%Sort 
sortVj=sort(Vj);
delVj=sortVj(round(0.9*m):end);
mcol=[]; mcct=1; %deleted 'm' collections 
for dt1=1:length(delVj)
    for dt2=1:length(Vj)
        if Vj(1,dt2)==delVj(dt1) 
            mcol(mcct)=dt2;

        mcct=mcct+1; 
        break;
        else
        continue;

        end
    end
end
 

deltaPHI2=[]; dPct=1; dPstart=1; 
for dt1=1:length(mcol)+1
    for dt2=dPstart:length(deltaPHI)
        if dt1<=length(mcol) & dt2==mcol(dt1) 
            dPstart=dt2+1;
        break; 
        else
        deltaPHI2(dPct)=deltaPHI(dt2); 
        dPct=dPct+1;


        end
    end
end
 

sortdeltaPHI2=sort(deltaPHI2);
delDPHI=sortdeltaPHI2(round(0.8*m):end);
for dt1=1:length(delDPHI)
    for dt2=1:length(deltaPHI) % Note the use of deltaPHI and not deltaPHI2 
        if deltaPHI(dt2)==delDPHI(dt1)
        mcol(mcct)=dt2; mcct=mcct+1; 
        break;
        else
        continue;
        end
    end
end
 

% Subset RR intervals 
rrmatsub(:,:)=0; rrsct3=1; 
for rrsct1=1:m
entr=0;
    for rrsct2=1:length(mcol) 
        if rrsct1==mcol(rrsct2)
        entr=1;
        break;
        
        else
        continue;
        end
        
    end
        if entr==1
        continue; 
        
        else
            for dmat=1:n
            rrmatsub(dmat,rrsct3)=rripolmat(dmat,rrsct1);
            end
        rrsct3=rrsct3+1;

    end
end 

 %########################################################################

[n1,m1]=size(rrmatsub); 
x1(n1,1)=0;
for i1=1:n1
x1ct=0;
    for i2=1:m1
        if rrmatsub(i1,i2)==0 
            continue;
        else
        x1(i1,1)=x1(i1,1)+rrmatsub(i1,i2); 
        x1ct=x1ct+1;
        end
    end
x1(i1,1)=x1(i1,1)/x1ct;
end

% RRI Deviation 
rrmatsubMean=[]; rridevmatsub(:,:)=0; 
for rrict1=1:m1
rrmatsubTemp=[]; rmct=1; 
    for rrict2=1:n1
        if rrmatsub(rrict2,rrict1)==0 
            continue;
        else
        rrmatsubTemp(rmct)=rrmatsub(rrict2,rrict1); rmct=rmct+1;
        end
    end
    
rrmatsubMean(rrict1)=mean(rrmatsubTemp);

    for rrict3=1:length(rrmatsubTemp) 
        rridevmatsub(rrict3,rrict1)=rrmatsubTemp(rrict3)-rrmatsubMean(rrict1);
    end
end


[n11,m11]=size(rridevmatsub); x1dev(n11,1)=0;
for rct1=1:n11 x1ct=0;
    for rct2=1:m11
        if rridevmatsub(rct1,rct2)==0 
            continue;
        else
        x1dev(rct1,1)=x1dev(rct1,1)+(rridevmatsub(rct1,rct2)); 
        x1ct=x1ct+1;
        end
    end
x1dev(rct1,1)=x1dev(rct1,1)/x1ct;
end

RSAamp=max(x1dev)-min(x1dev); % RSA amplitude 
phaseDom2=0:(2*pi)/(n11-1):2*pi;
phaseDelay=0;

for xdct=1:length(x1dev)
    if x1dev(xdct)==min(x1dev) 
        phaseDelay=phaseDom2(xdct); 
        break;
    end
end

phaseDelay2=0;
for xdct2=1:length(x1dev)
    if x1dev(xdct2)==max(x1dev) 
        phaseDelay2=phaseDom2(xdct2); 
        break;
    end
end


figure

% figure(11);
subplot(221)
plot(phaseDom2,rrmatsub); xlabel('Respiratory phase [rad]');... 
    ylabel('Avg RR Intervals (s)');


% figure(12); 
subplot(222)
plot(phaseDom2,rridevmatsub*1000);...
xlabel('Respiratory phase [rad]'); ylabel('RRI deviation (ms)');

% figure(13);
subplot(223)
plot(phaseDom2,x1); xlabel('Respiratory phase [rad]');... 
    ylabel('Avg RR Intervals (s)');

% figure(14);
subplot(224)
plot(phaseDom2,x1dev*1000); xlabel('Respiratory phase [rad]');... 
    ylabel('Avg RRI deviation (ms)');


end
end






end