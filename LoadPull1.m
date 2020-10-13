%% LoadPull
%DC 2020
clear
close all
clc
%% Set some paths
if ~sum(getenv('HPEESOF_DIR'))
    
    setenv('HPEESOF_DIR','C:\Program Files\Keysight\ADS2020');
    setenv('COMPL_DIR',getenv('HPEESOF_DIR'));
    setenv('SIMARCH','win32_64');
    
    setenv('PATH',[getenv('HPEESOF_DIR') '\bin\' ';' getenv('HPEESOF_DIR') ...
        '\lib\' getenv('SIMARCH') ';' getenv('HPEESOF_DIR') '\adsptolemy\lib.' ...
        getenv('SIMARCH') ';' getenv('HPEESOF_DIR') ...
        '\SystemVue\2014.10\win32_64\sveclient' ';' getenv('PATH')])
end
%% initial values
% inital magnitude and phase of the voltage sources
Z0 = 50;
iteracions = 0;
PavMax_dBm =25;
PavMin_dBm=1;
ztos = @(z,zref) (z-conj(zref))./(z+zref);
stoz = @(r,zref) (r*zref+conj(zref))./(1-r);
stos = @(r,z0,znew) ztos(stoz(r,zref),znew);

%% Generate the gamma
%parameters
Zcenter           = 31+ 16*i;  % center in which you want to create the gammas.
rhoMax            = 0.45;    % radius of the load-pull circle (in gamma amplitude)
nPts              = 100;      % number of loads

% gamma generation
rhoCenter         = ztos(Zcenter, Z0);
rhoMagMax         = min(rhoMax, 0.95-abs(rhoCenter));
rhoPh             = (1:nPts)*2.4;
rhoMag            = rhoMagMax*sqrt((1:nPts)/nPts);
rhoSwp            = rhoMag.*exp(1i*rhoPh) + rhoCenter;

GainDriver=[29.9 29.8 29.7 29.9 29.1 ];
%% run ADS from Matlab
% MenusADS > Simulate > GenerateNetlist > Save the netlist with the name of
% the schematic but with extension .log "for example: HB_ActiveLP.log"
FreqSwp= [1.8:0.1:2.2];
a_values = zeros(nPts,4,PavMax_dBm, 5);
Eff = zeros(PavMax_dBm,nPts,5);
Gain = zeros(PavMax_dBm,nPts,5);
Pavs_dBm = zeros(PavMax_dBm,nPts,5);
Pload_dBm = zeros(PavMax_dBm,nPts,5);
% Zin=zeros(nPts,1);
% Zout=zeros(nPts,1);
% gammaIn=zeros(nPts,1);
%%

for JJ=1:5
    
    Freq=FreqSwp(JJ);
    
    file_path = '..\FreqLP.mdf'; % name of the file (including the path)
    fid       = fopen(file_path,'w'); % open the file to Write
    
    Names    = [{'F'} {'idx'}]; % names of the parameters;
    
    Values = [Freq JJ]; % values for each parameter
    
    strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
    Nparam   = numel(Names);                    % Number of parameters to write
    
    % wirte the content
    fprintf(fid,'BEGIN block\n');
    fprintf(fid,['%% INDEX(int)\t' strnames]);
    fprintf(fid,['1\t' repmat('%f\t',1,Nparam) '\n'],Values);
    fprintf(fid,'END\n');
    fclose(fid);
    
    for KK=1:nPts
        
        P2m = 3e-9;
        P2p = 1;
        gammaMeas0 = -0.118 +1i*0.002;
        ZMeas0= stoz(gammaMeas0,Z0);
        P20= 0;
        gammaObj = rhoSwp(KK);
        ZObj= stoz(gammaObj,Z0);
        GainCheck = 0;
        GainMax =0;
        
        for QQ=1:PavMax_dBm
           
            if(QQ == 1 && KK > 1 )
                P2m = P2Maux;
                P2p = P2Paux;
            end
            if((GainMax-GainCheck) < 3)
                error = 10;
                P1m= QQ-GainDriver(JJ)+15;
                P1p = 0;
                
                if(QQ > 1)
                    P2m= P2m*1.12;
                    gammaMeas0 = -0.118 +1i*0.002;
                    ZMeas0= stoz(gammaMeas0,Z0);
                    P20=0;
                end
                P2mlog =10*log10(P2m)+30;
                file_path = '..\ParamLP_file.mdf'; % name of the file (including the path)
                fid       = fopen(file_path,'w'); % open the file to Write
                
                Names    = [{'P1m'} {'P1p'} {'P2m'} {'P2p'}]; % names of the parameters;
                
                Values = [P1m P1p P2mlog  P2p]; % values for each parameter
                
                strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
                Nparam   = numel(Names);                    % Number of parameters to write
                
                % wirte the content
                fprintf(fid,'BEGIN block\n');
                fprintf(fid,['%% INDEX(int)\t' strnames]);
                fprintf(fid,['1\t' repmat('%f\t',1,Nparam) '\n'],Values);
                fprintf(fid,'END\n');
                fclose(fid);
                
                while error> 0.0035
            
                    netPath = [fileparts(pwd) '\LoadPull.log']; % path where you saved the netlist
                    eval(['!hpeesofsim -q "' netPath '"'])
                    load ../MatlabScripts/LoadPull.mat
                    aux = who('LoadPull*');
                    eval(['RawData = ' aux{1} ';'])
                    clc
                    gammaMeas = RawData.dataBlocks.data.dependents(2);
                    P2= P2m*exp(1i*((pi*P2p)/180));
                    ZMeas= stoz(gammaMeas,Z0);
                    
                    figure(1)
                    smithchart;
                    hold on
                    plot(rhoSwp,'.k','MarkerSize',10)
                    plot(gammaMeas,'.r','MarkerSize',10)
                    title(['itr =',num2str(iteracions), '   F =', num2str(Freq),'   PdBm =', num2str(QQ), '  NPT = ', num2str(KK), '         ZObj = ', num2str(ZObj), '     error=', num2str(error)] )
                    if( gammaMeas >0.999)
                        P2m = 0.1 *P2m;
                    elseif(gammaMeas0 == gammaMeas)
                        P2m = 1.2*P2m;
                        
                    else
                        %NR
                        F0 = (gammaMeas0 - gammaObj);
                        F = (gammaMeas-gammaObj);
                        dF= (F-F0)/(P2 -P20);
                        P2New =  P2-(F/dF);
                        P2m= abs(P2New);
                        P2p = (angle(P2New)*180)/pi;
                    end
                    
                    %file
                    P2mlog =10*log10(P2m)+30;
                    file_path = '..\ParamLP_file.mdf'; % name of the file (including the path)
                    fid       = fopen(file_path,'w'); % open the file to Write
                    
                    Names    = [{'P1m'} {'P1p'} {'P2m'} {'P2p'}]; % names of the parameters;
                    
                    Values = [P1m P1p P2mlog P2p]; % values for each parameter
                    
                    strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
                    Nparam   = numel(Names);                    % Number of parameters to write
                    
                    % write the content
                    fprintf(fid,'BEGIN block\n');
                    fprintf(fid,['%% INDEX(int)\t' strnames]);
                    fprintf(fid,['1\t' repmat('%f\t',1,Nparam) '\n'],Values);
                    fprintf(fid,'END\n');
                    fclose(fid);
                    
                    %updates
                    P20 = P2;
                    gammaMeas0 = gammaMeas;
                    ZMeas0 = ZMeas;
                    error=abs(F);
                    iteracions =iteracions+1 ;
                    % dF0=dF; 
                   
                end
                
                load ../MatlabScripts/FOMs.mat
                aux = who('LoadPull*');
                eval(['RawData2 = ' aux{1} ';'])
                clc
                GainCheck = RawData2.dataBlocks(3).data.dependents(2);
                
                if(GainCheck > GainMax) % check MAX Gain
                    GainMax = GainCheck;
                end
                if(QQ == 1 )
                    P2Maux = P2m ;
                    P2Paux = P2p;
                end
                
                a_values(KK,1,QQ,JJ)= P1m;
                a_values(KK,2,QQ,JJ)= P1p;
                a_values(KK,3,QQ,JJ)= P2mlog;
                a_values(KK,4,QQ,JJ)= P2p;
                
                Eff(QQ:end,KK,JJ) = RawData2.dataBlocks(4).data.dependents(2);
                Gain(QQ:end,KK,JJ)= RawData2.dataBlocks(3).data.dependents(2);
                Pload_dBm(QQ:end,KK,JJ)= RawData2.dataBlocks(1).data.dependents(2);
                Pavs_dBm(QQ:end,KK,JJ)=  RawData2.dataBlocks(2).data.dependents(2);
                
%                 load ../MatlabScripts/Zs.mat
%                 aux = who('LoadPull*');
%                 eval(['RawData3 = ' aux{1} ';'])
%                 clc
%                 Zin(KK)=RawData3.dataBlocks(3).data.dependents(2);
%                 Zout(KK)=RawData3.dataBlocks(1).data.dependents(2);
%                 gammaIn(KK)=abs(RawData3.dataBlocks(2).data.dependents(2)); 
            end
        end
    end
 end

%%
Pav_dBm = Pavs_dBm;  
FcSwp = FreqSwp;
XdBcmp          = 3;
EffTarget       = 53:4:85;
PloadTarget     = 35:4:45;
GainRefOPBO     = 10;
gridPts = 100;
PavSwp_terp     = 15:0.1:40;

% Iload_dBn_terp
Pload_dBm_terp  = zeros(numel(PavSwp_terp),numel(rhoSwp),numel(FcSwp));
Gain_terp       = zeros(numel(PavSwp_terp),numel(rhoSwp),numel(FcSwp));
Eff_terp        = zeros(numel(PavSwp_terp),numel(rhoSwp),numel(FcSwp));
GainMax         = zeros(numel(rhoSwp),numel(FcSwp));
PavMax          = zeros(numel(rhoSwp),numel(FcSwp));
GainXdB         = zeros(numel(rhoSwp),numel(FcSwp));
PavXdB          = zeros(numel(rhoSwp),numel(FcSwp));
PloadXdB        = zeros(numel(rhoSwp),numel(FcSwp));
EffXdB          = zeros(numel(rhoSwp),numel(FcSwp));



for ff = 1:numel(FcSwp)
  for kk = 1:numel(rhoSwp)
    [Pav_ref,idx] = unique(Pav_dBm(~isnan(Pav_dBm(:,kk,ff)),kk,ff));
    
    
    Pload_ref = Pload_dBm(~isnan(Pload_dBm(:,kk,ff)),kk,ff);
    Gain_ref = Gain(~isnan(Gain(:,kk,ff)),kk,ff);
    Eff_ref = Eff(~isnan(Eff(:,kk,ff)),kk,ff);
    
    Pload_ref = Pload_ref(idx);
    Gain_ref = Gain_ref(idx);
    Eff_ref = Eff_ref(idx);
    
    Pload_dBm_terp(:,kk,ff) = interp1(Pav_ref,Pload_ref,PavSwp_terp,'pchip',NaN);
    Gain_terp(:,kk,ff)      = interp1(Pav_ref,Gain_ref,PavSwp_terp,'pchip',NaN);
    Eff_terp(:,kk,ff)       = interp1(Pav_ref,Eff_ref,PavSwp_terp,'pchip',NaN);
       %     ii                      = Pload_dBm_terp(:,kk,ff) > (Pload_ref(end)-GainRefOPBO);
    [GainMax(kk,ff),indmax] = max(Gain_terp(:,kk,ff));
    PavMax(kk,ff)           = PavSwp_terp(indmax);
    GainXdB(kk,ff)          = GainMax(kk,ff)-XdBcmp;
    idx                     = (Pload_dBm_terp(:,kk,ff) > (Pload_ref(end)-GainRefOPBO))&~isnan(Gain_terp(:,kk,ff));
    idx(1:indmax)           = false;
    if(sum(idx)==0)
        PloadXdB(kk,ff)         = NaN;
        EffXdB(kk,ff)           = NaN;
       
        
    else
        PavXdB(kk,ff)           = interp1(Gain_terp(idx,kk,ff),PavSwp_terp(idx),GainXdB(kk,ff),'pchip',NaN);
        PloadXdB(kk,ff)         = interp1(PavSwp_terp(idx),Pload_dBm_terp(idx,kk,ff),PavXdB(kk,ff),'pchip',NaN);
        EffXdB(kk,ff)           = interp1(PavSwp_terp(idx),Eff_terp(idx,kk,ff),PavXdB(kk,ff),'pchip',NaN);
              
    end
  end
end

% figure, plot(Pload_dBm_terp,Gain_terp)

rhoSwp_r                    = real(rhoSwp);
rhoSwp_i                    = imag(rhoSwp);
[rho_rgrid,rho_igrid]       = meshgrid( linspace(min(real(rhoSwp)),max(real(rhoSwp)),gridPts), ...
                                        linspace(min(imag(rhoSwp)),max(imag(rhoSwp)),gridPts));

Eff_interp                  = zeros(gridPts,gridPts,numel(FcSwp));
Pload_dBm_interp            = zeros(gridPts,gridPts,numel(FcSwp));
Gain_interp                 = zeros(gridPts,gridPts,numel(FcSwp));

FreqVct                     = FcSwp;
FreqIdxs                    = find(ismember(FcSwp,FreqVct));

% Plevels                     = repmat(PloadTarget,1,2);
Elevels                     = EffTarget; 
Plevels                     = PloadTarget;

figure(3), hold on
smithchart

for ff = 5
  figure(3), hold on
  F_Eff                     = scatteredInterpolant(rhoSwp_r',rhoSwp_i',EffXdB(:,ff));
  F_Pload_dBm               = scatteredInterpolant(rhoSwp_r',rhoSwp_i',PloadXdB(:,ff));
  F_Gain                    = scatteredInterpolant(rhoSwp_r',rhoSwp_i',GainXdB(:,ff));
  
  Eff_interp(:,:,ff)        = F_Eff(rho_rgrid,rho_igrid);
  Pload_dBm_interp(:,:,ff)  = F_Pload_dBm(rho_rgrid,rho_igrid);
  Gain_interp(:,:,ff)       = F_Gain(rho_rgrid,rho_igrid);
  
%   Ploadmax = ceil(max(max(max(PloadXdB))));
%   Effmax = ceil(max(max(max(EffXdB))));
%   [C1,h1]                   = contour(rho_rgrid,rho_igrid,Pload_dBm_interp(:,:,ff),Ploadmax-2.5:.5:Ploadmax);
%   [C2,h2]                   = contour(rho_rgrid,rho_igrid,Eff_interp(:,:,ff),Effmax-10:2:Effmax);
  [C1,h1]                   = contour(rho_rgrid,rho_igrid,Pload_dBm_interp(:,:,ff),Plevels);
  [C2,h2]                   = contour(rho_rgrid,rho_igrid,Eff_interp(:,:,ff),Elevels);
  plot(rhoSwp_r,rhoSwp_i,'k.')
  
  clabel(C1,h1,'FontSize',16,'Color','blue')
  clabel(C2,h2,'FontSize',16,'Color','red')
  h1.LineColor      = 'blue';
  h2.LineColor      = 'red';
  
figure(4)
hold on
yyaxis left
plot(Pload_dBm_terp(:,:,ff),Gain_terp(:,:,ff),'-r')
plot(PloadXdB(:,ff), GainXdB(:,ff),'ok')
axis([30 45 5 15])

yyaxis right
axis([30 45 0 90])
plot(Pload_dBm_terp(:,:,ff),Eff_terp(:,:,ff), '-b')
plot(PloadXdB(:,ff), EffXdB(:,ff),'ok')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
hold off 

  
  
end

%%
% figure(4)
% hold on
% plot(PavSwp_terp,Gain_terp,'r')
% plot(PavXdB, GainXdB,'ok')
% hold off 


%%
% index=string(1:1:100);
% file_path = '..\ZInZout.mdf'; % name of the file (including the path)
% fid       = fopen(file_path,'w'); % open the file to Write
% 
% Names    = [{'ZinR'} {'ZinC'} {'ZoutR'} {'ZoutC'}]; % names of the parameters;
% 
% Values = [real(Zin) imag(Zin) real(Zout) imag(Zout)]; % values for each parameter
% 
% strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
% Nparam   = numel(Names);                    % Number of parameters to write
% 
% 
% fprintf(fid,'BEGIN block\n');
% fprintf(fid,['%% INDEX(int)\t' strnames]);
% 
% for i=1:100
%     fprintf(fid,index(i));
%     fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
% end
% fprintf(fid,'END\n');
% fclose(fid);
%%
% Example how to write a Generic MDIF file (see ADS help of Generic MDIF and DataAccessComponent)
% index = string([1:1:KK]);
% index2 = [1:1:QQ];
% index3 =[1:1:JJ];
% file_path = '..\PARAMLP_test.mdf'; % name of the file (including the path)
% fid       = fopen(file_path,'w'); % open the file to Write
% 
% Names    = [{'P1m'} {'P1p'} {'P2m'} {'P2p'}  ]; % names of the parameters;
% 
% Values = a_values; % values for each parameter
% 
% strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n'];  %convert the names to a single string
% Nparam   = numel(Names);                    % Number of parameters to write
% 
% % wirte the content
% for j=1:JJ
%     
%     for q=1:QQ
%         fprintf(fid,'VAR Frq (int) = %d\n' ,index3(j));
%         fprintf(fid,'VAR Pin_dBm (int) = %d\n' ,index2(q));
%         fprintf(fid,'BEGIN block\n');
%         fprintf(fid,['%% INDEX(int)\t' strnames]);
%         for i=1:KK
%             fprintf(fid,index(i));
%             fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:,q,j));
%         end
%         fprintf(fid,'END\n');
%         fprintf(fid,' \n');
%     end
% end
% fclose(fid);

 %%
% index = string([1:1:KK]);
% 
% file_path = '..\gammasObj.mdf'; % name of the file (including the path)
% fid       = fopen(file_path,'w'); % open the file to Write
% 
% Names    = [{'gammar'} {'gammac'}]; % names of the parameters;
% 
% Values = [ real(rhoSwp); imag(rhoSwp)]'; % values for each parameter
% 
% strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n'];  %convert the names to a single string
% Nparam   = numel(Names);                    % Number of parameters to write
% 
% % wirte the content
% fprintf(fid,'BEGIN block\n');
% fprintf(fid,['%% INDEX(int)\t' strnames]);
% for i=1:KK
%     fprintf(fid,index(i));
%     fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
% end
% fprintf(fid,'END\n');
% fclose(fid);

