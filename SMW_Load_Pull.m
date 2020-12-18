%% SMW Load Pull
% Daniel

%%
clear;
close all;
close force all
clc
close_all_open_instruments;

%% Generate gamma

Z0 = 50;
iteracions = 0;
Pav_dBm = [0:1:15];
ztos = @(z,zref) (z-conj(zref))./(z+zref);
stoz = @(r,zref) (r*zref+conj(zref))./(1-r);
stos = @(r,z0,znew) ztos(stoz(r,zref),znew);

%parameters
Zcenter           = 25+1i*25;  % center in which you want to create the gammas.
rhoMax            = 0.4;       % radius of the load-pull circle (in gamma amplitude)
nPts              = 100;       % number of loads

% gamma generation
rhoCenter         = ztos(Zcenter, Z0);
rhoMagMax         = min(rhoMax, 0.95-abs(rhoCenter));
rhoPh             = (1:nPts)*2.4;
rhoMag            = rhoMagMax*sqrt((1:nPts)/nPts);
rhoSwp            = rhoMag.*exp(1i*rhoPh) + rhoCenter;

%% Calibration

load ../Calibration.mat
e00     = Calibration.e00;
e10e01  = Calibration.e10e01;
e11     = Calibration.e11;
e33r    = Calibration.e33r;
e23e32  = Calibration.e23e32;
e22r    = Calibration.e22r;
e22     = Calibration.e22;
e10e32  = Calibration.e10e32;
e11r    = Calibration.e11r;
e23e01  = Calibration.e23e01;

D    = @(S11Um, S22Um, S21Um, S12Um)     (1+((S11Um-e00)./(e10e01)).*e11).*(1+((S22Um-e33r)./(e23e32)).*e22r)-(S21Um./(e10e32)).*(S12Um./e23e01).*e22.*e11r;
S11C = @(S11Um, S22Um, S21Um, S12Um, D)  (((S11Um-e00)./e10e01).*(1+((S22Um-e33r)./e23e32).*e22r)-e22.*(S21Um./e10e32).*(S12Um./e23e01))./D;
S21C = @(S22Um, S21Um, D)                ((S21Um./e10e32).*(1+((S22Um-e33r)./e23e32).*(e22r-e22)))./D;
S12C = @(S11Um, S12Um, D)                ((S12Um/e23e01).*(1+((S11Um-e00)./e10e01).*(e11-e11r)))./D;
S22C = @(S11Um, S22Um, S21Um, S12Um, D)  (((S22Um-e33r)./e23e32).*(1+((S11Um-e00)./e10e01).*e11)-e11r.*(S21Um/e10e32).*(S12Um/e23e01))./D;


%% VSG Configuration

G.freqcStart    = 2.1e9;
G.powerCold     = -100;
G.Power1        = 0;
G.Power2        = 0;
G.phase         = 0;

[st, idn_SMW]   = SMW_connect('192.168.81.65');
[error1]        = SMW_preset(idn_SMW);
[error1]        = SMW_set_RF_ONOFF(idn_SMW,'OFF',1);
[error1]        = SMW_set_RF_ONOFF(idn_SMW,'OFF',2);
[error1]        = SMW_set_power(idn_SMW, G.powerCold, 1);
[error1]        = SMW_set_freq(idn_SMW,  G.freqcStart, 1);
[error1]        = SMW_set_power(idn_SMW, G.powerCold, 2);
[error1]        = SMW_set_freq(idn_SMW,  G.freqcStart, 2);
[error1]        = SMW_set_PhaseREF(idn_SMW, 1);
[error1]        = SMW_set_phase(idn_SMW, G.phase,2);

[error1]        = SMW_set_power(idn_SMW, G.Power, 1);
[error1]        = SMW_set_power(idn_SMW, G.Power, 2);

%% PNA Configuration

[~,idn_PNA]     = PNA_connect('192.168.81.201',5);

X.fileFormat    = 'RI';
X.byteOrder     = 'SWAP';
X.dataTypePNA   = 'REAL,64';
WaitForSystemReady(idn_PNA);

% Configure parameter to be acquired and initiate
ConfigurePNA(idn_PNA,X);

% Set up the measurements
PNA_CreateMeasRaw(idn_PNA)

% Select measurements
fprintf(idn_PNA,'CALC:PAR:SEL ''b_wave_1''');
fprintf(idn_PNA,'CALC:PAR:SEL ''a_wave_1''');
fprintf(idn_PNA,'CALC:PAR:SEL ''b_wave_2''');
fprintf(idn_PNA,'CALC:PAR:SEL ''a_wave_2''');

%% LOAD PULL

input('Start - Press Enter to continue');
[error1]        = SMW_set_RF_ONOFF(idn_SMW,'ON',1);
[error1]        = SMW_set_RF_ONOFF(idn_SMW,'ON',2);

gammaMeas_values    = zeros(nPts,size(Pav_dBm));

for KK=1:nPts
    
    G.Power2        = 0;
    G.phase         = 5;
    [error1]        = SMW_set_power(idn_SMW, G.Power2, 2);
    [error1]        = SMW_set_phase(idn_SMW, G.phase,2);        
    gammaMeas0      = 6.617e-22 +1i*1e-5;
    G2_0            = 0+1i*0;
    gammaObj        = rhoSwp(KK);
    ZObj            = stoz(gammaObj,Z0);
    
    for QQ=1:size(Pav_dBm)
        
       
       %   this is to not start the  NR method always in the center of the Smith Chart 
        if(QQ == 1 && KK > 1 )
            G.Power2        = Power_aux;
            G.phase         = Phase_aux;
            [error1]        = SMW_set_power(idn_SMW, G.Power2, 2);
            [error1]        = SMW_set_phase(idn_SMW, G.phase,2);
        end
        
        error = 10;
        
        G.Power1            = Pav_dBm(QQ);
        [error1]            = SMW_set_power(idn_SMW, G.Power1, 1);
        
        %this increases the power of g2 in 1bBm when g1 power increases 1dBm
        if(QQ > 1)
            G.Power2        = 1 + G.Power2;
            [error1]        = SMW_set_power(idn_SMW, G.Power2, 2);
            gammaMeas0      = 0;
            G2_0            = 0;
        end
        
        while error > 0.003
            
            % get Meas
             ConfigurePNA(idn_PNA,X);
             b1              = get_pna_data(idn_PNA,'"b_wave_1"');
             a1              = get_pna_data(idn_PNA,'"a_wave_1"');
             b2              = get_pna_data(idn_PNA,'"b_wave_2"');
             a2              = get_pna_data(idn_PNA,'"a_wave_2"');
             clrdevice(idn_PNA);
            % calib Meas
             Den             = D(b1/a1, b2/a2, b2/a1, b1/a2);
             S22             = S22C(b1/a1, b2/a2, b2/a1, b1/a2, Den);
             gammaMeas       = 1/S22;
             
%             Pin            = 0.5*abs(a1)^2  - 0.5*abs(b1)^2;
%             Pout           = 0.5*abs(b2)^2 - 0.5*abs(a2)^2;
            G2               = G.Power2*exp(1i*((pi*G.phase)/180));
            
            figure(1)
            smithchart;
            hold on
            plot(rhoSwp,'.k','MarkerSize',10)
            plot(gammaMeas,'.r','MarkerSize',10)
            title(['iteracion =',num2str(iteracions), '  PdBm =', num2str(QQ), '  NPT = ', num2str(KK), '   ZObj = ', num2str(ZObj), '    error=', num2str(error)] )
            
            if(gammaMeas0 == gammaMeas)
                % creates a disturbance
                G.Power2        = 1.05 *G.Power2;
                [error1]        = SMW_set_power(idn_SMW, G.Power2, 2);
            else
                %NR ALG
                F0              = gammaMeas0 - gammaObj;
                F               = gammaMeas-gammaObj;
                df              = (F-F0)/(G2 - G2_0);
                G2New           = G2-(F/(dF));
                G.Power2        = abs(G2New);
                G.phase         = (angle(G2New)*180)/pi;
                [error1]        = SMW_set_power(idn_SMW, G.Power2, 2);
                [error1]        = SMW_set_phase(idn_SMW, G.phase,2);
            end
            
            %updates
            G2_0                = G2;
            gammaMeas0          = gammaMeas;
            Zmeas               = 50*((1+gammaMeas)/(1-gammaMeas));
            error               = abs(F);
            iteracions          = iteracions+1 ;
        end
        
        gammaMeas_values(KK,QQ) = gammeMeas; 
        
        if(QQ == 1 )
            Power_aux       = G.Power2 ;
            Phase_aux       = G.phase;
        end
        
    end
end

%% VSG Shutdown

input('VSG Shutdown - Press Enter to continue');

[error1]             = SMW_set_RF_ONOFF(idn_SMW,'OFF',1);
[error1]             = SMW_set_power(idn_SMW, G.powerCold,1);
[error1]             = SMW_set_RF_ONOFF(idn_SMW,'OFF',2);
[error1]             = SMW_set_power(idn_SMW, G.powerCold,2);

%% Close Instruments
close_all_open_instruments;
clear idn*

%% Plot Results

figure(2)
smithchart;
hold on
plot(rhoSwp,'.k','MarkerSize',10)
plot(gammaMeas_values,'.r','MarkerSize',10)
