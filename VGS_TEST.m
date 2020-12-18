%% Workspace Cleanup

clear;
close all;
close force all
clc
close_all_open_instruments;

%% VSG Configuration
G.freqcStart    = 2.1e9;
G.powerCold     = -100;
G.Power         = -10;
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
%%
[error1]        = SMW_set_power(idn_SMW, G.Power, 1);
[error1]        = SMW_set_power(idn_SMW, G.Power, 2);

input('Press Enter to continue');


[error1]        = SMW_set_RF_ONOFF(idn_SMW,'ON',1);
[error1]        = SMW_set_RF_ONOFF(idn_SMW,'ON',2);


%% VSG Shutdown

input('Press Enter to continue');

[error1]             = SMW_set_RF_ONOFF(idn_SMW,'OFF',1);
[error1]             = SMW_set_power(idn_SMW, G.powerCold,1);
[error1]             = SMW_set_RF_ONOFF(idn_SMW,'OFF',2);
[error1]             = SMW_set_power(idn_SMW, G.powerCold,2);
%% Close Instruments
close_all_open_instruments;
clear idn*
