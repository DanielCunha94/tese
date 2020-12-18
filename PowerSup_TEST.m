%% Workspace Cleanup

clear;
close all;
close force all
clc
close_all_open_instruments;



%% Configure power supply
[idn_CPX400DP, error] = CPX400DP_Init('192.168.81.33');

tbias_set = 5;

% IGS and IDS current limits
[status] = CPX400DP_SetCurrentLimit(idn_CPX400DP, 1, 0.05);        pause(tbias_set);
[status] = CPX400DP_SetCurrentLimit(idn_CPX400DP, 2, 1.5);           pause(tbias_set);

% VGS and VDS equal to zero
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 2, 0);               pause(tbias_set);
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 1, 0);               pause(tbias_set);

% Power Supllies Turn-on

[status] =  CPX400DP_SetChannelON(idn_CPX400DP, 1);                pause(tbias_set);
[status] =  CPX400DP_SetChannelON(idn_CPX400DP, 2);                pause(tbias_set);

% VGS set
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 1, 5);              pause(tbias_set);

% VDS set
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 2, 28);              pause(tbias_set);

% VGS set
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 1, 3.150);


%% Read Current
[status, current] = CPX400DP_ReadCurrent(idn_CPX400DP, 1);

%%
input('Press Enter to continue');

%% Turn off power supply

% COLD operation
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 1, 5);               pause(tbias_set);

% VGS and VDS equal to zero
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 1, 0);               pause(tbias_set);
[status] =  CPX400DP_SetVoltage(idn_CPX400DP, 2, 0);               pause(tbias_set);


% Power Supllies Turn-off
[status] = CPX400DP_SetChannelOFF(idn_CPX400DP, 1);                pause(tbias_set);
[status] = CPX400DP_SetChannelOFF(idn_CPX400DP, 2);

%% Close Instruments
close_all_open_instruments;
clear idn*
