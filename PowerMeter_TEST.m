%% Workspace Cleanup

clear;
close all;
close force all
clc
close_all_open_instruments;

%% Configure Power Meter 

freq = 2.1e9; 

[idn_N1913] = N1913_Init('192.168.81.38');

[error1]    = N1913_Config(idn_N1913, freq); pause(1)

[error, aux] = N1913_ReadPower(idn_N1913);

%% Close Instruments
close_all_open_instruments;
clear idn*