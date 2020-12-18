%% Workspace Cleanup
clear all; close all force; clc
close_all_open_instruments;

%% Connect to interface object, idn_PNA
[~,idn_PNA] = PNA_connect('192.168.81.201',5);

X.fileFormat = 'RI';
X.byteOrder = 'SWAP';
X.dataTypePNA = 'REAL,64';
WaitForSystemReady(idn_PNA);
ConfigurePNA(idn_PNA,X);
numOfPoints = query(idn_PNA, 'SENS:SWE:POIN?','%s\n','%d');

%% GET S-Parameters

fprintf(idn_PNA, 'DISP:WIND ON');
fprintf(idn_PNA, 'CALC:PAR:DEF "CH1_S11_1", S11');
fprintf(idn_PNA, 'DISP:WIND:TRAC:FEED "CH1_S11_1"');
fprintf(idn_PNA, 'CALC:PAR:SEL "CH1_S11_1"');
fprintf(idn_PNA, 'CALC:FORM SMIT');
fprintf(idn_PNA, 'CALC:MARK1 ON');
fprintf(idn_PNA, 'CALC:MARK1:X 2.1e9');
fprintf(idn_PNA, 'DISP:WIND:ANN:MARK:SIZE LARG');

fprintf(idn_PNA, 'DISP:WIND2 ON');
fprintf(idn_PNA, 'CALC:PAR:DEF "CH1_S21_1", S21');
fprintf(idn_PNA, 'DISP:WIND2:TRAC2:FEED "CH1_S21_1"');
fprintf(idn_PNA, 'CALC:PAR:SEL "CH1_S21_1"');
fprintf(idn_PNA, 'CALC:FORM MLOG');
fprintf(idn_PNA, 'CALC:MARK2 ON');
fprintf(idn_PNA, 'CALC:MARK2:X 2.1e9');
fprintf(idn_PNA, 'DISP:WIND2:ANN:MARK:SIZE LARG');

fprintf(idn_PNA, 'CALC:PAR:DEF "CH1_S12_1", S12');
fprintf(idn_PNA, 'DISP:WIND2:TRAC3:FEED "CH1_S12_1"');
fprintf(idn_PNA, 'CALC:PAR:SEL "CH1_S12_1"');
fprintf(idn_PNA, 'CALC:FORM MLOG');
fprintf(idn_PNA, 'CALC:MARK2 ON');
fprintf(idn_PNA, 'CALC:MARK2:X 2.1e9');
fprintf(idn_PNA, 'DISP:WIND2:ANN:MARK:SIZE LARG');

fprintf(idn_PNA, 'CALC:PAR:DEF "CH1_S22_1", S22');
fprintf(idn_PNA, 'DISP:WIND:TRAC4:FEED "CH1_S22_1"');
fprintf(idn_PNA, 'CALC:PAR:SEL "CH1_S22_1"');
fprintf(idn_PNA, 'CALC:FORM SMIT');
fprintf(idn_PNA, 'CALC:MARK2 ON');
fprintf(idn_PNA, 'CALC:MARK2:X 2.1e9');
fprintf(idn_PNA, 'DISP:WIND:ANN:MARK:SIZE LARG');

% GET S11
ConfigurePNA(idn_PNA,X);
s11Data = PNA_FetchData(idn_PNA,'S11');
frequencies=s11Data(1,:);
s11 = complex(s11Data(2,:),s11Data(3,:));
clrdevice(idn_PNA);
% GET S12
ConfigurePNA(idn_PNA,X);
s12Data = PNA_FetchData(idn_PNA,'S12');
s12 = complex(s11Data(2,:),s11Data(3,:));
clrdevice(idn_PNA);
% GET S21
ConfigurePNA(idn_PNA,X);
s21Data = PNA_FetchData(idn_PNA,'S21');
s21 = complex(s11Data(2,:),s11Data(3,:));
clrdevice(idn_PNA);
% GET S22
ConfigurePNA(idn_PNA,X);
s22Data = PNA_FetchData(idn_PNA,'S22');
s22 = complex(s11Data(2,:),s11Data(3,:));
clrdevice(idn_PNA);

%% Close 

% fclose(idn_PNA);
% clear idn_PNA