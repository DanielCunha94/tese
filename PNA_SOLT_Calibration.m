%% PNA SOLT Calibration
%Daniel

%% Workspace Cleanup

clear all; close all force; clc
close_all_open_instruments;

%% Connect to interface object, idn_PNA
[~,idn_PNA] = PNA_connect('192.168.81.201',5);

X.fileFormat = 'RI';
X.byteOrder = 'SWAP';
X.dataTypePNA = 'REAL,64';
WaitForSystemReady(idn_PNA);

numOfPoints = query(idn_PNA, 'SENS:SWE:POIN?','%s\n','%d');

% Configure parameter to be acquired and initiate
ConfigurePNA(idn_PNA,X);

%  S Param in RI format.
PNA_CreateMeas(idn_PNA);

%% SOL PORT1

S11S = -1;
S11O = 1;
S11L = 0;

input('PORT1 OPEN - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S11OData = PNA_FetchData(idn_PNA,'S11');
Frequencies=S11Data(1,:);
S11Om = complex(S11OData(2,:),S11OData(3,:));
clrdevice(idn_PNA);

input('PORT1 SHORT - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S11SData = PNA_FetchData(idn_PNA,'S11');
S11Sm = complex(S11SData(2,:),S11SData(3,:));
clrdevice(idn_PNA);

input('PORT1 LOAD - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S11LData = PNA_FetchData(idn_PNA,'S11');
S11Lm = complex(S11LData(2,:),S11LData(3,:));
clrdevice(idn_PNA);

syms e00c e10c e11c
eqn1 = S11Sm == e00c + (e10c .* S11S)./(1- e11c.*S11S);
eqn2 = S11Om == e00c + (e10c .* S11O)./(1- e11c.*S11O);
eqn3 = S11Lm == e00c + (e10c .* S11L)./(1- e11c.*S11L);

%error model coefficients
sol = solve([eqn1, eqn2, eqn3],[e00c, e10c, e11c]);
e00 = double(sol.e00c);
e10e01 = double(sol.e10c);
e11 = double(sol.e11c);

%% SOL PORT 2

S22S = -1;
S22O = 1;
S22L = 0;

input('PORT2 OPEN - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S22OData = PNA_FetchData(idn_PNA,'S22');
S22Om = complex(S22OData(2,:),S22OData(3,:));
clrdevice(idn_PNA);

input('PORT2 SHORT - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S22SData = PNA_FetchData(idn_PNA,'S22');
S22Sm = complex(S22SData(2,:),S22SData(3,:));
clrdevice(idn_PNA);

input('PORT2 LOAD - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S11LData = PNA_FetchData(idn_PNA,'S22');
S11Lm = complex(S11LData(2,:),S11LData(3,:));
clrdevice(idn_PNA);

syms e33c e23e32c e22c
eqn1 = S22Sm == e33c + (e23e32c .* S22S)./(1- e22c.*S22S);
eqn2 = S22Om == e33c + (e23e32c .* S22O)./(1- e22c.*S22O);
eqn3 = S22Lm == e33c + (e23e32c .* S22L)./(1- e22c.*S22L);

%error model coefficients
sol = solve([eqn1, eqn2, eqn3],[e33c, e23e32c, e22c]);
e33r = double( sol.e33c);
e23e32 = double( sol.e23e32c);
e22r = double(sol.e22c);

%% thru

input('Thru - Press Enter to continue');
ConfigurePNA(idn_PNA,X);
S11TData = PNA_FetchData(idn_PNA,'S11');
S11Tm = complex(S11TData(2,:),S11TData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S12TData = PNA_FetchData(idn_PNA,'S12');
S12Tm = complex(S12TData(2,:),S12TData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S21TData = PNA_FetchData(idn_PNA,'S21');
S21Tm = complex(S21TData(2,:),S21TData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S22TData = PNA_FetchData(idn_PNA,'S22');
S22Tm = complex(S22TData(2,:),S22TData(3,:));
clrdevice(idn_PNA);

%error model coefficients
e22 = (S11Tm - e00)./(S11Tm.*e11-(e00.*e11-e10e01));
e10e32 = S21Tm.*(1-e11.*e22);
e11r =  (S22Tm - e33r)./(S22Tm.*e22r-(e33r.*e22r-e23e32));
e23e01 = S12Tm.*(1-e11r.*e22r);

%% Get Calibrated Meas

input('Get calibrated Meas - Press Enter to continue');

%Uncal Meas
ConfigurePNA(idn_PNA,X);
S11UData = PNA_FetchData(idn_PNA,'S11');
S11Um = complex(S11UData(2,:),S11UData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S12UData = PNA_FetchData(idn_PNA,'S12');
S12Um = complex(S12UData(2,:),S12UData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S21UData = PNA_FetchData(idn_PNA,'S21');
S21Um = complex(S21UData(2,:),S21UData(3,:));
clrdevice(idn_PNA);

ConfigurePNA(idn_PNA,X);
S22UData = PNA_FetchData(idn_PNA,'S22');
S22Um = complex(S22UData(2,:),S22UData(3,:));
clrdevice(idn_PNA);

%Cal Meas
D =(1+((S11Um-e00)./(e10e01)).*e11).*(1+((S22Um-e33r)./(e23e32)).*e22r)-(S21Um./(e10e32)).*(S12Um./e23e01).*e22.*e11r;
S11C = (((S11Um-e00)./e10e01).*(1+((S22Um-e33r)./e23e32).*e22r)-e22.*(S21Um./e10e32).*(S12Um./e23e01))./D;
S21C = ((S21Um./e10e32).*(1+((S22Um-e33r)./e23e32).*(e22r-e22)))./D;
S12C = ((S12Um/e23e01).*(1+((S11Um-e00)./e10e01).*(e11-e11r)))./D;
S22C = (((S22Um-e33r)./e23e32).*(1+((S11Um-e00)./e10e01).*e11)-e11r.*(S21Um/e10e32).*(S12Um/e23e01))./D;


