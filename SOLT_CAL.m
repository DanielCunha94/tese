%% LPNR
clear
close all
clc
%% Set some paths (change the ADS version for yours)
if ~sum(getenv('HPEESOF_DIR'))
    
    setenv('HPEESOF_DIR','C:\Program Files\Keysight\ADS2020');
    setenv('COMPL_DIR',getenv('HPEESOF_DIR'));
    setenv('SIMARCH','win32_64');
    
    setenv('PATH',[getenv('HPEESOF_DIR') '\bin\' ';' getenv('HPEESOF_DIR') ...
        '\lib\' getenv('SIMARCH') ';' getenv('HPEESOF_DIR') '\adsptolemy\lib.' ...
        getenv('SIMARCH') ';' getenv('HPEESOF_DIR') ...
        '\SystemVue\2014.10\win32_64\sveclient' ';' getenv('PATH')])
end
%% noise 

SNRdb = 500;
SNR = 10^(SNRdb/10);

syms Ps Pn
eqn1 = 1 == Ps+Pn;
eqn2 = SNR == Ps/Pn;

%error model coefficients
sol = solve([eqn1, eqn2],[Ps, Pn]);
Ps = double( sol.Ps);
Pn =double( sol.Pn);

errorMag =0.002;
errorPhase = 10;

% S11Ss = (1-errorMag*rand(1))*exp(1i*((180+errorPhase*rand(1)*randi([-1,1]))*pi)/180)
% S11Os = (1-errorMag*rand(1))*exp(1i*((360+errorPhase*rand(1)*randi([-1,1]))*pi)/180)
% S11Ls = (0+errorMag*rand(1)*randi([-1,1]))*exp(1i*((0+errorPhase*rand(1)*randi([-1,1]))*pi)/180)

S11Ss = -0.998-i*0.028;
S11Os = 0.9985;
S11Ls = -0.0004;

%% R write  Generic MDIF file
index = string([1:1:3]);
freq =[1:0.5:10];
file_path = '..\SOL_file.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'Sr'} {'Sc'}]; % names of the parameters;
Values = [real(S11Ss) imag(S11Ss); real(S11Os) imag(S11Os); real(S11Ls) imag(S11Ls)]; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% wirte the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);
for i=1:3
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);
%% Freq write a Generic MDIF file
index2 = string([1:1:5]);
file_path = '..\Freq_file.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'Freq'}]; % names of the parameters;
Values = [1.8 1.9 2.0 2.1 2.2 ]'; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% wirte the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);
for i=1:5
    fprintf(fid,index2(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);
%% %% SOL PORT1
% First you have to save the netlist of your schematic:
% MenusADS > Simulate > GenerateNetlist > Save the netlist with the name of
% the schematic but with extension .log "for example: HB_ActiveLP.log"
error_values = zeros(5,6);
e00 = zeros(5,1);
e10e01 = zeros(5,1);
e11 = zeros(5,1);

S11S = -1;
S11O = 1;
S11L = 0;

netPath = [fileparts(pwd) '\SOL1.log']; % path where you saved the netlist
eval(['!hpeesofsim -q "' netPath '"'])

load ../MatlabScripts/Port1SOLsimulation.mat
aux = who('SOL1*');

eval(['RawData = ' aux{1} ';'])
clc
for i=1:5
S11Sm = RawData.dataBlocks.data(1+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
S11Om = RawData.dataBlocks.data(2+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
S11Lm = RawData.dataBlocks.data(3+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));

syms e00c e10c e11c
eqn1 = S11Sm == e00c + (e10c * S11S)/(1- e11c*S11S);
eqn2 = S11Om == e00c + (e10c * S11O)/(1- e11c*S11O);
eqn3 = S11Lm == e00c + (e10c * S11L)/(1- e11c*S11L);

%error model coefficients
sol = solve([eqn1, eqn2, eqn3],[e00c, e10c, e11c]);
e00(i) = double(sol.e00c);
e10e01(i) = double(sol.e10c);
e11(i) = double(sol.e11c);

error_values(i,1)= real(e00(i));
error_values(i,2)= imag(e00(i)); 
error_values(i,3)= real(e10e01(i));
error_values(i,4)= imag(e10e01(i));
error_values(i,5)= real(e11(i));
error_values(i,6)= imag(e11(i));
end


%%
index=string(1:1:5);
file_path = '..\Port1Error_terms.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'e00R'} {'e00C'} {'e10e01R'} {'e10e01C'} {'e11R'} {'e11C'}]; % names of the parameters;

Values = error_values; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% wirte the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);
for i=1:5
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);

%% %% SOL PORT2
% First you have to save the netlist of your schematic:
% MenusADS > Simulate > GenerateNetlist > Save the netlist with the name of
% the schematic but with extension .log "for example: HB_ActiveLP.log"
error_values2 = zeros(5,6);
e33r = zeros(5,1);
e23e32 = zeros(5,1);
e22r = zeros(5,1);

S22S = S11S;
S22O = S11O;
S22L = S11L;

netPath = [fileparts(pwd) '\SOL2.log']; % path where you saved the netlist
eval(['!hpeesofsim -q "' netPath '"'])

load ../MatlabScripts/Port2SOLsimulation.mat
aux = who('SOL2*');

eval(['RawData2 = ' aux{1} ';'])
clc
for i=1:5
S22Sm = RawData2.dataBlocks.data(1+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
S22Om = RawData2.dataBlocks.data(2+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
S22Lm = RawData2.dataBlocks.data(3+3*(i-1)).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));

syms e33c e23e32c e22c
eqn1 = S22Sm == e33c + (e23e32c * S22S)/(1- e22c*S22S);
eqn2 = S22Om == e33c + (e23e32c * S22O)/(1- e22c*S22O);
eqn3 = S22Lm == e33c + (e23e32c * S22L)/(1- e22c*S22L);

%error model coefficients
sol = solve([eqn1, eqn2, eqn3],[e33c, e23e32c, e22c]);
e33r(i) = double( sol.e33c);
e23e32(i) =double( sol.e23e32c);
e22r(i) = double(sol.e22c);

error_values2(i,1)= real(e33r(i));
error_values2(i,2)= imag(e33r(i)); 
error_values2(i,3)= real(e23e32(i));
error_values2(i,4)= imag(e23e32(i));
error_values2(i,5)= real(e22r(i));
error_values2(i,6)= imag(e22r(i));
end
%%
index=string(1:1:5);
file_path = '..\Port2Error_terms.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'e33rR'} {'e33rC'} {'e23e32R'} {'e23e32C'} {'e22rR'} {'e22rC'}]; % names of the parameters;

Values = error_values2; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% wirte the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);
for i=1:5
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);

%% thru

index = string([1:1:2]);
a1M = [0.01,0];
a2M = [0,0.01];
file_path = '..\Thru.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'a1M(real)'} {'a2M(real)'}]; % names of the parameters;

Values = [ a1M; a2M]'; % values for each parameter

strnames = [strrep(strjoin(Names),' ','\t') '\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% write the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);
for i=1:2
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);
%% thru
error_values3 = zeros(5,8);
e22 = zeros(5,1);
e10e32= zeros(5,1);
e11r =  zeros(5,1);
e23e01 = zeros(5,1);

netPath = [fileparts(pwd) '\THRU.log']; % path where you saved the netlist
eval(['!hpeesofsim -q "' netPath '"'])

load ../MatlabScripts/ThruSimulation.mat
aux = who('THRU*');

eval(['RawData3 = ' aux{1} ';'])
clc
for i=1:5
a1u = RawData3.dataBlocks(1).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
b1u = RawData3.dataBlocks(2).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
b2u = RawData3.dataBlocks(3).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));

S11Tm(i) = b1u/a1u;
S21Tm(i) = b2u/a1u;

e22(i) = (S11Tm(i) - e00(i))/(S11Tm(i)*e11(i)-(e00(i)*e11(i)-e10e01(i)));
e10e32(i)= S21Tm(i)*(1-e11(i)*e22(i));

b1u = RawData3.dataBlocks(2).data(2*i).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
a2u = RawData3.dataBlocks(4).data(2*i).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
b2u = RawData3.dataBlocks(3).data(2*i).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));

S22Tm(i) = b2u/a2u;
S12Tm(i) = b1u/a2u;

e11r(i) =  (S22Tm(i) - e33r(i))/(S22Tm(i)*e22r(i)-(e33r(i)*e22r(i)-e23e32(i)));
e23e01(i) = S12Tm(i)*(1-e11r(i)*e22r(i));

error_values3(i,1)= real(e22(i));
error_values3(i,2)= imag(e22(i)); 
error_values3(i,3)= real(e10e32(i));
error_values3(i,4)= imag(e10e32(i));
error_values3(i,5)= real(e11r(i));
error_values3(i,6)= imag(e11r(i));
error_values3(i,7)= real(e23e01(i));
error_values3(i,8)= imag(e23e01(i));
end

%%
index=string(1:1:5);
file_path = '..\SOLTError_terms.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'e00R'} {'e00C'} {'e10e01R'} {'e10e01C'} {'e11R'} {'e11C'} {'e33rR'} {'e33rC'} {'e23e32R'} {'e23e32C'} {'e22rR'} {'e22rC'} {'e22R'} {'e22C'} {'e10e32R'} {'e10e32C'} {'e11rR'} {'e11rC'} {'e23e01R'} {'e23e01C'}]; % names of the parameters;
Values = [error_values error_values2 error_values3] ; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write

% wirte the content
fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);

for i=1:5
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);

for i=1:5
D(i) =(1+((S11Tm(i)-e00(i))/(e10e01(i)))*e11(i))*(1+((S22Tm(i)-e33r(i))/(e23e32(i)))*e22r(i))-(S21Tm(i)/(e10e32(i)))*(S12Tm(i)/e23e01(i))*e22(i)*e11r(i);
S11C(i) = (((S11Tm(i)-e00(i))/e10e01(i))*(1+((S22Tm(i)-e33r(i))/e23e32(i))*e22r(i))-e22(i)*(S21Tm(i)/e10e32(i))*(S12Tm(i)/e23e01(i)))/D(i);
S21C(i) = ((S21Tm(i) /e10e32(i) )*(1+((S22Tm(i) -e33r(i) )/e23e32(i) )*(e22r(i) -e22(i) )))/D(i) ;
S12C(i) = ((S12Tm(i) /e23e01(i) )*(1+((S11Tm(i) -e00(i) )/e10e01(i) )*(e11(i) -e11r(i) )))/D(i) ;
S22C(i) = (((S22Tm(i) -e33r(i) )/e23e32(i) )*(1+((S11Tm(i) -e00(i) )/e10e01(i) )*e11(i) )-e11r(i) *(S21Tm(i) /e10e32(i) )*(S12Tm(i) /e23e01(i) ))/D(i) ;
end

%% POWERCAL

netPath = [fileparts(pwd) '\PowerCal.log']; % path where you saved the netlist
eval(['!hpeesofsim -q "' netPath '"'])

load ../MatlabScripts/PowerCal.mat
aux = who('PowerCal*');

eval(['RawData4 = ' aux{1} ';'])
clc
for i=1:5
a1u = RawData4.dataBlocks(4).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
b1u = RawData4.dataBlocks(3).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
a2u = RawData4.dataBlocks(1).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));
b2u = RawData4.dataBlocks(2).data(2*i-1).dependents(2)*Ps+Pn*(rand(1)+i*rand(1));

S11un = b1u/a1u;
S22un = b2u/a2u;

S11C(i)= (S11un-e00(i))/(e10e01(i)+e11(i)*(S11un-e00(i)));
S22C(i)= (S22un-e33r(i))/(e23e32(i)+e22r(i)*(S22un-e33r(i)));

Pin(i)= RawData4.dataBlocks(5).data(2*i-1).dependents(2);

a1 = sqrt((2*Pin(i))/(1-abs(S11C(i))^2));
b1 = a1*S11C(i);
e10(i) = (a1-e11(i)*b1)/a1u;
e01(i)= e10e01(i)/e10(i);

a1cal(i) = ((e10e01(i) - e00(i)*e11(i))/(e01(i)))*a1u+(e11(i)/e01(i))*b1u;
b1cal(i) = (-e00(i)/e01(i))*a1u + (1/e01(i))*b1u;

e32f(i) = e10e32(i)/e10(i);
gammaSw2(i) =  (e22(i)-e22r(i))/(e23e32(i)+e33r(i)*(e22(i)-e22r(i)));
e23(i) = e23e32(i)/ (e32f(i)*(1-e33r(i)*gammaSw2(i)));
e32(i)= e23e32(i)/e23(i);
e01r(i) = e23e01(i)/e23(i); 


% a2 =b1;
% b2 = a1;
% e23(i) = (a2-e22r(i)*b2)/a2u;
% e32(i)= e23e32(i)/e23(i);

a2cal(i) = ((e23e32(i) - e33r(i)*e22r(i))/(e32(i)))*a2u+(e22r(i)/e32(i))*b2u;
b2cal(i) = (-e33r(i)/e32(i))*a2u + (1/e32(i))*b2u; 

error_values4(i,1)= real(e10(i));
error_values4(i,2)= imag(e10(i)); 
error_values4(i,3)= real(e01(i));
error_values4(i,4)= imag(e01(i));
error_values4(i,5)= real(e23(i));
error_values4(i,6)= imag(e23(i));
error_values4(i,7)= real(e32(i));
error_values4(i,8)= imag(e32(i));
end

S11P= b1cal./a1cal;
S22Pi= a2cal./b2cal;

%%
file_path = '..\POWER_CAL_terms.mdf'; % name of the file (including the path)
fid       = fopen(file_path,'w'); % open the file to Write

Names    = [{'e10R'} {'e10C'} {'e01R'} {'e01C'} {'e23R'} {'e23C'} {'e32R'} {'e32C'}]; % names of the parameters;

Values = error_values4; % values for each parameter

strnames = [strrep(strjoin(Names),' ','(real)\t') '(real)\n']; %convert the names to a single string
Nparam   = numel(Names);                    % Number of parameters to write


fprintf(fid,'BEGIN block\n');
fprintf(fid,['%% INDEX(int)\t' strnames]);

for i=1:5
    fprintf(fid,index(i));
    fprintf(fid,['\t' repmat('%f\t',1,Nparam) '\n'],Values(i,:));
end
fprintf(fid,'END\n');
fclose(fid);




