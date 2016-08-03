% Automated script to plot ATCA data
% William Peters June 22, 2005
clear
channelstr = { 'T1','T12','T20','M1','M20','B1','B12','B20' };
setstr = { 'Amador','Synth','Calaveras','Alt_Synth','Amaveras','Tuolmador','Synthest','peters010605'};
channels = [1:8]; % Select channels from channelstr
sets = [6,8]; % Select measurement set from setstr
% specify output plots
% Key     21 PR 11 22 N  F  MX ICR
plots = [ 0  0  0  0  0  1  0  0 ];
channelplots = 1; setplots = 1; plotlegend = 1;
file_prefix='20050623A'; % output plot file prefix
xlower = [   0 -0.5    0    0    0    0    0  0.5 ];
xupper = [  10    3   10   10   10   10   10    5 ];
ylower = [ -60 -Inf  -20  -20  -50  -50  -50  -10 ];
yupper = [   0  Inf    0    0  -20  -20  -20  -10 ];
fontsize = 16;
global ui; maxfreq=200; ui=0.1; % pulse response parameters

% Channel data specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amador
thru(1).dir = '\\kopernik\data1\Electrical Design Team\Lab\NEW_CTO_BACKPLANE\BackplaneMeasurements\NormalBoard2\10MHz_Steps\';
thru(1).str = { 'B2_A2B2_C2D2_L2_1p25inch_10MHz.s4p';
    'B2_A2B2_C2D2_L2_12inch_10MHz.s4p';
    'B2_A2B2_C2D2_L2_20inch_10MHz.s4p';
    'B2_G6H6_E6F6_L11_1p25inch_10MHz.s4p';
    'B2_G6H6_E6F6_L11_20inch_10MHz.s4p';
    'B2_A6B6_C6D6_L17_1p25inch_10MHz.s4p';
    'B2_A6B6_C6D6_L17_12inch_10MHz.s4p';
    'B2_A6B6_C6D6_L17_20inch_10MHz.s4p'};
next(1).dir = '\\kopernik\data1\Electrical Design Team\Lab\NEW_CTO_BACKPLANE\BackplaneMeasurements\NormalBoard2\New_NEXT\';
next(1).str = { 'B2_RC5_A1B1_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_A1B1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p';
    'B2_RC5_A1B1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p';
    'B2_RC5_G5H5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_G5H5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p'};
fext(1).dir = '\\kopernik\data1\Electrical Design Team\Lab\NEW_CTO_BACKPLANE\BackplaneMeasurements\NormalBoard2\New_FEXT\';
fext(1).str = { 'B2_A1B1_C2D2_L2_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_A1B1_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p';
    'B2_A1B1_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p';
    'B2_G5H5_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_G5H5_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p'};
% Improved 3/11
thru(2).dir = '\\kopernik\data2\Bill\ATCA_HVM\T1_replacement\replacementGen\';
thru(2).str = { 'HVM_ATCA_5_1p25_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_12_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_1p25_5_7_11_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_11_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_1p25_5_7_17_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_12_5_7_17_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_17_7_ty_ty_ty_AB_CD.s4p'};
next(2).dir = '\\kopernik\data1\Electrical Design Team\Lab\NEW_CTO_BACKPLANE\BackplaneMeasurements\NormalBoard2\New_NEXT\';
next(2).str = { 'B2_RC5_A1B1_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_A1B1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p';
    'B2_RC5_A1B1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A2B2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A3B3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E1F1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E2F2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E3F3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p';
    'B2_RC5_G5H5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_G5H5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p';
    'B2_RC5_A5B5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p'; };
fext(2).dir = '\\kopernik\data1\Electrical Design Team\Lab\NEW_CTO_BACKPLANE\BackplaneMeasurements\NormalBoard2\New_FEXT\';
fext(2).str = { 'B2_A1B1_C2D2_L2_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_A1B1_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p';
    'B2_A1B1_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p','B2_A3B3_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p';
    'B2_G5H5_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_G5H5_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p';
    'B2_A5B5_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p'; };
% Calaveras Alt1
thru(3).dir = '\\kopernik\data1\Electrical Design Team\Lab\Calaveras\measurements\';
thru(3).str = { 'C1_A2B2_C2D2_L2_THRU_1inch_10MHz_mag_ang.s4p';
    'C1_A2B2_C2D2_L2_THRU_12inch_10MHz_mag_ang.s4p';
    'C1_A2B2_C2D2_L2_THRU_20inch_10MHz_mag_ang.s4p';
    'C1_E6F6_G6H6_L11_THRU_1inch_10MHzs4p.s4p';
    'C1_E6F6_G6H6_L11_THRU_20inch_10MHz_mag_ang.s4p';
    'C1_A6B6_C6D6_L17_THRU_1inch_10MHz._mag_ang.s4p';
    'C1_A6B6_C6D6_L17_THRU_12inch_10MHz._mag_ang.s4p';
    'C1_A6B6_C6D6_L17_THRU_20inch_10MHz._mag_ang.s4p'};
next(3).dir = '\\kopernik\data1\Electrical Design Team\Lab\Calaveras\measurements\';
next(3).str = { ...
    'C1_A1B1_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_E1F1_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_E2F2_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_E3F3_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_E1F1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_E2F2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_E3F3_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_E1F1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_E2F2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_E3F3_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_G5H5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_G5H5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p'; ...
    };
fext(3).dir = '\\kopernik\data1\Electrical Design Team\Lab\Calaveras\measurements\';
fext(3).str = { ...
    'C1_A1B1_C2D2_L2_FEXT_1inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_1inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_G5H5_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_G5H5_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    };
% Improved 3/16
thru(4).dir = '\\kopernik\data2\Bill\ATCA_HVM\T1_replacement\replacementGen\';
thru(4).str = { 'HVM_ATCA_031605_5_1p25_5_2_2_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_12_5_2_2_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_20_5_2_2_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_1p25_5_2_11_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_20_5_2_11_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_1p25_5_2_17_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_12_5_2_17_2_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_031605_5_20_5_2_17_2_ty_ty_ty_AB_CD.s4p'};
% Amaveras (Amador test card with Calaveras backplane)
thru(5).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Amador_Line_Card\';
thru(5).str = { 'Ed_C1_A2B2_C2D2_L2_THRU_1inch_10MHz_mag_ang.s4p';
    'Ed_C1_A2B2_C2D2_L2_THRU_12inch_10MHz_mag_ang.s4p';
    'Ed_C1_A2B2_C2D2_L2_THRU_20inch_10MHz_mag_ang.s4p';
    'C1_E6F6_G6H6_L11_THRU_1inch_10MHzs4p.s4p';
    'C1_E6F6_G6H6_L11_THRU_20inch_10MHz_mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_1inch_10MHz._mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_12inch_10MHz._mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_20inch_10MHz._mag_ang.s4p'};
next(5).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Amador_Line_Card\';
next(5).str = { ...
    'C1_A1B1_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_1inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_1inch_10MHz_real_imag.s4p','C1_E1F1_C2D2_L2_NEXT_1inch_10MHz_real_imag.s4p','C1_E2F2_C2D2_L2_NEXT_1inch_10MHz_real_imag.s4p','C1_E3F3_C2D2_L2_NEXT_1inch_10MHz_real_imag.s4p'; ...
    'C1_A1B1_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_12inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_12inch_10MHz_real_imag.s4p','C1_E1F1_C2D2_L2_NEXT_12inch_10MHz_real_imag.s4p','C1_E2F2_C2D2_L2_NEXT_12inch_10MHz_real_imag.s4p','C1_E3F3_C2D2_L2_NEXT_12inch_10MHz_real_imag.s4p'; ...
    'C1_A1B1_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_A2B2_C2D2_L2_NEXT_20inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_NEXT_20inch_10MHz_real_imag.s4p','C1_E1F1_C2D2_L2_NEXT_20inch_10MHz_real_imag.s4p','C1_E2F2_C2D2_L2_NEXT_20inch_10MHz_real_imag.s4p','C1_E3F3_C2D2_L2_NEXT_20inch_10MHz_real_imag.s4p'; ...
    'B2_RC5_G5H5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_G5H5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G6H6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_G7H7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C5D5_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C6D6_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_C7D7_E6F6_L11_NEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_12inch_10MHz_mag_ang.s4p'; ...
    'B2_RC5_A5B5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A6B6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_A7B7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E5F5_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E6F6_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p','B2_RC5_E7F7_C6D6_L17_NEXT_20inch_10MHz_mag_ang.s4p'; ...
    };
fext(5).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Amador_Line_Card\';
fext(5).str = { ...
    'C1_A1B1_C2D2_L2_FEXT_1inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_1inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_12inch_10MHz_mag_ang.s4p'; ...
    'C1_A1B1_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p','C1_A3B3_C2D2_L2_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_G5H5_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_G5H5_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p','B2_G7H7_E6F6_L11_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_1p25inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_12inch_10MHz_mag_ang.s4p'; ...
    'B2_A5B5_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p','B2_A7B7_C6D6_L17_FEXT_20inch_10MHz_mag_ang.s4p'; ...
    };
% Tuolmador test card with Calaveras backplane
thru(6).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Tuolmador_LC\';
thru(6).str = { ...
    'Ed_C1_A2B2_C2D2_L2_THRU_1INCH_10MHz_MAG_ANG.s4p';
    'Ed_C1_A2B2_C2D2_L2_THRU_12INCH_10MHz_MAG_ANG.s4p';
    'Ed_C1_A2B2_C2D2_L2_THRU_20INCH_10MHz_MAG_ANG.s4p';
    'Ed_C1_E6F6_G6H6_L11_THRU_1INCH_10Mhz_MAG_ANG.s4p';
    'Ed_C1_E6F6_G6H6_L11_THRU_20inch_10MHz_mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_1inch_10MHz_mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_12inch_10MHz_mag_ang.s4p';
    'Ed_C1_A6B6_C6D6_L17_THRU_20inch_10MHz_mag_ang.s4p'};
next(6).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Tuolmador_LC\';
next(6).str = { ...
    'Ed_C1_A1B1_C2D2_L2_NEXT1_1inch_10MHz_mag_ang.s4p','Ed_C1_A2B2_C2D2_L2_NEXT2_1inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_NEXT3_1inch_10MHz_mag_ang.s4p','Ed_C1_E1F1_C2D2_L2_NEXT4_1inch_10MHz_mag_ang.s4p','Ed_C1_E2F2_C2D2_L2_NEXT5_1inch_10MHz_mag_ang.s4p','Ed_C1_E3F3_C2D2_L2_NEXT6_1INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_A1B1_C2D2_L2_NEXT1_12inch_10MHz_mag_ang.s4p','Ed_C1_A2B2_C2D2_L2_NEXT2_12inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_NEXT3_12inch_10MHz_mag_ang.s4p','Ed_C1_E1F1_C2D2_L2_NEXT4_12inch_10MHz_mag_ang.s4p','Ed_C1_E2F2_C2D2_L2_NEXT5_12inch_10MHz_mag_ang.s4p','Ed_C1_E3F3_C2D2_L2_NEXT6_12INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_A1B1_C2D2_L2_NEXT1_20inch_10MHz_mag_ang.s4p','Ed_C1_A2B2_C2D2_L2_NEXT2_20inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_NEXT3_20inch_10MHz_mag_ang.s4p','Ed_C1_E1F1_C2D2_L2_NEXT4_20inch_10MHz_mag_ang.s4p','Ed_C1_E2F2_C2D2_L2_NEXT5_20inch_10MHz_mag_ang.s4p','Ed_C1_E3F3_C2D2_L2_NEXT6_20INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_E5F5_G6H6_L11_NEXT1_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E6F6_G6H6_L11_NEXT2_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E7F7_G6H6_L11_NEXT3_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C5D5_G6H6_L11_NEXT4_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C6D6_G6H6_L11_NEXT5_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C7D7_G6H6_L11_NEXT6_1INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_E5F5_G6H6_L11_NEXT1_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E6F6_G6H6_L11_NEXT2_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E7F7_G6H6_L11_NEXT3_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C5D5_G6H6_L11_NEXT4_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C6D6_G6H6_L11_NEXT5_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_C7D7_G6H6_L11_NEXT6_20INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_NEXT1_1inch_10MHz_mag_ang.s4p','Ed_C1_A6B6_C6D6_L17_NEXT2_1inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_NEXT3_1inch_10MHz_mag_ang.s4p','Ed_C1_E5F5_C6D6_L17_NEXT4_1inch_10MHz_mag_ang.s4p','Ed_C1_E6F6_C6D6_L17_NEXT5_1inch_10MHz_mag_ang.s4p','Ed_C1_E7F7_C6D6_L17_NEXT6_1inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_NEXT1_12inch_10MHz_mag_ang.s4p','Ed_C1_A6B6_C6D6_L17_NEXT2_12inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_NEXT3_12inch_10MHz_mag_ang.s4p','Ed_C1_E5F5_C6D6_L17_NEXT4_12inch_10MHz_mag_ang.s4p','Ed_C1_E6F6_C6D6_L17_NEXT5_12inch_10MHz_mag_ang.s4p','Ed_C1_E7F7_C6D6_L17_NEXT6_12inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_NEXT1_20inch_10MHz_mag_ang.s4p','Ed_C1_A6B6_C6D6_L17_NEXT2_20inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_NEXT3_20inch_10MHz_mag_ang.s4p','Ed_C1_E5F5_C6D6_L17_NEXT4_20inch_10MHz_mag_ang.s4p','Ed_C1_E6F6_C6D6_L17_NEXT5_20inch_10MHz_mag_ang.s4p','Ed_C1_E7F7_C6D6_L17_NEXT6_20inch_10MHz_mag_ang.s4p'; ...
    };
fext(6).dir = '\\Kopernik\data1\Electrical Design Team\Lab\Measurement\Calaveras_with_Tuolmador_LC\';
fext(6).str = { ...
    'Ed_C1_A1B1_C2D2_L2_FEXT1_1inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_FEXT2_1inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A1B1_C2D2_L2_FEXT1_12inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_FEXT2_12inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A1B1_C2D2_L2_FEXT1_20inch_10MHz_mag_ang.s4p','Ed_C1_A3B3_C2D2_L2_FEXT2_20inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_E5F5_G6H6_L11_FEXT1_1INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E7F7_G6H6_L11_FEXT2_1INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_E5F5_G6H6_L11_FEXT1_20INCH_10Mhz_MAG_ANG.s4p','Ed_C1_E7F7_G6H6_L11_FEXT2_20INCH_10Mhz_MAG_ANG.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_FEXT1_1inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_FEXT2_1inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_FEXT1_12inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_FEXT2_12inch_10MHz_mag_ang.s4p'; ...
    'Ed_C1_A5B5_C6D6_L17_FEXT1_20inch_10MHz_mag_ang.s4p','Ed_C1_A7B7_C6D6_L17_FEXT2_20inch_10MHz_mag_ang.s4p'; ...
    };
% Improved 3/11 with crosstalk estimates
thru(7).dir = '\\kopernik\data2\Bill\ATCA_HVM\T1_replacement\replacementGen\';
thru(7).str = { 'HVM_ATCA_5_1p25_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_12_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_2_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_1p25_5_7_11_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_11_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_1p25_5_7_17_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_12_5_7_17_7_ty_ty_ty_AB_CD.s4p';
    'HVM_ATCA_5_20_5_7_17_7_ty_ty_ty_AB_CD.s4p'};
next(7).dir = '\\kopernik\data2\Bill\test\';
next(7).str = { ...
    'T1_NEXT1_estimated.s4p','T1_NEXT2_estimated.s4p','T1_NEXT3_estimated.s4p','T1_NEXT4_estimated.s4p','T1_NEXT5_estimated.s4p','T1_NEXT6_estimated.s4p'; ...
    'T12_NEXT1_estimated.s4p','T12_NEXT2_estimated.s4p','T12_NEXT3_estimated.s4p','T12_NEXT4_estimated.s4p','T12_NEXT5_estimated.s4p','T12_NEXT6_estimated.s4p'; ...
    'T20_NEXT1_estimated.s4p','T20_NEXT2_estimated.s4p','T20_NEXT3_estimated.s4p','T20_NEXT4_estimated.s4p','T20_NEXT5_estimated.s4p','T20_NEXT6_estimated.s4p'; ...
    'M1_NEXT1_estimated.s4p','M1_NEXT2_estimated.s4p','M1_NEXT3_estimated.s4p','M1_NEXT4_estimated.s4p','M1_NEXT5_estimated.s4p','M1_NEXT6_estimated.s4p'; ...
    'M20_NEXT1_estimated.s4p','M20_NEXT2_estimated.s4p','M20_NEXT3_estimated.s4p','M20_NEXT4_estimated.s4p','M20_NEXT5_estimated.s4p','M20_NEXT6_estimated.s4p'; ...
    'B1_NEXT1_estimated.s4p','B1_NEXT2_estimated.s4p','B1_NEXT3_estimated.s4p','B1_NEXT4_estimated.s4p','B1_NEXT5_estimated.s4p','B1_NEXT6_estimated.s4p'; ...
    'B12_NEXT1_estimated.s4p','B12_NEXT2_estimated.s4p','B12_NEXT3_estimated.s4p','B12_NEXT4_estimated.s4p','B12_NEXT5_estimated.s4p','B12_NEXT6_estimated.s4p'; ...
    'B20_NEXT1_estimated.s4p','B20_NEXT2_estimated.s4p','B20_NEXT3_estimated.s4p','B20_NEXT4_estimated.s4p','B20_NEXT5_estimated.s4p','B20_NEXT6_estimated.s4p'; ...
    };
fext(7).dir = '\\kopernik\data2\Bill\test\';
fext(7).str = { ...
    'T1_FEXT1_estimated.s4p','T1_FEXT2_estimated.s4p'; ...
    'T12_FEXT1_estimated.s4p','T12_FEXT2_estimated.s4p'; ...
    'T20_FEXT1_estimated.s4p','T20_FEXT2_estimated.s4p'; ...
    'M1_FEXT1_estimated.s4p','M1_FEXT2_estimated.s4p'; ...
    'M20_FEXT1_estimated.s4p','M20_FEXT2_estimated.s4p'; ...
    'B1_FEXT1_estimated.s4p','B1_FEXT2_estimated.s4p'; ...
    'B12_FEXT1_estimated.s4p','B12_FEXT2_estimated.s4p'; ...
    'B20_FEXT1_estimated.s4p','B20_FEXT2_estimated.s4p'; ...
    };
% Tuolmador test card with Calaveras backplane (IEEE version)
thru(8).dir = '\\kopernik\data2\BackplaneEthernet\Measurement\peters_01_0605\';
thru(8).str = { ...
    'peters_01_0605_T1_thru.s4p';
    'peters_01_0605_T12_thru.s4p';
    'peters_01_0605_T20_thru.s4p';
    'peters_01_0605_M1_thru.s4p';   
    'peters_01_0605_M20_thru.s4p';
    'peters_01_0605_B1_thru.s4p';
    'peters_01_0605_B12_thru.s4p';
    'peters_01_0605_B20_thru.s4p'};
next(8).dir = '\\kopernik\data2\BackplaneEthernet\Measurement\peters_01_0605\';
next(8).str = { ...
    'peters_01_0605_T1_next1.s4p','peters_01_0605_T1_next2.s4p','peters_01_0605_T1_next3.s4p','peters_01_0605_T1_next4.s4p','peters_01_0605_T1_next5.s4p','peters_01_0605_T1_next6.s4p'; ...
    'peters_01_0605_T12_next1.s4p','peters_01_0605_T12_next2.s4p','peters_01_0605_T12_next3.s4p','peters_01_0605_T12_next4.s4p','peters_01_0605_T12_next5.s4p','peters_01_0605_T12_next6.s4p'; ...
    'peters_01_0605_T20_next1.s4p','peters_01_0605_T20_next2.s4p','peters_01_0605_T20_next3.s4p','peters_01_0605_T20_next4.s4p','peters_01_0605_T20_next5.s4p','peters_01_0605_T20_next6.s4p'; ...
    'peters_01_0605_M1_next1.s4p','peters_01_0605_M1_next2.s4p','peters_01_0605_M1_next3.s4p','peters_01_0605_M1_next4.s4p','peters_01_0605_M1_next5.s4p','peters_01_0605_M1_next6.s4p'; ...
    'peters_01_0605_M20_next1.s4p','peters_01_0605_M20_next2.s4p','peters_01_0605_M20_next3.s4p','peters_01_0605_M20_next4.s4p','peters_01_0605_M20_next5.s4p','peters_01_0605_M20_next6.s4p'; ...
    'peters_01_0605_B1_next1.s4p','peters_01_0605_B1_next2.s4p','peters_01_0605_B1_next3.s4p','peters_01_0605_B1_next4.s4p','peters_01_0605_B1_next5.s4p','peters_01_0605_B1_next6.s4p'; ...
    'peters_01_0605_B12_next1.s4p','peters_01_0605_B12_next2.s4p','peters_01_0605_B12_next3.s4p','peters_01_0605_B12_next4.s4p','peters_01_0605_B12_next5.s4p','peters_01_0605_B12_next6.s4p'; ...
    'peters_01_0605_B20_next1.s4p','peters_01_0605_B20_next2.s4p','peters_01_0605_B20_next3.s4p','peters_01_0605_B20_next4.s4p','peters_01_0605_B20_next5.s4p','peters_01_0605_B20_next6.s4p'; ...
    };
fext(8).dir = '\\kopernik\data2\BackplaneEthernet\Measurement\peters_01_0605\';
fext(8).str = { ...
    'peters_01_0605_T1_fext1.s4p','peters_01_0605_T1_fext2.s4p'; ...
    'peters_01_0605_T12_fext1.s4p','peters_01_0605_T12_fext2.s4p'; ...
    'peters_01_0605_T20_fext1.s4p','peters_01_0605_T20_fext2.s4p'; ...
    'peters_01_0605_M1_fext1.s4p','peters_01_0605_M1_fext2.s4p'; ...
    'peters_01_0605_M20_fext1.s4p','peters_01_0605_M20_fext2.s4p'; ...
    'peters_01_0605_B1_fext1.s4p','peters_01_0605_B1_fext2.s4p'; ...
    'peters_01_0605_B12_fext1.s4p','peters_01_0605_B12_fext2.s4p'; ...
    'peters_01_0605_B20_fext1.s4p','peters_01_0605_B20_fext2.s4p'; ...
    };

% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(channels)
    for n = 1:length(sets)
        % THRU
        if any(plots([1:4,8]))
            [checker(channels(m),sets(n)).thru.S,checker(channels(m),sets(n)).thru.f] = get_snp([thru(sets(n)).dir,thru(sets(n)).str{channels(m)}]);
            checker(channels(m),sets(n)).thru.SDD = mixedmode(checker(channels(m),sets(n)).thru.S,1);
        end
        % Pulse Response
        if plots(2)
            if checker(channels(m),sets(n)).thru.f(1) == 0 % check for DC point, add if necessary
                checker(channels(m),sets(n)).thru.SDD21td = squeeze(checker(channels(m),sets(n)).thru.SDD(2,1,:));
                checker(channels(m),sets(n)).thru.ftd = checker(channels(m),sets(n)).thru.f;
            elseif checker(channels(m),sets(n)).thru.f(1) < 1
                checker(channels(m),sets(n)).thru.f(1) = 0;
                checker(channels(m),sets(n)).thru.SDD21td = squeeze(checker(channels(m),sets(n)).thru.SDD(2,1,:));
                checker(channels(m),sets(n)).thru.ftd = checker(channels(m),sets(n)).thru.f;
            else
                checker(channels(m),sets(n)).thru.SDD21td = [1;squeeze(checker(channels(m),sets(n)).thru.SDD(2,1,:))];
                checker(channels(m),sets(n)).thru.ftd = [0;checker(channels(m),sets(n)).thru.f];
            end
            [checker(channels(m),sets(n)).thru.PR,checker(channels(m),sets(n)).thru.t,checker(channels(m),sets(n)).thru.tstep] = gen_pulse_response(checker(channels(m),sets(n)).thru.SDD21td,checker(channels(m),sets(n)).thru.ftd,maxfreq,ui);
            checker(channels(m),sets(n)).thru.t0 = checker(channels(m),sets(n)).thru.t - checker(channels(m),sets(n)).thru.t(checker(channels(m),sets(n)).thru.PR==max(checker(channels(m),sets(n)).thru.PR));
        end
        % NEXT
        if any(plots([5,7,8]))
            % NEXT Aggressors
            for p=1:size(next(sets(n)).str,2)
                [checker(channels(m),sets(n)).next(p).S,checker(channels(m),sets(n)).next(p).f] = get_snp([next(sets(n)).dir,next(sets(n)).str{channels(m),p}]);
                checker(channels(m),sets(n)).next(p).SDD = mixedmode(checker(channels(m),sets(n)).next(p).S,1);
            end
            % Calculate power sum of multiple NEXT aggressors
            clear mnext
            for ii = 1:length(checker(channels(m),sets(n)).next)
                mnext(:,ii) = abs( squeeze( checker(channels(m),sets(n)).next(ii).SDD(2,1,:) ) ).^2;
            end
            checker(channels(m),sets(n)).mnext = sqrt(sum(mnext,2));
        end
        % FEXT
        if any(plots([6:8]))
            % FEXT Aggressors
            for p=1:size(fext(sets(n)).str,2)
                [checker(channels(m),sets(n)).fext(p).S,checker(channels(m),sets(n)).fext(p).f] = get_snp([fext(sets(n)).dir,fext(sets(n)).str{channels(m),p}]);
                checker(channels(m),sets(n)).fext(p).SDD = mixedmode(checker(channels(m),sets(n)).fext(p).S,1);
            end
            % Calculate power sum of multiple FEXT aggressors
            clear mfext
            for ii = 1:length(checker(channels(m),sets(n)).fext)
                mfext(:,ii) = abs( squeeze( checker(channels(m),sets(n)).fext(ii).SDD(2,1,:) ) ).^2;
            end
            checker(channels(m),sets(n)).mfext = sqrt(sum(mfext,2));
        end
        % MXT
        if plots(7)
            clear mxt
            for ii = 1:length(checker(channels(m),sets(n)).fext)
                mxt(:,ii) = abs( squeeze( checker(channels(m),sets(n)).fext(ii).SDD(2,1,:) ) ).^2;
            end
            for jj = (ii+1):length(checker(channels(m),sets(n)).next)
                mxt(:,jj) = abs( squeeze( checker(channels(m),sets(n)).next(jj-ii).SDD(2,1,:) ) ).^2;
            end
            checker(channels(m),sets(n)).mxt = sqrt(sum(mxt,2));
        end
        % ICR
        if plots(8)
           icr_il = checker(channels(m),sets(n)).thru.SDD(2,1,checker(channels(m),sets(n)).thru.f >= checker(channels(m),sets(n)).fext(1).f(1));
           checker(channels(m),sets(n)).icr = dB(squeeze(icr_il))-dB(checker(channels(m),sets(n)).mxt);
end,end,end

% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=1; % insertion loss
if plots(p)
    f00 = [5e7:1e7:15e9];
    ILlim = specline(f00,'IL');
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).thru.f*1e-9,dB(squeeze(checker(channels(m),sets(%d)).thru.SDD(2,1,:))),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = sprintf('%sf00*1e-9,dB(ILlim),''k'',',plotstr);
            legendstr = sprintf('%s''limit'',',legendstr);
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthWest'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD21: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('SDD21 (dB)')
            print('-dpng',sprintf('%s_%s_SDD21.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).thru.f*1e-9,dB(squeeze(checker(channels(%d),sets(n)).thru.SDD(2,1,:))),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = sprintf('%sf00*1e-9,dB(ILlim),''k'',',plotstr);
            legendstr = sprintf('%s''limit'',',legendstr);
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthWest'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD21: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('SDD21 (dB)')
            print('-dpng',sprintf('%s_%s_SDD21.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=2; % pulse response
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).thru.t0.*checker(channels(m),sets(%d)).thru.tstep,checker(channels(m),sets(%d)).thru.PR,',plotstr,n,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''NorthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[-5:0.5:20])
            title(['1V 100ps Pulse Response: ',channelstr{channels(m)}])
            xlabel('time (ns)')
            ylabel('pulse response (V)')
            print('-dpng',sprintf('%s_%s_PR.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).thru.t0.*checker(channels(%d),sets(n)).thru.tstep,checker(channels(%d),sets(n)).thru.PR,',plotstr,m,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''NorthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[-5:0.5:20])
            title(['1V 100ps Pulse Response: ',setstr{sets(n)}])
            xlabel('time (ns)')
            ylabel('pulse response (V)')
            print('-dpng',sprintf('%s_%s_PR.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=3; % SDD11
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).thru.f*1e-9,dB(squeeze(checker(channels(m),sets(%d)).thru.SDD(1,1,:))),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD11: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('SDD11 (dB)')
            print('-dpng',sprintf('%s_%s_SDD11.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).thru.f*1e-9,dB(squeeze(checker(channels(%d),sets(n)).thru.SDD(1,1,:))),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD11: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('SDD11 (dB)')
            print('-dpng',sprintf('%s_%s_SDD11.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=4; % SDD22
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).thru.f*1e-9,dB(squeeze(checker(channels(m),sets(%d)).thru.SDD(2,2,:))),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD22: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('SDD22 (dB)')
            print('-dpng',sprintf('%s_%s_SDD22.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).thru.f*1e-9,dB(squeeze(checker(channels(%d),sets(n)).thru.SDD(2,2,:))),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['SDD22: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('SDD22 (dB)')
            print('-dpng',sprintf('%s_%s_SDD22.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=5; % MNEXT
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).next(1).f*1e-9,dB(checker(channels(m),sets(%d)).mnext),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MNEXT: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('MNEXT (dB)')
            print('-dpng',sprintf('%s_%s_MNEXT.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).next(1).f*1e-9,dB(checker(channels(%d),sets(n)).mnext),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MNEXT: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('MNEXT (dB)')
            print('-dpng',sprintf('%s_%s_MNEXT.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=6; % MFEXT
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).fext(1).f*1e-9,dB(checker(channels(m),sets(%d)).mfext),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MFEXT: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('MFEXT (dB)')
            print('-dpng',sprintf('%s_%s_MFEXT.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).fext(1).f*1e-9,dB(checker(channels(%d),sets(n)).mfext),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MFEXT: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('MFEXT (dB)')
            print('-dpng',sprintf('%s_%s_MFEXT.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=7; % MXT
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).next(1).f*1e-9,dB(checker(channels(m),sets(%d)).mxt),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MXT: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('MXT (dB)')
            print('-dpng',sprintf('%s_%s_MXT.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).next(1).f*1e-9,dB(checker(channels(%d),sets(n)).mxt),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MXT: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('MXT (dB)')
            print('-dpng',sprintf('%s_%s_MXT.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end
p=8; % ICR
if plots(p)
    if channelplots
        for m = 1:length(channels)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for n = 1:length(sets)
                plotstr = sprintf('%schecker(channels(m),sets(%d)).next(1).f*1e-9,dB(checker(channels(m),sets(%d)).mxt),',plotstr,n,n);
                legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['ICR: ',channelstr{channels(m)}])
            xlabel('frequency (GHz)')
            ylabel('ICR (dB)')
            print('-dpng',sprintf('%s_%s_ICR.png',file_prefix,channelstr{channels(m)}));
            close(gcf)
    end,end
    if setplots
        for n = 1:length(sets)
            % generate plot commands
            plotstr='plot(';
            legendstr='legend(';
            for m = 1:length(channels)
                plotstr = sprintf('%schecker(channels(%d),sets(n)).next(1).f*1e-9,dB(checker(channels(%d),sets(n)).mxt),',plotstr,m,m);
                legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
            end
            plotstr = [plotstr(1:end-1),');'];
            legendstr = [legendstr,'''Location'',''SouthEast'');'];
            ATCA_plotscript;
            set(gca,'XTick',[0:15])
            title(['MXT: ',setstr{sets(n)}])
            xlabel('frequency (GHz)')
            ylabel('MXT (dB)')
            print('-dpng',sprintf('%s_%s_MXT.png',file_prefix,setstr{sets(n)}));
            close(gcf)
end,end,end

if plots(p)
    for m = 1:length(channels)
        % Plot ICR %%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1000*channels(m)+10*p);
        % set(gcf,'Position',[100,100,900,600],'Color',[1,1,1]);
        set(gcf,'WindowStyle','docked','Color',[1,1,1]);
        % generate plot commands
        plotstr='semilogx(';
        legendstr='legend(';
        for n = 1:length(sets)
            plotstr = sprintf('%schecker(channels(m),sets(%d)).fext(1).f*1e-9,checker(channels(m),sets(%d)).icr,',plotstr,n,n);
            legendstr = sprintf('%ssetstr{sets(%d)},',legendstr,n);
        end
        plotstr = [plotstr(1:end-1),');'];
        legendstr = [legendstr,'''Location'',''SouthWest'');'];
        eval(plotstr);
        hold on
        f00=[5e7:1e7:5e9];
        ICRlim=12.5-20*log10(f00/5e9);
        semilogx(f00*1e-9,ICRlim,'k')
        hold off
        if plot_legend(1), eval(legendstr); end
        axis([.05,5,ploty1{5},ploty2{5}])
        set(get(gca,'Children'),'LineWidth',2)
        set(gca,'XTick',[.05,.1,.5,1,2,3,4,5],'LineWidth',2)
        set(gca,'FontSize',fontsize);
        grid on
        title(['ICR: ',channelstr{channels(m)}])
        xlabel('frequency (GHz)')
        ylabel('ICR (dB)')
        print('-dpng',sprintf('%s_%s_ICR.png',file_prefix,channelstr{channels(m)}));
        close(gcf)
    end
    % Plot ICR %%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(10*p+1);
    % set(gcf,'Position',[100,100,900,600],'Color',[1,1,1]);
    set(gcf,'WindowStyle','docked','Color',[1,1,1]);
    % generate plot commands
    plotstr='semilogx(';
    legendstr='legend(';
    for m = 1:length(channels)
        plotstr = sprintf('%schecker(channels(%d),sets(1)).fext(1).f*1e-9,checker(channels(%d),sets(1)).icr,',plotstr,m,m);
        legendstr = sprintf('%schannelstr{channels(%d)},',legendstr,m);
    end
    plotstr = [plotstr(1:end-1),');'];
    legendstr = [legendstr,'''Location'',''SouthWest'');'];
    eval(plotstr);
    hold on
    f00=[5e7:1e7:5e9];
    ICRlim=12.5-20*log10(f00/5e9);
    semilogx(f00*1e-9,ICRlim,'k')
    hold off
    if plot_legend(2), eval(legendstr); end
    axis([.05,5,ploty1{5},ploty2{5}])
    set(get(gca,'Children'),'LineWidth',2)
    set(gca,'XTick',[.05,.1,.5,1,2,3,4,5],'LineWidth',2)
    set(gca,'FontSize',fontsize);
    grid on
    title(['ICR: ',setstr{sets(1)}])
    xlabel('frequency (GHz)')
    ylabel('ICR (dB)')
    print('-dpng',sprintf('%s_%s_ICR.png',file_prefix,setstr{sets(1)}));
    close(gcf)
end

