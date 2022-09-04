%% Get MTR and MTsat percentage difference varying F and T1
% Simulate complete qMT-SPGR protocol
Model_qmtSPGR = qmt_spgr;
Model_qmtSPGR.Prot.MTdata.Mat = [710, 1200]; % 710 deg, 1.2 kHz

% Timing - TR = 30.3 ms
Model_qmtSPGR.Prot.TimingTable.Mat(1) = 0.0103;
Model_qmtSPGR.Prot.TimingTable.Mat(2) = 0.0032;
Model_qmtSPGR.Prot.TimingTable.Mat(3) = 0.0018;
Model_qmtSPGR.Prot.TimingTable.Mat(4) = 0.0150;
Model_qmtSPGR.Prot.TimingTable.Mat(5) = 0.3030; % 30.3 ms

Model_qmtSPGR.options.Readpulsealpha = 7;

% Set simulation options
Opt.Method = 'Analytical equation';
Opt.ResetMz = false;
Opt.SNR = 1000;

% Input parameters
% MTon
x = struct;
x.F = 0.16; % Vary F
x.kr = 30;
x.R1f = 1; % Vary R1f
x.R1r = 1;
x.T2f = 0.03;
x.T2r = 1.3e-05;

% Timing parameters for vfa_t1.analytical_solution function
% PDw - MToff
PDw_Model = vfa_t1;
paramsPDw.EXC_FA = 7;
paramsPDw.TR = 30.3; % ms

% T1w
T1w_Model = vfa_t1; 
paramsT1w.EXC_FA = 20;
paramsT1w.TR = 18; % ms

% Timing parameters for MTSAT_exec function
PDParams(1) = paramsPDw.EXC_FA;
PDParams(2) = paramsPDw.TR;
T1Params(1) = paramsT1w.EXC_FA;
T1Params(2) = paramsT1w.TR;
MTParams(1) = 6;
MTParams(2) = 32;
B1Params(1) = 1;

%%%%%%% VARY F AND R1F %%%%%%%
F = 0.1:0.01:0.30;
R1f = [0.35:0.01:0.50, 0.6:0.1:1.5];

diffMTR = zeros(length(F),length(R1f));
diffMTsat = zeros(length(F),length(R1f));

% Get signal using different F and T1w
for ii=1:length(F)
    % Update input params
    x.F = F(ii);
    R1obs = 1/convert_T1f_T1meas([x.F, x.kr*x.F, 1, 1/x.R1r],'t1f_2_t1meas');
    x.R1f = R1obs;
    [refMTR,refMTsat] = getMTR_MTsat(Model_qmtSPGR,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params);
    for jj=1:length(R1f)
        % Update input params
        x.R1f = R1f(jj);
        
        [MTR,MTsat] = getMTR_MTsat(Model_qmtSPGR,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params);
        diffMTR(ii,jj) = 100*(refMTR - MTR)/refMTR;
        diffMTsat(ii,jj) = 100*(refMTsat - MTsat)/refMTsat;
    end
end

function [MTR,MTsat] = getMTR_MTsat(Model,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params)
% Signal
Signal_qmtSPGR = equation(Model, x, Opt);
% PDw - MToff
paramsPDw.T1 = 1/x.R1f*1000; % ms
PDw = vfa_t1.analytical_solution(paramsPDw);
% T1w
paramsT1w.T1 = 1/x.R1f*1000; % ms
T1w = vfa_t1.analytical_solution(paramsT1w);
% MTon
MT = Signal_qmtSPGR*PDw;
        
% MTR calculation
MTR = 100*(PDw - MT)/PDw;
        
% MTsat calculation
dataMTsat.PDw = PDw;
dataMTsat.T1w = T1w;
dataMTsat.MTw = MT;
[MTsaturation,~] = MTSAT_exec(dataMTsat, MTParams, PDParams, T1Params, B1Params);
MTsat = MTsaturation;
end