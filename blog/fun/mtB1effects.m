%% Get MTR and MTsat percentage difference varying the B1 map quality
% Simulate complete qMT-SPGR protocol
Model_qmtSPGR = qmt_spgr;
Model_qmtSPGR.Prot.MTdata.Mat = [490, 1200]; % 490�, 1.2 kHz

% Timing - TR = 32 ms
Model_qmtSPGR.Prot.TimingTable.Mat(1) = 0.0120;
Model_qmtSPGR.Prot.TimingTable.Mat(2) = 0.0032;
Model_qmtSPGR.Prot.TimingTable.Mat(3) = 0.0018;
Model_qmtSPGR.Prot.TimingTable.Mat(4) = 0.0150;
Model_qmtSPGR.Prot.TimingTable.Mat(5) = 0.0320;

Model_qmtSPGR.options.Readpulsealpha = 6;

% Set simulation options
Opt.Method = 'Analytical equation';
Opt.ResetMz = false;
Opt.SNR = 1000;

% Input parameters
% MTon
x = struct;
x.F = 0.16;
x.kr = 30;
x.R1f = 1;
x.R1r = 1;
x.T2f = 0.03;
x.T2r = 1.3e-05;

% Timing parameters for vfa_t1.analytical_solution function
% PDw - MToff
PDw_Model = vfa_t1;
paramsPDw.EXC_FA = 6;
paramsPDw.TR = 32; % ms

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
B1Params(1) = 0.4;

%%%%%% VARY F AND B1 %%%%%%%
F = 0.1:0.01:0.30;
B1map = 0.7:0.05:1.3;

diffMTR = zeros(length(F),length(B1map));
diffMTsat = zeros(length(F),length(B1map));

% Get MTR and MTsat using different B1 maps
for ii=1:length(F)
    % Update input parameteres
    x.F = F(ii);
    inputB1 = 1;
    [refMTR,refMTsat] = getMTR_MTsat(Model_qmtSPGR,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params,inputB1);
    for jj=1:length(B1map)
        inputB1 = B1map(jj);
        [MTR,MTsat] = getMTR_MTsat(Model_qmtSPGR,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params,inputB1);
        diffMTR(ii,jj) = 100*(refMTR - MTR)/refMTR;
        diffMTsat(ii,jj) = 100*(refMTsat - MTsat)/refMTsat;
    end
end

function [MTR,MTsat] = getMTR_MTsat(Model,x,Opt,paramsPDw,paramsT1w,MTParams,PDParams,T1Params,B1Params,inputB1)
% Signal
Model.Prot.MTdata.Mat(1) = Model.Prot.MTdata.Mat(1)*inputB1;
Signal_qmtSPGR = equation(Model, x, Opt);
% PDw - MToff
paramsPDw.T1 = 1/x.R1f*1000; % ms
paramsPDw.EXC_FA = paramsPDw.EXC_FA*inputB1;
PDw = vfa_t1.analytical_solution(paramsPDw);
% T1w
paramsT1w.T1 = 1/x.R1f*1000; % ms
paramsT1w.EXC_FA = paramsT1w.EXC_FA*inputB1;
T1w = vfa_t1.analytical_solution(paramsT1w);
% MTon
MT = Signal_qmtSPGR*PDw;
        
% MTR calculation
MTR = 100*(PDw - MT)/PDw;
        
% MTsat calculation
dataMTsat.PDw = PDw;
dataMTsat.T1w = T1w;
dataMTsat.MTw = MT;
dataMTsat.B1map = inputB1;
[MTsaturation,~] = MTSAT_exec(dataMTsat, MTParams, PDParams, T1Params, B1Params);
MTsat = MTsaturation;
end