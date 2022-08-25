%% Get MTR, MTsat and qMT using Bloch simulations
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
Opt.Method = 'Block equation';
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
B1Params(1) = 1;

% PDw - MToff
paramsPDw.T1 = 1/x.R1f*1000; % ms
PDw = vfa_t1.analytical_solution(paramsPDw);

% T1w
paramsT1w.T1 = 1/x.R1f*1000; % ms
T1w = vfa_t1.analytical_solution(paramsT1w);

%%%%%%% NUMBER OF REPETITIONS: BLOCH SIMULATIONS %%%%%%%
numPulses = [1:1:10, 20:10:90, 100:100:400];
MTR = zeros(1,length(numPulses));
MTsat = zeros(1,length(numPulses));
qMT = zeros(1,length(numPulses));
% Get signal using different F and T1w
for ii=1:length(numPulses)
    % Set number of pulses
    Model_qmtSPGR.options.MT_Pulse_NofMTpulses = numPulses(ii);
    
    % Signal (MTon is the only changing
    Signal_qmtSPGR = equation(Model_qmtSPGR, x, Opt);
    
    % MTon
    MT = Signal_qmtSPGR*PDw;
    qMT(ii) = MT;
    
    % MTR calculation
    MTR(ii) = 100*(PDw - MT)/PDw;
    
    % MTsat calculation
    dataMTsat.PDw = PDw;
    dataMTsat.T1w = T1w;
    dataMTsat.MTw = MT;
    [MTsaturation,~] = MTSAT_exec(dataMTsat, MTParams, PDParams, T1Params, B1Params);
    MTsat(ii) = MTsaturation;
end
SmodelBloch(1,:) = MTR;
SmodelBloch(2,:) = MTsat;
SmodelBloch(3,:) = qMT;