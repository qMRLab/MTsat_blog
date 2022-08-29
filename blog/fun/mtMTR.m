%% Get parameters that give the same MTR
% Simulate complete qMT-SPGR protocol
Model_qmtSPGR = qmt_spgr;
Model_qmtSPGR.Prot.MTdata.Mat = [490, 1200]; % 490°, 1.2 kHz

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

%%%%%%% VARY F AND R1F %%%%%%%
F = 0.05:0.02:0.40;
T1 = 0.5:0.05:1.9;
R1f = 1./T1;

% Get MTR
for ii=1:length(F)
    for jj=1:length(R1f)
        % Signal
        x.F = F(ii);
        x.R1f = R1f(jj);
        Signal_qmtSPGR = equation(Model_qmtSPGR, x, Opt);
        % PDw - MToff
        paramsPDw.T1 = 1/x.R1f*1000; % ms
        PDw = vfa_t1.analytical_solution(paramsPDw);
        % MTon
        MT = Signal_qmtSPGR*PDw;
        
        % MTR calculation
        MTR(ii,jj) = 100*(PDw - MT)/PDw;
        
        % Matrix of parameters
        paramMatrix(ii,jj) = struct('F', x.F, 'T1', 1/x.R1f*1000);
    end
end

% Round MTR matrix
rMTR = round(MTR);

% Get indices with equal MTR values
cnt = 1;
for i=1:100
    idx{i} = find(rMTR==i);
    if ~isempty(idx{i})
        repeatedMTR = length(idx{i});
        for j=1:repeatedMTR
            outMTR(cnt,1) = rMTR(idx{i}(j));
            outMTR(cnt,2) = paramMatrix(idx{i}(j)).F;
            outMTR(cnt,3) = paramMatrix(idx{i}(j)).T1;
            cnt = cnt + 1;
        end
    end
end