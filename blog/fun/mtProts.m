%% Simulate spectrum 
% Define protocol (Vendor, Field strength, MT pulse shape, MT pulse alpha, off-resonance)
customProt = struct();

%%% GE %%%
% 1
customProt.Vendor{1} = 'GE';
customProt.FieldStrength{1} = 1.5; % Tesla
customProt.excFA{1} = 7; % Excitation FA (deg)
customProt.MTshape{1} = 'sinc';
customProt.MTalpha{1} = 710; % deg
customProt.MTdur{1} = 10.3; % ms
customProt.offset{1} = 1200; % Hz

% 2
customProt.Vendor{2} = 'GE';
customProt.FieldStrength{2} = 1.5; % Tesla
customProt.MTshape{2} = 'sinc';
customProt.excFA{2} = 10; % Excitation FA (deg)
customProt.MTalpha{2} = 1388; % deg
customProt.MTdur{2} = 30.7; % ms
customProt.offset{2} = 1200; % Hz

% 3
customProt.Vendor{3} = 'GE';
customProt.FieldStrength{3} = 1.5; % Tesla
customProt.MTshape{3} = 'sinc';
customProt.excFA{3} = 7; % Excitation FA (deg)
customProt.MTalpha{3} = 1430; % deg
customProt.MTdur{3} = 64; % ms
customProt.offset{3} = 2000; % Hz

% 4
customProt.Vendor{4} = 'GE';
customProt.FieldStrength{4} = 1.5; % Tesla
customProt.MTshape{4} = 'fermi';
customProt.excFA{4} = 5; % Excitation FA (deg)
customProt.MTalpha{4} = 670; % deg
customProt.MTdur{4} = 8; % ms
customProt.offset{4} = 1600; % Hz

% 5
customProt.Vendor{5} = 'Siemens';
customProt.FieldStrength{5} = 1.5; % Tesla
customProt.MTshape{5} = 'gaussian';
customProt.excFA{5} = 30; % Excitation FA (deg)
customProt.MTalpha{5} = 500; % deg
customProt.MTdur{5} = 7.68; % ms
customProt.offset{5} = 1500; % Hz

% 6
customProt.Vendor{6} = 'Siemens';
customProt.FieldStrength{6} = 3; % Tesla
customProt.MTshape{6} = 'gaussian';
customProt.excFA{6} = 15; % Excitation FA (deg)
customProt.MTalpha{6} = 1000; % deg
customProt.MTdur{6} = 9.5; % ms
customProt.offset{6} = 1500; % Hz

% 7
customProt.Vendor{7} = 'Siemens';
customProt.FieldStrength{7} = 3; % Tesla
customProt.MTshape{7} = 'gaussian';
customProt.excFA{7} = 7; % Excitation FA (deg)
customProt.MTalpha{7} = 630; % deg
customProt.MTdur{7} = 15; % ms
customProt.offset{7} = 1000; % Hz

% 8
customProt.Vendor{8} = 'Philips';
customProt.FieldStrength{8} = 3; % Tesla
customProt.MTshape{8} = 'sincgauss';
customProt.excFA{8} = 18; % Excitation FA (deg)
customProt.MTalpha{8} = 700; % deg
customProt.MTdur{8} = 15; % ms
customProt.offset{8} = 1100; % Hz

% Simulate all protocols
for ii=1:length(customProt.Vendor)
    % Load qMT model
    Model_qmtSPGR = qmt_spgr;
    Model_qmtSPGR.Prot.MTdata.Mat = [customProt.MTalpha{ii}, customProt.offset{ii}]; % 710 deg, 1.2 kHz
    
    % Timing
    Model_qmtSPGR.Prot.TimingTable.Mat(1) = customProt.MTdur{ii}/1000; % MT pulse duration (s)
    Model_qmtSPGR.Prot.TimingTable.Mat(2) = 0.0032;
    Model_qmtSPGR.Prot.TimingTable.Mat(3) = 0.0015;
    Model_qmtSPGR.Prot.TimingTable.Mat(4) = 0.0100;
    Model_qmtSPGR.Prot.TimingTable.Mat(5) = Model_qmtSPGR.Prot.TimingTable.Mat(1) + Model_qmtSPGR.Prot.TimingTable.Mat(2) + Model_qmtSPGR.Prot.TimingTable.Mat(3) + Model_qmtSPGR.Prot.TimingTable.Mat(4);
    
    % MT Shape
    Model_qmtSPGR.options.MT_Pulse_Shape = customProt.MTshape{ii};
    % Fermi pulse
    if strcmp(customProt.MTshape{ii},'fermi')
        % MT pulse fermi transition Tmt/33.81
        Model_qmtSPGR.options.MT_Pulse_Fermitransitiona = Model_qmtSPGR.Prot.TimingTable.Mat(1)/33.81;
        %Model_qmtSPGR.options.MT_Pulse_Fermitransitiona = 0.35;
        Model_qmtSPGR.options.MT_Pulse_Shape = customProt.MTshape{ii};
    end
    
    % Excitation FA
    Model_qmtSPGR.options.Readpulsealpha = customProt.excFA{ii}; % deg 
    
    % Set simulation options
    Opt.Method = 'Analytical equation';
    Opt.ResetMz = false;
    Opt.SNR = 1000;

    % MT parameters for GM, WM, NAGM, NAWM, Lesion 1, Lesion 2 (Sled and Pike, 2001)
    % Pool size ratio
    F = [0.072, 0.161, 0.068, 0.150, 0.120, 0.094]; % [GM, WM, NAGM, NAWM, L1, L2]
    % Exchange rate (free pool)
    kf = [2.4, 4.3, 2.6, 4.9, 3.6, 2.7];
    kr = kf./F;
    % Relaxation rate (free pool)
    R1f = [0.93, 1.8, 0.89, 1.78, 1.52, 1.26]; % s
    % Tranverse relaxation time (free pool)
    T2f = [0.056, 0.037, 0.062, 0.038, 0.043, 0.052];
    % Tranverse relaxation time (restricted pool)
    T2r = [11.1e-6, 12.3e-6, 9.6e-6, 11.4e-6, 10.3e-6, 10.9e-6];
    
    % x gets different tissue params
    x = struct;
    x.R1r = 1;
    
    % Timing parameters for vfa_t1.analytical_solution function
    % PDw - MToff
    PDw_Model = vfa_t1;
    paramsPDw.EXC_FA = customProt.excFA{ii}; % deg
    paramsPDw.TR = Model_qmtSPGR.Prot.TimingTable.Mat(5); % ms
    
    % Calculate MTR
    for jj=1:length(F)
        % Update input params for qmt_spgr model
        x.F = F(jj);
        x.kr = kr(jj);
        x.R1f = R1f(jj);
        x.T2f = T2f(jj);
        x.T2r = T2r(jj);
        
        % Signal
        Signal_qmtSPGR = equation(Model_qmtSPGR, x, Opt);
        % PDw - MToff
        paramsPDw.T1 = 1/x.R1f*1000; % ms
        PDw = vfa_t1.analytical_solution(paramsPDw);
        % MTon
        MT = Signal_qmtSPGR*PDw;
        
        % MTR calculation
        MTR(ii,jj) = 100*(PDw - MT)/PDw;
    end
end