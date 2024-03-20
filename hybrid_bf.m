clearvars;
close all;
clc



N_tx = 64; % N_tx antenna
N_tx_RF = 4; % N_rf chain

N_rx = 16; % N_rx antenna
N_rx_RF = 4; % N_rf chain on rx side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numBits = 1e6; % Total number of bits
modOrder = 4;
bitsPerSymbol = log2(modOrder);

% Generate random bits
dataBits = randi([0 1], numBits, 1);

% Assume that each antenna is connected to all RF chains (4 phase shifters)

rng(4096);
light_spd = 3e8;
fc = 28e9;
lambda = light_spd/fc;
tx_array = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(N_tx) sqrt(N_tx)],lambda/2),...
    'SubarraySelection',ones(N_tx_RF,N_tx),'SubarraySteering','Custom');
rx_array = phased.PartitionedArray(...
    'Array',phased.URA([sqrt(N_rx) sqrt(N_rx)],lambda/2),...
    'SubarraySelection',ones(N_rx_RF,N_rx),'SubarraySteering','Custom');


% RF chain can be used to send an independent data stream. 
% In this case, the system can support up to 4 streams.

N_scattering_clus = 6;
N_scatt = 8;
N_tot_scatter = N_scatt*N_scattering_clus;
angle_spread = 5;
% compute randomly placed scatterer clusters
txclang = [rand(1,N_scattering_clus)*120-60;rand(1,N_scattering_clus)*60-30];
rxclang = [rand(1,N_scattering_clus)*120-60;rand(1,N_scattering_clus)*60-30];
txang = zeros(2,N_tot_scatter);
rxang = zeros(2,N_tot_scatter);
% compute the rays within each cluster
for m = 1:N_scattering_clus
    txang(:,(m-1)*N_scatt+(1:N_scatt)) = randn(2,N_scatt)*sqrt(angle_spread)+txclang(:,m);
    rxang(:,(m-1)*N_scatt+(1:N_scatt)) = randn(2,N_scatt)*sqrt(angle_spread)+rxclang(:,m);
end

g = (randn(1,N_tot_scatter)+1i*randn(1,N_tot_scatter))/sqrt(N_tot_scatter);


% The channel matrix can be formed as

txpos = getElementPosition(tx_array)/lambda;
rxpos = getElementPosition(rx_array)/lambda;
H = scatteringchanmtx(txpos,rxpos,txang,rxang,g);

% Hybrid Weights Computation

figure;
F = diagbfweights(H);
F = F(1:N_tx_RF,:);
pattern(tx_array,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',F','PropagationSpeed',light_spd);



% The hybrid weights, on the other hand, can be computed as
At = steervec(txpos,txang);
Ar = steervec(rxpos,rxang);

Ns = N_tx_RF;
[Fbb,Frf] = omphybweights(H,Ns,N_tx_RF,At);

figure;
pattern(tx_array,fc,-90:90,-90:90,'Type','efield',...
    'ElementWeights',Frf'*Fbb','PropagationSpeed',light_spd);

% Spectral Efficiency Comparison

snr_param = -40:2:-5;
Nsnr = numel(snr_param);
Ns_param = [1 2];
NNs = numel(Ns_param);

N_tx_RF = 4;
N_rx_RF = 4;

Ropt = zeros(Nsnr,NNs);
Rhyb = zeros(Nsnr,NNs);
Niter = 50;

BER_hyb = cell(1, 2); BER_fulldig = cell(1, 2);

for m = 1:Nsnr
    snr = db2pow(snr_param(m));
    for n = 1:Niter
        % Channel realization
        txang = [rand(1,N_tot_scatter)*60-30;rand(1,N_tot_scatter)*20-10];
        rxang = [rand(1,N_tot_scatter)*180-90;rand(1,N_tot_scatter)*90-45];
        At = steervec(txpos,txang);
        Ar = steervec(rxpos,rxang);
        g = (randn(1,N_tot_scatter)+1i*randn(1,N_tot_scatter))/sqrt(N_tot_scatter);
        H = scatteringchanmtx(txpos,rxpos,txang,rxang,g);
        
        for k = 1:NNs
            Ns = Ns_param(k);
            % Compute optimal weights and its spectral efficiency
            [Fopt,Wopt] = helperOptimalHybridWeights(H,Ns,1/snr);
            Ropt(m,k) = Ropt(m,k)+helperComputeSpectralEfficiency(H,Fopt,Wopt,Ns,snr);

            % Compute hybrid weights and its spectral efficiency
            [Fbb,Frf,Wbb,Wrf] = omphybweights(H,Ns,N_tx_RF,At,N_rx_RF,Ar,1/snr);
            Rhyb(m,k) = Rhyb(m,k)+helperComputeSpectralEfficiency(H,Fbb*Frf,Wrf*Wbb,Ns,snr);
            
            if m == 1 && n == 1
                [BER_bf_out, BER_opt_out] = simulateBeamformingBER(dataBits, bitsPerSymbol, Ns, modOrder, snr_param, Frf, Fbb, Wrf, Wbb, Fopt, Wopt, H);
                BER_hyb{Ns} =  BER_bf_out; BER_fulldig{Ns} = BER_opt_out;
                
            end

        end
    end
end

Ropt = Ropt/Niter;
Rhyb = Rhyb/Niter;

figure;
plot(snr_param, Ropt(:,1), '-.o', 'Color', 'red', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(snr_param, Ropt(:,2), '-.d', 'Color', 'green', 'LineWidth', 1.5, 'MarkerSize', 8);

% For the second pair, using black for Rhyb(:,1) and blue for Rhyb(:,2)
plot(snr_param, Rhyb(:,1), '--*', 'Color', 'black', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(snr_param, Rhyb(:,2), '--x', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 8);

xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)');
title('Spectral Efficiency VS SNR (1 & 2 Data Streams)')
legend({'Ns=1 fully digital', 'Ns=2 fully digital', 'Ns=1 hybrid', 'Ns=2 hybrid'}, 'Location', 'best');
grid on; hold off;



figure;
semilogy(snr_param, BER_fulldig{1}, '-.o', 'Color', 'red', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(snr_param, BER_hyb{1}, '--*', 'Color', 'green', 'LineWidth', 1.5, 'MarkerSize', 8);

% For the second pair, using black for full digital and blue for hybrid
semilogy(snr_param, BER_fulldig{2}, '-.d', 'Color', 'black', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(snr_param, BER_hyb{2}, '--x', 'Color', 'blue', 'LineWidth', 1.5, 'MarkerSize', 8);

xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER VS SNR (1 & 2 Data Streams)')
legend({'Ns=1 Fully digital', 'Ns=2 Hybrid BF', 'Ns=2 Fully digital', 'Ns=2 Hybrid BF'}, 'Location', 'best');
grid on; hold off;
