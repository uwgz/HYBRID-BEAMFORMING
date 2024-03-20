function [BER_bf, BER_opt] = simulateBeamformingBER(dataBits, bitsPerSymbol, Ns, modOrder, snr_param, Frf, Fbb, Wrf, Wbb, Fopt, Wopt, H)
    % Reshape bits into numStreams rows
    reshapedBits = reshape(dataBits, bitsPerSymbol * Ns, []);

    % Modulate the bit stream (QPSK)
    dataSyms = pskmod(reshapedBits, modOrder, pi/4, 'gray', 'InputType', 'bit');

    % Initialize BER storage
    BER_bf = zeros(1, length(snr_param));
    BER_opt = zeros(1, length(snr_param));

    for m = 1:length(snr_param)
        snr = snr_param(m);

        %% For hybrid beamforming
        % Apply transmit beamforming weights
        txSignal = (Frf.' * Fbb.') * dataSyms;

        % Pass the signal through the channel (with noise)
        rxSignal = awgn(H.' * txSignal, snr, 'measured');

        % Apply receive beamforming weights
        rxBeamformed = (Wrf * Wbb).' * rxSignal;

        % Demodulate
        rxDataBits = pskdemod(rxBeamformed, modOrder, pi/4, 'gray', 'OutputType', 'bit');

        %% For optimal beamforming
        % Pass the signal through the channel (with noise) for optimal beamforming
        txSignal_opt = Fopt.' * dataSyms;
        rxSignal_opt = awgn(H.' * txSignal_opt, snr, 'measured');

        % Apply receive beamforming weights (optimal)
        rxBeamformed_opt = Wopt.' * rxSignal_opt;

        % Demodulate
        rxDataBits_opt = pskdemod(rxBeamformed_opt, modOrder, pi/4, 'gray', 'OutputType', 'bit');

        % Calculate and store BER for both beamforming techniques
        [~, BER_opt(m)] = biterr(dataBits, rxDataBits_opt(:));
        [~, BER_bf(m)] = biterr(dataBits, rxDataBits(:));
    end
end
