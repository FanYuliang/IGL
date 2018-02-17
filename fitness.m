function output = fitness( A, tMax, per, varargin )
%A - signal vector
%tMax - length of timespan of signal
%P - period being selected for
N = length(A);
Ak = fft(A);
pow = sqrt(Ak.*conj(Ak));

if mod(N,2) == 0
    P = pow(2:N/2); %exclude 0 mode
    f = (1:N/2 - 1)/tMax; %frequency
else
    P = pow(2:(N-1)/2); %exclude 0 mode
    f = (1:(N-1)/2 - 1)/tMax; %frequency
end
T = 1./f; %period

if ~isempty(varargin)
    plotPower = varargin{1};
    if plotPower == true
        stem(f,P)
        vline(1/per)
        xlabel('frequency')
        ylabel('power')
    end
end
    
[m,i] = min(sqrt((T - per).^2));
if m < per/5
%     output = P(i)/sum(P); %return the fraction of total power in the selected period.
    try
        output = sum(P(i-2:i+2))/sum(P);
    catch
        disp('Error: signal not well enough resolved for selected period.')
        output = NaN;
    end
else
    disp('Error: signal not well enough resolved for selected period.')
    output = NaN;
end

end

