function [ eqns, N, A0, C, K0, K1, K2, K3, K4 ] = removeProt( modProt, eqns, N, A0, C, K0, K1, K2, K3, K4 )
% keyboard
N = N-1; %reduce variable count
A0(modProt) = []; %remove initial condition
C(modProt,:)  = []; %remove complex index entry
K0(modProt) = []; %remove transcription rate entry
K1(modProt,:) = []; %remove modProt row
K1(:,modProt) = []; %remove modprot column
K2(modProt,:) = [];
K2(:,modProt) = [];

%empty out any promoter interactions associated with modProt
for jj = 1:2
    if ~isempty(K4)
        rows = find(K4(:,jj)==modProt);
        K4(rows,:) = [];
    end
end
K4(K4>modProt) = K4(K4>modProt) - 1;

%remove all complexes that involve modProt - recurses to include complexes
%of complexes
% keyboard
if ~isempty(K3)
    rows = find(K3(:,1)==modProt);
    K3(rows,:) = [];
    rows = unique([find(K3(:,2)==modProt);find(K3(:,3)==modProt)]);
    eqns = sort([eqns(2:end); K3(rows,1)]);
    eqns(eqns>modProt) = eqns(eqns>modProt)-1;
    K3(K3>modProt) = K3(K3>modProt) - 1; %because i removed modProt
    K3(rows,:) = [];
else
    eqns = eqns(2:end);
    eqns(eqns>modProt) = eqns(eqns>modProt) - 1;
end
% keyboard
if ~isempty(eqns)
    [eqns, N,A0,C,K0,K1,K2,K3,K4] = removeProt(eqns(1),eqns,N,A0,C,K0,K1,K2,K3,K4);
    if ~isempty(eqns)
        eqns(eqns>eqns(1)) = eqns(eqns>eqns(1)) - 1;
    end
end
% keyboard
end

