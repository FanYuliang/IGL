function cell = initialize( nProt )

%create first protein
tau = rand(1);
del = rand(1);

K0 = tau;
K1 = -del;
K2 = 0;
K3 = [];
K4 = [];

A0 = 1;
C  = [0 0]; %initial variable is not a protein complex

%create additional proteins
for i = 2:nProt
    %add a new protein
    N = length(A0) + 1;
    A0 = [A0; 1]; %set new initial condition
    C  = [C; 0 0]; %label new variable as not protein complex or phosphorylated

    tau = rand(1); %add new transcription rate
    K0(N) = tau;

    del = rand(1); %add new degradation rate
    K1(N,N) = -del;
    K2(N,N) = 0; %add new empty row/column to K2
end

%store variables in cell structure
cell.A0 = A0;
cell.K0 = K0;
cell.K1 = K1;
cell.K2 = K2;
cell.K3 = K3;
cell.K4 = K4;
cell.C  = C;

end

