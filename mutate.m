function cellOut = mutate(cellIn, nMut )

A0 = cellIn.A0;
K0 = cellIn.K0;
K1 = cellIn.K1;
K2 = cellIn.K2;
K3 = cellIn.K3; 
K4 = cellIn.K4;
C  = cellIn.C;

for ii = 1:nMut
    mut = randi(6); %select mutation type
    N = length(A0); %number of equations
    while N <= 10 && mut == 6
        mut = randi(6);
    end
    switch mut
        case 1 %modify degradation rate      
            modProt = randi(N); %choose which protein to degrade
            del = 2*rand(1)*K1(modProt,modProt);
            K1(modProt,modProt) = del; %change the rate of the selected protein
        case 2 %modify random kinetic constant
            %list of unique numbers
            K0temp = K0(:);
            K1temp = K1(:);
            K2temp = K2(:);
<===============why????
            if ~isempty(K3) && ~isempty(K4)
                K3temp = K3(:,4);
                K3temp = K3temp(:);
                K4temp = K4(:,[3,4]);
                K4temp = K4temp(:);
                list = unique(abs([unique(K0temp); unique(K1temp); unique(K2temp); unique(K3temp); unique(K4temp)]));
            elseif ~isempty(K3) && isempty(K4)
                K3temp = K3(:,4);
                K3temp = K3temp(:);
                list = unique(abs([unique(K0temp); unique(K1temp); unique(K2temp); unique(K3temp)]));
            elseif isempty(K3) && ~isempty(K4)
                K4temp = K4(:,[3,4]);
                K4temp = K4temp(:);
                list = unique(abs([unique(K0temp); unique(K1temp); unique(K2temp); unique(K4temp)]));
            else
                list = unique(abs([unique(K0)'; unique(K1); unique(K2)]));
            end
            param = datasample(list(list~=0),1);
            modParam = 2*rand(1)*param;
            
            K0(K0==param) = modParam;
            K1(K1==param) = modParam;
            K2(K2==param) = modParam;
            K3(K3==param) = modParam;
            K4(K4==param) = modParam;
            
            K0(K0==-param) = -modParam;
            K1(K1==-param) = -modParam;
            K2(K2==-param) = -modParam;
            K3(K3==-param) = -modParam;
            K4(K4==-param) = -modParam;
            
%             %verify that promoter repression does not outweigh transcription
%             modProt = find(K0(:)==abs(modParam));
%             if ~isempty(modProt) && ~isempty(K4)
% % %                 keyboard
%                 rows = find(K4(:,1)==modProt);
%                 for jj = 1:length(rows)
%                     if K4(rows(jj),3) < 0 && K0(modProt) < -K4(rows(jj),3)
%                         K4(rows(jj),3) = -K0(modProt);
%                     end
%                 end
%             end
        case 3 %add a new protein
            N = N + 1;
            A0 = [A0; 1]; %set new initial condition
            C  = [C; 0 0]; %label new variable as not protein complex

            tau = rand(1); %add new transcription rate
            K0(N) = tau;

            del = rand(1); %add new degradation rate
            K1(N,N) = -del;

            K2(N,N) = 0; %add new empty row/column to K2
        case 4
            modPromot = randi(N); %choose protein whose promoter will be modified
            while sum(C(modPromot,:)) > 0 %make sure variable is actually not a complex or phosphorylated
                modPromot = randi(N);
            end
            modProt = randi(N); %choose binding protein
<============why?where are the 3 reactions?????
            tauP = 2*rand(1)-1;%%?????why
%             tau  = K0(modPromot);
%             if tauP < 0 %prevent negative transcription, which leads to negative concentration
%                 tauP = max(tauP,-tau);
%             end
            Kd = rand(1);

            K4 = [K4;modPromot modProt tauP Kd];
        case 5
            if randi(2) == 1 %select one protein
                modProt = randi(N); %choose which protein to modify
                if C(modProt,1) == 0%if not a protein complex
                    N = N + 1;
                    A0 = [A0; 1]; %set new initial condition
                    C  = [C; 0 1]; %label new variable as not a protein complex and as phosphorylated

                    tau = rand(1); %add new phosphorylation rate

                    K0(N) = 0;
                    K1(N,N) = 2*rand(1)*K1(modProt,modProt); %give phosphorylated protein a modified degradation rate
                    K1(N,modProt) = tau;
                    K1(modProt,modProt) = K1(modProt,modProt) - tau;

                    K2(N,N) = 0; %add new empty row/column to K2
                else %if a protein complex, partial degradation
                    %figure out what variables are in the complex
                    try
                        row = datasample(find(K3(:,1)==modProt),1);
                    catch
                        keyboard
                    end
                    var1 = K3(row,2);
                    var2 = K3(row,3);

                    del = rand(1);
                    K1(modProt,modProt) = -del;
                    if randi(2) == 1
                        %degrade var2 (only var1 gets replenished)
                        K1(var1,modProt) = del;
                    else
                        %degrade var1 (only var2 gets replenished)
                        K1(var2,modProt) = del;
                    end
                end   
            else %select two proteins
                modProt1 = randi(N);
                modProt2 = randi(N);
                if randi(2) == 1 %dimerize, add new variable for complex
                    N = N + 1;
                    A0 = [A0; 1]; %set new initial condition
                    C  = [C; 1 0]; %label new variable as a protein complex

                    gamma = 2*rand(1);
                    del   = 2*rand(1);

                    K0(N) = 0;
                    K1(N,N) = -del;

                    K2(N,N) = 0; %add new empty row/column to K2
                    K2(modProt1,modProt2) = -gamma; %create a term in equation modprot1 that is -gamma*modProt1*modProt2
                    K2(modProt2,modProt1) = -gamma; 

                    K3 = [K3; N modProt1 modProt2 gamma];
                else %degradation
                    nComplex = C(modProt1,1) + C(modProt2,1);
                    if nComplex == 1 %if only one is a complex
                        %partial catalytic degradation
                        if C(modProt1) == 1 %figure out which protein is the complex
                            compProt = modProt1;
                            enzProt  = modProt2;
                        else
                            compProt = modProt2;
                            enzProt  = modProt1;
                        end
                        row = datasample(find(K3(:,1)==compProt),1); %figure out which row in K3 tells us what the complex is made of
                        if randi(2) == 1 %var 2 within complex is the thing that gets degraded, and var1 gets returned to the pool
                            var1 = K3(row,2);
                            var2 = K3(row,3);
                        else
                            var1 = K3(row,3);
                            var2 = K3(row,2);
                        end
                        del = rand(1); %generate degradation rate
                        K2(compProt,enzProt) = -del;
                        K2(enzProt,compProt) = -del;
                        K3 = [K3;var1 compProt var2 del];
                    else %catalytic degradation if neither are complexes or if both are complexes
                        del = rand(1);
                        if randi(2) == 1 %degrade modProt1
                            K2(modProt1,modProt2) = -del;
                        else %degrade modProt2
                            K2(modProt2,modProt1) = -del;
                        end
                    end
                end
            end
        case 6 %loss of function mutation (gene deletion)
            modProt = randsample(N-1,1)+1; %cannot remove the first protein (which we are tracking)
            [~,N,A0,C,K0,K1,K2,K3,K4] = removeProt(modProt,[],N,A0,C,K0,K1,K2,K3,K4);
    end
end

if ~isempty(K4)
    for i = 1:size(K4,1)
        modProt = K4(i,1);
        if K4(i,3) < 0 && K0(modProt) < -K4(i,3)
            K4(i,3) = -K0(modProt);
        end
    end
end

cellOut.A0 = A0;
cellOut.K0 = K0;
cellOut.K1 = K1;
cellOut.K2 = K2;
cellOut.K3 = K3; 
cellOut.K4 = K4;
cellOut.C  = C;

end

