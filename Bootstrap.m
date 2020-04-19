
% code for Parametric Bootstrapping for the calculations reported in Purkayastha et al, Cellular
% Evolution on a Biomaterial, Under review
clear all
close all
T=readtable('Gene_expression_FPKM.csv');
genenames=T.geneID;
A1(:,1)=T.A_1_R1; A1(:,2)=T.A_1_R2; A1(:,3)=T.A_1_R3;
A1(:,4)=T.A_1_R4; A1(:,5)=T.A_1_R5;

A308(:,1)=T.A_308_R1;A308(:,2)=T.A_308_R2;
A308(:,3)=T.A_308_R3;A308(:,4)=T.A_308_R4;A308(:,5)=T.A_308_R5;

trunc=0;

% Both North and South cluster
E1(:,1)=T.E_1_R1; E1(:,2)=T.E_1_R2;E1(:,3)=T.E_1_R4;E1(:,4)=T.E_1_R5;
E1(:,5)=T.E_1_R6;E1(:,5)=T.E_1_R7;E1(:,5)=T.E_1_R8;E1(:,5)=T.E_1_R9;
[PC, GC, RV, RI, CRI, CRV, delta,Lo,Lp,La]= bootstrap(A1,A308,E1,trunc);
%

%North cluster bootstrap is run on four samples; see paper for definition
%of North cluster
clear E1
E1(:,1)=T.E_1_R6; E1(:,2)=T.E_1_R7;E1(:,3)=T.E_1_R8;E1(:,4)=T.E_1_R2; % north cluster
[PCnorth, GCnorth, RVnorth, RInorth, CRInorth, CRVnorth, deltanorth,Lonorth,Lpnorth,Lanorth]= bootstrap(A1,A308,E1,trunc);
%

% South cluster bootstrap is run on four samples; see paper for definition
% of North cluster
clear E1
E1(:,1)=T.E_1_R1;E1(:,2)=T.E_1_R4;E1(:,3)=T.E_1_R5;E1(:,4)=T.E_1_R9; %
[PCsouth, GCsouth, RVsouth, RIsouth, CRIsouth, CRVsouth, deltasouth,Losouth,Lpsouth,Lasouth]= bootstrap(A1,A308,E1,trunc);
%


GenesRevertinall=genenames(find(RV>950));
GenesRevertinnorth=genenames(find(RVnorth>950));

figure(1)
h1=histogram(deltadistr_all,'Normalization','probability');
hold on
h2=histogram(deltadistr_north,'Normalization','probability');
hold on
h3=histogram(deltadistr_south,'Normalization','probability');
legend('Full','North','South')
xlabel('delta')
ylabel('Probability')
figure(2)

h2= histogram(mean(GCnorth,1)./mean(PCnorth,1),1000,'Normalization','probability');
hold on
h3= histogram(mean(GCsouth,1)./mean(PCsouth,1),1000,'Normalization','probability');
legend('North','South')
xlim([-30 30])
xlabel('GC/PC')
ylabel('Probability')



function [PCtrunc, GCtrunc, RV, RI, CRI, CRV, delta,Lo,Lp,La]= bootstrap(A1,A308,E1,trunc)

Lp=mean(A1,2); % mean values and standard deviations for each gene for a given condiition are calculated
SE_Lp=std(A1,0,2);
Lo=mean(A308,2);
SE_Lo=std(A308,0,2);
La=mean(E1,2);

SE_La=std(E1,0,2);

counter1=0;counter2=0;k = 1;

for j = 1:length(Lp) % gene of interest
for i=1:1000 % number of sampling for each gene

    PC = normrnd(Lp(j),SE_Lp(j))-normrnd(Lo(j),SE_Lo(j)); % we assume that the FPKM values follow a normal distribution
    GC = normrnd(La(j),SE_La(j))-normrnd(Lp(j),SE_Lp(j));
       if abs(PC)>=trunc % truncating PC
           if abs(GC)>=trunc % truncating GC
               
               PCtrunc(k,j)=PC;
               GCtrunc(k,j)=GC;
               if sign(PCtrunc(k,j)*GCtrunc(k,j))==-1
                   counter1=counter1+1; % this is an event of reversion
                   
               else
                   counter2=counter2+1; %this is an event of reinforcement
               end
               
                k = k+1;
           end
       end

end
k=1;
RI(j)=counter2;
RV(j)=counter1;
counter1=0;
counter2=0;
end
CRI=length(find(RI>0.95*1000))/length(Lp)*100;
CRV=length(find(RV>0.95*1000))/length(Lp)*100;
delta=CRV-CRI;

end

function [deltadistr,meanGC,stdevGC,distPC,meanPC,stdevPC] = sampleGCPC(GC,PC)

k=1;
for j=1:size(GC,2)
    distGC=fitdist(GC(:,j),'Normal');
    meanGC(j)=distGC.mu;
    stdevGC(j)=distGC.sigma;
    distPC=fitdist(PC(:,j),'Normal');
    meanPC(j)=distPC.mu;
    stdevPC(j)=distPC.sigma;
end

for i=1:5000
    revert=0;
    reinforce=0;
for j=1:size(GC,2);
    
          GC1(j)=mean(mean(normrnd(meanGC(j),stdevGC(j),2))); % mean across four samples
        
          PC1(j)=mean(mean(normrnd(meanPC(j),stdevPC(j),2))); %mean across four samples
          
           if sign(PC1(j)*GC1(j))==-1
                   revert=revert+1; % this is an event of reversion         
               else
                   reinforce=reinforce+1; %this is an event of reinforcement
           end
                       
          
          
end
CRI=reinforce/size(GC,2)*100;
CRV=revert/size(GC,2)*100;
deltadistr(i)=CRV-CRI;

end
end