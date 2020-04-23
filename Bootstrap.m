% This program runs a parametric bootstrap as described in Purkayastha et
% al, April 23, 2020
clear all
close all
T=readtable('A3_1_Diff_fpkm.csv'); % FPKM values from mRNA seq experiments should be stored in this table
A1(:,1)=T.A_1_R1; A1(:,2)=T.A_1_R2; A1(:,3)=T.A_1_R3; % Ancestral cells on 1 kPa substrate
A1(:,4)=T.A_1_R4; A1(:,5)=T.A_1_R5; 
A308(:,1)=T.A_308_R1;A308(:,2)=T.A_308_R2; % Ancestral cells on 308 kPa substrate
A308(:,3)=T.A_308_R3;A308(:,4)=T.A_308_R4;A308(:,5)=T.A_308_R5;
trunc=0; % This is optional and set to zero here. One could truncate the PC or GC as required.
% Both North and South cluster % See paper for definition of North and
% South cluster
E1(:,1)=T.E_1_R1; E1(:,2)=T.E_1_R2;E1(:,3)=T.E_1_R4;E1(:,4)=T.E_1_R5; % Cells evolved on 1 kPa substrate for 90 d
E1(:,5)=T.E_1_R6;E1(:,5)=T.E_1_R7;E1(:,5)=T.E_1_R8;E1(:,5)=T.E_1_R9;
[PC, GC, RV, RI, CRI, CRV, delta,Lo,Lp,La]= bootstrap(A1,A308,E1,trunc); % Bootstrap runs the parametric bootstrap and returns 
% the quantities specified. See paper for definitions of these parameters. 

%
%North cluster bootstrap
clear E1
E1(:,1)=T.E_1_R6; E1(:,2)=T.E_1_R7;E1(:,3)=T.E_1_R8;E1(:,4)=T.E_1_R2; % north cluster
[PCnorth, GCnorth, RVnorth, RInorth, CRInorth, CRVnorth, deltanorth,Lonorth,Lpnorth,Lanorth]= bootstrap(A1,A308,E1,trunc);
%
% South cluster bootstrap
clear E1
E1(:,1)=T.E_1_R1;E1(:,2)=T.E_1_R4;E1(:,3)=T.E_1_R5;E1(:,4)=T.E_1_R9; %
[PCsouth, GCsouth, RVsouth, RIsouth, CRIsouth, CRVsouth, deltasouth,Losouth,Lpsouth,Lasouth]= bootstrap(A1,A308,E1,trunc);
%
[delta1,delta2,meanGC_all,stdevGC_all,distPC_all,meanPC_all,stdevPC_all]=sampleGCPC(GC,PC); % This is sampling two groups of four samples
% to calculate delta1 and delta2; see description in paper 

[deltadistr_north,meanGC_north,stdevGC_north,distPC_north,meanPC_north,stdevPC_north]=sampleGCPC(GCnorth,PCnorth);
[deltadistr_south,meanGC_south,stdevGC_south,distPC_south,meanPC_south,stdevPC_south]=sampleGCPC(GCsouth,PCsouth);

figure(1)
histogram(delta1-delta2,'Normalization','probability')
xlabel('{\delta_1}-{\delta_2}')
ylabel('Probability')

figure(2)
h2= histogram(mean(GCnorth,1)./mean(PCnorth,1),1000,'Normalization','probability');
hold on
h3= histogram(mean(GCsouth,1)./mean(PCsouth,1),1000,'Normalization','probability');
legend('North','South')
xlim([-30 30])
xlabel('GC/PC')
ylabel('Probability')
 
CV_Lo=std(A308,0,2)./Lo;
CV_Lp=std(A1,0,2)./Lp;
CV_La=std(E1,0,2)./La;
 
function [PCtrunc, GCtrunc, RV, RI, CRI, CRV, delta,Lo,Lp,La]= bootstrap(A1,A308,E1,trunc)
Lp=mean(A1,2);
SD_Lp=std(A1,0,2);
Lo=mean(A308,2);
SD_Lo=std(A308,0,2);
La=mean(E1,2);
SD_La=std(E1,0,2);
counter1=0;counter2=0;k = 1;
for j = 1:length(Lp) % gene of interest
for i=1:10000 % number of samples for each gene
    PC = normrnd(Lp(j),SD_Lp(j))-normrnd(Lo(j),SD_Lo(j));
    GC = normrnd(La(j),SD_La(j))-normrnd(Lp(j),SD_Lp(j));
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
CRI=length(find(RI>0.95*10000))/length(Lp)*100;
CRV=length(find(RV>0.95*10000))/length(Lp)*100;
delta=CRV-CRI;
end

function [delta1,delta2,meanGC,stdevGC,distPC,meanPC,stdevPC] = sampleGCPC(GC,PC)
k=1;
for j=1:size(GC,2)
    distGC=fitdist(GC(:,j),'Normal');
    meanGC(j)=distGC.mu;
    stdevGC(j)=distGC.sigma;
    distPC=fitdist(PC(:,j),'Normal');
    meanPC(j)=distPC.mu;
    stdevPC(j)=distPC.sigma;
end
for i=1:10000
    revert1=0;revert2=0;
    reinforce1=0;reinforce2=0;
for j=1:size(GC,2);
   
          GC1(j)=mean(mean(normrnd(meanGC(j),stdevGC(j),2))); % mean across four samples
          GC2(j)=mean(mean(normrnd(meanGC(j),stdevGC(j),2))); % mean across four samples
          PC1(j)=mean(mean(normrnd(meanPC(j),stdevPC(j),2))); %mean across four samples
          PC2(j)=mean(mean(normrnd(meanPC(j),stdevPC(j),2))); %mean across four samples
           if sign(PC1(j)*GC1(j))==-1
                   revert1=revert1+1; % this is an event of reversion        
               else
                   reinforce1=reinforce1+1; %this is an event of reinforcement
           end
                      
           if sign(PC2(j)*GC2(j))==-1
                   revert2=revert2+1; % this is an event of reversion        
               else
                   reinforce2=reinforce2+1; %this is an event of reinforcement
           end
         
end
CRI1=reinforce1/size(GC,2)*100;
CRV1=revert1/size(GC,2)*100;
delta1(i)=CRV1-CRI1;
CRI2=reinforce2/size(GC,2)*100;
CRV2=revert2/size(GC,2)*100;
delta2(i)=CRV2-CRI2;
end
end
