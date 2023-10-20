clear all
close all
clc

load InPutData

Pattern = InPutSequences;
TimePattern = InPutTimeStamps;

Gender = [zeros(length(Pattern),1), [ones(length(Pattern)/2,1); ones(length(Pattern)/2,1)*2]];

clear InPutSequences InPutTimeStamps

%% Write sequences for SPMF file (Elapsed Time: 2.24 sec)

names={'input.txt'};

str=[];

fid = fopen(string(names(1 )),'w');
for i=1:size(Pattern,1)
        str=Pattern{i,:};
    if i==length(Pattern)
        str=join([string(str(1:end-1)) string([char(string(str(end))) char('')])]);
        fprintf(fid,str);
        continue
    end
    str=join([string(str(1:end-1)) string([char(string(str(end))) char('\n')])]);
    fprintf(fid,str);    
end

fclose(fid);

clear str fid


%% Exectuting Sequence Mining

minsup=0.003;
MinL=1;
MaxL=12;
Gap=200; % to infinity

runner=['!java -jar spmf.jar run CM-SPAM ','input.txt ','output.txt ', num2str(minsup),' ',num2str(MinL),' ',num2str(MaxL),' ','""',' ',num2str(Gap), ' true'];
eval(runner)

%% Read the outputs of Sequence Analysis


Rules=[];
Support=[];
RuleLength=[];

fid = fopen('output.txt');
tline = fgetl(fid);
k=0;

while ischar(tline)
    k=k+1;
    parts=strsplit(tline, { '#SUP: ', '#SID: '});
    dum=str2num(parts{1});
    Rules{k,1}=dum; % Rules ID
    RuleLength(k,1)=length(dum)/2;
    Support(k,1)=str2num(parts{2});
    TID{k,1}=str2num(parts{3})+1;
    tline = fgetl(fid);
end
fclose(fid);

clear parts dum dumstr tline

%% Matching Rules

%Match=30;
MatchNum = [15];

MatchedRules=zeros(size(Rules,1),1);


for i=1:length(MatchedRules)
    dum=Rules{i,1};

    if sum(ismember(MatchNum,dum))==length(MatchNum)
        if length(dum(dum~=-1))==length(MatchNum)
        MatchedRules(i,1)=5;    
        end
        dummem=ismember(dum,MatchNum);
        if dummem(end-1)~=1
        MatchedRules(i,1)=1;
        end
    end
end


MatchedRules=find(MatchedRules); % Saving results
MatchedRules(:,2)=zeros(length(MatchedRules),1); % Prepearing the next Consistation vector


%% Finding Matchies and their extensions
% Match BNO-s must be on first places, it checked here
for i=1:size(MatchedRules,1)
    
    dum=Rules{MatchedRules(i),1};
    dum=dum(dum~=-1);
    if isequal(dum(1:length(MatchNum)),MatchNum) & length(dum)>length(MatchNum) % find Extension
        MatchedRules(i,2)=1;
    elseif isequal(dum(1:length(MatchNum)),MatchNum) & length(dum)==length(MatchNum) % find the Base
        MatchedRules(i,2)=2;
    end   
end

BaseRules=MatchedRules(MatchedRules(:,2)==2,1);
ExtendedRules=MatchedRules(MatchedRules(:,2)==1,1);

%% Gathering times

for i=1:size(ExtendedRules,1)
    
    dumID=ExtendedRules(i,1);
    dumRule=Rules{dumID,:};
    dumTID=TID{dumID,:};
    MatchNumCont=dumRule(find(dumRule==MatchNum(end))+2);
    
    dumSurT=[];
    dumGender=[];
    counter=0;
    male=0;
    female=0;
    dumCorrTID=[];
    for k=1:length(dumTID)
    
        dumBNO=Pattern{dumTID(k),:};
        dumTime=TimePattern{dumTID(k),:};
        dumTimeL=dumTime(find(dumBNO==MatchNum(end)));
        dumTimeH=dumTime(find(dumBNO==MatchNumCont));
        
        if isempty(dumTimeL)
            dumSurT(k)=NaN;
            dumGender(k)=NaN;
            continue
        end
        
        if (dumTimeH)-(dumTimeL)<0
            dumSurT(k)=NaN;
            dumGender(k)=NaN;
            continue
        end 
        
        
        dumSurT(k)=[dumTimeH-dumTimeL]/365;
        counter=counter+1;
        dumCorrTID(counter)=dumTID(k);
        dumGender(k)=Gender(dumTID(k),2);
        
        if Gender(dumTID(k),2)==1
            male=male+1;
        end
        
        if Gender(dumTID(k),2)==2
            female=female+1;
        end
    end
    
    CorrTID{i,1}=dumCorrTID;
    CorrSupport(i,1)=counter;    
    SurTime{i,1}=dumSurT;
    [dumy dumx]=ecdf(dumSurT,'function','survivor');
    SAFunc{i,1}=[dumx dumy];
    
    
    CorrGenSupport(i,1)=male;
    CorrGenSupport(i,2)=female; 
    
    if isempty(dumSurT(dumGender==1))
        SAFuncMale{i,1}=[NaN NaN];
    else
        
        [dumy dumx]=ecdf(dumSurT(dumGender==1),'function','survivor');
        SAFuncMale{i,1}=[dumx dumy];
        Time1=[dumSurT(dumGender==1)' ones(sum(dumGender==1),1)];
    end
    
    if isempty(dumSurT(dumGender==2))
        SAFuncFemale{i,1}=[NaN NaN];
    else
        
        [dumy dumx]=ecdf(dumSurT(dumGender==2),'function','survivor');
        SAFuncFemale{i,1}=[dumx dumy];
        Time2=[dumSurT(dumGender==2)' 2*ones(sum(dumGender==2),1)];
    end
    Time=[Time1;Time2];
    Time(isnan(Time(:,1)),:)=[];
    GenTime{i,1}=Time;
    Ttest(i,1)=kruskalwallis(Time(:,1),Time(:,2),'off');

    
end

clear dumBNO dumCorrTID dumID dummem dumRule dumSurT dumTimeH dumTimeL dumx dumy Time1 Time2 Time male female


%% Eliminating tricky base sickness based in diagnose day and the last date of the dataset

Jelol={'b-','bo-','bx-','b+-','bd-','r-','ro-','rx-','r+-','rd-','g-','go-','gx-','g+-','gd-','k-','ko-','kx-','k+-','kd-','b--','bo--','bx--','b+--','bd--','r--','ro--','rx--','r+--','rd--','g--','go--','gx--','g+--','gd--','k--','ko--','kx--','k+--','kd--'};

initL = 732313;
initU = 733408;


EndDate=initU; % Last record of data in time
StartDate=initL; % First record of data in time

YearBackTime=1;

dumTID=TID{BaseRules(1),1};
counter=0;
count=1;
male=0;
female=0;
info=[];
for i=1:length(dumTID)
    dumRule=Pattern{dumTID(i),1};
    dumindex=find(dumRule==MatchNum);
    dumTime=TimePattern{dumTID(i),1}(dumindex);
    dumdiff=EndDate-dumTime;
    
    if isempty(dumdiff)
        counter=counter+1;
        continue
    end

    if (dumdiff) < YearBackTime*365
        counter=counter+1;
    end
    CorrTIDBase(1,count)=dumTID(i);
    count=count+1;
    
    if Gender(dumTID(i),2)==1
        male=male+1;
    end
    
    if Gender(dumTID(i),2)==2
        female=female+1;
    end
end


BaseGenSupport(1,1)=male;
BaseGenSupport(1,2)=female;
BaseSupport=Support(BaseRules(1),1)-counter;

clear dumRule dumindex dumTime dumdiff counter female male

%% Calculating confidence


for j=1:size(SAFunc,1)
    Conf(j,1)=CorrSupport(j,1)/BaseSupport;
    ConfMale(j,1)=CorrGenSupport(j,1)/BaseGenSupport(1,1);
    ConfFemale(j,1)=CorrGenSupport(j,2)/BaseGenSupport(1,2);
end

clear dumRule dumindex dumTime dumdiff counter female male

%% Confidence Function

PlotVar=2; 

figure(3)
hold on

Test=[];

Legen=[];
for j=PlotVar

    
    SAVal=SAFuncMale{j,1};
    
    Legen=[Legen; string('Male')];
    stairs(SAVal(:,1),(1-SAVal(:,2))*ConfMale(j,1),'b','LineWidth',2)
end

Legen=string(Legen);
legend(Legen)

k=1;

for j=PlotVar

    
    SAVal=SAFuncFemale{j,1};

    Legen=[Legen; string('Female')];
    stairs(SAVal(:,1),(1-SAVal(:,2))*ConfFemale(j,1),'r','LineWidth',2)
end


Legen=string(Legen);
legend(Legen)
xlabel('Time [Years]')
ylabel('Probability')
disp('The result of statistical test')

Time=GenTime{j,1};
%logrank(Time(Time(:,2)==1,1)',Time(Time(:,2)==2,1)',0)
ax=gca;
ax.FontSize=12;



%% Plot functions
% Bootstrapping
clear suppB BSSAx BSSAy BSSurT BSConf BSSurY

it=150;
BSn=length(Pattern)*0.4;
suppB=[];
alpha=0.05;

for ik=1:it
    
    
    BSId=randperm(length(Pattern),BSn);
    suppA=sum(ismember(CorrTIDBase,BSId));
    
    for i=PlotVar
        dumTID=ismember(CorrTID{i,1},BSId);
        suppB(ik,i)=sum(dumTID);
        dumTID=CorrTID{i,1}(dumTID);
        
        
        
        dumID=ExtendedRules(i,1);
        dumRule=Rules{dumID,:};
        
        
        
        dumF=find(dumRule==-1);
        
        if length(dumF)<2
            continue
        end
        

        indexN=dumRule(dumF(end-1)+1);
        indexO=dumRule(dumF(end-1)-1);
        dumSurT=[];
        
        for k=1:length(dumTID)
            
            dumBNO=Pattern{dumTID(k),:};
            dumTime=TimePattern{dumTID(k),:};
            dumTimeL=dumTime(dumBNO==indexO);
            dumTimeH=dumTime(dumBNO==indexN);
            
            if isempty(dumTimeL)|isempty(dumTimeH)
                dumSurT(k)=NaN;
                continue
            end
            
            if (dumTimeH)-(dumTimeL)<0
                dumSurT(k)=NaN;
                continue
            end
            
            dumSurT(k)=[dumTime(dumBNO==indexN)-dumTime(dumBNO==indexO)]/365;
            
        end
        BSSurT{ik,i}=dumSurT;
        if (~isempty(dumSurT)) & (sum(isnan(dumSurT))<length(dumSurT))
            [dumy dumx]=ecdf(dumSurT,'Function','survivor');
            BSSAx{ik,i}=dumx;
            BSSAy{ik,i}=(1-dumy)*(suppB(ik,i)/suppA);
        end
        
    end
    
    
end

BSConf=suppB/suppA;

BSLC = prctile(BSConf,(alpha/2)*100);
BSHC = prctile(BSConf,(1-(alpha/2))*100);

index=1;

domf=linspace(0,max(SAFunc{PlotVar,1}(:,1)),1000);

figure(10)
hold on
for i=1:size(BSSAx,1)
    
    dumx=BSSAx{i,PlotVar};
    if isempty(dumx)
        continue
    end
    dumy=BSSAy{i,PlotVar};
    index=1;
    indexF=2;
    for j=domf
        if j>dumx(indexF)
            indexF=indexF+1;
            if indexF>length(dumx)
                indexF=indexF-1;%length(dumx);
                BSSurY(i,index)=dumy(indexF);
                index=index+1;
                continue
            end
        end

        BSSurY(i,index)=dumy(indexF-1);
        index=index+1;
    end
    if i==1
       p(1) = stairs(dumx,dumy,'g');
    end
    stairs(dumx,dumy,'g')
    
end

p(2) = stairs(SAFunc{PlotVar,1}(:,1),(1-SAFunc{PlotVar,1}(:,2))*Conf(PlotVar)/3,'k','linewidth',2);

stairs(domf,prctile(BSSurY,(alpha/2)*100),'r','linewidth',2)
p(3) = stairs(domf,prctile(BSSurY,(1-(alpha/2))*100),'r','linewidth',2);


xlabel('Time [Years]')
ylabel('Probability')

legend([p(1) p(2) p(3)],{'Confidence Function of bootstraps','Confidence Function of total data','Confidence Bounds on 5% significance level'}) 
ax=gca;
ax.FontSize=12;
ylim([0 .25])
% saveas(gcf, 'confidence_func_v1.eps','epsc')
clear dumBNO dumF dumID dumRule dumSurT dumTID dumTime dumTimeH dumTimeL dumx dumy


