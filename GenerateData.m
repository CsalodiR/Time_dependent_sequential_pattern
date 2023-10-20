clear all
close all
clc


% Data is synthetized with Markov models

NoEvents = 100; % No of unique events
SeqEndProb = 0.2; % Probability of ending the sequence
StampEndProb = 1; % Probability of event is handled together
NoInPuts = 15000; % No of samples


StateTransition = rand(NoEvents, NoEvents);
StateTransition = StateTransition ./ sum(StateTransition, 2); % Making state transition matrix

InPutSequences = cell(NoInPuts,1);

for i = 1:NoInPuts % Generate events with roulette wheel selection
    actseq = [];
    firstEvent = floor(rand(1) * NoEvents) + 1;
    actseq = [actseq, firstEvent];
    actseq = [actseq, -1];
    
    makeMore = true;
    j = 3;
    
    while makeMore
        
        Probs = StateTransition(:, actseq(j-2));
        Probs = [0; Probs];
        Probs(end) = [];
        
        Probs = cumsum(Probs);
        
        moreEventA = rand(1);
        moreEvent = find(Probs < moreEventA);
        moreEvent = moreEvent(end);
        if ismember(moreEvent,actseq)
            continue
        end
        actseq = [actseq, moreEvent];
        
        
        
        if rand(1) < SeqEndProb
            actseq = [actseq, -2];
            makeMore = false;
        else
            actseq = [actseq, -1];
        end
        
        j = j + 2;
    end
    
    
    InPutSequences{i,1} = actseq;
end

% Parameters of event time distribution (exponential)
% Timestamps are provided in datenum format
% Initial date of sequence is between 2005 (732313) and 2008 (733408)

initL = 732313;
initU = 733408;

% Parameters of exponential distribution of elapsed times
ExponentialParameters = randi([50, 100], 1, NoEvents);


InPutTimeStamps = cell(NoInPuts,1);

for i = 1:NoInPuts
    dumSeq = InPutSequences{i};
    times = [];
    
    startDate = randi([initL, initU], 1);
    times = [times, startDate];
    actDate = startDate;
    
    for j = 2:length(dumSeq)
        if dumSeq(j) == -1
            times = [times, -1];
            actDate = actDate + ceil(random('Exponential', ExponentialParameters(dumSeq(j-1)), 1, 1));
        elseif dumSeq(j) == -2
            times = [times, -2];
        else
            times = [times, actDate];
        end
    end
    
    
    InPutTimeStamps{i,1} = times;
end

clear actDate actseq dumSeq ExponentialParameters firstEvent
clear i initL initU j makeMore moreEvent moreEventA NoEvents NoInPuts
clear Probs SeqEndProb StampEndProb startDate StateTransition times

save InPutData