% Author: Chris Fietkiewicz
% Reference: Christopher Fietkiewicz and Kenneth A. Loparo, “Analysis and Enhancements of a Prolific Macroscopic Model of Epilepsy,” Scientifica, vol. 2016, Article ID 3628247, 10 pages, 2016. doi:10.1155/2016/3628247
% Description: Original model used for Figure 2 in Fietkiewicz and Loparo 2016. It is based on the original model specified in Wendling F, % Bartolomei F, Bellanger JJ, Chauvel P. Epileptic fast activity can be
% explained by a model of impaired GABAergic dendritic inhibition. Eur J Neurosci, 2002;15(9):1499-1508.
function OriginalModel
stepSize = 0.001;
totalDuration = 13 + 1 + 10 + 1 + 2 + 1 + 9 + 1 + 11;
totalSteps = totalDuration / stepSize;
duration = [13 10 2 9 11];
bVals = [45 38 37 8 11];
gVals = [20 20 20 20 2];

A = 5 * ones(totalSteps, 1);
B = [];
G = [];
for i = 1:5
    steps = duration(i) / stepSize;
    bSegment = bVals(i) * ones(steps, 1);
    B = [B; bSegment];
    gSegment = gVals(i) * ones(steps, 1);
    G = [G; gSegment];
    if (i < 5)
        bStep = (bVals(i + 1) - bVals(i)) * stepSize; % Transitions are all 1 second
        bTransition = bVals(i) : bStep : bVals(i + 1);
        if (i < 4)
            gTransition = gVals(i) * ones(1, 1000);
        else
            gStep = (gVals(i + 1) - gVals(i)) * stepSize; % Transitions are all 1 second
            gTransition = gVals(i) : gStep : gVals(i + 1);
        end
        B = [B; bTransition(1:1000)'];
        G = [G; gTransition(1:1000)'];
    end
end

eeg = simulate(0.001, totalDuration, A, B, G, 1, [0 0 0 0 0 0 0 0 0 0], 0);
t = stepSize * [0:length(eeg)-1];
plot(t, eeg)
xlabel('Time (sec)')
ylabel('y-out (mV)')
title('Original Model')

% Main simulation
function eeg = simulate(stepSize, endTime, A, B, G, transientTime, yInit, plotFlag)
steps = int32(endTime / stepSize);
transientSteps = int32(transientTime / stepSize);
eeg = zeros(steps, 1);
y = yInit;

for i = 1:transientSteps
	y = takeStep(y, stepSize, A(1), B(1), G(1));
end

for i = 1:steps
	y = takeStep(y, stepSize, A(i), B(i), G(i));
	eeg(i) = y(2) - y(3) - y(4);
end

if (plotFlag == 1)
	plot(eeg)
	set(gcf, 'Position',[0 260 790 260])
end

% ***************************************************************************************************************
% Take a step
function yNew = takeStep(y, stepSize, A, B, G)

a = 100;
b = 50;
g = 350;
C = 135;
C1 = C;
C2 = 0.8 * C;
C3 = 0.25 * C;
C4 = 0.25 * C;
C5 = 0.3 * C;
C6 = 0.1 * C;
C7 = 0.8 * C;
MEAN = 90;
SIGMA = 30;

yNew(1) = y(1) + y(6) * stepSize;
yNew(6) = y(6) + (A * a * S(y(2)-y(3)-y(4)) - 2 * a * y(6) - a * a * y(1)) * stepSize;
yNew(2) = y(2) + y(7) * stepSize;
yNew(7) = y(7) + (A * a * (randn * SIGMA + MEAN + C2 * S(C1* y(1))) - 2 * a * y(7) - a * a * y(2)) * stepSize;
yNew(3) = y(3) + y(8) * stepSize;
yNew(8) = y(8) + (B * b * (C4 * S(C3 * y(1)) ) - 2 * b * y(8) - b * b * y(3)) * stepSize;
yNew(4) = y(4) + y(9) * stepSize;
yNew(9) = y(9) + (G * g * (C7 * S((C5 * y(1) - C6 * y(5))) ) - 2 * g * y(9) - g * g * y(4)) * stepSize;
yNew(5) = y(5) + y(10) * stepSize;
yNew(10) = y(10) + (B * b * ( S(C3 * y(1)) ) - 2 * b * y(10) - b * b * y(5)) * stepSize;

% ***************************************************************************************************************
% Sigmoid function
function r = S(v)
r = 5.0 / (1.0 + exp(0.56 * (6.0 - v)));