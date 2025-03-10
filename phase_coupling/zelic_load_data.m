% Assuming your table is already loaded into MATLAB and is named 'dataTable'
% If it's not loaded, you need to use readtable or a similar function to load your data.
% Set the path to the CSV file
gait_low  = '/Users/simonknogler/Desktop/MS1/Low/H/StepEventResults.csv';
gait_high = '/Users/simonknogler/Desktop/MS1/High/H/StepEventResults.csv';
ecg_low   = '/Users/simonknogler/Desktop/MS1/Low/H/ECG_RPeaks.csv';
ecg_high  = '/Users/simonknogler/Desktop/MS1/High/H/ECG_RPeaks.csv';

% Load the data into a MATLAB table
gaitEventsTable_low = readtable(gait_low);
gaitEventsTable_high = readtable(gait_high);
ecgEventsTable_low = readtable(ecg_low);
ecgEventsTable_high = readtable(ecg_high);


% Find rows where the EventType column contains 'TO'
rowsWithTO_low = contains(gaitEventsTable_low.EventType, 'TO');
rowsWithTO_high = contains(gaitEventsTable_high.EventType, 'TO');


% Delete these rows from the table
gaitEventsTable_low(rowsWithTO_low, :) = [];
gaitEventsTable_high(rowsWithTO_high, :) = [];


% Assuming gaitEventsTable_high is already loaded in the workspace

% Extract all RPeakTime for participant 31
x = ecgEventsTable_low.RPeakTime(ecgEventsTable_low.ParticipantID == 31);
y = gaitEventsTable_low.Time(gaitEventsTable_low.ParticipantID == 31);