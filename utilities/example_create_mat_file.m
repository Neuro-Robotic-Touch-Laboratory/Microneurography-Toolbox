% Settings
start_time = 55; % seconds

% Signal 1
fs1 = 10000; % Sampling frequency [Hz]
f1 = 0.01;   % Signal frequency [Hz]
duration1 = 100; % Let's say 100 seconds duration
t1 = start_time + (0:1/fs1:duration1-1/fs1); % timestamps
signal1 = sin(2*pi*f1*(t1 - start_time)); % remove start offset for pure sine

% Signal 2
fs2 = 5000;  % Sampling frequency [Hz]
f2 = 0.05;   % Signal frequency [Hz]
duration2 = 100;
t2 = start_time + (0:1/fs2:duration2-1/fs2);
signal2 = sin(2*pi*f2*(t2 - start_time));

% Signal 3
fs3 = 1000;  % Sampling frequency [Hz]
f3 = 0.07;   % Signal frequency [Hz]
duration3 = 100;
t3 = start_time + (0:1/fs3:duration3-1/fs3);
signal3 = sin(2*pi*f3*(t3 - start_time));

create_mat_file('testfile',...
                 'channel', signal1, t1 ,'MSNA', 'mV',...
                 'channel', signal2, [fs2,55],'ECG','mV',...
                 'channel', signal3, t3 ,'fingerpressure', 'mmHg',...
                 'comments', {'bla1','blubb2','blablubb3','blubbbla4','blabla5','blubblubb6','blablablub7','blublubbla8'},[55,60,65,70,75,80,85,90])