# -*- coding: utf-8 -*-


import numpy as np
from matfile_creator import create_mat_file

# --- Settings ---
start_time = 55  # seconds

# Signal 1
fs1 = 10000  # Sampling frequency [Hz]
f1 = 0.01    # Signal frequency [Hz]
duration1 = 100  # seconds
t1 = start_time + np.arange(0, duration1, 1/fs1)  # timestamps
signal1 = np.sin(2 * np.pi * f1 * (t1 - start_time))  # sine wave, no start offset

# Signal 2
fs2 = 5000   # Sampling frequency [Hz]
f2 = 0.05    # Signal frequency [Hz]
duration2 = 100
t2 = start_time + np.arange(0, duration2, 1/fs2)
signal2 = np.sin(2 * np.pi * f2 * (t2 - start_time))

# Signal 3
fs3 = 1000   # Sampling frequency [Hz]
f3 = 0.07    # Signal frequency [Hz]
duration3 = 100
t3 = start_time + np.arange(0, duration3, 1/fs3)
signal3 = np.sin(2 * np.pi * f3 * (t3 - start_time))

# Comments and timestamps
comments = ['bla1', 'blubb2', 'blablubb3', 'blubbbla4', 'blabla5', 'blubblubb6', 'blablablub7', 'blublubbla8']
comment_times = [55, 60, 65, 70, 75, 80, 85, 90]

# --- Call your mat file creation function ---
create_mat_file(
    'testfile',
    'channel', signal1, t1, 'MSNA', 'mV',
    'channel', signal2, [fs2, 55], 'ECG', 'mV',
    'channel', signal3, t3, 'fingerpressure', 'mmHg',
    'comments', comments, comment_times
)

