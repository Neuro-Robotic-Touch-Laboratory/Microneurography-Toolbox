import numpy as np
from scipy.io import savemat

def create_mat_file(filename, *args):

    temp_data = []
    temp_comm = []
    args = list(args)

    while args:
        if args[0] == 'channel':
            data = np.array(args[1])
            ts_fs = np.array(args[2])
            name = args[3]
            unit = args[4]

            # Ensure column vector
            data = np.atleast_2d(data)
            if data.shape[0] < data.shape[1]:
                data = data.T

            entry = {}
            entry['data'] = data
            entry['n_samples'] = len(data)

            if len(ts_fs) == 2:
                sampling_freq, start_time = ts_fs
                entry['ts'] = np.arange(len(data)) / sampling_freq + start_time
            else:
                ts_fs = np.atleast_2d(ts_fs)
                if ts_fs.shape[0] < ts_fs.shape[1]:
                    ts_fs = ts_fs.T
                entry['ts'] = ts_fs.flatten()

            entry['dt'] = round(np.mean(np.diff(entry['ts'])), 6)
            entry['name'] = name
            entry['units'] = unit

            temp_data.append(entry)
            args = args[5:]

        elif args[0] == 'comments':
            comms = args[1]
            ts = args[2]
            for i in range(len(comms)):
                temp_comm.append({'com_str': comms[i], 'ts': ts[i]})
            args = args[3:]

        else:
            args = args[1:]

    # Sort by dt
    idx = sorted(range(len(temp_data)), key=lambda i: temp_data[i]['dt'])
    
    len_n = max(len(td['name']) for td in temp_data)
    len_u = max(len(td['units']) for td in temp_data)

    data = np.array([], dtype=float)
    datastart = []
    dataend = []
    titles = []
    unittext = []
    unittextmap = []
    samplerate = []

    tickrate = 1 / temp_data[idx[0]]['dt']
    blocktimes = temp_data[idx[0]]['ts'][0] * tickrate

    for i in range(len(idx)):
        chan = temp_data[idx[i]]
        datastart.append(len(data)+1)
        data = np.hstack((data, chan['data'].flatten()))
        dataend.append(len(data))

        # Fill titles and units padded
        titles.append(chan['name'].ljust(len_n))
        unittext.append(chan['units'].ljust(len_u))
        unittextmap.append(i+1)
        samplerate.append(1 / chan['dt'])

    # Prepare comments
    com = []
    if temp_comm:
        maxlen = max(len(c['com_str']) for c in temp_comm)
        comtext = [''.ljust(maxlen) for _ in temp_comm]
        for i, comm in enumerate(temp_comm):
            pos = np.searchsorted(temp_data[idx[0]]['ts'], comm['ts'], side='left')
            com.append([-1, 1, pos+1, 1, i+1])
            comtext[i] = comm['com_str'].ljust(maxlen)
    else:
        comtext = []

    # Assemble final dict to save
    mat_struct = {
        'data': data,
        'datastart': np.array(datastart).reshape(-1, 1),  # Ensure column vector
        'dataend': np.array(dataend).reshape(-1, 1),      # Ensure column vector
        'titles': np.array(titles, dtype='U'),
        'unittext': np.array(unittext, dtype='U'),
        'unittextmap': np.array(unittextmap).reshape(-1, 1),  # Ensure column vector
        'tickrate': tickrate,
        'blocktimes': blocktimes,
        'samplerate': np.array(samplerate).reshape(-1, 1),  # Ensure column vector
        'com': np.array(com) if com else np.zeros((0, 5)),
        'comtext': np.array(comtext, dtype='U') if comtext else np.array([], dtype='U'),
    }

    # Ensure the filename includes the .mat extension
    if not filename.endswith('.mat'):
        filename += '.mat'

    # Save the .mat file
    savemat(filename, mat_struct)

    # --- Print information like MATLAB ---
    print(f'\n✅ {filename} was successfully created')
    print(f'Start time: {blocktimes / tickrate:.6f} s, data length: {dataend[0] / tickrate:.6f} s')
    print(f'The file contains {len(samplerate)} channels:')

    for i in range(len(samplerate)):
        print(f"Channel {i+1} name: {titles[i]}, units: {unittext[i]}, "
              f"length: {dataend[i] - datastart[i] + 1} samples, sampling frequency: {samplerate[i]:.2f} Hz")

    if len(comtext) > 0:
        print(f'\nThe file contains {len(comtext)} comments:')
        for i, text in enumerate(comtext):
            ts = (blocktimes + com[i][2] - 1) / tickrate
            print(f'Comment {i+1}: {text.strip()}, at: {ts:.6f} s')
    else:
        print('\nNo comments.')