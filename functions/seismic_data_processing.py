#!/usr/bin/python
# -*- coding: UTF-8 -*-

#__modification time__ = 2024-02-23
#__author__ = Qi Zhou, Helmholtz Centre Potsdam - GFZ German Research Centre for Geosciences
#__find me__ = qi.zhou@gfz-potsdam.de, qi.zhou.geo@gmail.com, https://github.com/Nedasd
# Please do not distribute this code without the author's permission

import os
from obspy import read, Stream, read_inventory, signal
from obspy.core import UTCDateTime # default is UTC+0 time zone


def load_seismic_signal(seismic_network, station, component, data_start, data_end, remove_sensor_response=True):
    sac_path = f"/storage/vast-gfz-hpc-01/project/seismic_data_qi/seismic/EU/Illgraben/"

    d1 = UTCDateTime(data_start)
    d2 = UTCDateTime(data_end)

    file_dir = f"{sac_path}{d1.year}/{station}/{component}/"

    if d1.julday == d2.julday:
        data_name = f"{seismic_network}.{station}.{component}.{d1.year}.{str(d1.julday).zfill(3)}.mseed"
        st = read(file_dir + data_name)
    else:
        st = Stream()
        for n in np.arange(d1.julday-1, d2.julday+1):
            data_name = f"{seismic_network}.{station}.{component}.{d1.year}.{str(n).zfill(3)}.mseed"
            st += read(file_dir + data_name)

    st.merge(method=1, fill_value='latest', interpolation_samples=0)
    st._cleanup()
    st.detrend('linear')
    st.detrend('demean')

    if remove_sensor_response is True:
        meta_file = [f for f in os.listdir(directory) if f.startswith(seismic_network)][0]
        inv = read_inventory(f"{sac_path}meta_data/{meta_file}")
        st.remove_response(inventory=inv)

    st.filter("bandpass", freqmin=1, freqmax=45)
    st = st.trim(starttime=d1, endtime=d2, nearest_sample=False)
    st.detrend('linear')
    st.detrend('demean')

    return st
