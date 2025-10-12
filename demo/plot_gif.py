#!/usr/bin/python
# -*- coding: UTF-8 -*-

# __modification time__ = 2025-01-20
# __author__ = Qi Zhou and Sibashish Dash, GFZ Helmholtz Centre for Geosciences
# __find me__ = qi.zhou@gfz.de, qi.zhou.geo@gmail.com, https://github.com/Qi-Zhou-Geo
# Please do not distribute this code without the author's permission

import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import matplotlib.ticker as ticker

from tqdm import tqdm

from obspy import Stream, Trace, read
from obspy.core import UTCDateTime # default is UTC+0 time zone


# <editor-fold desc="add the sys.path to search for custom modules">
from pathlib import Path
current_dir = Path(__file__).resolve().parent
# using ".parent" on a "pathlib.Path" object moves one level up the directory hierarchy
project_root = current_dir.parent
import sys
sys.path.append(str(project_root))
# </editor-fold>


# import the custom functions
from functions.Type_A_features import calBL_feature
from functions.welch_spectrum import welch_psd


plt.rcParams.update( {'font.size':7,
                      'axes.formatter.limits': (-3, 6),
                      'axes.formatter.use_mathtext': True} )

# prepare the feature
st = read(f"{project_root}/data/Illgraben_9J_IGB02_HHZ_2014-07-12T00:00:00_2014-07-13T00:00:00.mseed")
amp = st[0].data

sampling_rate = st[0].stats.sampling_rate
sub_window_size = 60 # unit is second
scaling = 1e9
num_window = int((st[0].stats.endtime - st[0].stats.starttime)/sub_window_size)

BL_feature = np.empty((num_window, 19), dtype=object) # time_float, time_str, 17 BL features
freq_psd = np.empty((num_window, 491), dtype=object)

for i in tqdm(range(num_window)):

    t = st[0].stats.starttime + i * sub_window_size
    data_array = amp[int(i * sampling_rate * sub_window_size): int((i+1) * sampling_rate * sub_window_size)]
    data_array_nm = data_array * scaling # converty m/s to nm/s
    output = calBL_feature(data=data_array_nm, ruler=300, epsilon=1e-8)

    temp = [t.strftime("%Y-%m-%dT%H:%M:%S")] + [float(t)] + output.tolist()

    BL_feature[i, :] = temp

    freq, psd, psd_unit = welch_psd(data=data_array, sampling_freq=sampling_rate,
                                    f_min=1, f_max=50, segment_window=10,
                                    scaling="density", unit_dB=True)
    freq_psd[i, :] = psd


# prepare data
selecte_time = ['2014-07-12T13:00:00', '2014-07-12T21:00:00']
x_interval = 2 # hours
id1 = np.where(BL_feature[:, 0] == selecte_time[0])[0][0]
id2 = np.where(BL_feature[:, 0] == selecte_time[1])[0][0] + 1
feature = BL_feature[id1:id2, :]

amp =  st.trim(UTCDateTime(selecte_time[0]), UTCDateTime(selecte_time[1]))[0].data[:-1] # sps by 12000 data per second
amp_t = np.arange(len(amp))
digit_frequency = feature[:, 2:11].astype(float) # sps by 1 data per minute
goodness = feature[:, 12].astype(float) # sps by 1 data per minute
ks = feature[:, 16].astype(float) # sps by 1 data per minute
MannWhitneyU =  feature[:, 17].astype(float) # sps by 1 data per minute
time_index = np.arange(len(goodness))
sps_data = 1/60 # Hz



# <editor-fold desc="hide the following codes">
fig = plt.figure(figsize=(6, 5.5))
gs = gridspec.GridSpec(3, 2, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])  # top-left
ax2 = fig.add_subplot(gs[1, 0])  # middle-left
ax3 = fig.add_subplot(gs[2, 0])  # bottom-left
ax4 = fig.add_subplot(gs[:, 1])  # entire right column
axes = [ax1, ax2, ax3, ax4]


amp_line, = ax1.plot([], [], color="black", lw=1, zorder=2)
ax1.set_ylabel("Amplitude [m/s]", fontweight="bold")
ax1.grid(axis='y', ls="--", lw=0.5, zorder=1)
ax1.set_xlim(0, amp_t[-1])
ax1.set_ylim(-2.5e-4, 2.5e-4)
ax1.axes.xaxis.set_ticklabels([])
ax1.xaxis.set_major_locator(ticker.MultipleLocator(3600 *  x_interval * sampling_rate))


ks_scatter = ax2.scatter([], [], color="#3C405B", s=20, alpha=0.5, zorder=2, label="Kolmogorov–Smirnov")
mwu_scatter = ax2.scatter([], [], color="#82B29A", s=20, alpha=0.5, zorder=2, label="Mann–Whitney U")
ax2.axhline(y=0.95, ls="--", lw=1, color="red", alpha=0.7, label="95% Confidence Level", zorder=1)
ax2.legend(loc="upper right", fontsize=6)
ax2.set_xlim(0, time_index[-1])
ax2.set_ylim(-0.05, 1.05)
ax2.set_ylabel("p-value", fontweight="bold")
ax2.grid(axis='y', ls="--", lw=0.5, zorder=1)
ax2.axes.xaxis.set_ticklabels([])
ax2.xaxis.set_major_locator(ticker.MultipleLocator(60 * x_interval))


goodness_scatter = ax3.scatter([], [], color="black", s=20, alpha=0.6, zorder=2)
ax3.set_xlim(0, time_index[-1])
ax3.set_ylim(-200, 100)
ax3.set_ylabel("Goodness of Fit [%]", fontweight="bold")
ax3.grid(axis='y', ls="--", lw=0.5, zorder=1)
ax3.axes.xaxis.set_ticklabels([])
ax3.xaxis.set_major_locator(ticker.MultipleLocator(60 * x_interval))
duration = int((UTCDateTime(selecte_time[1]) - UTCDateTime(selecte_time[0])) / 3600)
xLocation = np.arange(0, sps_data * 3600 * (duration + x_interval), sps_data * 3600 * x_interval)
xTicks = []
for idx, i in enumerate(xLocation):
    if idx == 0:
        xTicks.append((UTCDateTime(selecte_time[0]) + i * 1 / sps_data).strftime('%Y-%m-%d' + '\n' + '%H:%M:%S'))
    else:
        xTicks.append((UTCDateTime(selecte_time[0]) + i * 1 / sps_data).strftime('%H:%M:%S'))

ax3.set_xticks(xLocation, xTicks)
ax3.set_xlabel(f"Time [UTC+0]", fontweight='bold')


data0 = [0.301, 0.176, 0.125, 0.097, 0.079, 0.067, 0.058, 0.051, 0.046]
ax4.plot(data0, color="#EA6B66", marker="o", ls="-", lw=1, markersize=4.5, label="Theoretical value", zorder=3)
digit_frequency_line, = ax4.plot([], color="black", marker="o", ls="-", lw=1, markersize=4.5, label="Observed value", zorder=3)
text_plot = ax4.text(x=0, y=1.0, s="", fontsize=7, fontweight="bold", horizontalalignment="left")

ax4.legend(loc="upper right", fontsize=6)
ax4.set_xticks(np.arange(0, 9), np.arange(1, 10))
ax4.set_xlim(-0.5, 8.5)
ax4.set_ylim(-0.05, 1.05)
ax4.grid(axis='y', ls="--", lw=0.5, zorder=1)
ax4.grid(axis='x', ls="--", lw=0.5, zorder=1)
ax4.set_ylabel("Frequency", fontweight="bold")
ax4.set_xlabel("First-Digit", fontweight="bold")
ax4.yaxis.tick_right()
ax4.yaxis.set_label_position("right")


plt.tight_layout()


def init():
    amp_line.set_data([], [])
    ks_scatter.set_offsets(np.empty((0, 2)))
    mwu_scatter.set_offsets(np.empty((0, 2)))
    goodness_scatter.set_offsets(np.empty((0, 2)))
    digit_frequency_line.set_data([], [])
    text_plot.set_text("")

    # FIXED: Added text_plot to return
    return [amp_line, ks_scatter, mwu_scatter, goodness_scatter, digit_frequency_line, text_plot]


def update(frame):
    # Skip frame 0 to avoid empty data issues
    if frame == 0:
        return [amp_line, ks_scatter, mwu_scatter, goodness_scatter, digit_frequency_line, text_plot]

    # Amplitude: 12,000 samples per minute
    amp_samples = frame * 12000
    amp_line.set_data(amp_t[:amp_samples], amp[:amp_samples])

    # Scatter plots
    ks_scatter.set_offsets(np.column_stack((time_index[:frame], ks[:frame])))
    mwu_scatter.set_offsets(np.column_stack((time_index[:frame], MannWhitneyU[:frame])))
    goodness_scatter.set_offsets(np.column_stack((time_index[:frame], goodness[:frame])))

    # First-digit frequency line
    digit = np.arange(0, 9)
    digit_frequency_line.set_data(digit, digit_frequency[frame - 1, :])
    text_plot.set_text(f" Goodness of Fit:\n {goodness[frame - 1]:.2f}%")  # FIXED: frame-1

    return [amp_line, ks_scatter, mwu_scatter, goodness_scatter, digit_frequency_line, text_plot]


ani = animation.FuncAnimation(
    fig,
    update,
    frames=len(time_index),
    init_func=init,
    interval=50, # delay between frames (milliseconds)
    blit=False # only redraw changed parts for faster rendering
)

ani.save("BL_2014-07-12.gif", writer="pillow", fps=20, dpi=300)

# </editor-fold>
