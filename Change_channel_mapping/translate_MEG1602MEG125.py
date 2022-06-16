#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 12:11:15 2022

@author: mq20096022
"""

import scipy.io
import mne
import numpy as np




data_dir = "E:\Downloads\Change_channel_mapping\\"
# load up the mapping lists - stored as .mat files
MEG125 = scipy.io.loadmat(
    data_dir + "labels_systemmatch_child_sensors.mat", simplify_cells=True
)
MEG160 = scipy.io.loadmat(
    data_dir + "labels_systemmatch_adult_sensors.mat", simplify_cells=True
)

MEG125_picks = MEG125["lab125"]
MEG160_picks = MEG160["lab160"]

for idx, ch in enumerate(MEG160_picks):
    MEG160_picks[idx] = ch.replace("AG", "MEG ")

for idx, ch in enumerate(MEG125_picks):
    MEG125_picks[idx] = ch.replace("AG", "MEG ")


# load up the MEG160 data we want to change channels for
MEG160_raw = mne.io.read_raw_kit(
    data_dir + "MEG160_sample.con",  # change depending on file i want
    stim=[*range(182, 190)],  # Just the trigger chans
    slope="+",
    stim_code="channel",
    stimthresh=1,  # 2 for adults
    preload=True,
    allow_unknown_format=False,
    verbose=True,
)

#%% Finding events
events = mne.find_events(
    MEG160_raw,
    output="onset",
    consecutive=False,
    min_duration=0,
    shortest_event=1,  # 5 for adults
    mask=None,
    uint_cast=False,
    mask_type="and",
    initial_event=False,
    verbose=None,
)

MEG160_epochs = mne.Epochs(MEG160_raw, events, tmin=-0.1, tmax=0.4, preload=True)
MEG160_epochs.pick_types("mag")

MEG160_epochs_trans = MEG160_epochs.copy().pick_channels(MEG160_picks)
MEG160_epochs_trans.reorder_channels(MEG160_picks)

MEG125_raw = mne.io.read_raw_kit(
    data_dir + "MEG125_sample.con",
    allow_unknown_format=False,
    verbose=True,
)

MEG125_raw.pick_channels(MEG125_picks)

MEG160_epochs_trans.info["chs"] = MEG125_raw.info["chs"]
MEG160_epochs_trans.info["ch_names"] = MEG125_raw.info["ch_names"]

mne.viz.plot_evoked_topomap(
    MEG160_epochs_trans.average(), times=np.arange(-0.1, 0.41, 0.05)
)
mne.viz.plot_evoked_topomap(MEG160_epochs.average(), times=np.arange(-0.1, 0.41, 0.05))

mne.viz.plot_sensors(
    MEG160_epochs.info,
    kind="topomap",
    ch_type=None,
    title=None,
    show_names=False,
    ch_groups=None,
    to_sphere=True,
    axes=None,
    block=False,
    show=True,
    sphere=None,
    #pointsize=None,
    #linewidth=2,
    verbose=None,
)

mne.viz.plot_sensors(
    MEG160_epochs_trans.info,
    kind="topomap",
    ch_type=None,
    title=None,
    show_names=False,
    ch_groups=None,
    to_sphere=True,
    axes=None,
    block=False,
    show=True,
    sphere=None,
    #pointsize=None,
    #linewidth=2,
    verbose=None,
)
