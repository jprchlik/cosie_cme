# Code Suite for CME Analysis

This directory contains the code used to make plots used in the COSIE proposal and analyze data for the proposal.

cosie_sun.py
------------
This script takes the CACTUS HALO CME parameters from 2011 and matches them up with a ISS orbit in the same year.
Assuming COSIE was on the ISS that year and had a 30 cadence, the output plot (COSE_CME_plot.png) shows at which points in space
COSIE would have observed the CME.

cosie_sun_rad_sc.py
-------------------
This script takes and analyzes 11 years of CATCUS HALO CME parameters (2004-10-31 to 2015-10-31).
The CME start times are varied with a uniform random distribution in the year 2011,
so that the ISS orbit for that year may be used to statistically understand the observing probabilities of COSIE.
This script creates the file simulated_cosie_cme.csv, which contains the observational infromation for each simulated CME.

cosie_sun_rad_sc_static_time.py
-------------------------------
This script takes and analyzes 11 years of CATCUS HALO CME parameters (2004-10-31 to 2015-10-31).
The CME start times are varied with a uniform random distribution in a year with a fix observing window for the ISS (30 mintues).
This script creates the file simulated_cosie_cme.csv, which contains the observational infromation for each simulated CME.

fancy_plot.py
-------------
A small module which formats plots in a way I find appealing.

plot_simulated_cmes.py
----------------------
This script creates plots that characterize the simulated observations from the cosie_sun_rad_sc.py script.
The three plots created by this program are as follows:

cactus_durat_obs_vs_vel_reverse_2Dhist.png
=========
A three panel plot showing the detection fraction a function of speed, the number of CMEs observed for a given duration as a function of speed, and the cumulative fraction of observed CME durations.

cactus_cme_two_year.png
========
The number of CMEs observed as a function of Speed for a 2 year long mission, assuming an average year.

cactus_cme_orbits_obs.png
=========================
The fraction of observed CMEs that are visible in a single orbit.

plot_simulated_cmes_static_time.py
----------------------
This script creates plots that characterize the simulated observations from the cosie_sun_rad_sc_static_time.py script. The last part of the file name signifies the assumed COSIE observing window in minutes.

cactus_durat_obs_vs_vel_reverse_2Dhist_XX.png
=================================================

A two panel plot showing the detection fraction a function of speed and number of CMEs observed per year and the observed CME duration as a function of speed. This plot was used in the proposal.