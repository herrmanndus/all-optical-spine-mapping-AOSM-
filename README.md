# all-optical-spine-mapping-AOSM-

Dustin Herrmann 2021/2022
Note: 2022-09-16, this repository is under construction. Dependencies may be incomplete, please contact me if any issues arise. 

This collection of scripts allows reading in 2-photon FOV and widefield retinotopic mapping data to semi-automatically select suitable target locations to be used for all-optical connectivity mapping. Generates SLM phase masks, galvo positions (.gpl) and position list (.xml) files to be read by Bruker's Prairie View software. After stimulating target grid locations, the function "" 

Note: SLM Phase mask maker script, .gpl maker and .xml maker and Trigger Builder are based on Lloyd Russell's NAPARM script (https://github.com/llerussell/Naparm) and were modified by Dustin Herrmann

Instructions: 
- Start with AOSM_grid_maker.m --> read 1-photon and 2-photon imaging data. line up and generate transform between the two. select soma location from a large 2-photon FOV. semi-automatically select area to target for holographic optogenetic stimulation. make target grid to probe these locations. 
- run_make_exptFiles --> from the grid positions generate SLM phase masks and galvo position lists to interface with Bruker's Prairie view software. Also automatically computes command to send to laser for power control. Note: reads in settings_github.ylm file. This script is based on https://github.com/llerussell/Naparm 
- online_spine_detection --> load in extracted spine fluorescence traces and quickly find candidates of target location - responding spine pairs 
