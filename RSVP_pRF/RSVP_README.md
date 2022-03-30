Overview:
This version is a branch of iEye that has been configured to work with eye tracking collected with the RSVP pRF program.

Using this software to process RSVP data:
Everything specific to RSVP data will be located under the RSVP_pRF/ directory. When first starting, be sure to include a copy of rsvp_params.txt (the parameters file that configures how the RSVP program will be run), as the software will use these details to configure the data and resulting plots. 
To add a new subject to run, add a new directory under RSVP_pRF/Subjects/, and place all EDFs in that directory. Then open rsvp_prf_analysis.m, and add a line that will direct iEye where to look for the new subject:  edf_prefix = 'example_sub_1/';
Configure your analysis parameters as you normally would with iEye, and that's it!
Currently the only crietria used to flag TRs is presence of a saccade, but this can be expanded.
For all other specifications, default to general iEye instructions.


Additions to iEye for RSVP:
- rsvp_prf_loadparams.m
	- Load parameters from rsvp_params.txt to properly configure analysis
- rsvp_prf_proc.m
	- Process RSVP parameters, select TRs to flag for further inspection
- rsvp_prf_plotQC_exclusions.m
	- Plot matrix of flagged TRs to evaluate parameter changes and identify patterns

Major edits of iEye functions ("ii_" replaced with "rsvp_prf_" unless noted):
- rsvp_prf_driftcorrect.m
	- Added a drift correction mode that does not rely on ITI fixation, as some subjects may purposefully blink during ITI. New mode uses average X,Y fixation of other ITI fixations when none is available for a particular ITI
- rsvp_prf_plottimeseries.m
	- Added TR structure, trace of X,Y location of stimuli on screen

Minor edits of iEye function ("ii_" replaced with "rsvp_prf_" unless noted):
- rsvp_prf_collectTrials (original = ii_scoreMGS.m)
	- Collects and combines trials
- rsvp_prf_combineruns.m
	- Changed path to function
- rsvp_prf_import_edf.m
	- Added checks to make sure eye tracker messages have printed when they should (this slows the software and is not required in most cases)
- rsvp_prf_calibratebytrial.m
	- Prevents bad calibrations from clearing eye traces
- rsvp_prf_prerpoc.m
	- Calls scripts unique to RSVP, calls ii_extractsaccades.m earlier in stream


General troubleshooting:
If you experience an issue using edf2asc to translate files, you may need to tell matlab the location of your edf2asc function. You can find this in ii_import_edf.m (should be on lines 99 & 125). To find your edf2asc, on your bash terminal, type:
which edf2asc
Then on the first call (should be line 99) for example, change:
[status,result] = system(['edf2asc -t -c -s -miss 0 ' edf_file]);
to the location of your edf2asc:
[status,result] = system(['/usr/local/bin/edf2asc -t -c -s -miss 0 ' edf_file]);
and the same for the second call to the function.
