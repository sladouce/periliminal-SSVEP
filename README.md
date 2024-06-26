# periliminal-SSVEP
Dataset and code accompanying the paper: Frequency-tagging of spatial attention using periliminal flickers

The published manuscript can be found at the following address: [https://doi.org/10.1162/imag_a_00223]

  preprocess.m applies the EEG signal processing pipeline used in this study (see Methods section)


  RESS.m computes the SSVEP responses to both the target and non-target stimuli using the Rhythmic Entrainment Source Separation (RESS) Cohen & Gulbinaite 2017 analysis [https://doi.org/10.1016/j.neuroimage.2016.11.036].

The preprocessed datasets of the 24 subjects that took part in this experiment are epoched around the experimental events of the left and right cues (-4 to 7seconds) to capture the preceding 3 seconds of fixation cross, and the following 3 seconds cueing, and later target detection phases of the trials.
Processed data can be downloaded from the following osf repository [https://osf.io/rx2fs/?view_only=77bc7357608a44639ea5b0ec5e1fccd0]

Example video can be downloaded from this repo or directly streamed from: [https://www.youtube.com/watch?v=knQHjKJv14g]



