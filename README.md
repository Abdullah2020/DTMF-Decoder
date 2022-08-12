# DTMF-Decoder

The target of this task is to develop a DTMF decoder to decode the signals present in an audio file by returning the phone numbers pressed in the audio tones signals. The step by step approach adopted is explained below:

+ Read the audio file to see the frequency components
+ Window size of 65ms was created which generated 520 samples per window.
+ The power of the signal was calculated and converted to dB
+ The windows with signal power greater than 20dB were selected across the entire duration of the signal.
+ For windows with power greater than 20dB, for every window the highest frequency obtained via findpeak() function. This frequency corresponds to one of the high frequency of the DTMF frequency table.
+ A cheby1 band-stop filter, cheby1() was designed to remove this high frequency
+ The second highest frequency was obtained using findpeak() function. This frequency corresponds to one of the low frequency of the DTMF frequency table
+ The obtained frequencies (high and low after factoring the 1.5% frequency tolerance) was compared to the DTMF frequency table.
+ The corresponding numbers were returned.
