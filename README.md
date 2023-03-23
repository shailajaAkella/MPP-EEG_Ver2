# MPP-EEG_Ver2
Matlab code that decomposes EEG, ECoG or LFP data exploiting a Marked Point Process (MPP) framework. Improvements in implementation includes a robust Correntropy based denoising approach.  

Running instructions: 

Inputs
X - banpassed, multi - trial, single channel EEG signal in cell format or row format. (Dimension: number of trials X Time)
M - Maximum length of event in samples. For ex., sleep spindles are approximately 0.5s in length. In that case, M = 0.5 X sampling frequency
K - Total number of dictionary atoms/ templates/ clusters of events.

Outputs
MPP - structure comprising of trial wise detections (as a marked point process (MPP)): 
       PhEv - normalized snippet of the exracted event,
       tau - time point of occurrence (mid - point),
       alph - maximum absolute amplitude of the event,
       D_idx - dictionary index that maximally correlates
       pow - power of the extracted snippet
D - Dictionary atoms/ centers of the cluster/ templates of the events
th - minimum norm of burst in training data
ar, bw - properties of event norm distribition.  

For training, call: 

```
[D,MPP,th,ar,bw] = PhEv_Learn_fast_2(X, M, K); 
```

For testing, call:  

```
th = GetThreshold(X,M,ar,bw); # to adjust threshold on unseen data
MPP = Decomp_EEG(X,D,th,1);  
```
