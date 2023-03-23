# MPP-EEG_Ver2
Matlab code that decomposes EEG, ECoG or LFP data exploiting a Marked Point Process (MPP) framework. Improvements in implementation includes a robust Correntropy based denoising approach.  

For training, call: 

```
[D,MPP,th,ar,bw] = PhEv_Learn_fast_2(X, M, K); 
```

For testing, call:  

```
th = GetThreshold(X,M,ar,bw); # to adjust threshold on unseen data
MPP = Decomp_EEG(X,D,th,1); 
```
