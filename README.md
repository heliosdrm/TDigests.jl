# TDigests

t-digest data structure in Julia

*Note*:

This is an attempt to implement t-digests in order to calculate extreme quantiles of streamed data efficiently, without having to store large histograms. However, this implementation does not seem to be faster or consume less memory than the `quantile` function of Julia's stdlib. So it is not currently maintained.

(Anyone who feels like taking the challeng is invited to take over the code and improve it.) 
