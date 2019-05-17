# core

# core.hydromet_conv

## main
```python
main(binData:list, incr_excess:pandas.core.frame.DataFrame, tempE:float, convE:float, volE:float, tsthresh:float, display_print:bool=True) -> dict
```
Function for grouping incremental excess rainfall curves using a novel
test statistic that quantifies the incremental and cumulative
volumentric differences between two curves. The mean of each group of
curves is calculated and used in place of the original set of curves
in order to improve modeling efficiency by reducing redundancy.

