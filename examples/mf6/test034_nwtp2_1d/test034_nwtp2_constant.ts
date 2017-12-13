# The ATTRIBUTES block is required.
# Number of names on NAME line indicates 
#     number of time series to be created.
# NAME must precede METHOD (or METHODS).
# Either METHOD or METHODS must be specified, but not both.
# If METHOD is specified, all time series in file
#     share the same method (either LINEAR or STEPWISE).
# IF METHODS is specified, a method is specified for each time series.
#
BEGIN ATTRIBUTES
  NAME    L01HEAD  L10HEAD
  METHODS stepwise stepwise
END ATTRIBUTES

BEGIN TIMESERIES
  0.000000000000  37.63686464  28.44462746
  19.00000000000  37.63686464  28.44462746
  38.00000000000  39.51799627  30.55890831
  57.00000000000  41.62772905  32.73199899
  76.00000000000  42.71697008  33.97117642
  95.00000000000  43.90046995  35.13032014
  114.0000000000  44.90059323  36.06101855
  133.0000000000  45.75376623  36.87839325
  152.0000000000  46.27354747  37.49465934
  171.0000000000  46.70308883  38.01006194
  190.0000000000  47.26094677  38.57244534
  449.0000000000  50.60481107  41.97585373
  708.0000000000  52.50088247  44.01869733
  821.0000000000  53.20954030  44.75986558
  934.0000000000  53.78208746  45.36381486
  1047.000000000  54.22216834  45.87527002
  1160.000000000  54.60210681  46.32624919
  1273.000000000  54.96386883  46.73354163
  1386.000000000  55.35611763  47.13251346
  1499.000000000  55.66545298  47.46545300
  1612.000000000  55.92173546  47.74877509
  1725.000000000  56.16184531  48.01652110
  1838.000000000  56.38765322  48.26771419
  1951.000000000  56.59626192  48.49846716
  2064.000000000  56.81870870  48.74058283
  2177.000000000  57.02131198  48.96065078
  2290.000000000  57.20695055  49.16275529
  2403.000000000  57.38191508  49.35396568
  2516.000000000  57.53877477  49.52423535
  2629.000000000  57.69639210  49.69799744
  2630.000000000  63.56956148  56.34615365
END TIMESERIES
