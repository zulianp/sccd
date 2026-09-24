# TO BE FIXED

- AxisSum3 class is unnecessary
- Make Cell2DGrid a data only struct, and the instance method should lean inline functions taking all arguments outside of the class
- sccd_broadphase_celld2d.hpp containts code that is too object oriented (and clever c++), make it leaner and cleaner, avoid modernisms, inline things where possible and try to promote ILP and DLP
- Aabbs -> AABBs
- In the paper for the broadphase plots and tables make sure that is always clear what phases are included in the measurement
    - Tag: meaning
	-  BP full: prep + queries
	-  BP prep: prep only
	-  BP queries: queries only
  also make sure that is clear in the captions everywhere for every method that we mention. Use NP for narrow phase

# VERIFY

- Use compiler option to see if loops that are meant to be vectorized actually are