#include "jupiter.h"


inline real TVDslope (slope1, slope2)
     real slope1, slope2;
{
  real slopec;
  slopec = .5*(slope1+slope2);	/* slope1 on the right side,  */
  if (slope1 * slope2 <= 0.0)	/* slope2 on the left side */
    slopec = 0.0;
  else {
    if (slopec < 0.0) {
      if (slopec < slope1) /* check no overshoot condition on right side */
	slopec = slope1;
      if (slopec < slope2) /* check no overshoot condition on left side */
	slopec = slope2;
    }
    if (slopec > 0.0) {
      if (slopec > slope1) /* check no overshoot condition on right side */
	slopec = slope1;
      if (slopec > slope2) /* check no overshoot condition on left side */
	slopec = slope2;
    }
  }
  return slopec;
}


inline real minmod (slope1, slope2)
     real slope1, slope2;
{
  real slopec;
  slopec = .5*(slope1+slope2);	/* slope1 on the right side,  */
  if (slope1 * slope2 <= 0.0)	/* slope2 on the left side */
    slopec = 0.0;
  else {
    if (slopec < 0.0) {
      if (slopec < slope1) /* check no overshoot condition on right side */
	slopec = slope1;
      if (slopec < slope2) /* check no overshoot condition on left side */
	slopec = slope2;
    }
    if (slopec > 0.0) {
      if (slopec > slope1) /* check no overshoot condition on right side */
	slopec = slope1;
      if (slopec > slope2) /* check no overshoot condition on left side */
	slopec = slope2;
    }
  }
  return slopec;
}


inline real superbee (slope1, slope2)
     real slope1, slope2;
{
  real slopec, sigma1, sigma2;

  sigma1 = minmod(slope1,2.0 * slope2);
  sigma2 = minmod(2.0 * slope1,slope2);

  if (sigma1 * sigma2 <= 0.0){
    slopec = 0.0;
  }else{

    if( sigma1*sigma1 >= sigma2*sigma2){
      slopec = sigma1;
    }else{
      slopec = sigma2;
    }
  }

  return slopec;
}
