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