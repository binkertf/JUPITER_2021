#include "jupiter.h"

void SetCommSource () {
  long gncell[3], stride[3];
  long i[3],k,m, size;
  char *CommSource, *CommDest, *SrcCommPerp, *DestCommPerp;
  Communicator *com;
  tGrid_CPU *cpug;
  boolean srcsame, srcup, srcsame_ext, srcup_ext;
  cpug=Grid_CPU_list;
  while (cpug != NULL) {
    if (cpug->cpu == CPU_Rank) {
      getgridsize (cpug, gncell, stride);
      com=ComListGhost;
      size = gncell[0]*gncell[1]*gncell[2];
      CommSource = (char *)prs_malloc(size*sizeof(char));
      CommDest   = (char *)prs_malloc(size*sizeof(char));
      SrcCommPerp = (char *)prs_malloc(size*sizeof(char));
      DestCommPerp = (char *)prs_malloc(size*sizeof(char));
      for (k = 0; k < size; k++) {
	CommSource[k] = EVERYWHERE;
	CommDest[k] = EVERYWHERE;
	SrcCommPerp[k] = 0;
	DestCommPerp[k] = 0;
      }
      cpug->CommSource = CommSource;
      cpug->CommDest   = CommDest;
      cpug->DestCommPerp = DestCommPerp;
      cpug->SrcCommPerp = SrcCommPerp;
      while (com != NULL) {
	if (com->srcg == cpug) {
	  for (i[0] = 0; i[0] < gncell[0]; i[0]++) {
	    for (i[1] = 0; i[1] < gncell[1]; i[1]++) {
	      for (i[2] = 0; i[2] < gncell[2]; i[2]++) {
		m = i[0]*stride[0]+i[1]*stride[1]+i[2]*stride[2];
		srcsame = srcup = srcsame_ext = srcup_ext = 0;
		if (com->dest_level == com->src_level) {
		  srcsame = SRCGHOSTSAME;
		  srcsame_ext = SRCGHOSTSAMEEXTEND;
		}
		else {
		  srcup = SRCGHOSTUP;
		  srcup_ext = SRCGHOSTUPEXTEND;
		}
		for (k = 0; k < 3; k++) {
		  if ((i[k] < com->imin_src[k]) || (i[k] >= com->imax_src[k]))
		    srcsame = srcup = 0;
		  if ((i[k] < com->imin_src[k]-1) || (i[k] >= com->imax_src[k]+1))
		    srcsame_ext = srcup_ext = 0;
		}
		CommSource[m] |= (srcsame | srcsame_ext | srcup | srcup_ext);
		if (srcup)
		  SrcCommPerp[m] |= 1<<(com->facedim);
	      }
	    }
	  }
	}
	if (com->destg == cpug) {
	  for (i[0] = 0; i[0] < gncell[0]; i[0]++) {
	    for (i[1] = 0; i[1] < gncell[1]; i[1]++) {
	      for (i[2] = 0; i[2] < gncell[2]; i[2]++) {
		m = i[0]*stride[0]+i[1]*stride[1]+i[2]*stride[2];
		srcsame = srcup = srcsame_ext = srcup_ext = 0;
		/* src.. in order to use same variables as before, but */
		/* these have meaning of 'dest' instead of 'src'. */
		if (com->dest_level == com->src_level) {
		  srcsame = DESTGHOSTSAME;
		  srcsame_ext = DESTGHOSTSAMEEXTEND;
		}
		else {
		  srcup = DESTGHOSTUP;
		  srcup_ext = DESTGHOSTUPEXTEND;
		}
		for (k = 0; k < 3; k++) {
		  if ((i[k] < com->imin_dest[k]) || (i[k] >= com->imax_dest[k]))
		    srcsame = srcup = 0;
		  if ((i[k] < com->imin_dest[k]-1) || (i[k] >= com->imax_dest[k]+1))
		    srcsame_ext = srcup_ext = 0;
		}
		CommDest[m] |= (srcsame | srcsame_ext | srcup | srcup_ext);
		if (srcup)
		  DestCommPerp[m] |= 1<<(com->facedim);
	      }
	    }
	  }
	}
	com = com->next;
      }
    }
    cpug = cpug->next;
  }
}
