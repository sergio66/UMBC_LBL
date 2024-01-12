#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "hdefs.h"


main (int argc, char *argv[]) {

  struct HVAL *p;
  extern struct HVAL * read_hitran();

  p = read_hitran(700.0, 800.0, 0.0, 3, "hitran98.by.gas/g3.dat");

  while (p != NULL) {
    write_hval(p, stdout);
    p = p->next;
  }
}

