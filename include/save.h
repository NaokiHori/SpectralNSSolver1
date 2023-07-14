#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"

typedef struct save_t_ {
  int (* const init)(
    const domain_t * domain,
    const double time
  );
  int (* const prepare)(
      const domain_t * domain,
      const size_t step,
      char ** dirname
  );
  double (* const get_next_time)(
      void
  );
} save_t;

extern const save_t save;

#endif // SAVE_H
