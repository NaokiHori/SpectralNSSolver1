#if !defined(LOGGING_H)
#define LOGGING_H

#include "domain.h"
#include "fluid.h"

typedef struct {
  int (* const init)(
      const domain_t * domain,
      const double time
  );
  void (* const check_and_output)(
      const domain_t * domain,
      const size_t step,
      const double time,
      const double dt,
      const double wtime,
      const fluid_t * fluid
  );
  double (* const get_next_time)(
      void
  );
} logging_t;

extern const logging_t logging;

#endif // LOGGING_H
