#ifndef TT_SPECTRAL_VALIDATION_H
#define TT_SPECTRAL_VALIDATION_H
#include "turbotuv.h"
/* Shared structural checks for direct evaluation and prepared constructors. */
tt_status tt_cross_section_validate(tt_cross_section_kind kind, size_t wavelengths, size_t count,
                                    const tt_spectral_table *tables, size_t *required);
tt_status tt_quantum_yield_validate(tt_quantum_yield_kind kind, size_t wavelengths, size_t count,
                                    const tt_spectral_table *tables,
                                    const tt_yield_options *options, int *needs_data);
#endif
