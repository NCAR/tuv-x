/* Internal readers for compiled text tables; no filesystem access. */
#ifndef TT_TABLE_INTERNAL_H
#define TT_TABLE_INTERNAL_H
#include "turbotuv.h"
/* Return an owned column-major array. Skip the first line for grid files,
 * whose first record is a count or header. Ignore blank/comment lines. */
tt_status tt_text_columns(const tt_dataset *data, size_t columns, int skip_header, size_t *rows,
                          double **values);
#endif
