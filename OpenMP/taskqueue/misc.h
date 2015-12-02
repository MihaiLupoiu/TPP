
/* Copyright (c) 1995 Gary W. Flake -- Permission granted for any use
 * provied that the author's comments are neither modified nor removed.
 * No warranty is given or implied.
 *
 * NAME
 *   misc.h
 * PURPOSE
 *   Miscellaneous routines for program examples.  Most of the programs
 *   will need to be linked to this code.
 */

/* Insures that this file is only loaded once. */

#ifndef __MISC_H__
#define __MISC_H__

#ifdef  __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Command line option routines... */

/* Possible command line option types.  OPT_NULL is a tag to indicate
   the end of options. */

typedef enum OPT_TYPE {
  OPT_INT, OPT_DOUBLE, OPT_STRING, OPT_SWITCH, OPT_OTHER, OPT_NULL
} OPT_TYPE;


/* An option consists of a name, a type enum, and a ptr to either the
   location of the storage or a ptr to a special user supplied function
   that partially parses the line for that one type. */

typedef struct OPTION {
  const char *name;
  OPT_TYPE type;
  void *ptr;
  const char *help;
} OPTION;


/* Function to automagically parse command line arguments. */

extern void get_options(int argc, char **argv, OPTION *options, char *help);

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Function to get memory with check for failure built in. */

void *xmalloc(size_t bytes);


/* Miscelaneous macros. */

#define MIN(x, y)     ((x) < (y) ? (x) : (y))
#define MAX(x, y)     ((x) > (y) ? (x) : (y))
#define ABS(x)        ((x) >= 0 ? (x) : -(x))
#define SIGN(x)       ((x) >= 0 ? 1 : -1)
#define SQR(x)        ((x)*(x))
#define GETBIT(c,i)   (((1 << (i)) & (c)) ? 1 : 0)
#define SETBIT(c,i,b) ((b) ? ((c) | (1 << (i))) : ((c) & ~(1 << (i))))

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef  __cplusplus
}
#endif


#endif /* __MISC_H__ */
