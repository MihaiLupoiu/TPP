
/* Copyright (c) 1995 Gary W. Flake -- Permission granted for any use
 * provied that the author's comments are neither modified nor removed.
 * No warranty is given or implied.
 *
 * NAME
 *   misc.c
 * PURPOSE
 *   Miscellaneous routines for program examples.  Most of the programs
 *   will need to be linked to this code.
 */

#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "misc.h"

/* A little buffer for formatting option help string entries. */

#define BUFFERLEN 4096

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static void format_string(char *str, int leftind, int rightind, int firstind)
{
  int i, start, stop, len, save, first = 1;
  char *indent;
  
  indent = xmalloc(sizeof(char) * (MAX(leftind, firstind) + 1));
  for(i = 0; i < firstind; i++)
    indent[i] = ' ';
  indent[firstind] = 0;
  fprintf(stderr, "%s", indent);
  for(i = 0; i < leftind; i++)
    indent[i] = ' ';
  indent[leftind] = 0;
  len = strlen(str);
  start = 0;
  while(1) {
    if(start >= len) break;
    if(!first) {
      while(str[start] == ' ') start++;
      fprintf(stderr, "%s", indent);
      stop = start + (80 - leftind - rightind);
    }
    else
      stop = start + (80 - firstind - rightind);
    if(stop > len)
      stop = len - 1;
    else
      while(!isspace(str[stop]) && stop > start) stop--;
    if(stop <= start) break;
    save = str[stop]; str[stop] = 0;
    fprintf(stderr, "%s", str + start);
    str[stop] = save;
    if(str[stop] != ' ')
      fputc(str[stop], stderr);
    start = stop + 1;
    fputc('\n', stderr);
    first = 0;
  }
  free(indent);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define HELP_INDENT 4

static void display_options(char **argv, OPTION *options, char *help)
{
  int i, j, maxlen, len;
  char buffer[BUFFERLEN], temp[256];
  
  fprintf(stderr, "Usage: %s [ options ]\n\n", argv[0]);
  format_string(help, HELP_INDENT, HELP_INDENT, HELP_INDENT);
  fprintf(stderr, "\nOptions with defaults in parentheses are:\n\n");
  
  i = maxlen = 0;
  while(options[i].name != NULL) {
    if(maxlen < (len = strlen(options[i].name)))
      maxlen = len;
    i++;
  }
  
  i = 0;
  while(options[i].name != NULL) {
    buffer[0] = 0;
    strcat(buffer, options[i].name);    
    len = strlen(options[i].name);
    for(j = 0; j < maxlen - len + 2; j++)
      buffer[len + j] = ' ';
    buffer[maxlen + 2] = 0;
    strcat(buffer, options[i].help);
    switch(options[i].type) {
      case OPT_INT:
        sprintf(temp, " (%d)", *(int *)options[i].ptr);
        strcat(buffer, temp);
        format_string(buffer, maxlen + 6, 4, 4);
        break;
      case OPT_DOUBLE:
        sprintf(temp, " (%g)", *(double *)options[i].ptr);
        strcat(buffer, temp);
        format_string(buffer, maxlen + 6, 4, 4);
        break;
      case OPT_STRING:
        sprintf(temp, " (\"%s\")", *(char **)options[i].ptr);
        strcat(buffer, temp);
        format_string(buffer, maxlen + 6, 4, 4);
        break;
      case OPT_SWITCH:
        sprintf(temp, (*(int *)options[i].ptr) ? " (ON)" : " (OFF)");
        strcat(buffer, temp);
        format_string(buffer, maxlen + 6, 4, 4);
        break;
      case OPT_OTHER:
        sprintf(temp, " (no default)");
        strcat(buffer, temp);
        format_string(buffer, maxlen + 6, 4, 4);
        break;
      default:
        break;
    }
    i++;
  }
  fprintf(stderr, "\n");
  exit(1);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void bad_option(char **argv, int badopt)
{
  fprintf(stderr, "%s: unknown or incorrectly used option \"%s\".\n",
          argv[0], argv[badopt]);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void get_options(int argc, char **argv, OPTION *options, char *help)
{
  int i, j, found;
  
  /* For each argument in the command line.
   */
  i = 1;
  while(i < argc) {
    found = 0, j = 0;
    /*
     * Search for the current word in the option list.
     */
    while(!found) {
      /*
       * If the option was not found.
       */
      if(strcmp(argv[i], "-help") == 0) {
        display_options(argv, options, help);
      }
      else if(options[j].name == NULL) {
        bad_option(argv, i);
        display_options(argv, options, help);
      }
      /*
       * If the word is found in the option list.
       */ 
      else if(strcmp(argv[i], options[j].name) == 0) {
        /*
         * Check the type, make sure we have another arg
         * coming if needed, eat up the next arg, set
         * the value of the parameter, and display help
         * if needed.
         */
        found = 1;
        switch(options[j].type) {
          case OPT_INT:
            if(i + 1 >= argc) {
              bad_option(argv, i);
              display_options(argv, options, help);
            }
            *(int *)options[j].ptr = atoi(argv[i + 1]);
            i += 2;
            break;
          case OPT_DOUBLE:
            if(i + 1 >= argc) {
              bad_option(argv, i);
              display_options(argv, options, help);
            }
            *(double *)options[j].ptr = atof(argv[i + 1]);
            i += 2;
            break;
          case OPT_STRING:
            if(i + 1 >= argc) {
              bad_option(argv, i);
              display_options(argv, options, help);
            }
            *(char **)options[j].ptr = argv[i + 1];
            i += 2;
            break;
          case OPT_SWITCH:
            *(int *)options[j].ptr = !*(int *)options[j].ptr;
            i += 1;
            break;
          case OPT_OTHER:
            if(((int (*)())options[j].ptr)(argv, argc, options, &i, j)) {
              bad_option(argv, i);
              display_options(argv, options, help);
            }
            break;
          default:
            break;
        }
      }
      j++;
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void *xmalloc(size_t bytes)
{
  void *ptr;

  if((ptr = malloc(bytes)) == NULL) {
    fprintf(stderr, "malloc failed to get %d bytes.\n", (int)bytes);
    exit(1);
  }
  return(ptr);
}

