#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*********************************************************************/
/*FUNCTION to parse float parameters from input line*/
int PARSE_PAR_LONG(int argc, char *argv[], long *output_long, int *i_arg)
{
    int err_status = EXIT_FAILURE;

    if (argc >= (*i_arg)) {
       *output_long = atoi(argv[*i_arg - 1]);
       err_status = EXIT_SUCCESS;
       ++(*i_arg);
    }
    return(err_status);
}/*END PARSE_PAR_INT*/
/*********************************************************************/
/*FUNCTION to parse float parameters from input line*/
int PARSE_PAR_FLOAT(int argc, char *argv[], float *output_flt, int *i_arg)
{
    int err_status = EXIT_FAILURE;
    if (argc >= (*i_arg)) {
       *output_flt = atof(argv[*i_arg - 1]);
       err_status = EXIT_SUCCESS;
       ++(*i_arg);
    }
    return(err_status);
}/*END PARSE_PAR_FLOAT*/
/*************************************************************************/
/*FUNCTION to parse char parameters from input line*/
int PARSE_PAR_CHAR(int argc, char *argv[], char *output_str, int *i_arg)
{
    int err_status = EXIT_FAILURE;
    if (argc >= (*i_arg)) {
       strcpy(output_str, argv[*i_arg - 1]);
       err_status = EXIT_SUCCESS;
       ++(*i_arg);
    }
    return(err_status);
}/*END PARSE_PAR_CHAR*/
/*********************************************************************/
/*FUNCTION to parse int parameters from input line*/
int PARSE_PAR_INT(int argc, char *argv[], int *output_int, int *i_arg)
{
    int err_status = EXIT_FAILURE;
    if (argc >= (*i_arg)) {
       *output_int = atoi(argv[*i_arg - 1]);
       err_status = EXIT_SUCCESS;
       ++(*i_arg);
    }
    return(err_status);
}/*END PARSE_PAR_INT*/
/*********************************************************************/
/*FUNCTION to parse int parameters from input line*/
int PARSE_PAR_DOUBLE(int argc, char *argv[], double *output_dbl, int *i_arg)
{
    int err_status = EXIT_FAILURE;
    if (argc >= (*i_arg)) {
       *output_dbl = atof(argv[*i_arg - 1]);
       err_status = EXIT_SUCCESS;
       ++(*i_arg);
    }
    return(err_status);
}/*END PARSE_PAR_DOUBLE*/
/*********************************************************************/

