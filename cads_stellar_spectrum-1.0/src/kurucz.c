/***************************************************************************
 * This File is a part of the CADS/UV Stellar Spectrum model software      *
 *   Copyright (C) 2012 by CADS/UV Software Team,                          *
 *                         Indian Institute of Astrophysics                *
 *                         Bangalore 560034                                *
 *                         cads_AT_iiap.res.in                             *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*-------------------------------------------------------------------------*
File Description: The following functions are defined in this file
1. READ_TABLE() Selects are reads input file MKSpTypewwww.in
2. RCross_Sec() Estimates the Absorption coefficient using "crossec0.dat"
3. GET_REF1()   Gets the flux in various bands
4. ALPHA()      Calculate Scale value for Kurucz flux]
5. WRITE_LINE() Write spectrum to output file (one line per call)
*-------------------------------------------------------------------------*/

#include "stellar_spectrum.h"

/*---------------- Select & Reads input file MKSpTypewwww.in--------------*/
int READ_TABLE(struct S_TABLE *mkstTG)
{
    FILE           *fp =NULL;
    int             j;
    int             Nmax;
    char            s_tmp[MAX_TEXT], TG_tmp[MAX_TEXT], tok[MAX_TEXT];
    char            KFdir[MAX_TEXT];
    char            tab_line[TAB_LEN], *token, filename[MAX_TEXT];

    strcpy(filename,mkstTG->w_filename);

    if ((fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open file %s\n",filename);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "%i", &Nmax);
    mkstTG->N_tab = Nmax;

    fgets(tab_line, TAB_LEN, fp);

    fgets(tab_line, TAB_LEN, fp);
    sscanf(tab_line, "%s", KFdir);
    strcpy(mkstTG->KFluxdir0, KFdir);

    fgets(tab_line, TAB_LEN, fp);
    fgets(tab_line, TAB_LEN, fp);

    /*LOOP BEGINS HERE*/
    for (j=1; j <= Nmax; j++)
    {

        fgets(tab_line, TAB_LEN, fp);

        token = strtok(tab_line, TAB_DELIM);
        /* Reads Sp_Type */
        strcpy(tok, token);
        sscanf(tok, "%s", s_tmp);
        strcpy(mkstTG->star_name[j], s_tmp);

        token = strtok(NULL, TAB_DELIM);
        /* Reads Temp. Gravity Filename */
        strcpy(tok, token);
        sscanf(tok, "%s", TG_tmp);
        strcpy(mkstTG->TG_name[j], TG_tmp);

        token = strtok(NULL, TAB_DELIM);
        /* Reads Bflux */
        mkstTG->Bflux[j] = atof(token);
        token = strtok(NULL, TAB_DELIM);

        /* Reads Vflux */
        mkstTG->Vflux[j] = atof(token);
        token = strtok(NULL, TAB_DELIM);

        /* Reads Observed flux */
        mkstTG->Oflux[j] = atof(token);
    }/*End for */

    fclose(fp);
    return (EXIT_SUCCESS);

}/* End READ_TABLE */

/*---------------------------- RCross_Sec -------------------------------*/
int RCross_Sec( struct STARS *line, char filename[MAX_TEXT])
{
    FILE           *fp = NULL;
    int             n1;
    float           Wo, WW1, WW2, C_Sectiono;

    if ((fp = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open file %s \n",filename);
        exit(EXIT_FAILURE);
    }
    else
    {
        fscanf(fp,"%i%f", &n1, &line->Wlength[1]);
    }

    fclose(fp);
    fp = NULL;

    /* The following steps are to find out the Absorption CrossSection from
           "crossec0.dat"  corresponding to our wavelength of interest*/

    if ((fp = fopen(line->SigmaFile, "r")) == NULL)
    {
        fprintf(stderr, "ERROR  : Unable to open Cross Section File %s\n",
                line->SigmaFile);
        exit(EXIT_FAILURE);
    }

    WW1 = pow(line->Wlength[1],2);
    fscanf(fp, "%f %e", &Wo, &C_Sectiono);
    while (feof(fp) == 0)
    {
        WW2 = pow((Wo - line->Wlength[1]),2);
        if (WW2 <= WW1)
        {
            WW1 = WW2;
            line->C_Section[1] = C_Sectiono;
        }
        fscanf(fp, "%f %e", &Wo, &C_Sectiono);
    }    /* End while feof */
    fclose(fp);

    return (EXIT_SUCCESS);
}


/*----------------------------- GET_REF1 --------------------------------*/
int GET_REF1(struct STARS * line, struct S_TABLE * mkstTG,
             char input_file[MAX_TEXT])
{
    int    i = 1, gr = 1;

    while ((strcmp(line->sp_type, mkstTG->star_name[i]) != 0)
            && (i<=mkstTG->N_tab))
        i++;
    if (i <= mkstTG->N_tab)
    {
        strcpy(line->tg_type, mkstTG->TG_name[i]);
        line->Bflux = mkstTG->Bflux[i];
        line->Vflux = mkstTG->Vflux[i];
        line->Oflux = mkstTG->Oflux[i];
        gr = 0;
    }
    else
    {
        fprintf(stderr, "ERROR  : Unknown Spectral-Type : '%s'\n",
                line->sp_type);
        exit(EXIT_FAILURE);
    }
    return (gr);
}

/*------------------ Calculate Scale value for Kurucz flux ----------*/
int ALPHA(struct STARS *line)
{
    double     Tau, Fmv;
    double flux_zero = 3.64e-9;

    if (line->Ebv >=0.0)
    {
        Tau = ( line->Ebv * 3.1) / 1.0863 ;
        Fmv = flux_zero * pow(10.0,(-0.4*line->mv));
        if (line->Vflux > 0.0) line->alpha = ( Fmv * exp(Tau))/line->Vflux;
    } /* end if */
    else
        line->alpha = 0.0 ;
    return (EXIT_SUCCESS);
}

/*---------------- Write output file "hip_wwww.dat" ---------------------*/

int WRITE_LINE(FILE *staro_file, struct STARS *line)
{
    double           f_tmp,fstar,tau;
    double   gas_to_dust = 5.8e21;

    f_tmp = (line->alpha * line->Oflux);
    tau = line->C_Section[1] * line->Ebv * gas_to_dust;
    fstar = f_tmp * exp(-tau);
    fprintf(staro_file, " %5.8f  %5.8e \n", line->Wlength[1],fstar);
    return (EXIT_SUCCESS);
}