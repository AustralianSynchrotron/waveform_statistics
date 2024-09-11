/* $File: //ASP/tec/epics/wfs/trunk/Waveform_StatisticsSup/Waveform_Statistics_Subroutines.c $
 * $Revision: #8 $
 * $DateTime: 2023/09/15 15:31:35 $
 * Last checked in by: $Author: starritt $
 *
 * Description
 * Waveform_Statistics calculation sub-routine, for use with the aSubRecord.
 * This module allows the calculation of statistics of values extracted from
 * multi-valued records such as a waveform, subArray or concat record.
 *
 * Field usage:
 * INPA: epicsFloat64 array - input values
 * INPB: long - input data size - must be > 0
 * INPC: long - data offset - must be >= 0 and < NOA.
 * INPD: epicsFloat64 - sample interval - defaults to 1 if input not DOUBLE.
 * INPE: long array - input mask. LSB defines mask value - other bits ignored.
 *
 * Notes:
 * Number of elements processed is minimum of (NOA - INPC, INPB).
 * The data offset (INPC) applied to both INPA and INPE.
 * For INPE - if input type is not LONG, then all mask elements are deemed to be true.
 * For INPE - if element index exceeds NOE, then the mask element is deemed to be true.
 * Thus, as the default INPE type is DOUBLE, an unspecified mask implies all elements used.
 *
 *
 * OUTA: epicsFloat64 - mean value
 * OUTB: epicsFloat64 - minimum value
 * OUTC: epicsFloat64 - maximum value
 * OUTD: epicsFloat64 - standard deviation (based on sample variance, i.e. / (N-1))
 * OUTE: epicsFloat64 - total
 * OUTF: epicsFloat64 - median value
 * OUTG: epicsFloat64 - least squares fit slope (m)     as in y = m.x + c
 * OUTH: epicsFloat64 - least squares fit interset (c)  as in y = m.x + c
 * OUTI: epicsFloat64 - maximum absolute value
 * OUTJ: epicsFloat64 - root mean square (RMS) value
 * OUTK: epicsFloat64 - standard deviation (based on population variance, i.e. / N)
 * OUTL: epicsInt32   - acutal number of elements used to calculate the statistics
 *                      taking into account INPB, INPC and INPE.
 *
 * Note: G/H assume equi-spaced samples.
 *
 * Others outputs may be added as and when required.
 *
 *
 * Source code formatting:  indent -kr -pcs -i3 -cli3 -nbbo -nut
 *
 * Copyright (c) 2010-2024  Australian Synchrotron
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * Licence as published by the Free Software Foundation; either
 * version 2.1 of the Licence, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public Licence for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * Licence along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Contact details:
 * andrew.starritt@synchrotron.org.au
 * 800 Blackburn Road, Clayton, Victoria 3168, Australia.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <alarm.h>
#include <cantProceed.h>
#include <dbDefs.h>
#include <dbStaticLib.h>
#include <epicsExport.h>
#include <epicsTypes.h>
#include <errlog.h>
#include <aSubRecord.h>
#include <menuFtype.h>
#include <recGbl.h>
#include <registryFunction.h>


/*-------------------------------------------------------------------------------
 * Common macro functions.
 *-------------------------------------------------------------------------------
 */
#define ABS(a)       ((a) >= 0  ? (a) : -(a))
#define MIN(a, b)    ((a) <= (b) ? (a) : (b))
#define MAX(a, b)    ((a) >= (b) ? (a) : (b))


/*-------------------------------------------------------------------------------
 * Local utilities
 *-------------------------------------------------------------------------------
 */
static epicsInt32 long_at (const void* item)
{
   return *((epicsInt32*) item);
}

/* -----------------------------------------------------------------------------
 */
static epicsFloat64 double_at (const void* item)
{
   return *((epicsFloat64*) item);
}

/* -----------------------------------------------------------------------------
 */
static int cmp_double_p (const void* p1, const void* p2)
{
   epicsFloat64 d1, d2;
   int result;

   /* The actual arguments to this function are "pointers to
    * epicsFloat64", hence the following cast plus dereference.
    */
   d1 = *(epicsFloat64* const) p1;
   d2 = *(epicsFloat64* const) p2;

   if (d1 < d2) {
      result = -1;
   } else if (d1 > d2) {
      result = +1;
   } else {
      result = 0;
   }
   return result;
}


/*------------------------------------------------------------------------------
 *
 */
static epicsFloat64 calcMedian (const int n, const epicsFloat64* data_set)
{
   epicsFloat64 result = 0.0;
   epicsFloat64* work_list = NULL;
   int j;
   int mid_index;

   /* Make working copy
    */
   work_list = (epicsFloat64*) calloc ((size_t) n, sizeof (epicsFloat64));

   for (j = 0; j < n; j++) {
      work_list[j] = data_set [j];
   }

   /* Sort - use quick sort see "man 3 qsort" for details.
    *
    *   void qsort (void *base, size_t nmemb, size_t size,
    *               int(*compar)(const void *, const void *));
    */
   qsort (work_list, (size_t) n, sizeof (epicsFloat64), cmp_double_p);

   /* Choose middle item (round up if even number)
    */
   mid_index = n / 2;
   result = work_list [mid_index];

   free (work_list);  /* free working copy */

   return result;
}


/*------------------------------------------------------------------------------
 * Record functions
 *------------------------------------------------------------------------------
 */
static long Waveform_Statistics_Init (aSubRecord* precord)
{
   /* This function is just a place holder.
    */
   return 0;
}


/*------------------------------------------------------------------------------
 */
static long Waveform_Statistics_Process (aSubRecord* precord)
{
   int number_elements;
   int requested_size;
   int offset;
   int size;
   int last_point;
   int actual_size;
   int j;
   epicsInt32* m_ptr;
   epicsFloat64* a_ptr;
   epicsFloat64 sample_interval;
   epicsFloat64 a;
   epicsFloat64 a_sum;
   epicsFloat64 a_squared_sum;
   epicsFloat64 a_mean;
   epicsFloat64 a_median;
   epicsFloat64 a_squared_mean;
   epicsFloat64 a_min;
   epicsFloat64 a_max;
   epicsFloat64 a_max_abs;
   epicsFloat64 a_rms;
   epicsFloat64 pop_variance;
   epicsFloat64 sam_variance;
   epicsFloat64 pop_std_dev;
   epicsFloat64 sam_std_dev;
   epicsFloat64 x;
   epicsFloat64 x_sum;
   epicsFloat64 xx_sum;
   epicsFloat64 xa_sum;
   epicsFloat64 delta;
   epicsFloat64 slope;
   epicsFloat64 intersect;

   /* Verify that field types are as expected.
    */
   if ((precord->fta != menuFtypeDOUBLE) ||
       (precord->ftb != menuFtypeLONG) ||
       (precord->ftc != menuFtypeLONG)) {
      errlogPrintf
          ("WFS: (%s) incorrect FTA, FTB and/or FTC type specified.\n",
           precord->name);
      return -1;
   }

   /* Read available data size and specified data size.
    */
   number_elements = precord->noa;
   requested_size = (int) long_at (precord->b);
   offset = (int) long_at (precord->c);

   if (precord->ftd == menuFtypeDOUBLE) {
      sample_interval = double_at (precord->d);
      if (sample_interval == 0.0) sample_interval = 1.0;
   } else {
      sample_interval = 1.0;
   }

   /* Find lesser of two and verify sensible.
    */
   size = MIN (number_elements - offset, requested_size);
   if (size < 1) {
      errlogPrintf
          ("WFS: (%s) size, min of (noa=%d - inpc=%d, inpb=%d), must be at least 1\n",
           precord->name, number_elements, offset, requested_size);
      return -1;
   }

   /* Process input values.
    */
   a_ptr = (epicsFloat64 *) precord->a;
   m_ptr = (epicsInt32 *) precord->e;
   a_sum = 0.0;
   a_squared_sum = 0.0;
   a_min = +INFINITY;
   a_max = -INFINITY;
   a_max_abs = 0.0;
   a_median = 0.0;

   x_sum = 0.0;
   xx_sum = 0.0;
   xa_sum = 0.0;

   last_point = offset + size - 1;
   actual_size = 0;

   /* Macro to check is this element in use.
    */
   #define IN_USE(x) ( (precord->fte != menuFtypeLONG) ||  \
                       (x >= precord->noe)             ||  \
                       ((m_ptr [x] & 1) == 1) )


   for (j = offset; j < offset + size; j++) {
      /* INPE must be long and define as false to mask out a value.
       */
      int inUse = IN_USE (j);
      if (!inUse) continue;
      a = a_ptr[j];

      a_sum += a;
      a_squared_sum += (a * a);

      a_min = MIN (a_min, a);
      a_max = MAX (a_max, a);
      a_max_abs = MAX(a_max_abs, ABS(a));

      /* least squares - use last point as origin
       * a plays the role of y.
       */
      x = (epicsFloat64) (j - last_point) * sample_interval;
      x_sum  += x;
      xx_sum += (x * x);
      xa_sum += (x * a);

      actual_size++;
   }

   if (actual_size < 1) {
      errlogPrintf
          ("WFS: (%s) at least one element must be included\n",
           precord->name);
      return -1;
   }

   /* Perform actual statistical calculations.
    */
   a_mean = a_sum / (epicsFloat64) actual_size;
   a_squared_mean = a_squared_sum / (epicsFloat64) actual_size;


   /* Calculate population variance, and then sample variance.
    */
   pop_variance = a_squared_mean - (a_mean * a_mean);
   if (actual_size >= 2) {
      sam_variance = (actual_size * pop_variance) / (actual_size - 1.0);
   } else {
      sam_variance = 0.0;
   }

   pop_std_dev = sqrt (pop_variance);
   sam_std_dev = sqrt (sam_variance);

   /* Calc root mean square.
    */
   a_rms = sqrt (a_squared_mean);

   /* Least squares fit algo.
    */
   delta = (actual_size * xx_sum) - (x_sum * x_sum);
   if (delta == 0.0) delta = 1.0;
   slope = ((actual_size * xa_sum) - (x_sum * a_sum)) / delta;
   intersect = ((a_sum * xx_sum) - (x_sum * xa_sum)) / delta;

   if (actual_size < size) {
      /* shuffle a to calculate median - fill the gaps - this is safe */
      int t = 0;
      for (j = offset; j < offset + size; j++) {
         /* INPE must be long and define as false to mask out a value.
          */
         int inUse = IN_USE (j);
         if (!inUse) continue;
         a_ptr [t] = a_ptr[j];
         t++;
      }
      a_median = calcMedian (actual_size, &a_ptr [0]);
   } else {
      a_median = calcMedian (actual_size, &a_ptr [offset]);
   }

   /* Write to output value fields iff defined, i.e. of type epicsFloat64.
    */
   if (precord->ftva == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->vala;
      *out_ptr = a_mean;
   }

   if (precord->ftvb == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valb;
      *out_ptr = a_min;
   }

   if (precord->ftvc == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valc;
      *out_ptr = a_max;
   }

   if (precord->ftvd == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->vald;
      *out_ptr = sam_std_dev;
   }

   if (precord->ftve == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->vale;
      *out_ptr = a_sum;
   }

   if (precord->ftvf == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valf;
      *out_ptr = a_median;
   }

   if (precord->ftvg == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valg;
      *out_ptr = slope;
   }

   if (precord->ftvh == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valh;
      *out_ptr = intersect;
   }

   if (precord->ftvi == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->vali;
      *out_ptr = a_max_abs;
   }

   if (precord->ftvj == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valj;
      *out_ptr = a_rms;
   }

   if (precord->ftvk == menuFtypeDOUBLE) {
      epicsFloat64* out_ptr = (epicsFloat64*) precord->valk;
      *out_ptr = pop_std_dev;
   }

   if (precord->ftvl == menuFtypeLONG) {
      epicsInt32 *out_ptr = (epicsInt32*) precord->vall;
      *out_ptr = actual_size;
   }

   return 0;
}


/* -----------------------------------------------------------------------------
 */
epicsRegisterFunction (Waveform_Statistics_Init);
epicsRegisterFunction (Waveform_Statistics_Process);

/* end */
