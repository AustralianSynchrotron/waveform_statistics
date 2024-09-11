# waveform_statistics

## Introduction

Waveform_Statistics calculation sub-routine, for use with the EPICS aSubRecord.
This module allows the calculation of statistics of values extracted from
multi-valued records such as a waveform, subArray or concat records.

The  available values are

- mean
- minimum value
- maximum value
- standard deviation (based on sample variance)
- total
- median value
- least squares fit slope (m)     as in y = m.x + c
- least squares fit intersect (c)  as in y = m.x + c
- maximum absolute value
- standard deviation (based on population variance)
- total number of samples (integer)



## aSub record field usage

The comments/documentation within the Waveform_Statistics_Subroutines.c file is extensive and should be read.

### processing

| processing | fucntion name               |
|:-----------|:----------------------------|
| INAM       | Waveform_Statistics_Init    |
| SNAM       | Waveform_Statistics_Process |


### inputs

| field  |  FTx   |   NOx      | Comment                                               |
|:------|:-------|:------------|:------------------------------------------------------|
| INPA   | DOUBLE | array size | Input values                                          |
| INPB   | LONG   | 1          | Number of elements to be processes (must be > 0)      |
| INPC   | LONG   | 1          | Data offset - must be >= 0 and < NOA.                 |
| INPD   | DOUBLE | 1          | Sample interval - defaults to 1.0 if input not DOUBLE |
| INPE   | LONG   | 1          | Input mask. LSB defines mask value - other bits ignored (optional) |


### outputs

| field | FTVx   | NOVx | Comment                                           |
|:------|:-------|:-----|:--------------------------------------------------|
| OUTA  | DOUBLE | 1    | mean                                              |
| OUTB  | DOUBLE | 1    | minimum                                           |
| OUTC  | DOUBLE | 1    | maximum                                           |
| OUTD  | DOUBLE | 1    | standard deviation (based on sample variance)     |
| OUTE  | DOUBLE | 1    | total                                             |
| OUTF  | DOUBLE | 1    | median                                            |
| OUTG  | DOUBLE | 1    | least squares fit slope                           |
| OUTH  | DOUBLE | 1    | least squares fit intersect                       |
| OUTI  | DOUBLE | 1    | maximum absolute value                            |
| OUTJ  | DOUBLE | 1    | root mean square (RMS) value                      |
| OUTK  | DOUBLE | 1    | standard deviation (based on population variance) |
| OUTL  | LONG   | 1    | number of elements used                           |

Notes:

- acutal number of elements used to calculate the statistics taking into account INPB, INPC and INPE.



<font size="-1">Last updated: Wed Sep 11 13:45:03 AEST 2024</font>
<br>
