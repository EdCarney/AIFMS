/*
 * File: A_Star_S_Func_2.cpp
 *
 *
 *   --- THIS FILE GENERATED BY S-FUNCTION BUILDER: 3.0 ---
 *
 *   This file is an S-function produced by the S-Function
 *   Builder which only recognizes certain fields.  Changes made
 *   outside these fields will be lost the next time the block is
 *   used to load, edit, and resave this file. This file will be overwritten
 *   by the S-function Builder block. If you want to edit this file by hand, 
 *   you must change it only in the area defined as:  
 *
 *        %%%-SFUNWIZ_defines_Changes_BEGIN
 *        #define NAME 'replacement text' 
 *        %%% SFUNWIZ_defines_Changes_END
 *
 *   DO NOT change NAME--Change the 'replacement text' only.
 *
 *   For better compatibility with the Simulink Coder, the
 *   "wrapper" S-function technique is used.  This is discussed
 *   in the Simulink Coder's Manual in the Chapter titled,
 *   "Wrapper S-functions".
 *
 *  -------------------------------------------------------------------------
 * | See matlabroot/simulink/src/sfuntmpl_doc.c for a more detailed template |
 *  ------------------------------------------------------------------------- 
 *
 * Created: Thu Jun 13 17:57:42 2019
 */

#define S_FUNCTION_LEVEL 2
#define S_FUNCTION_NAME A_Star_S_Func_2
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/* %%%-SFUNWIZ_defines_Changes_BEGIN --- EDIT HERE TO _END */
#define NUM_INPUTS            6
/* Input Port  0 */
#define IN_PORT_0_NAME        nodeArray
#define INPUT_0_WIDTH         1000000
#define INPUT_DIMS_0_COL      4
#define INPUT_0_DTYPE         real_T
#define INPUT_0_COMPLEX       COMPLEX_NO
#define IN_0_FRAME_BASED      FRAME_NO
#define IN_0_BUS_BASED        0
#define IN_0_BUS_NAME         
#define IN_0_DIMS             2-D
#define INPUT_0_FEEDTHROUGH   1
#define IN_0_ISSIGNED         0
#define IN_0_WORDLENGTH       8
#define IN_0_FIXPOINTSCALING  1
#define IN_0_FRACTIONLENGTH   9
#define IN_0_BIAS             0
#define IN_0_SLOPE            0.125
/* Input Port  1 */
#define IN_PORT_1_NAME        edgeArray
#define INPUT_1_WIDTH         1000000
#define INPUT_DIMS_1_COL      3
#define INPUT_1_DTYPE         real_T
#define INPUT_1_COMPLEX       COMPLEX_NO
#define IN_1_FRAME_BASED      FRAME_NO
#define IN_1_BUS_BASED        0
#define IN_1_BUS_NAME         
#define IN_1_DIMS             2-D
#define INPUT_1_FEEDTHROUGH   1
#define IN_1_ISSIGNED         0
#define IN_1_WORDLENGTH       8
#define IN_1_FIXPOINTSCALING  1
#define IN_1_FRACTIONLENGTH   9
#define IN_1_BIAS             0
#define IN_1_SLOPE            0.125
/* Input Port  2 */
#define IN_PORT_2_NAME        numNodesPtr
#define INPUT_2_WIDTH         1
#define INPUT_DIMS_2_COL      1
#define INPUT_2_DTYPE         real_T
#define INPUT_2_COMPLEX       COMPLEX_NO
#define IN_2_FRAME_BASED      FRAME_NO
#define IN_2_BUS_BASED        0
#define IN_2_BUS_NAME         
#define IN_2_DIMS             1-D
#define INPUT_2_FEEDTHROUGH   1
#define IN_2_ISSIGNED         0
#define IN_2_WORDLENGTH       8
#define IN_2_FIXPOINTSCALING  1
#define IN_2_FRACTIONLENGTH   9
#define IN_2_BIAS             0
#define IN_2_SLOPE            0.125
/* Input Port  3 */
#define IN_PORT_3_NAME        numEdgesPtr
#define INPUT_3_WIDTH         1
#define INPUT_DIMS_3_COL      1
#define INPUT_3_DTYPE         real_T
#define INPUT_3_COMPLEX       COMPLEX_NO
#define IN_3_FRAME_BASED      FRAME_NO
#define IN_3_BUS_BASED        0
#define IN_3_BUS_NAME         
#define IN_3_DIMS             1-D
#define INPUT_3_FEEDTHROUGH   1
#define IN_3_ISSIGNED         0
#define IN_3_WORDLENGTH       8
#define IN_3_FIXPOINTSCALING  1
#define IN_3_FRACTIONLENGTH   9
#define IN_3_BIAS             0
#define IN_3_SLOPE            0.125
/* Input Port  4 */
#define IN_PORT_4_NAME        startNodePtr
#define INPUT_4_WIDTH         1
#define INPUT_DIMS_4_COL      1
#define INPUT_4_DTYPE         real_T
#define INPUT_4_COMPLEX       COMPLEX_NO
#define IN_4_FRAME_BASED      FRAME_NO
#define IN_4_BUS_BASED        0
#define IN_4_BUS_NAME         
#define IN_4_DIMS             1-D
#define INPUT_4_FEEDTHROUGH   1
#define IN_4_ISSIGNED         0
#define IN_4_WORDLENGTH       8
#define IN_4_FIXPOINTSCALING  1
#define IN_4_FRACTIONLENGTH   9
#define IN_4_BIAS             0
#define IN_4_SLOPE            0.125
/* Input Port  5 */
#define IN_PORT_5_NAME        goalNodePtr
#define INPUT_5_WIDTH         1
#define INPUT_DIMS_5_COL      1
#define INPUT_5_DTYPE         real_T
#define INPUT_5_COMPLEX       COMPLEX_NO
#define IN_5_FRAME_BASED      FRAME_NO
#define IN_5_BUS_BASED        0
#define IN_5_BUS_NAME         
#define IN_5_DIMS             1-D
#define INPUT_5_FEEDTHROUGH   1
#define IN_5_ISSIGNED         0
#define IN_5_WORDLENGTH       8
#define IN_5_FIXPOINTSCALING  1
#define IN_5_FRACTIONLENGTH   9
#define IN_5_BIAS             0
#define IN_5_SLOPE            0.125

#define NUM_OUTPUTS           4
/* Output Port  0 */
#define OUT_PORT_0_NAME       finalPath
#define OUTPUT_0_WIDTH        1
#define OUTPUT_DIMS_0_COL     1000000
#define OUTPUT_0_DTYPE        real_T
#define OUTPUT_0_COMPLEX      COMPLEX_NO
#define OUT_0_FRAME_BASED     FRAME_NO
#define OUT_0_BUS_BASED       0
#define OUT_0_BUS_NAME        
#define OUT_0_DIMS            2-D
#define OUT_0_ISSIGNED        1
#define OUT_0_WORDLENGTH      8
#define OUT_0_FIXPOINTSCALING 1
#define OUT_0_FRACTIONLENGTH  3
#define OUT_0_BIAS            0
#define OUT_0_SLOPE           0.125
/* Output Port  1 */
#define OUT_PORT_1_NAME       finalPathCost
#define OUTPUT_1_WIDTH        1
#define OUTPUT_DIMS_1_COL     1
#define OUTPUT_1_DTYPE        real_T
#define OUTPUT_1_COMPLEX      COMPLEX_NO
#define OUT_1_FRAME_BASED     FRAME_NO
#define OUT_1_BUS_BASED       0
#define OUT_1_BUS_NAME        
#define OUT_1_DIMS            1-D
#define OUT_1_ISSIGNED        1
#define OUT_1_WORDLENGTH      8
#define OUT_1_FIXPOINTSCALING 1
#define OUT_1_FRACTIONLENGTH  3
#define OUT_1_BIAS            0
#define OUT_1_SLOPE           0.125
/* Output Port  2 */
#define OUT_PORT_2_NAME       numNodesCheck
#define OUTPUT_2_WIDTH        1
#define OUTPUT_DIMS_2_COL     1
#define OUTPUT_2_DTYPE        real_T
#define OUTPUT_2_COMPLEX      COMPLEX_NO
#define OUT_2_FRAME_BASED     FRAME_NO
#define OUT_2_BUS_BASED       0
#define OUT_2_BUS_NAME        
#define OUT_2_DIMS            1-D
#define OUT_2_ISSIGNED        1
#define OUT_2_WORDLENGTH      8
#define OUT_2_FIXPOINTSCALING 1
#define OUT_2_FRACTIONLENGTH  3
#define OUT_2_BIAS            0
#define OUT_2_SLOPE           0.125
/* Output Port  3 */
#define OUT_PORT_3_NAME       numEdgesCheck
#define OUTPUT_3_WIDTH        1
#define OUTPUT_DIMS_3_COL     1
#define OUTPUT_3_DTYPE        real_T
#define OUTPUT_3_COMPLEX      COMPLEX_NO
#define OUT_3_FRAME_BASED     FRAME_NO
#define OUT_3_BUS_BASED       0
#define OUT_3_BUS_NAME        
#define OUT_3_DIMS            1-D
#define OUT_3_ISSIGNED        1
#define OUT_3_WORDLENGTH      8
#define OUT_3_FIXPOINTSCALING 1
#define OUT_3_FRACTIONLENGTH  3
#define OUT_3_BIAS            0
#define OUT_3_SLOPE           0.125

#define NPARAMS               0

#define SAMPLE_TIME_0         INHERITED_SAMPLE_TIME
#define NUM_DISC_STATES       0
#define DISC_STATES_IC        [0]
#define NUM_CONT_STATES       0
#define CONT_STATES_IC        [0]

#define SFUNWIZ_GENERATE_TLC  1
#define SOURCEFILES           "__SFB__"
#define PANELINDEX            8
#define USE_SIMSTRUCT         0
#define SHOW_COMPILE_STEPS    0
#define CREATE_DEBUG_MEXFILE  0
#define SAVE_CODE_ONLY        0
#define SFUNWIZ_REVISION      3.0
/* %%%-SFUNWIZ_defines_Changes_END --- EDIT HERE TO _BEGIN */
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
#include "simstruc.h"


extern void A_Star_S_Func_2_Outputs_wrapper(const real_T *nodeArray,
			const real_T *edgeArray,
			const real_T *numNodesPtr,
			const real_T *numEdgesPtr,
			const real_T *startNodePtr,
			const real_T *goalNodePtr,
			real_T *finalPath,
			real_T *finalPathCost,
			real_T *numNodesCheck,
			real_T *numEdgesCheck);
/*=============================*
 * Data Transposition Routines *
 *=============================*/
int linear_idx(const int srcdims[], const int numdims, const int dstk)
{
    int idxCurDim;
    int i, j;
    int srck = 0, dstk_remainingDims = dstk;
    for(i=0,j=numdims-1; i<numdims; i++, j--){
        idxCurDim = dstk_remainingDims % srcdims[j];
        srck = srck * srcdims[j] + idxCurDim;
        dstk_remainingDims = dstk_remainingDims / srcdims[j];
    }
    return srck;
}

/*
  A is [N1xN2xN3] and T is [N3xN2xN1]
  A[n1,n2,n3]=T[n3,n2,n1]
  
  1) Given linear index k = n1+N1*(n2+N2*n3), find [n1,n2,n3]
  2) Compute corresponding linear index k' = n3+N3*(n2+N2*n1)
  3) Use linear index to assign A[k] = T[k']

  Example input:
  src:     A matrix memory pointer
  srcdims: [N1,N2, N3]
  numdims: 3
  dst:     T matrix memory pointer
  elsize:  sizeof(double)/sizeof(char) 
 */
void NDTransposeBySrcSpecs(void *dst, const void *src, const int srcdims[], const int numdims, const int elsize)
{
    int w = srcdims[0];
    int k;
    for (k = 1; k < numdims; k ++) {
        w *= srcdims[k];
    }
   
    for (k = 0; k < w; k ++) {
        int sk = linear_idx(srcdims, numdims, k);
        int offset = k * elsize;
        memcpy((char*)dst + k * elsize, (const char*)src + sk * elsize, elsize);
    }
}

void NDTransposeByDstSpecs(void *dst, const void *src, const int dstdims[], const int numdims, const int elsize)
{
    int w = dstdims[0];
    int k;
    for (k = 1; k < numdims; k ++) {
        w *= dstdims[k];
    }
   
    for (k = 0; k < w; k ++) {
        int dk = linear_idx(dstdims, numdims, k);
        int offset = k * elsize;
        memcpy((char*)dst + dk * elsize, (const char*)src + k * elsize, elsize);
    }
}
/*====================*
 * S-function methods *
 *====================*/
/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *   Setup sizes of the various vectors.
 */
static void mdlInitializeSizes(SimStruct *S)
{

    DECL_AND_INIT_DIMSINFO(inputDimsInfo);
    DECL_AND_INIT_DIMSINFO(outputDimsInfo);
    ssSetNumSFcnParams(S, NPARAMS);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    ssSetArrayLayoutForCodeGen(S, SS_ROW_MAJOR);

    ssSetOperatingPointCompliance(S, USE_DEFAULT_OPERATING_POINT);

    ssSetNumContStates(S, NUM_CONT_STATES);
    ssSetNumDiscStates(S, NUM_DISC_STATES);


    if (!ssSetNumInputPorts(S, NUM_INPUTS)) return;
    /* Input Port 0 */
    inputDimsInfo.width = INPUT_0_WIDTH;
    ssSetInputPortDimensionInfo(S, 0, &inputDimsInfo);
    ssSetInputPortMatrixDimensions(S, 0, INPUT_0_WIDTH, INPUT_DIMS_0_COL);
    ssSetInputPortFrameData(S, 0, IN_0_FRAME_BASED);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 0, INPUT_0_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 0, INPUT_0_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/

    /* Input Port 1 */
    inputDimsInfo.width = INPUT_1_WIDTH;
    ssSetInputPortDimensionInfo(S, 1, &inputDimsInfo);
    ssSetInputPortMatrixDimensions(S, 1, INPUT_1_WIDTH, INPUT_DIMS_1_COL);
    ssSetInputPortFrameData(S, 1, IN_1_FRAME_BASED);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, INPUT_1_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 1, INPUT_1_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/

    /* Input Port 2 */
    ssSetInputPortWidth(S, 2, INPUT_2_WIDTH);
    ssSetInputPortDataType(S, 2, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 2, INPUT_2_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 2, INPUT_2_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 2, 1); /*direct input signal access*/

    /* Input Port 3 */
    ssSetInputPortWidth(S, 3, INPUT_3_WIDTH);
    ssSetInputPortDataType(S, 3, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 3, INPUT_3_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 3, INPUT_3_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 3, 1); /*direct input signal access*/

    /* Input Port 4 */
    ssSetInputPortWidth(S, 4, INPUT_4_WIDTH);
    ssSetInputPortDataType(S, 4, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 4, INPUT_4_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 4, INPUT_4_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 4, 1); /*direct input signal access*/

    /* Input Port 5 */
    ssSetInputPortWidth(S, 5, INPUT_5_WIDTH);
    ssSetInputPortDataType(S, 5, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 5, INPUT_5_COMPLEX);
    ssSetInputPortDirectFeedThrough(S, 5, INPUT_5_FEEDTHROUGH);
    ssSetInputPortRequiredContiguous(S, 5, 1); /*direct input signal access*/


    if (!ssSetNumOutputPorts(S, NUM_OUTPUTS)) return;
    /* Output Port 0 */
    outputDimsInfo.width = OUTPUT_0_WIDTH;
    ssSetOutputPortDimensionInfo(S, 0, &outputDimsInfo);
    ssSetOutputPortMatrixDimensions(S, 0, OUTPUT_0_WIDTH, OUTPUT_DIMS_0_COL);
    ssSetOutputPortFrameData(S, 0, OUT_0_FRAME_BASED);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, OUTPUT_0_COMPLEX);
    /* Output Port 1 */
    ssSetOutputPortWidth(S, 1, OUTPUT_1_WIDTH);
    ssSetOutputPortDataType(S, 1, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 1, OUTPUT_1_COMPLEX);
    /* Output Port 2 */
    ssSetOutputPortWidth(S, 2, OUTPUT_2_WIDTH);
    ssSetOutputPortDataType(S, 2, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 2, OUTPUT_2_COMPLEX);
    /* Output Port 3 */
    ssSetOutputPortWidth(S, 3, OUTPUT_3_WIDTH);
    ssSetOutputPortDataType(S, 3, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 3, OUTPUT_3_COMPLEX);
    if (!ssSetNumDWork(S, 10)) return;

    /*
     * Configure the dwork 0 (nodeArray_t)
     */
    ssSetDWorkDataType(S, 0, ssGetInputPortDataType(S, 0));
    ssSetDWorkUsageType(S, 0, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 0, "nodeArray_t");
    ssSetDWorkWidth(S, 0, ssGetInputPortWidth(S, 0));
    ssSetDWorkComplexSignal(S, 0, COMPLEX_NO);

    /*
     * Configure the dwork 1 (edgeArray_t)
     */
    ssSetDWorkDataType(S, 1, ssGetInputPortDataType(S, 1));
    ssSetDWorkUsageType(S, 1, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 1, "edgeArray_t");
    ssSetDWorkWidth(S, 1, ssGetInputPortWidth(S, 1));
    ssSetDWorkComplexSignal(S, 1, COMPLEX_NO);

    /*
     * Configure the dwork 2 (numNodesPtr_t)
     */
    ssSetDWorkDataType(S, 2, ssGetInputPortDataType(S, 2));
    ssSetDWorkUsageType(S, 2, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 2, "numNodesPtr_t");
    ssSetDWorkWidth(S, 2, ssGetInputPortWidth(S, 2));
    ssSetDWorkComplexSignal(S, 2, COMPLEX_NO);

    /*
     * Configure the dwork 3 (numEdgesPtr_t)
     */
    ssSetDWorkDataType(S, 3, ssGetInputPortDataType(S, 3));
    ssSetDWorkUsageType(S, 3, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 3, "numEdgesPtr_t");
    ssSetDWorkWidth(S, 3, ssGetInputPortWidth(S, 3));
    ssSetDWorkComplexSignal(S, 3, COMPLEX_NO);

    /*
     * Configure the dwork 4 (startNodePtr_t)
     */
    ssSetDWorkDataType(S, 4, ssGetInputPortDataType(S, 4));
    ssSetDWorkUsageType(S, 4, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 4, "startNodePtr_t");
    ssSetDWorkWidth(S, 4, ssGetInputPortWidth(S, 4));
    ssSetDWorkComplexSignal(S, 4, COMPLEX_NO);

    /*
     * Configure the dwork 5 (goalNodePtr_t)
     */
    ssSetDWorkDataType(S, 5, ssGetInputPortDataType(S, 5));
    ssSetDWorkUsageType(S, 5, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 5, "goalNodePtr_t");
    ssSetDWorkWidth(S, 5, ssGetInputPortWidth(S, 5));
    ssSetDWorkComplexSignal(S, 5, COMPLEX_NO);

    /*
     * Configure the dwork 6 (finalPath_t)
     */
    ssSetDWorkDataType(S, 6, ssGetOutputPortDataType(S, 0));
    ssSetDWorkUsageType(S, 6, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 6, "finalPath_t");
    ssSetDWorkWidth(S, 6, ssGetOutputPortWidth(S, 0));
    ssSetDWorkComplexSignal(S, 6, COMPLEX_NO);

    /*
     * Configure the dwork 7 (finalPathCost_t)
     */
    ssSetDWorkDataType(S, 7, ssGetOutputPortDataType(S, 1));
    ssSetDWorkUsageType(S, 7, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 7, "finalPathCost_t");
    ssSetDWorkWidth(S, 7, ssGetOutputPortWidth(S, 1));
    ssSetDWorkComplexSignal(S, 7, COMPLEX_NO);

    /*
     * Configure the dwork 8 (numNodesCheck_t)
     */
    ssSetDWorkDataType(S, 8, ssGetOutputPortDataType(S, 2));
    ssSetDWorkUsageType(S, 8, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 8, "numNodesCheck_t");
    ssSetDWorkWidth(S, 8, ssGetOutputPortWidth(S, 2));
    ssSetDWorkComplexSignal(S, 8, COMPLEX_NO);

    /*
     * Configure the dwork 9 (numEdgesCheck_t)
     */
    ssSetDWorkDataType(S, 9, ssGetOutputPortDataType(S, 3));
    ssSetDWorkUsageType(S, 9, SS_DWORK_USED_AS_SCRATCH);
    ssSetDWorkName(S, 9, "numEdgesCheck_t");
    ssSetDWorkWidth(S, 9, ssGetOutputPortWidth(S, 3));
    ssSetDWorkComplexSignal(S, 9, COMPLEX_NO);
    ssSetNumPWork(S, 0);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    ssSetSimulinkVersionGeneratedIn(S, "9.3");

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}


#define MDL_SET_INPUT_PORT_DIMENSION_INFO
void mdlSetInputPortDimensionInfo(SimStruct        *S, 
                                  int              portIndex, 
                                  const DimsInfo_T *dimsInfo)
{
    DECL_AND_INIT_DIMSINFO(portDimsInfo);
    int_T dims[2] = { OUTPUT_0_WIDTH, 1 };
    bool  frameIn = (ssGetInputPortFrameData(S, 0) == FRAME_YES);

    ssSetInputPortDimensionInfo(S, 0, dimsInfo);

    if (ssGetOutputPortNumDimensions(S, 0) == (-1)) {
        /* the output port has not been set */

        portDimsInfo.width   = OUTPUT_0_WIDTH;
        portDimsInfo.numDims = frameIn ? 2 : 1;
        portDimsInfo.dims    = frameIn ? dims : &portDimsInfo.width;

        ssSetOutputPortDimensionInfo(S, 0, &portDimsInfo);
    }
}


#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO
void mdlSetOutputPortDimensionInfo(SimStruct        *S,         
                                   int_T            portIndex,
                                   const DimsInfo_T *dimsInfo)
{
    DECL_AND_INIT_DIMSINFO(portDimsInfo);
    int_T dims[2] = { OUTPUT_0_WIDTH, 1 };
    bool  frameOut = (ssGetOutputPortFrameData(S, 0) == FRAME_YES);

    ssSetOutputPortDimensionInfo(S, 0, dimsInfo);

    if (ssGetInputPortNumDimensions(S, 0) == (-1)) {
      /* the input port has not been set */

        portDimsInfo.width   = INPUT_0_WIDTH;
        portDimsInfo.numDims = frameOut ? 2 : 1;
        portDimsInfo.dims    = frameOut ? dims : &portDimsInfo.width;

        ssSetInputPortDimensionInfo(S, 0, &portDimsInfo);
    }
}


#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
    DECL_AND_INIT_DIMSINFO(portDimsInfo);
    int_T dims[2] = { INPUT_0_WIDTH, 1 };
    bool  frameIn = ssGetInputPortFrameData(S, 0) == FRAME_YES;

    /* Neither the input nor the output ports have been set */

    portDimsInfo.width   = INPUT_0_WIDTH;
    portDimsInfo.numDims = frameIn ? 2 : 1;
    portDimsInfo.dims    = frameIn ? dims : &portDimsInfo.width;
    if (ssGetInputPortNumDimensions(S, 0) == (-1)) {  
        ssSetInputPortDimensionInfo(S, 0, &portDimsInfo);
    }
    portDimsInfo.width   = OUTPUT_0_WIDTH;
    dims[0]              = OUTPUT_0_WIDTH;
    if (ssGetOutputPortNumDimensions(S, 0) == (-1)) {  
        ssSetOutputPortDimensionInfo(S, 0, &portDimsInfo);
    }
    return;
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy  the sample time.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLE_TIME_0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_SET_INPUT_PORT_DATA_TYPE
static void mdlSetInputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetInputPortDataType(S, 0, dType);
}

#define MDL_SET_OUTPUT_PORT_DATA_TYPE
static void mdlSetOutputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetOutputPortDataType(S, 0, dType);
}

#define MDL_SET_DEFAULT_PORT_DATA_TYPES
static void mdlSetDefaultPortDataTypes(SimStruct *S)
{
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
}

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START)
/* Function: mdlStart =======================================================
 * Abstract:
 *    This function is called once at start of model execution. If you
 *    have states that should be initialized once, this is the place
 *    to do it.
 */
static void mdlStart(SimStruct *S)
{
}
#endif /*  MDL_START */

/* Function: mdlOutputs =======================================================
 *
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T *nodeArray = (real_T *) ssGetInputPortRealSignal(S, 0);
    const real_T *edgeArray = (real_T *) ssGetInputPortRealSignal(S, 1);
    const real_T *numNodesPtr = (real_T *) ssGetInputPortRealSignal(S, 2);
    const real_T *numEdgesPtr = (real_T *) ssGetInputPortRealSignal(S, 3);
    const real_T *startNodePtr = (real_T *) ssGetInputPortRealSignal(S, 4);
    const real_T *goalNodePtr = (real_T *) ssGetInputPortRealSignal(S, 5);
    real_T *finalPath = (real_T *) ssGetOutputPortRealSignal(S, 0);
    real_T *finalPathCost = (real_T *) ssGetOutputPortRealSignal(S, 1);
    real_T *numNodesCheck = (real_T *) ssGetOutputPortRealSignal(S, 2);
    real_T *numEdgesCheck = (real_T *) ssGetOutputPortRealSignal(S, 3);

    /* S-Function Builder Row Major Support has been enabled for custom
     * code, a transposed copy will be created for any array signals.
     */
    real_T *nodeArray_t = (real_T *)ssGetDWork(S, 0);
    real_T *edgeArray_t = (real_T *)ssGetDWork(S, 1);
    real_T *numNodesPtr_t = (real_T *)ssGetDWork(S, 2);
    real_T *numEdgesPtr_t = (real_T *)ssGetDWork(S, 3);
    real_T *startNodePtr_t = (real_T *)ssGetDWork(S, 4);
    real_T *goalNodePtr_t = (real_T *)ssGetDWork(S, 5);
    real_T *finalPath_t = (real_T *)ssGetDWork(S, 6);
    real_T *finalPathCost_t = (real_T *)ssGetDWork(S, 7);
    real_T *numNodesCheck_t = (real_T *)ssGetDWork(S, 8);
    real_T *numEdgesCheck_t = (real_T *)ssGetDWork(S, 9);

    NDTransposeBySrcSpecs((void*)nodeArray_t, (const void*)nodeArray, ssGetInputPortDimensions(S, 0), ssGetInputPortNumDimensions(S, 0), sizeof(real_T));
    NDTransposeBySrcSpecs((void*)edgeArray_t, (const void*)edgeArray, ssGetInputPortDimensions(S, 1), ssGetInputPortNumDimensions(S, 1), sizeof(real_T));
    NDTransposeBySrcSpecs((void*)numNodesPtr_t, (const void*)numNodesPtr, ssGetInputPortDimensions(S, 2), ssGetInputPortNumDimensions(S, 2), sizeof(real_T));
    NDTransposeBySrcSpecs((void*)numEdgesPtr_t, (const void*)numEdgesPtr, ssGetInputPortDimensions(S, 3), ssGetInputPortNumDimensions(S, 3), sizeof(real_T));
    NDTransposeBySrcSpecs((void*)startNodePtr_t, (const void*)startNodePtr, ssGetInputPortDimensions(S, 4), ssGetInputPortNumDimensions(S, 4), sizeof(real_T));
    NDTransposeBySrcSpecs((void*)goalNodePtr_t, (const void*)goalNodePtr, ssGetInputPortDimensions(S, 5), ssGetInputPortNumDimensions(S, 5), sizeof(real_T));
    A_Star_S_Func_2_Outputs_wrapper(nodeArray_t, edgeArray_t, numNodesPtr_t, numEdgesPtr_t, startNodePtr_t, goalNodePtr_t, finalPath_t, finalPathCost_t, numNodesCheck_t, numEdgesCheck_t);
    NDTransposeByDstSpecs((void*)finalPath, (const void*)finalPath_t, ssGetOutputPortDimensions(S, 0), ssGetOutputPortNumDimensions(S, 0), sizeof(real_T));
    NDTransposeByDstSpecs((void*)finalPathCost, (const void*)finalPathCost_t, ssGetOutputPortDimensions(S, 1), ssGetOutputPortNumDimensions(S, 1), sizeof(real_T));
    NDTransposeByDstSpecs((void*)numNodesCheck, (const void*)numNodesCheck_t, ssGetOutputPortDimensions(S, 2), ssGetOutputPortNumDimensions(S, 2), sizeof(real_T));
    NDTransposeByDstSpecs((void*)numEdgesCheck, (const void*)numEdgesCheck_t, ssGetOutputPortDimensions(S, 3), ssGetOutputPortNumDimensions(S, 3), sizeof(real_T));

}

/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{

}


#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif



