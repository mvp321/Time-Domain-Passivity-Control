/*
 * MTDPC.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "MTDPC".
 *
 * Model version              : 1.16
 * Simulink Coder version : 9.3 (R2020a) 18-Nov-2019
 * C source code generated on : Sun Oct 11 18:37:30 2020
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_MTDPC_h_
#define RTW_HEADER_MTDPC_h_
#include <float.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#ifndef MTDPC_COMMON_INCLUDES_
# define MTDPC_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#endif                                 /* MTDPC_COMMON_INCLUDES_ */

#include "MTDPC_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T UnitDelay;                    /* '<Root>/Unit Delay' */
  real_T UnitDelay3;                   /* '<Root>/Unit Delay3' */
  real_T SineWave3;                    /* '<Root>/Sine Wave3' */
  real_T Fth;                          /* '<Root>/Variable Transport Delay1' */
  real_T Vp;                           /* '<Root>/1//Zm' */
  real_T Transpose1;                   /* '<Root>/Transpose1' */
  real_T UnitDelay1;                   /* '<Root>/Unit Delay1' */
  real_T UnitDelay2;                   /* '<Root>/Unit Delay2' */
  real_T E1;                           /* '<Root>/Integrator2' */
  real_T E2;                           /* '<Root>/Integrator3' */
  real_T UnitDelay5;                   /* '<Root>/Unit Delay5' */
  real_T UnitDelay4;                   /* '<Root>/Unit Delay4' */
  real_T Fth_b;                        /* '<Root>/Sum3' */
  real_T UnitDelay7;                   /* '<Root>/Unit Delay7' */
  real_T UnitDelay6;                   /* '<Root>/Unit Delay6' */
  real_T SineWave;                     /* '<Root>/Sine Wave' */
  real_T Fp;                           /* '<Root>/Sum' */
  real_T Vth;                          /* '<Root>/1//Zs' */
  real_T Pp;                           /* '<Root>/Integrator' */
  real_T Pth;                          /* '<Root>/Integrator1' */
  real_T FilterCoefficient;            /* '<S44>/Filter Coefficient' */
  real_T FilterCoefficient_j;          /* '<S92>/Filter Coefficient' */
  real_T Product;                      /* '<Root>/Product' */
  real_T Product1;                     /* '<Root>/Product1' */
  real_T Sum1;                         /* '<Root>/Sum1' */
  real_T Fth_l;                        /* '<Root>/Sum2' */
  real_T TransportDelay;               /* '<Root>/Transport Delay' */
  real_T TransportDelay1;              /* '<Root>/Transport Delay1' */
  real_T LOP;                          /* '<Root>/Energy Domain TDPC' */
  real_T alpha;                        /* '<Root>/Energy Domain TDPC' */
  real_T F_mod;                        /* '<Root>/Energy Domain TDPC' */
  real_T SOP;                       /* '<Root>/EOP of environment Calculator' */
  real_T EOP;                       /* '<Root>/EOP of environment Calculator' */
  real_T D;                         /* '<Root>/EOP of environment Calculator' */
  real_T N;                         /* '<Root>/EOP of environment Calculator' */
  real_T SOP_h;     /* '<Root>/EOP of Communication + Environment Calculator' */
  real_T EOP_k;     /* '<Root>/EOP of Communication + Environment Calculator' */
  real_T B;         /* '<Root>/EOP of Communication + Environment Calculator' */
  real_T A;         /* '<Root>/EOP of Communication + Environment Calculator' */
} B_MTDPC_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<Root>/Unit Delay' */
  real_T UnitDelay3_DSTATE;            /* '<Root>/Unit Delay3' */
  real_T UnitDelay1_DSTATE;            /* '<Root>/Unit Delay1' */
  real_T UnitDelay2_DSTATE;            /* '<Root>/Unit Delay2' */
  real_T UnitDelay5_DSTATE;            /* '<Root>/Unit Delay5' */
  real_T UnitDelay4_DSTATE;            /* '<Root>/Unit Delay4' */
  real_T UnitDelay7_DSTATE;            /* '<Root>/Unit Delay7' */
  real_T UnitDelay6_DSTATE;            /* '<Root>/Unit Delay6' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } VariableTransportDelay1_RWORK;     /* '<Root>/Variable Transport Delay1' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } VariableTransportDelay_RWORK;      /* '<Root>/Variable Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK;              /* '<Root>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay1_RWORK;             /* '<Root>/Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } VariableTransportDelay1_PWORK;     /* '<Root>/Variable Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } VariableTransportDelay_PWORK;      /* '<Root>/Variable Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<Root>/Transport Delay' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay1_PWORK;             /* '<Root>/Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } VariableTransportDelay1_IWORK;     /* '<Root>/Variable Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } VariableTransportDelay_IWORK;      /* '<Root>/Variable Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<Root>/Transport Delay' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay1_IWORK;             /* '<Root>/Transport Delay1' */
} DW_MTDPC_T;

/* Continuous states (default storage) */
typedef struct {
  real_T uZm_CSTATE;                   /* '<Root>/1//Zm' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T Ze_CSTATE;                    /* '<Root>/Ze' */
  real_T TransferFcn2_CSTATE;          /* '<S8>/Transfer Fcn2' */
  real_T TransferFcn_CSTATE;           /* '<S8>/Transfer Fcn' */
  real_T uZs_CSTATE;                   /* '<Root>/1//Zs' */
  real_T Integrator_CSTATE;            /* '<Root>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Filter_CSTATE;                /* '<S36>/Filter' */
  real_T TransferFcn_CSTATE_d;         /* '<S1>/Transfer Fcn' */
  real_T Cm_CSTATE;                    /* '<Root>/Cm' */
  real_T Filter_CSTATE_i;              /* '<S84>/Filter' */
} X_MTDPC_T;

/* State derivatives (default storage) */
typedef struct {
  real_T uZm_CSTATE;                   /* '<Root>/1//Zm' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T Ze_CSTATE;                    /* '<Root>/Ze' */
  real_T TransferFcn2_CSTATE;          /* '<S8>/Transfer Fcn2' */
  real_T TransferFcn_CSTATE;           /* '<S8>/Transfer Fcn' */
  real_T uZs_CSTATE;                   /* '<Root>/1//Zs' */
  real_T Integrator_CSTATE;            /* '<Root>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Filter_CSTATE;                /* '<S36>/Filter' */
  real_T TransferFcn_CSTATE_d;         /* '<S1>/Transfer Fcn' */
  real_T Cm_CSTATE;                    /* '<Root>/Cm' */
  real_T Filter_CSTATE_i;              /* '<S84>/Filter' */
} XDot_MTDPC_T;

/* State disabled  */
typedef struct {
  boolean_T uZm_CSTATE;                /* '<Root>/1//Zm' */
  boolean_T Integrator2_CSTATE;        /* '<Root>/Integrator2' */
  boolean_T Integrator3_CSTATE;        /* '<Root>/Integrator3' */
  boolean_T Ze_CSTATE;                 /* '<Root>/Ze' */
  boolean_T TransferFcn2_CSTATE;       /* '<S8>/Transfer Fcn2' */
  boolean_T TransferFcn_CSTATE;        /* '<S8>/Transfer Fcn' */
  boolean_T uZs_CSTATE;                /* '<Root>/1//Zs' */
  boolean_T Integrator_CSTATE;         /* '<Root>/Integrator' */
  boolean_T Integrator1_CSTATE;        /* '<Root>/Integrator1' */
  boolean_T Filter_CSTATE;             /* '<S36>/Filter' */
  boolean_T TransferFcn_CSTATE_d;      /* '<S1>/Transfer Fcn' */
  boolean_T Cm_CSTATE;                 /* '<Root>/Cm' */
  boolean_T Filter_CSTATE_i;           /* '<S84>/Filter' */
} XDis_MTDPC_T;

#ifndef ODE4_INTG
#define ODE4_INTG

/* ODE4 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[4];                        /* derivatives */
} ODE4_IntgData;

#endif

/* Parameters (default storage) */
struct P_MTDPC_T_ {
  real_T PIDController_D;              /* Mask Parameter: PIDController_D
                                        * Referenced by: '<S35>/Derivative Gain'
                                        */
  real_T Cs_D;                         /* Mask Parameter: Cs_D
                                        * Referenced by: '<S83>/Derivative Gain'
                                        */
  real_T PIDController_InitialConditionF;
                              /* Mask Parameter: PIDController_InitialConditionF
                               * Referenced by: '<S36>/Filter'
                               */
  real_T Cs_InitialConditionForFilter;
                                 /* Mask Parameter: Cs_InitialConditionForFilter
                                  * Referenced by: '<S84>/Filter'
                                  */
  real_T PIDController_N;              /* Mask Parameter: PIDController_N
                                        * Referenced by: '<S44>/Filter Coefficient'
                                        */
  real_T Cs_N;                         /* Mask Parameter: Cs_N
                                        * Referenced by: '<S92>/Filter Coefficient'
                                        */
  real_T PIDController_P;              /* Mask Parameter: PIDController_P
                                        * Referenced by: '<S46>/Proportional Gain'
                                        */
  real_T Cs_P;                         /* Mask Parameter: Cs_P
                                        * Referenced by: '<S94>/Proportional Gain'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay3'
                                        */
  real_T SineWave3_Amp;                /* Expression: 0.02
                                        * Referenced by: '<Root>/Sine Wave3'
                                        */
  real_T SineWave3_Bias;               /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave3'
                                        */
  real_T SineWave3_Freq;               /* Expression: 12
                                        * Referenced by: '<Root>/Sine Wave3'
                                        */
  real_T SineWave3_Phase;              /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave3'
                                        */
  real_T VariableTransportDelay1_MaxDela;/* Expression: 1
                                          * Referenced by: '<Root>/Variable Transport Delay1'
                                          */
  real_T VariableTransportDelay1_InitOut;/* Expression: 0
                                          * Referenced by: '<Root>/Variable Transport Delay1'
                                          */
  real_T uZm_A;                        /* Computed Parameter: uZm_A
                                        * Referenced by: '<Root>/1//Zm'
                                        */
  real_T uZm_C;                        /* Computed Parameter: uZm_C
                                        * Referenced by: '<Root>/1//Zm'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay1'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay2'
                                        */
  real_T Ep_Value;                     /* Expression: 3
                                        * Referenced by: '<Root>/Ep'
                                        */
  real_T Er_Value;                     /* Expression: 0
                                        * Referenced by: '<Root>/Er'
                                        */
  real_T Integrator2_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator2'
                                        */
  real_T Integrator3_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator3'
                                        */
  real_T UnitDelay5_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay5'
                                        */
  real_T UnitDelay4_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay4'
                                        */
  real_T Ze_A;                         /* Computed Parameter: Ze_A
                                        * Referenced by: '<Root>/Ze'
                                        */
  real_T Ze_C;                         /* Computed Parameter: Ze_C
                                        * Referenced by: '<Root>/Ze'
                                        */
  real_T Fth_Value;                    /* Expression: 0
                                        * Referenced by: '<Root>/Fth+'
                                        */
  real_T UnitDelay7_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay7'
                                        */
  real_T UnitDelay6_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<Root>/Unit Delay6'
                                        */
  real_T SineWave_Amp;                 /* Expression: 0.02
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Bias;                /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Freq;                /* Expression: 12
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T SineWave_Phase;               /* Expression: 0
                                        * Referenced by: '<Root>/Sine Wave'
                                        */
  real_T VariableTransportDelay_MaxDelay;/* Expression: 1
                                          * Referenced by: '<Root>/Variable Transport Delay'
                                          */
  real_T VariableTransportDelay_InitOutp;/* Expression: 0
                                          * Referenced by: '<Root>/Variable Transport Delay'
                                          */
  real_T SineWave1_Amp;                /* Expression: 1
                                        * Referenced by: '<S7>/Sine Wave1'
                                        */
  real_T SineWave1_Bias;               /* Expression: 0
                                        * Referenced by: '<S7>/Sine Wave1'
                                        */
  real_T SineWave1_Freq;               /* Expression: 5
                                        * Referenced by: '<S7>/Sine Wave1'
                                        */
  real_T SineWave1_Phase;              /* Expression: 0
                                        * Referenced by: '<S7>/Sine Wave1'
                                        */
  real_T SineWave2_Amp;                /* Expression: 1
                                        * Referenced by: '<S7>/Sine Wave2'
                                        */
  real_T SineWave2_Bias;               /* Expression: 0
                                        * Referenced by: '<S7>/Sine Wave2'
                                        */
  real_T SineWave2_Freq;               /* Expression: 8
                                        * Referenced by: '<S7>/Sine Wave2'
                                        */
  real_T SineWave2_Phase;              /* Expression: pi/2
                                        * Referenced by: '<S7>/Sine Wave2'
                                        */
  real_T TransferFcn2_A;               /* Computed Parameter: TransferFcn2_A
                                        * Referenced by: '<S8>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_C;               /* Computed Parameter: TransferFcn2_C
                                        * Referenced by: '<S8>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_D;               /* Computed Parameter: TransferFcn2_D
                                        * Referenced by: '<S8>/Transfer Fcn2'
                                        */
  real_T TransferFcn_A;                /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S8>/Transfer Fcn'
                                        */
  real_T TransferFcn_C;                /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S8>/Transfer Fcn'
                                        */
  real_T Gain1_Gain;                   /* Expression: 10
                                        * Referenced by: '<S8>/Gain1'
                                        */
  real_T uZs_A;                        /* Computed Parameter: uZs_A
                                        * Referenced by: '<Root>/1//Zs'
                                        */
  real_T uZs_C;                        /* Computed Parameter: uZs_C
                                        * Referenced by: '<Root>/1//Zs'
                                        */
  real_T Integrator_IC;                /* Expression: 0
                                        * Referenced by: '<Root>/Integrator'
                                        */
  real_T Integrator1_IC;               /* Expression: 0
                                        * Referenced by: '<Root>/Integrator1'
                                        */
  real_T TransferFcn_A_e;              /* Computed Parameter: TransferFcn_A_e
                                        * Referenced by: '<S1>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_e;              /* Computed Parameter: TransferFcn_C_e
                                        * Referenced by: '<S1>/Transfer Fcn'
                                        */
  real_T TransferFcn_D;                /* Computed Parameter: TransferFcn_D
                                        * Referenced by: '<S1>/Transfer Fcn'
                                        */
  real_T C2_Gain;                      /* Expression: 1
                                        * Referenced by: '<S2>/C2'
                                        */
  real_T C5_Gain;                      /* Expression: -1
                                        * Referenced by: '<Root>/C5 '
                                        */
  real_T C6_Gain;                      /* Expression: 0
                                        * Referenced by: '<Root>/C6 '
                                        */
  real_T Cm_A;                         /* Computed Parameter: Cm_A
                                        * Referenced by: '<Root>/Cm'
                                        */
  real_T Cm_C;                         /* Computed Parameter: Cm_C
                                        * Referenced by: '<Root>/Cm'
                                        */
  real_T Cm_D;                         /* Computed Parameter: Cm_D
                                        * Referenced by: '<Root>/Cm'
                                        */
  real_T TransportDelay_Delay;         /* Expression: 0.1
                                        * Referenced by: '<Root>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<Root>/Transport Delay'
                                        */
  real_T TransportDelay1_Delay;        /* Expression: 0.1
                                        * Referenced by: '<Root>/Transport Delay1'
                                        */
  real_T TransportDelay1_InitOutput;   /* Expression: 0
                                        * Referenced by: '<Root>/Transport Delay1'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_MTDPC_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;
  RTWSolverInfo solverInfo;
  X_MTDPC_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[13];
  real_T odeF[4][13];
  ODE4_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T tFinal;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block parameters (default storage) */
extern P_MTDPC_T MTDPC_P;

/* Block signals (default storage) */
extern B_MTDPC_T MTDPC_B;

/* Continuous states (default storage) */
extern X_MTDPC_T MTDPC_X;

/* Block states (default storage) */
extern DW_MTDPC_T MTDPC_DW;

/* Model entry point functions */
extern void MTDPC_initialize(void);
extern void MTDPC_step(void);
extern void MTDPC_terminate(void);

/* Real-time Model object */
extern RT_MODEL_MTDPC_T *const MTDPC_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'MTDPC'
 * '<S1>'   : 'MTDPC/C1'
 * '<S2>'   : 'MTDPC/C2'
 * '<S3>'   : 'MTDPC/Cs'
 * '<S4>'   : 'MTDPC/EOP of Communication + Environment Calculator'
 * '<S5>'   : 'MTDPC/EOP of environment Calculator'
 * '<S6>'   : 'MTDPC/Energy Domain TDPC'
 * '<S7>'   : 'MTDPC/Fh+'
 * '<S8>'   : 'MTDPC/Zh'
 * '<S9>'   : 'MTDPC/C1/PID Controller'
 * '<S10>'  : 'MTDPC/C1/PID Controller/Anti-windup'
 * '<S11>'  : 'MTDPC/C1/PID Controller/D Gain'
 * '<S12>'  : 'MTDPC/C1/PID Controller/Filter'
 * '<S13>'  : 'MTDPC/C1/PID Controller/Filter ICs'
 * '<S14>'  : 'MTDPC/C1/PID Controller/I Gain'
 * '<S15>'  : 'MTDPC/C1/PID Controller/Ideal P Gain'
 * '<S16>'  : 'MTDPC/C1/PID Controller/Ideal P Gain Fdbk'
 * '<S17>'  : 'MTDPC/C1/PID Controller/Integrator'
 * '<S18>'  : 'MTDPC/C1/PID Controller/Integrator ICs'
 * '<S19>'  : 'MTDPC/C1/PID Controller/N Copy'
 * '<S20>'  : 'MTDPC/C1/PID Controller/N Gain'
 * '<S21>'  : 'MTDPC/C1/PID Controller/P Copy'
 * '<S22>'  : 'MTDPC/C1/PID Controller/Parallel P Gain'
 * '<S23>'  : 'MTDPC/C1/PID Controller/Reset Signal'
 * '<S24>'  : 'MTDPC/C1/PID Controller/Saturation'
 * '<S25>'  : 'MTDPC/C1/PID Controller/Saturation Fdbk'
 * '<S26>'  : 'MTDPC/C1/PID Controller/Sum'
 * '<S27>'  : 'MTDPC/C1/PID Controller/Sum Fdbk'
 * '<S28>'  : 'MTDPC/C1/PID Controller/Tracking Mode'
 * '<S29>'  : 'MTDPC/C1/PID Controller/Tracking Mode Sum'
 * '<S30>'  : 'MTDPC/C1/PID Controller/Tsamp - Integral'
 * '<S31>'  : 'MTDPC/C1/PID Controller/Tsamp - Ngain'
 * '<S32>'  : 'MTDPC/C1/PID Controller/postSat Signal'
 * '<S33>'  : 'MTDPC/C1/PID Controller/preSat Signal'
 * '<S34>'  : 'MTDPC/C1/PID Controller/Anti-windup/Disabled'
 * '<S35>'  : 'MTDPC/C1/PID Controller/D Gain/Internal Parameters'
 * '<S36>'  : 'MTDPC/C1/PID Controller/Filter/Cont. Filter'
 * '<S37>'  : 'MTDPC/C1/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S38>'  : 'MTDPC/C1/PID Controller/I Gain/Disabled'
 * '<S39>'  : 'MTDPC/C1/PID Controller/Ideal P Gain/Passthrough'
 * '<S40>'  : 'MTDPC/C1/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S41>'  : 'MTDPC/C1/PID Controller/Integrator/Disabled'
 * '<S42>'  : 'MTDPC/C1/PID Controller/Integrator ICs/Disabled'
 * '<S43>'  : 'MTDPC/C1/PID Controller/N Copy/Disabled'
 * '<S44>'  : 'MTDPC/C1/PID Controller/N Gain/Internal Parameters'
 * '<S45>'  : 'MTDPC/C1/PID Controller/P Copy/Disabled'
 * '<S46>'  : 'MTDPC/C1/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S47>'  : 'MTDPC/C1/PID Controller/Reset Signal/Disabled'
 * '<S48>'  : 'MTDPC/C1/PID Controller/Saturation/Passthrough'
 * '<S49>'  : 'MTDPC/C1/PID Controller/Saturation Fdbk/Disabled'
 * '<S50>'  : 'MTDPC/C1/PID Controller/Sum/Sum_PD'
 * '<S51>'  : 'MTDPC/C1/PID Controller/Sum Fdbk/Disabled'
 * '<S52>'  : 'MTDPC/C1/PID Controller/Tracking Mode/Disabled'
 * '<S53>'  : 'MTDPC/C1/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S54>'  : 'MTDPC/C1/PID Controller/Tsamp - Integral/Disabled wSignal Specification'
 * '<S55>'  : 'MTDPC/C1/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S56>'  : 'MTDPC/C1/PID Controller/postSat Signal/Forward_Path'
 * '<S57>'  : 'MTDPC/C1/PID Controller/preSat Signal/Forward_Path'
 * '<S58>'  : 'MTDPC/Cs/Anti-windup'
 * '<S59>'  : 'MTDPC/Cs/D Gain'
 * '<S60>'  : 'MTDPC/Cs/Filter'
 * '<S61>'  : 'MTDPC/Cs/Filter ICs'
 * '<S62>'  : 'MTDPC/Cs/I Gain'
 * '<S63>'  : 'MTDPC/Cs/Ideal P Gain'
 * '<S64>'  : 'MTDPC/Cs/Ideal P Gain Fdbk'
 * '<S65>'  : 'MTDPC/Cs/Integrator'
 * '<S66>'  : 'MTDPC/Cs/Integrator ICs'
 * '<S67>'  : 'MTDPC/Cs/N Copy'
 * '<S68>'  : 'MTDPC/Cs/N Gain'
 * '<S69>'  : 'MTDPC/Cs/P Copy'
 * '<S70>'  : 'MTDPC/Cs/Parallel P Gain'
 * '<S71>'  : 'MTDPC/Cs/Reset Signal'
 * '<S72>'  : 'MTDPC/Cs/Saturation'
 * '<S73>'  : 'MTDPC/Cs/Saturation Fdbk'
 * '<S74>'  : 'MTDPC/Cs/Sum'
 * '<S75>'  : 'MTDPC/Cs/Sum Fdbk'
 * '<S76>'  : 'MTDPC/Cs/Tracking Mode'
 * '<S77>'  : 'MTDPC/Cs/Tracking Mode Sum'
 * '<S78>'  : 'MTDPC/Cs/Tsamp - Integral'
 * '<S79>'  : 'MTDPC/Cs/Tsamp - Ngain'
 * '<S80>'  : 'MTDPC/Cs/postSat Signal'
 * '<S81>'  : 'MTDPC/Cs/preSat Signal'
 * '<S82>'  : 'MTDPC/Cs/Anti-windup/Disabled'
 * '<S83>'  : 'MTDPC/Cs/D Gain/Internal Parameters'
 * '<S84>'  : 'MTDPC/Cs/Filter/Cont. Filter'
 * '<S85>'  : 'MTDPC/Cs/Filter ICs/Internal IC - Filter'
 * '<S86>'  : 'MTDPC/Cs/I Gain/Disabled'
 * '<S87>'  : 'MTDPC/Cs/Ideal P Gain/Passthrough'
 * '<S88>'  : 'MTDPC/Cs/Ideal P Gain Fdbk/Disabled'
 * '<S89>'  : 'MTDPC/Cs/Integrator/Disabled'
 * '<S90>'  : 'MTDPC/Cs/Integrator ICs/Disabled'
 * '<S91>'  : 'MTDPC/Cs/N Copy/Disabled'
 * '<S92>'  : 'MTDPC/Cs/N Gain/Internal Parameters'
 * '<S93>'  : 'MTDPC/Cs/P Copy/Disabled'
 * '<S94>'  : 'MTDPC/Cs/Parallel P Gain/Internal Parameters'
 * '<S95>'  : 'MTDPC/Cs/Reset Signal/Disabled'
 * '<S96>'  : 'MTDPC/Cs/Saturation/Passthrough'
 * '<S97>'  : 'MTDPC/Cs/Saturation Fdbk/Disabled'
 * '<S98>'  : 'MTDPC/Cs/Sum/Sum_PD'
 * '<S99>'  : 'MTDPC/Cs/Sum Fdbk/Disabled'
 * '<S100>' : 'MTDPC/Cs/Tracking Mode/Disabled'
 * '<S101>' : 'MTDPC/Cs/Tracking Mode Sum/Passthrough'
 * '<S102>' : 'MTDPC/Cs/Tsamp - Integral/Disabled wSignal Specification'
 * '<S103>' : 'MTDPC/Cs/Tsamp - Ngain/Passthrough'
 * '<S104>' : 'MTDPC/Cs/postSat Signal/Forward_Path'
 * '<S105>' : 'MTDPC/Cs/preSat Signal/Forward_Path'
 */
#endif                                 /* RTW_HEADER_MTDPC_h_ */
