/*
 * MTDPC_data.c
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

#include "MTDPC.h"
#include "MTDPC_private.h"

/* Block parameters (default storage) */
P_MTDPC_T MTDPC_P = {
  /* Mask Parameter: PIDController_D
   * Referenced by: '<S35>/Derivative Gain'
   */
  1.0,

  /* Mask Parameter: Cs_D
   * Referenced by: '<S83>/Derivative Gain'
   */
  1.0,

  /* Mask Parameter: PIDController_InitialConditionF
   * Referenced by: '<S36>/Filter'
   */
  0.0,

  /* Mask Parameter: Cs_InitialConditionForFilter
   * Referenced by: '<S84>/Filter'
   */
  0.0,

  /* Mask Parameter: PIDController_N
   * Referenced by: '<S44>/Filter Coefficient'
   */
  1000.0,

  /* Mask Parameter: Cs_N
   * Referenced by: '<S92>/Filter Coefficient'
   */
  1000.0,

  /* Mask Parameter: PIDController_P
   * Referenced by: '<S46>/Proportional Gain'
   */
  1000.0,

  /* Mask Parameter: Cs_P
   * Referenced by: '<S94>/Proportional Gain'
   */
  1000.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay3'
   */
  0.0,

  /* Expression: 0.02
   * Referenced by: '<Root>/Sine Wave3'
   */
  0.02,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave3'
   */
  0.0,

  /* Expression: 12
   * Referenced by: '<Root>/Sine Wave3'
   */
  12.0,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave3'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/Variable Transport Delay1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Variable Transport Delay1'
   */
  0.0,

  /* Computed Parameter: uZm_A
   * Referenced by: '<Root>/1//Zm'
   */
  -10.0,

  /* Computed Parameter: uZm_C
   * Referenced by: '<Root>/1//Zm'
   */
  10.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay1'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay2'
   */
  0.0,

  /* Expression: 3
   * Referenced by: '<Root>/Ep'
   */
  3.0,

  /* Expression: 0
   * Referenced by: '<Root>/Er'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator2'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator3'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay4'
   */
  0.0,

  /* Computed Parameter: Ze_A
   * Referenced by: '<Root>/Ze'
   */
  -1000.0,

  /* Computed Parameter: Ze_C
   * Referenced by: '<Root>/Ze'
   */
  100000.0,

  /* Expression: 0
   * Referenced by: '<Root>/Fth+'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay7'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Unit Delay6'
   */
  0.0,

  /* Expression: 0.02
   * Referenced by: '<Root>/Sine Wave'
   */
  0.02,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave'
   */
  0.0,

  /* Expression: 12
   * Referenced by: '<Root>/Sine Wave'
   */
  12.0,

  /* Expression: 0
   * Referenced by: '<Root>/Sine Wave'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<Root>/Variable Transport Delay'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Variable Transport Delay'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S7>/Sine Wave1'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S7>/Sine Wave1'
   */
  0.0,

  /* Expression: 5
   * Referenced by: '<S7>/Sine Wave1'
   */
  5.0,

  /* Expression: 0
   * Referenced by: '<S7>/Sine Wave1'
   */
  0.0,

  /* Expression: 1
   * Referenced by: '<S7>/Sine Wave2'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<S7>/Sine Wave2'
   */
  0.0,

  /* Expression: 8
   * Referenced by: '<S7>/Sine Wave2'
   */
  8.0,

  /* Expression: pi/2
   * Referenced by: '<S7>/Sine Wave2'
   */
  1.5707963267948966,

  /* Computed Parameter: TransferFcn2_A
   * Referenced by: '<S8>/Transfer Fcn2'
   */
  -1000.0,

  /* Computed Parameter: TransferFcn2_C
   * Referenced by: '<S8>/Transfer Fcn2'
   */
  -1.0E+6,

  /* Computed Parameter: TransferFcn2_D
   * Referenced by: '<S8>/Transfer Fcn2'
   */
  1000.0,

  /* Computed Parameter: TransferFcn_A
   * Referenced by: '<S8>/Transfer Fcn'
   */
  -0.0,

  /* Computed Parameter: TransferFcn_C
   * Referenced by: '<S8>/Transfer Fcn'
   */
  5.0,

  /* Expression: 10
   * Referenced by: '<S8>/Gain1'
   */
  10.0,

  /* Computed Parameter: uZs_A
   * Referenced by: '<Root>/1//Zs'
   */
  -10.0,

  /* Computed Parameter: uZs_C
   * Referenced by: '<Root>/1//Zs'
   */
  10.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<Root>/Integrator1'
   */
  0.0,

  /* Computed Parameter: TransferFcn_A_e
   * Referenced by: '<S1>/Transfer Fcn'
   */
  -1000.0,

  /* Computed Parameter: TransferFcn_C_e
   * Referenced by: '<S1>/Transfer Fcn'
   */
  -99000.0,

  /* Computed Parameter: TransferFcn_D
   * Referenced by: '<S1>/Transfer Fcn'
   */
  100.0,

  /* Expression: 1
   * Referenced by: '<S2>/C2'
   */
  1.0,

  /* Expression: -1
   * Referenced by: '<Root>/C5 '
   */
  -1.0,

  /* Expression: 0
   * Referenced by: '<Root>/C6 '
   */
  0.0,

  /* Computed Parameter: Cm_A
   * Referenced by: '<Root>/Cm'
   */
  -1000.0,

  /* Computed Parameter: Cm_C
   * Referenced by: '<Root>/Cm'
   */
  99000.0,

  /* Computed Parameter: Cm_D
   * Referenced by: '<Root>/Cm'
   */
  -100.0,

  /* Expression: 0.1
   * Referenced by: '<Root>/Transport Delay'
   */
  0.1,

  /* Expression: 0
   * Referenced by: '<Root>/Transport Delay'
   */
  0.0,

  /* Expression: 0.1
   * Referenced by: '<Root>/Transport Delay1'
   */
  0.1,

  /* Expression: 0
   * Referenced by: '<Root>/Transport Delay1'
   */
  0.0
};
