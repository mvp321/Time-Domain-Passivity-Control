/*
 * MTDPC.c
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

/* Block signals (default storage) */
B_MTDPC_T MTDPC_B;

/* Continuous states */
X_MTDPC_T MTDPC_X;

/* Block states (default storage) */
DW_MTDPC_T MTDPC_DW;

/* Real-time model */
RT_MODEL_MTDPC_T MTDPC_M_;
RT_MODEL_MTDPC_T *const MTDPC_M = &MTDPC_M_;

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                 /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * fabs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

/*
 * This function updates continuous states using the ODE4 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE4_IntgData *id = (ODE4_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T temp;
  int_T i;
  int_T nXc = 13;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  MTDPC_derivatives();

  /* f1 = f(t + (h/2), y + (h/2)*f0) */
  temp = 0.5 * h;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f0[i]);
  }

  rtsiSetT(si, t + temp);
  rtsiSetdX(si, f1);
  MTDPC_step();
  MTDPC_derivatives();

  /* f2 = f(t + (h/2), y + (h/2)*f1) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (temp*f1[i]);
  }

  rtsiSetdX(si, f2);
  MTDPC_step();
  MTDPC_derivatives();

  /* f3 = f(t + h, y + h*f2) */
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (h*f2[i]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f3);
  MTDPC_step();
  MTDPC_derivatives();

  /* tnew = t + h
     ynew = y + (h/6)*(f0 + 2*f1 + 2*f2 + 2*f3) */
  temp = h / 6.0;
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + temp*(f0[i] + 2.0*f1[i] + 2.0*f2[i] + f3[i]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model step function */
void MTDPC_step(void)
{
  /* local block i/o variables */
  real_T rtb_VariableTransportDelay;
  real_T rtb_Sum_o;
  real_T rtb_C2;
  real_T LOP;
  real_T alpha;
  real_T rtb_Fp;
  real_T SineWave3_tmp;
  real_T LOP_tmp;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* set solver stop time */
    if (!(MTDPC_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&MTDPC_M->solverInfo, ((MTDPC_M->Timing.clockTickH0
        + 1) * MTDPC_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&MTDPC_M->solverInfo, ((MTDPC_M->Timing.clockTick0 +
        1) * MTDPC_M->Timing.stepSize0 + MTDPC_M->Timing.clockTickH0 *
        MTDPC_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(MTDPC_M)) {
    MTDPC_M->Timing.t[0] = rtsiGetT(&MTDPC_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* UnitDelay: '<Root>/Unit Delay' */
    MTDPC_B.UnitDelay = MTDPC_DW.UnitDelay_DSTATE;

    /* UnitDelay: '<Root>/Unit Delay3' */
    MTDPC_B.UnitDelay3 = MTDPC_DW.UnitDelay3_DSTATE;
  }

  /* Sin: '<Root>/Sine Wave3' incorporates:
   *  Sin: '<Root>/Sine Wave'
   *  Sin: '<S7>/Sine Wave1'
   *  Sin: '<S7>/Sine Wave2'
   */
  SineWave3_tmp = MTDPC_M->Timing.t[0];
  MTDPC_B.SineWave3 = sin(MTDPC_P.SineWave3_Freq * SineWave3_tmp +
    MTDPC_P.SineWave3_Phase) * MTDPC_P.SineWave3_Amp + MTDPC_P.SineWave3_Bias;

  /* VariableTransportDelay: '<Root>/Variable Transport Delay1' */
  {
    real_T **uBuffer = (real_T**)
      &MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[1];
    real_T simTime = MTDPC_M->Timing.t[0];
    real_T appliedDelay;
    appliedDelay = MTDPC_B.SineWave3;

    /* For variable time delay, output here */
    if (appliedDelay > MTDPC_P.VariableTransportDelay1_MaxDela) {
      appliedDelay = MTDPC_P.VariableTransportDelay1_MaxDela;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    MTDPC_B.Fth = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      MTDPC_DW.VariableTransportDelay1_IWORK.CircularBufSize,
      &MTDPC_DW.VariableTransportDelay1_IWORK.Last,
      MTDPC_DW.VariableTransportDelay1_IWORK.Tail,
      MTDPC_DW.VariableTransportDelay1_IWORK.Head,
      MTDPC_P.VariableTransportDelay1_InitOut,
      0,
      0);
  }

  /* TransferFcn: '<Root>/1//Zm' */
  MTDPC_B.Vp = 0.0;
  MTDPC_B.Vp += MTDPC_P.uZm_C * MTDPC_X.uZm_CSTATE;

  /* Math: '<Root>/Transpose1' */
  MTDPC_B.Transpose1 = MTDPC_B.Vp;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* UnitDelay: '<Root>/Unit Delay1' */
    MTDPC_B.UnitDelay1 = MTDPC_DW.UnitDelay1_DSTATE;

    /* UnitDelay: '<Root>/Unit Delay2' */
    MTDPC_B.UnitDelay2 = MTDPC_DW.UnitDelay2_DSTATE;
  }

  /* MATLAB Function: '<Root>/Energy Domain TDPC' incorporates:
   *  Constant: '<Root>/Ep'
   *  Constant: '<Root>/Er'
   *  Math: '<Root>/Transpose2'
   *  Product: '<Root>/Product'
   */
  LOP_tmp = MTDPC_B.Fth * MTDPC_B.Vp;
  LOP = (((MTDPC_P.Ep_Value - MTDPC_P.Er_Value) * MTDPC_B.Transpose1 *
          MTDPC_B.Vp + LOP_tmp) + MTDPC_B.UnitDelay3 * MTDPC_B.UnitDelay2 *
         MTDPC_B.UnitDelay1) * 0.001 + MTDPC_B.UnitDelay;
  if (LOP >= 0.0) {
    alpha = 0.0;
  } else {
    alpha = -LOP / (0.001 * MTDPC_B.Transpose1 * MTDPC_B.Vp + 1.0E-5);
  }

  MTDPC_B.F_mod = (alpha * MTDPC_B.Vp + MTDPC_B.Fth) + 1.0E-5;
  MTDPC_B.LOP = LOP;
  MTDPC_B.alpha = alpha;

  /* End of MATLAB Function: '<Root>/Energy Domain TDPC' */
  if (rtmIsMajorTimeStep(MTDPC_M)) {
  }

  /* Integrator: '<Root>/Integrator2' */
  MTDPC_B.E1 = MTDPC_X.Integrator2_CSTATE;

  /* Integrator: '<Root>/Integrator3' */
  MTDPC_B.E2 = MTDPC_X.Integrator3_CSTATE;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* UnitDelay: '<Root>/Unit Delay5' */
    MTDPC_B.UnitDelay5 = MTDPC_DW.UnitDelay5_DSTATE;

    /* UnitDelay: '<Root>/Unit Delay4' */
    MTDPC_B.UnitDelay4 = MTDPC_DW.UnitDelay4_DSTATE;
  }

  /* Sum: '<Root>/Sum3' incorporates:
   *  Constant: '<Root>/Fth+'
   *  TransferFcn: '<Root>/Ze'
   */
  MTDPC_B.Fth_b = MTDPC_P.Ze_C * MTDPC_X.Ze_CSTATE + MTDPC_P.Fth_Value;

  /* MATLAB Function: '<Root>/EOP of environment Calculator' incorporates:
   *  Math: '<Root>/Transpose3'
   *  Math: '<Root>/Transpose4'
   */
  alpha = MTDPC_B.Fth_b * MTDPC_B.Vp * 0.001 + MTDPC_B.UnitDelay4;
  rtb_Fp = MTDPC_B.Vp * MTDPC_B.Vp * 0.001 + MTDPC_B.UnitDelay5;
  LOP = alpha / rtb_Fp;
  MTDPC_B.EOP = -LOP;
  MTDPC_B.SOP = LOP;
  MTDPC_B.D = rtb_Fp;
  MTDPC_B.N = alpha;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* UnitDelay: '<Root>/Unit Delay7' */
    MTDPC_B.UnitDelay7 = MTDPC_DW.UnitDelay7_DSTATE;

    /* UnitDelay: '<Root>/Unit Delay6' */
    MTDPC_B.UnitDelay6 = MTDPC_DW.UnitDelay6_DSTATE;
  }

  /* Sin: '<Root>/Sine Wave' */
  MTDPC_B.SineWave = sin(MTDPC_P.SineWave_Freq * SineWave3_tmp +
    MTDPC_P.SineWave_Phase) * MTDPC_P.SineWave_Amp + MTDPC_P.SineWave_Bias;

  /* VariableTransportDelay: '<Root>/Variable Transport Delay' */
  {
    real_T **uBuffer = (real_T**)
      &MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)
      &MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = MTDPC_M->Timing.t[0];
    real_T appliedDelay;
    appliedDelay = MTDPC_B.SineWave;

    /* For variable time delay, output here */
    if (appliedDelay > MTDPC_P.VariableTransportDelay_MaxDelay) {
      appliedDelay = MTDPC_P.VariableTransportDelay_MaxDelay;
    }

    if (appliedDelay < 0.0) {
      /* negative delay is not supported
       *  set delay to 0
       */
      appliedDelay = 0.0;
    }

    rtb_VariableTransportDelay = rt_TDelayInterpolate(
      simTime - appliedDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      MTDPC_DW.VariableTransportDelay_IWORK.CircularBufSize,
      &MTDPC_DW.VariableTransportDelay_IWORK.Last,
      MTDPC_DW.VariableTransportDelay_IWORK.Tail,
      MTDPC_DW.VariableTransportDelay_IWORK.Head,
      MTDPC_P.VariableTransportDelay_InitOutp,
      0,
      0);
  }

  /* MATLAB Function: '<Root>/EOP of Communication + Environment Calculator' incorporates:
   *  Math: '<Root>/Transpose5'
   *  Math: '<Root>/Transpose6'
   */
  alpha = MTDPC_B.Fth * rtb_VariableTransportDelay * 0.001 + MTDPC_B.UnitDelay6;
  rtb_Fp = rtb_VariableTransportDelay * rtb_VariableTransportDelay * 0.001 +
    MTDPC_B.UnitDelay7;
  LOP = alpha / rtb_Fp;
  MTDPC_B.EOP_k = -LOP;
  MTDPC_B.SOP_h = LOP;
  MTDPC_B.B = rtb_Fp;
  MTDPC_B.A = alpha;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
  }

  /* Sum: '<Root>/Sum' incorporates:
   *  Gain: '<S8>/Gain1'
   *  Sin: '<S7>/Sine Wave1'
   *  Sin: '<S7>/Sine Wave2'
   *  Sum: '<S8>/Sum4'
   *  TransferFcn: '<S8>/Transfer Fcn'
   *  TransferFcn: '<S8>/Transfer Fcn2'
   */
  MTDPC_B.Fp = ((sin(MTDPC_P.SineWave1_Freq * SineWave3_tmp +
                     MTDPC_P.SineWave1_Phase) * MTDPC_P.SineWave1_Amp +
                 MTDPC_P.SineWave1_Bias) + (sin(MTDPC_P.SineWave2_Freq *
    SineWave3_tmp + MTDPC_P.SineWave2_Phase) * MTDPC_P.SineWave2_Amp +
    MTDPC_P.SineWave2_Bias)) - (((MTDPC_P.TransferFcn2_C *
    MTDPC_X.TransferFcn2_CSTATE + MTDPC_P.TransferFcn2_D * MTDPC_B.Vp) +
    MTDPC_P.TransferFcn_C * MTDPC_X.TransferFcn_CSTATE) + MTDPC_P.Gain1_Gain *
    MTDPC_B.Vp);

  /* TransferFcn: '<Root>/1//Zs' */
  MTDPC_B.Vth = 0.0;
  MTDPC_B.Vth += MTDPC_P.uZs_C * MTDPC_X.uZs_CSTATE;

  /* Integrator: '<Root>/Integrator' */
  MTDPC_B.Pp = MTDPC_X.Integrator_CSTATE;

  /* Integrator: '<Root>/Integrator1' */
  MTDPC_B.Pth = MTDPC_X.Integrator1_CSTATE;
  if (rtmIsMajorTimeStep(MTDPC_M)) {
  }

  /* Gain: '<S44>/Filter Coefficient' incorporates:
   *  Gain: '<S35>/Derivative Gain'
   *  Integrator: '<S36>/Filter'
   *  Sum: '<S36>/SumD'
   */
  MTDPC_B.FilterCoefficient = (MTDPC_P.PIDController_D * MTDPC_B.Vp -
    MTDPC_X.Filter_CSTATE) * MTDPC_P.PIDController_N;

  /* Sum: '<S1>/Sum' incorporates:
   *  Gain: '<S46>/Proportional Gain'
   *  Sum: '<S50>/Sum'
   *  TransferFcn: '<S1>/Transfer Fcn'
   */
  rtb_Sum_o = (MTDPC_P.TransferFcn_C_e * MTDPC_X.TransferFcn_CSTATE_d +
               MTDPC_P.TransferFcn_D * MTDPC_B.Vp) + (MTDPC_P.PIDController_P *
    MTDPC_B.Vp + MTDPC_B.FilterCoefficient);

  /* Gain: '<S2>/C2' */
  rtb_C2 = MTDPC_P.C2_Gain * MTDPC_B.Fth_b;

  /* Gain: '<S92>/Filter Coefficient' incorporates:
   *  Gain: '<S83>/Derivative Gain'
   *  Integrator: '<S84>/Filter'
   *  Sum: '<S84>/SumD'
   */
  MTDPC_B.FilterCoefficient_j = (MTDPC_P.Cs_D * MTDPC_B.Vth -
    MTDPC_X.Filter_CSTATE_i) * MTDPC_P.Cs_N;

  /* Product: '<Root>/Product' */
  MTDPC_B.Product = LOP_tmp;

  /* Product: '<Root>/Product1' */
  MTDPC_B.Product1 = MTDPC_B.F_mod * MTDPC_B.Vp;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Gain: '<Root>/C6 '
   *  TransferFcn: '<Root>/Cm'
   */
  MTDPC_B.Sum1 = ((MTDPC_P.C6_Gain * MTDPC_B.Fp - (MTDPC_P.Cm_C *
    MTDPC_X.Cm_CSTATE + MTDPC_P.Cm_D * MTDPC_B.Vp)) + MTDPC_B.Fp) -
    MTDPC_B.F_mod;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Gain: '<Root>/C5 '
   *  Gain: '<S94>/Proportional Gain'
   *  Sum: '<S98>/Sum'
   */
  MTDPC_B.Fth_l = ((rtb_VariableTransportDelay - MTDPC_B.Fth_b) -
                   MTDPC_P.C5_Gain * MTDPC_B.Fth_b) - (MTDPC_P.Cs_P *
    MTDPC_B.Vth + MTDPC_B.FilterCoefficient_j);

  /* TransportDelay: '<Root>/Transport Delay' */
  {
    real_T **uBuffer = (real_T**)&MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[1];
    real_T simTime = MTDPC_M->Timing.t[0];
    real_T tMinusDelay = simTime - MTDPC_P.TransportDelay_Delay;
    MTDPC_B.TransportDelay = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      MTDPC_DW.TransportDelay_IWORK.CircularBufSize,
      &MTDPC_DW.TransportDelay_IWORK.Last,
      MTDPC_DW.TransportDelay_IWORK.Tail,
      MTDPC_DW.TransportDelay_IWORK.Head,
      MTDPC_P.TransportDelay_InitOutput,
      0,
      0);
  }

  /* TransportDelay: '<Root>/Transport Delay1' */
  {
    real_T **uBuffer = (real_T**)&MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs[0];
    real_T **tBuffer = (real_T**)&MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs[1];
    real_T simTime = MTDPC_M->Timing.t[0];
    real_T tMinusDelay = simTime - MTDPC_P.TransportDelay1_Delay;
    MTDPC_B.TransportDelay1 = rt_TDelayInterpolate(
      tMinusDelay,
      0.0,
      *tBuffer,
      *uBuffer,
      MTDPC_DW.TransportDelay1_IWORK.CircularBufSize,
      &MTDPC_DW.TransportDelay1_IWORK.Last,
      MTDPC_DW.TransportDelay1_IWORK.Tail,
      MTDPC_DW.TransportDelay1_IWORK.Head,
      MTDPC_P.TransportDelay1_InitOutput,
      0,
      0);
  }

  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* Matfile logging */
    rt_UpdateTXYLogVars(MTDPC_M->rtwLogInfo, (MTDPC_M->Timing.t));
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(MTDPC_M)) {
    if (rtmIsMajorTimeStep(MTDPC_M)) {
      /* Update for UnitDelay: '<Root>/Unit Delay' */
      MTDPC_DW.UnitDelay_DSTATE = MTDPC_B.LOP;

      /* Update for UnitDelay: '<Root>/Unit Delay3' */
      MTDPC_DW.UnitDelay3_DSTATE = MTDPC_B.alpha;
    }

    /* Update for VariableTransportDelay: '<Root>/Variable Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)
        &MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[1];
      real_T simTime = MTDPC_M->Timing.t[0];
      MTDPC_DW.VariableTransportDelay1_IWORK.Head =
        ((MTDPC_DW.VariableTransportDelay1_IWORK.Head <
          (MTDPC_DW.VariableTransportDelay1_IWORK.CircularBufSize-1)) ?
         (MTDPC_DW.VariableTransportDelay1_IWORK.Head+1) : 0);
      if (MTDPC_DW.VariableTransportDelay1_IWORK.Head ==
          MTDPC_DW.VariableTransportDelay1_IWORK.Tail) {
        MTDPC_DW.VariableTransportDelay1_IWORK.Tail =
          ((MTDPC_DW.VariableTransportDelay1_IWORK.Tail <
            (MTDPC_DW.VariableTransportDelay1_IWORK.CircularBufSize-1)) ?
           (MTDPC_DW.VariableTransportDelay1_IWORK.Tail+1) : 0);
      }

      (*tBuffer)[MTDPC_DW.VariableTransportDelay1_IWORK.Head] = simTime;
      (*uBuffer)[MTDPC_DW.VariableTransportDelay1_IWORK.Head] =
        MTDPC_B.TransportDelay1;

      /* when use fixed buffer, reset solver at when buffer is updated
       * to avoid output consistency fail.
       */
    }

    if (rtmIsMajorTimeStep(MTDPC_M)) {
      /* Update for UnitDelay: '<Root>/Unit Delay1' */
      MTDPC_DW.UnitDelay1_DSTATE = MTDPC_B.Vp;

      /* Update for UnitDelay: '<Root>/Unit Delay2' */
      MTDPC_DW.UnitDelay2_DSTATE = MTDPC_B.Transpose1;

      /* Update for UnitDelay: '<Root>/Unit Delay5' */
      MTDPC_DW.UnitDelay5_DSTATE = MTDPC_B.D;

      /* Update for UnitDelay: '<Root>/Unit Delay4' */
      MTDPC_DW.UnitDelay4_DSTATE = MTDPC_B.N;

      /* Update for UnitDelay: '<Root>/Unit Delay7' */
      MTDPC_DW.UnitDelay7_DSTATE = MTDPC_B.B;

      /* Update for UnitDelay: '<Root>/Unit Delay6' */
      MTDPC_DW.UnitDelay6_DSTATE = MTDPC_B.A;
    }

    /* Update for VariableTransportDelay: '<Root>/Variable Transport Delay' */
    {
      real_T **uBuffer = (real_T**)
        &MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)
        &MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1];
      real_T simTime = MTDPC_M->Timing.t[0];
      MTDPC_DW.VariableTransportDelay_IWORK.Head =
        ((MTDPC_DW.VariableTransportDelay_IWORK.Head <
          (MTDPC_DW.VariableTransportDelay_IWORK.CircularBufSize-1)) ?
         (MTDPC_DW.VariableTransportDelay_IWORK.Head+1) : 0);
      if (MTDPC_DW.VariableTransportDelay_IWORK.Head ==
          MTDPC_DW.VariableTransportDelay_IWORK.Tail) {
        MTDPC_DW.VariableTransportDelay_IWORK.Tail =
          ((MTDPC_DW.VariableTransportDelay_IWORK.Tail <
            (MTDPC_DW.VariableTransportDelay_IWORK.CircularBufSize-1)) ?
           (MTDPC_DW.VariableTransportDelay_IWORK.Tail+1) : 0);
      }

      (*tBuffer)[MTDPC_DW.VariableTransportDelay_IWORK.Head] = simTime;
      (*uBuffer)[MTDPC_DW.VariableTransportDelay_IWORK.Head] =
        MTDPC_B.TransportDelay;

      /* when use fixed buffer, reset solver at when buffer is updated
       * to avoid output consistency fail.
       */
    }

    /* Update for TransportDelay: '<Root>/Transport Delay' */
    {
      real_T **uBuffer = (real_T**)&MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)&MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[1];
      real_T simTime = MTDPC_M->Timing.t[0];
      MTDPC_DW.TransportDelay_IWORK.Head = ((MTDPC_DW.TransportDelay_IWORK.Head <
        (MTDPC_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
        (MTDPC_DW.TransportDelay_IWORK.Head+1) : 0);
      if (MTDPC_DW.TransportDelay_IWORK.Head ==
          MTDPC_DW.TransportDelay_IWORK.Tail) {
        MTDPC_DW.TransportDelay_IWORK.Tail =
          ((MTDPC_DW.TransportDelay_IWORK.Tail <
            (MTDPC_DW.TransportDelay_IWORK.CircularBufSize-1)) ?
           (MTDPC_DW.TransportDelay_IWORK.Tail+1) : 0);
      }

      (*tBuffer)[MTDPC_DW.TransportDelay_IWORK.Head] = simTime;
      (*uBuffer)[MTDPC_DW.TransportDelay_IWORK.Head] = rtb_Sum_o;
    }

    /* Update for TransportDelay: '<Root>/Transport Delay1' */
    {
      real_T **uBuffer = (real_T**)&MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs
        [0];
      real_T **tBuffer = (real_T**)&MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs
        [1];
      real_T simTime = MTDPC_M->Timing.t[0];
      MTDPC_DW.TransportDelay1_IWORK.Head =
        ((MTDPC_DW.TransportDelay1_IWORK.Head <
          (MTDPC_DW.TransportDelay1_IWORK.CircularBufSize-1)) ?
         (MTDPC_DW.TransportDelay1_IWORK.Head+1) : 0);
      if (MTDPC_DW.TransportDelay1_IWORK.Head ==
          MTDPC_DW.TransportDelay1_IWORK.Tail) {
        MTDPC_DW.TransportDelay1_IWORK.Tail =
          ((MTDPC_DW.TransportDelay1_IWORK.Tail <
            (MTDPC_DW.TransportDelay1_IWORK.CircularBufSize-1)) ?
           (MTDPC_DW.TransportDelay1_IWORK.Tail+1) : 0);
      }

      (*tBuffer)[MTDPC_DW.TransportDelay1_IWORK.Head] = simTime;
      (*uBuffer)[MTDPC_DW.TransportDelay1_IWORK.Head] = rtb_C2;
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(MTDPC_M)) {
    /* signal main to stop simulation */
    {                                  /* Sample time: [0.0s, 0.0s] */
      if ((rtmGetTFinal(MTDPC_M)!=-1) &&
          !((rtmGetTFinal(MTDPC_M)-(((MTDPC_M->Timing.clockTick1+
               MTDPC_M->Timing.clockTickH1* 4294967296.0)) * 0.0001)) >
            (((MTDPC_M->Timing.clockTick1+MTDPC_M->Timing.clockTickH1*
               4294967296.0)) * 0.0001) * (DBL_EPSILON))) {
        rtmSetErrorStatus(MTDPC_M, "Simulation finished");
      }
    }

    rt_ertODEUpdateContinuousStates(&MTDPC_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick0 and the high bits
     * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++MTDPC_M->Timing.clockTick0)) {
      ++MTDPC_M->Timing.clockTickH0;
    }

    MTDPC_M->Timing.t[0] = rtsiGetSolverStopTime(&MTDPC_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.0001s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.0001, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       * Timer of this task consists of two 32 bit unsigned integers.
       * The two integers represent the low bits Timing.clockTick1 and the high bits
       * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
       */
      MTDPC_M->Timing.clockTick1++;
      if (!MTDPC_M->Timing.clockTick1) {
        MTDPC_M->Timing.clockTickH1++;
      }
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void MTDPC_derivatives(void)
{
  XDot_MTDPC_T *_rtXdot;
  _rtXdot = ((XDot_MTDPC_T *) MTDPC_M->derivs);

  /* Derivatives for VariableTransportDelay: '<Root>/Variable Transport Delay1' */
  {
  }

  /* Derivatives for TransferFcn: '<Root>/1//Zm' */
  _rtXdot->uZm_CSTATE = 0.0;
  _rtXdot->uZm_CSTATE += MTDPC_P.uZm_A * MTDPC_X.uZm_CSTATE;
  _rtXdot->uZm_CSTATE += MTDPC_B.Sum1;

  /* Derivatives for Integrator: '<Root>/Integrator2' */
  _rtXdot->Integrator2_CSTATE = MTDPC_B.Product;

  /* Derivatives for Integrator: '<Root>/Integrator3' */
  _rtXdot->Integrator3_CSTATE = MTDPC_B.Product1;

  /* Derivatives for TransferFcn: '<Root>/Ze' */
  _rtXdot->Ze_CSTATE = 0.0;
  _rtXdot->Ze_CSTATE += MTDPC_P.Ze_A * MTDPC_X.Ze_CSTATE;
  _rtXdot->Ze_CSTATE += MTDPC_B.Vth;

  /* Derivatives for VariableTransportDelay: '<Root>/Variable Transport Delay' */
  {
  }

  /* Derivatives for TransferFcn: '<S8>/Transfer Fcn2' */
  _rtXdot->TransferFcn2_CSTATE = 0.0;
  _rtXdot->TransferFcn2_CSTATE += MTDPC_P.TransferFcn2_A *
    MTDPC_X.TransferFcn2_CSTATE;
  _rtXdot->TransferFcn2_CSTATE += MTDPC_B.Vp;

  /* Derivatives for TransferFcn: '<S8>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += MTDPC_P.TransferFcn_A *
    MTDPC_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += MTDPC_B.Vp;

  /* Derivatives for TransferFcn: '<Root>/1//Zs' */
  _rtXdot->uZs_CSTATE = 0.0;
  _rtXdot->uZs_CSTATE += MTDPC_P.uZs_A * MTDPC_X.uZs_CSTATE;
  _rtXdot->uZs_CSTATE += MTDPC_B.Fth_l;

  /* Derivatives for Integrator: '<Root>/Integrator' */
  _rtXdot->Integrator_CSTATE = MTDPC_B.Vp;

  /* Derivatives for Integrator: '<Root>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = MTDPC_B.Vth;

  /* Derivatives for Integrator: '<S36>/Filter' */
  _rtXdot->Filter_CSTATE = MTDPC_B.FilterCoefficient;

  /* Derivatives for TransferFcn: '<S1>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE_d = 0.0;
  _rtXdot->TransferFcn_CSTATE_d += MTDPC_P.TransferFcn_A_e *
    MTDPC_X.TransferFcn_CSTATE_d;
  _rtXdot->TransferFcn_CSTATE_d += MTDPC_B.Vp;

  /* Derivatives for TransferFcn: '<Root>/Cm' */
  _rtXdot->Cm_CSTATE = 0.0;
  _rtXdot->Cm_CSTATE += MTDPC_P.Cm_A * MTDPC_X.Cm_CSTATE;
  _rtXdot->Cm_CSTATE += MTDPC_B.Vp;

  /* Derivatives for Integrator: '<S84>/Filter' */
  _rtXdot->Filter_CSTATE_i = MTDPC_B.FilterCoefficient_j;
}

/* Model initialize function */
void MTDPC_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)MTDPC_M, 0,
                sizeof(RT_MODEL_MTDPC_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&MTDPC_M->solverInfo, &MTDPC_M->Timing.simTimeStep);
    rtsiSetTPtr(&MTDPC_M->solverInfo, &rtmGetTPtr(MTDPC_M));
    rtsiSetStepSizePtr(&MTDPC_M->solverInfo, &MTDPC_M->Timing.stepSize0);
    rtsiSetdXPtr(&MTDPC_M->solverInfo, &MTDPC_M->derivs);
    rtsiSetContStatesPtr(&MTDPC_M->solverInfo, (real_T **) &MTDPC_M->contStates);
    rtsiSetNumContStatesPtr(&MTDPC_M->solverInfo, &MTDPC_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&MTDPC_M->solverInfo,
      &MTDPC_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&MTDPC_M->solverInfo,
      &MTDPC_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&MTDPC_M->solverInfo,
      &MTDPC_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&MTDPC_M->solverInfo, (&rtmGetErrorStatus(MTDPC_M)));
    rtsiSetRTModelPtr(&MTDPC_M->solverInfo, MTDPC_M);
  }

  rtsiSetSimTimeStep(&MTDPC_M->solverInfo, MAJOR_TIME_STEP);
  MTDPC_M->intgData.y = MTDPC_M->odeY;
  MTDPC_M->intgData.f[0] = MTDPC_M->odeF[0];
  MTDPC_M->intgData.f[1] = MTDPC_M->odeF[1];
  MTDPC_M->intgData.f[2] = MTDPC_M->odeF[2];
  MTDPC_M->intgData.f[3] = MTDPC_M->odeF[3];
  MTDPC_M->contStates = ((X_MTDPC_T *) &MTDPC_X);
  rtsiSetSolverData(&MTDPC_M->solverInfo, (void *)&MTDPC_M->intgData);
  rtsiSetSolverName(&MTDPC_M->solverInfo,"ode4");
  rtmSetTPtr(MTDPC_M, &MTDPC_M->Timing.tArray[0]);
  rtmSetTFinal(MTDPC_M, 20.0);
  MTDPC_M->Timing.stepSize0 = 0.0001;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = NULL;
    MTDPC_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(MTDPC_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(MTDPC_M->rtwLogInfo, (NULL));
    rtliSetLogT(MTDPC_M->rtwLogInfo, "tout");
    rtliSetLogX(MTDPC_M->rtwLogInfo, "");
    rtliSetLogXFinal(MTDPC_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(MTDPC_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(MTDPC_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(MTDPC_M->rtwLogInfo, 0);
    rtliSetLogDecimation(MTDPC_M->rtwLogInfo, 1);
    rtliSetLogY(MTDPC_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(MTDPC_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(MTDPC_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &MTDPC_B), 0,
                sizeof(B_MTDPC_T));

  /* states (continuous) */
  {
    (void) memset((void *)&MTDPC_X, 0,
                  sizeof(X_MTDPC_T));
  }

  /* states (dwork) */
  (void) memset((void *)&MTDPC_DW, 0,
                sizeof(DW_MTDPC_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(MTDPC_M->rtwLogInfo, 0.0, rtmGetTFinal
    (MTDPC_M), MTDPC_M->Timing.stepSize0, (&rtmGetErrorStatus(MTDPC_M)));

  /* Start for VariableTransportDelay: '<Root>/Variable Transport Delay1' */
  {
    real_T *pBuffer = &MTDPC_DW.VariableTransportDelay1_RWORK.TUbufferArea[0];
    int_T j;
    MTDPC_DW.VariableTransportDelay1_IWORK.Tail = 0;
    MTDPC_DW.VariableTransportDelay1_IWORK.Head = 0;
    MTDPC_DW.VariableTransportDelay1_IWORK.Last = 0;
    MTDPC_DW.VariableTransportDelay1_IWORK.CircularBufSize = 1024;
    for (j=0; j < 1024; j++) {
      pBuffer[j] = MTDPC_P.VariableTransportDelay1_InitOut;
      pBuffer[1024 + j] = MTDPC_M->Timing.t[0];
    }

    MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
    MTDPC_DW.VariableTransportDelay1_PWORK.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* Start for VariableTransportDelay: '<Root>/Variable Transport Delay' */
  {
    real_T *pBuffer = &MTDPC_DW.VariableTransportDelay_RWORK.TUbufferArea[0];
    int_T j;
    MTDPC_DW.VariableTransportDelay_IWORK.Tail = 0;
    MTDPC_DW.VariableTransportDelay_IWORK.Head = 0;
    MTDPC_DW.VariableTransportDelay_IWORK.Last = 0;
    MTDPC_DW.VariableTransportDelay_IWORK.CircularBufSize = 1024;
    for (j=0; j < 1024; j++) {
      pBuffer[j] = MTDPC_P.VariableTransportDelay_InitOutp;
      pBuffer[1024 + j] = MTDPC_M->Timing.t[0];
    }

    MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
    MTDPC_DW.VariableTransportDelay_PWORK.TUbufferPtrs[1] = (void *) &pBuffer
      [1024];
  }

  /* Start for TransportDelay: '<Root>/Transport Delay' */
  {
    real_T *pBuffer = &MTDPC_DW.TransportDelay_RWORK.TUbufferArea[0];
    MTDPC_DW.TransportDelay_IWORK.Tail = 0;
    MTDPC_DW.TransportDelay_IWORK.Head = 0;
    MTDPC_DW.TransportDelay_IWORK.Last = 0;
    MTDPC_DW.TransportDelay_IWORK.CircularBufSize = 1024;
    pBuffer[0] = MTDPC_P.TransportDelay_InitOutput;
    pBuffer[1024] = MTDPC_M->Timing.t[0];
    MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
    MTDPC_DW.TransportDelay_PWORK.TUbufferPtrs[1] = (void *) &pBuffer[1024];
  }

  /* Start for TransportDelay: '<Root>/Transport Delay1' */
  {
    real_T *pBuffer = &MTDPC_DW.TransportDelay1_RWORK.TUbufferArea[0];
    MTDPC_DW.TransportDelay1_IWORK.Tail = 0;
    MTDPC_DW.TransportDelay1_IWORK.Head = 0;
    MTDPC_DW.TransportDelay1_IWORK.Last = 0;
    MTDPC_DW.TransportDelay1_IWORK.CircularBufSize = 1024;
    pBuffer[0] = MTDPC_P.TransportDelay1_InitOutput;
    pBuffer[1024] = MTDPC_M->Timing.t[0];
    MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
    MTDPC_DW.TransportDelay1_PWORK.TUbufferPtrs[1] = (void *) &pBuffer[1024];
  }

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay' */
  MTDPC_DW.UnitDelay_DSTATE = MTDPC_P.UnitDelay_InitialCondition;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay3' */
  MTDPC_DW.UnitDelay3_DSTATE = MTDPC_P.UnitDelay3_InitialCondition;

  /* InitializeConditions for TransferFcn: '<Root>/1//Zm' */
  MTDPC_X.uZm_CSTATE = 0.0;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay1' */
  MTDPC_DW.UnitDelay1_DSTATE = MTDPC_P.UnitDelay1_InitialCondition;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay2' */
  MTDPC_DW.UnitDelay2_DSTATE = MTDPC_P.UnitDelay2_InitialCondition;

  /* InitializeConditions for Integrator: '<Root>/Integrator2' */
  MTDPC_X.Integrator2_CSTATE = MTDPC_P.Integrator2_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator3' */
  MTDPC_X.Integrator3_CSTATE = MTDPC_P.Integrator3_IC;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay5' */
  MTDPC_DW.UnitDelay5_DSTATE = MTDPC_P.UnitDelay5_InitialCondition;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay4' */
  MTDPC_DW.UnitDelay4_DSTATE = MTDPC_P.UnitDelay4_InitialCondition;

  /* InitializeConditions for TransferFcn: '<Root>/Ze' */
  MTDPC_X.Ze_CSTATE = 0.0;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay7' */
  MTDPC_DW.UnitDelay7_DSTATE = MTDPC_P.UnitDelay7_InitialCondition;

  /* InitializeConditions for UnitDelay: '<Root>/Unit Delay6' */
  MTDPC_DW.UnitDelay6_DSTATE = MTDPC_P.UnitDelay6_InitialCondition;

  /* InitializeConditions for TransferFcn: '<S8>/Transfer Fcn2' */
  MTDPC_X.TransferFcn2_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S8>/Transfer Fcn' */
  MTDPC_X.TransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<Root>/1//Zs' */
  MTDPC_X.uZs_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<Root>/Integrator' */
  MTDPC_X.Integrator_CSTATE = MTDPC_P.Integrator_IC;

  /* InitializeConditions for Integrator: '<Root>/Integrator1' */
  MTDPC_X.Integrator1_CSTATE = MTDPC_P.Integrator1_IC;

  /* InitializeConditions for Integrator: '<S36>/Filter' */
  MTDPC_X.Filter_CSTATE = MTDPC_P.PIDController_InitialConditionF;

  /* InitializeConditions for TransferFcn: '<S1>/Transfer Fcn' */
  MTDPC_X.TransferFcn_CSTATE_d = 0.0;

  /* InitializeConditions for TransferFcn: '<Root>/Cm' */
  MTDPC_X.Cm_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S84>/Filter' */
  MTDPC_X.Filter_CSTATE_i = MTDPC_P.Cs_InitialConditionForFilter;
}

/* Model terminate function */
void MTDPC_terminate(void)
{
  /* (no terminate code required) */
}
