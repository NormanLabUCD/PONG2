#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* PONG2 C routines */
extern void PONG2_AlleleStrand(void);
extern void PONG2_BEDFlag(void);
extern void PONG2_Classifier_GetHaplos(void);
extern void PONG2_Close(void);
extern void PONG2_Confusion(void);
extern void PONG2_ConvBED(void);
extern void PONG2_Done(void);
extern void PONG2_ErrMsg(void);
extern void PONG2_GetNumClassifiers(void);
extern void PONG2_GetParam(void);
extern void PONG2_GPU_Init(void);
extern void PONG2_Idv_GetNumHaplo(void);
extern void PONG2_Init(void);
extern void PONG2_New(void);
extern void PONG2_NewClassifierHaplo(void);
extern void PONG2_NewClassifiers(void);
extern void PONG2_Predict_Prob(void);
extern void PONG2_Predict_Resp(void);
extern void PONG2_Predict_Resp_Prob(void);
extern void PONG2_SetParam(void);
extern void PONG2_SortAlleleStr(void);
extern void PONG2_Training(void);

static const R_CMethodDef CEntries[] = {
    {"PONG2_AlleleStrand",       (DL_FUNC) &PONG2_AlleleStrand,       0},
    {"PONG2_BEDFlag",            (DL_FUNC) &PONG2_BEDFlag,            0},
    {"PONG2_Classifier_GetHaplos",(DL_FUNC) &PONG2_Classifier_GetHaplos, 0},
    {"PONG2_Close",              (DL_FUNC) &PONG2_Close,              0},
    {"PONG2_Confusion",          (DL_FUNC) &PONG2_Confusion,          0},
    {"PONG2_ConvBED",            (DL_FUNC) &PONG2_ConvBED,            0},
    {"PONG2_Done",               (DL_FUNC) &PONG2_Done,               0},
    {"PONG2_ErrMsg",             (DL_FUNC) &PONG2_ErrMsg,             0},
    {"PONG2_GetNumClassifiers",  (DL_FUNC) &PONG2_GetNumClassifiers,  0},
    {"PONG2_GetParam",           (DL_FUNC) &PONG2_GetParam,           0},
    {"PONG2_GPU_Init",           (DL_FUNC) &PONG2_GPU_Init,           0},
    {"PONG2_Idv_GetNumHaplo",    (DL_FUNC) &PONG2_Idv_GetNumHaplo,   0},
    {"PONG2_Init",               (DL_FUNC) &PONG2_Init,               0},
    {"PONG2_New",                (DL_FUNC) &PONG2_New,                0},
    {"PONG2_NewClassifierHaplo", (DL_FUNC) &PONG2_NewClassifierHaplo, 0},
    {"PONG2_NewClassifiers",     (DL_FUNC) &PONG2_NewClassifiers,     0},
    {"PONG2_Predict_Prob",       (DL_FUNC) &PONG2_Predict_Prob,       0},
    {"PONG2_Predict_Resp",       (DL_FUNC) &PONG2_Predict_Resp,       0},
    {"PONG2_Predict_Resp_Prob",  (DL_FUNC) &PONG2_Predict_Resp_Prob,  0},
    {"PONG2_SetParam",           (DL_FUNC) &PONG2_SetParam,           0},
    {"PONG2_SortAlleleStr",      (DL_FUNC) &PONG2_SortAlleleStr,      0},
    {"PONG2_Training",           (DL_FUNC) &PONG2_Training,           0},
    {NULL, NULL, 0}
};

void R_init_PONG2(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);  /* keep TRUE for HIBAG compatibility */
}
