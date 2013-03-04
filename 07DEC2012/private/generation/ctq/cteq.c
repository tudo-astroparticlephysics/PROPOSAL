#include <jni.h>

void setctq6_(int *);
double ctq6pdf_(int *, double *, double *);

JNIEXPORT void JNICALL Java_gen_CteqPDF_SetCtq6(JNIEnv *e, jobject a, jint set){
  int iset=set;
  setctq6_(&iset);
};

JNIEXPORT jdouble JNICALL Java_gen_CteqPDF_Ctq6Pdf(JNIEnv *e, jobject a, jint p, jdouble x, jdouble q){
  int ip=p;
  double ix=x, iq=q;
  return ctq6pdf_(&ip, &ix, &iq);
};
