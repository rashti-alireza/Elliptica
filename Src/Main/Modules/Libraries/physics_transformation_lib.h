/* struct for transformation boost, rotation etc. */
typedef struct TRANSFORMATION_T
{
  /* Lorentz Boost */
  struct
  {
    double Bx;/* v_x/c in vector B */
    double By;/* v_y/c in vector B*/
    double Bz;/* v_z/c in vector B*/
    double B2;/* B.B, this needs to be given by user b/c of dot product */
    int inverse;/* if 1 it uses the INVERT of the Lorentz matrix transformation, 
                // 0 uses the Lorentz matrix transformation*/
  }boost[1];
}Transformation_T;

void Lorentz_boost(Transformation_T *const t,const double *const in,double *const out);
