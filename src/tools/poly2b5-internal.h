#ifndef POLY2B5_INTERNAL_H
#define POLY2B5_INTERNAL_H

// Implementation details of cpt_poly2b5, compiled as separate translation units.
namespace poly::detail {

void eval_ha1(const double* xH, double* pHa1);
void eval_ha2(const double* xH, double* pHa2);
void eval_hb1(const double* xH, double* pHb1);
void eval_hb2(const double* xH, double* pHb2);
void eval_oa(const double* xO, double* pOa);
void eval_ob(const double* xO, double* pOb);

} // namespace poly::detail

#endif // POLY2B5_INTERNAL_H
