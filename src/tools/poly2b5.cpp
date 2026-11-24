#include "poly2b5.h"
#include "poly2b5-internal.h"

namespace poly {

void cpt_poly2b5(const double* xH, const double* xO,
                 double* pHa1, double* pHa2, double* pHb1, double* pHb2,
                 double* pOa, double* pOb)
{
    // Preserve the original assignment order, including for overlapping arrays.
    detail::eval_ha1(xH, pHa1);
    detail::eval_ha2(xH, pHa2);
    detail::eval_hb1(xH, pHb1);
    detail::eval_hb2(xH, pHb2);
    detail::eval_oa(xO, pOa);
    detail::eval_ob(xO, pOb);
}

} // namespace poly
