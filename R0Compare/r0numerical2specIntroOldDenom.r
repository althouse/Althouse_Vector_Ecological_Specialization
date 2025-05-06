
################################################################################
#' r0old — Basic reproduction number (R0) for the *opportunistic-biting*
#'          (unweighted) formulation in Althouse *et al.* (2012)
#'
#' Implements the closed-form expression for the spectral radius of the
#' next-generation matrix given in Equation (6) (Appendix S1) of:
#'
#'   Althouse BM, Lessler J, Sall AA, Diallo M, Hanley KA, **et al.** (2012)
#'   *Synchrony of Sylvatic Dengue Isolations: A Multi-Host, Multi-Vector
#'   SIR Model of Dengue Virus Transmission in Senegal.*
#'   \emph{PLOS Neglected Tropical Diseases} 6 (11): e1928.  
#'   \doi{10.1371/journal.pntd.0001928}
#'
#' Parameter names follow the symbol set in the paper:
#'
#'   • b[m?]  : mosquito-to-host transmission probabilities  
#'   • b[h?]  : host-to-mosquito transmission probabilities  
#'   • r[m?h?]: relative feeding preferences  
#'   • V?, N? : vector and host abundances  
#'   • g?, m? : mortality rates  
#'   • intro  : seeding term (set 0 for analytic R0)  
#'   • …      : full list exactly matches Althouse *et al.* (2012)
#'
#' All arguments are numeric; default values are deliberately
#' \code{NULL} so the user must supply model-specific numbers.
#'
#' @author Benjamin M. Althouse
#' @copyright
#'   © 2012 Althouse *et al.* — analytic derivation  
#'   © 2025 Benjamin M. Althouse — this R implementation  
#'   Distributed under the Creative Commons Attribution 4.0 Licence (CC BY 4.0).
#'
#' @return A single numeric value: the basic reproduction number \(R_0\).
#'
#' @examples
#' ## Minimal illustrative call (toy numbers only):
#' # r0old(bm1h = 0.15, bm1p1 = 0.15, bm2h = 0.10, bm2p1 = 0.10,
#' #       bhm1 = 0.25, bp1m1 = 0.25, bhm2 = 0.20, bp1m2 = 0.20,
#' #       rm1h = 0.6, rm1p1 = 0.4, rm2h = 0.5, rm2p1 = 0.5,
#' #       Vm = 2e4, Vm2 = 1e4,
#' #       gh = 0.02, mh = 0.02, gp1 = 0.02, mp1 = 0.02,
#' #       Nh = 2500, Np1 = 1500, Nm1 = 0, Nm2 = 0, NN = 1)
#'
#' @note
#'   The algebraic expression is extremely long; for readability it is wrapped
#'   onto multiple lines but is otherwise verbatim from Appendix S1 of the
#'   source paper.
################################################################################

r0old <- function(
    bm1h, bm1p1, bm1p2 = 0,
    bm2h, bm2p1, bm2p2 = 0,
    bhm1, bp1m1, bp2m1 = 0,
    bhm2, bp1m2, bp2m2 = 0,
    rm1h, rm1p1, rm1p2 = 0,
    rm2h, rm2p1, rm2p2 = 0,
    Vm,  Vm2,
    gh, mh, gp1, mp1, mp2 = 0, gp2 = 0,
    Nh, Np1, Np2 = 0,
    Nm1, Nm2,
    NN,
    intro = 0)
{
(1/sqrt(2))*(sqrt((bhm2 * gp1 * intro * Nh^2 * NN^3 * rm2h * Vm+bhm2 * intro * mp1 * Nh^2 * NN^3 * rm2h * Vm+bhm2 * intro^2 * Nh^2 * NN^3 * Np1 * rm2h * Vm+bhm2 * bm2h * gp1 * Nh * Nm2 * NN^2 * rm2h^2 * Vm+bhm2 * bm2h * mp1 * Nh * Nm2 * NN^2 * rm2h^2 * Vm+bhm2 * bm2h * intro * Nh * Nm2 * NN^2 * Np1 * rm2h^2 * Vm+bp1m2 * gh * intro * NN^3 * Np1 * Np1 * rm2p1 * Vm+bp1m2 * intro * mh * NN^3 * Np1 * Np1 * rm2p1 * Vm+bp1m2 * intro^2 * Nh * NN^3 * Np1 * Np1 * rm2p1 * Vm+bm2p1 * bp1m2 * gh * Nm2 * NN^2 * Np1 * rm2p1^2 * Vm+bm2p1 * bp1m2 * mh * Nm2 * NN^2 * Np1 * rm2p1^2 * Vm+bm2p1 * bp1m2 * intro * Nh * Nm2 * NN^2 * Np1 * rm2p1^2 * Vm+bhm1 * gp1 * intro * Nh^2 * NN^3 * rm1h * Vm2+bhm1 * intro * mp1 * Nh^2 * NN^3 * rm1h * Vm2+bhm1 * intro^2 * Nh^2 * NN^3 * Np1 * rm1h * Vm2+bhm1 * bm1h * gp1 * Nh * Nm1 * NN^2 * rm1h^2 * Vm2+bhm1 * bm1h * mp1 * Nh * Nm1 * NN^2 * rm1h^2 * Vm2+bhm1 * bm1h * intro * Nh * Nm1 * NN^2 * Np1 * rm1h^2 * Vm2+bp1m1 * gh * intro * NN^3 * Np1 * Np1 * rm1p1 * Vm2+bp1m1 * intro * mh * NN^3 * Np1 * Np1 * rm1p1 * Vm2+bp1m1 * intro^2 * Nh * NN^3 * Np1 * Np1 * rm1p1 * Vm2+bm1p1 * bp1m1 * gh * Nm1 * NN^2 * Np1 * rm1p1^2 * Vm2+bm1p1 * bp1m1 * mh * Nm1 * NN^2 * Np1 * rm1p1^2 * Vm2+bm1p1 * bp1m1 * intro * Nh * Nm1 * NN^2 * Np1 * rm1p1^2 * Vm2+sqrt(NN^4 * (4 * Nh * (gh+mh+intro * Nh) * (gp1+mp1+intro * Np1) * Np1 * (bhm2 * bp1m1 * rm1p1 * rm2h-bhm1 * bp1m2 * rm1h * rm2p1) * (-bm1p1 * Nm1 * rm1p1 * (intro * Nh * NN+bm2h * Nm2 * rm2h)+intro * Nm2 * NN * (-bm2h * Np1 * rm2h+bm2p1 * Nh * rm2p1)+bm1h * Nm1 * rm1h * (intro * NN * Np1+bm2p1 * Nm2 * rm2p1)) * Vm * Vm2+(bhm2 * Nh * (gp1+mp1+intro * Np1) * rm2h * (intro * Nh * NN+bm2h * Nm2 * rm2h) * Vm+bp1m2 * (gh+mh+intro * Nh) * Np1 * rm2p1 * (intro * NN * Np1+bm2p1 * Nm2 * rm2p1) * Vm+(bhm1 * Nh * (gp1+mp1+intro * Np1) * rm1h * (intro * Nh * NN+bm1h * Nm1 * rm1h)+bp1m1 * (gh+mh+intro * Nh) * Np1 * rm1p1 * (intro * NN * Np1+bm1p1 * Nm1 * rm1p1)) * Vm2)^2)))/((gh+mh+intro * Nh) * NN^4 * (gp1+mp1+intro * Np1) * Vm * Vm2)))
}