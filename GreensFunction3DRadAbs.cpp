#include <stdexcept>
#include <vector>
#include <sstream>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include "factorial.hpp"
#include "funcSum.hpp"
#include "findRoot.hpp"
#include "freeFunctions.hpp"
#include "SphericalBesselGenerator.hpp"
#include "GreensFunction3DRadAbs.hpp"

const Real GreensFunction3DRadAbs::TOLERANCE = 1e-8;
const Real GreensFunction3DRadAbs::THETA_TOLERANCE = 1e-5;
const Real GreensFunction3DRadAbs::MIN_T_FACTOR = 1e-8;
const uint GreensFunction3DRadAbs::MAX_ORDER = GF_MAX_ORDER;
const uint GreensFunction3DRadAbs::MAX_ALPHA_SEQ = 2000;

Logger& GreensFunction3DRadAbs::log_(Logger::get_logger("GreensFunction3DRadAbs"));

GreensFunction3DRadAbs::GreensFunction3DRadAbs(Real D, Real kf, Real r0, Real Sigma, Real a)
    : PairGreensFunction(D, kf, r0, Sigma), a(a), h(kf / (4.0 * M_PI * Sigma * Sigma * D)), hsigma_p_1(1.0 + h * Sigma)
{
    if (a < sigma)
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: a >= sigma : a=%.16g, sigma=%.16g") % a % sigma).str());
    }
    clearAlphaTable();
}

//
// Alpha-related methods
//

void GreensFunction3DRadAbs::clearAlphaTable() const
{
    std::for_each(alphaTable.begin(), alphaTable.end(), boost::mem_fn(&RealVector::clear));
    alphaOffsetTable[0] = 0;
    std::fill(alphaOffsetTable.begin() + 1, alphaOffsetTable.end(), -1);
}

Real GreensFunction3DRadAbs::f_alpha0(Real alpha) const
{
    const Real alpha_a_m_sigma(alpha * (a - sigma));
    const Real term1(alpha * sigma * cos(alpha_a_m_sigma));
    const Real term2(hsigma_p_1 * sin(alpha_a_m_sigma));
    return term1 + term2;
}

Real GreensFunction3DRadAbs::f_alpha0_aux(Real alpha) const
{
    const Real term1((a - sigma) * alpha);
    const Real angle(hsigma_p_1 / (sigma * alpha));
    const Real term2(std::atan(angle));
    return term1 - term2;
}

struct f_alpha0_aux_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real value;
};

static Real f_alpha0_aux_F(Real alpha, f_alpha0_aux_params const* params)
{
    return params->gf->f_alpha0_aux(alpha) - params->value;
}

Real GreensFunction3DRadAbs::alpha0_i(int i) const
{
    if (!(i >= 0))
    {
        throw std::out_of_range((boost::format("GreensFunction3DRadAbs: i >= 0 : i=%.16g") % i).str());
    }

    const Real target(i * M_PI + M_PI_2);
    f_alpha0_aux_params params = { this, target };

    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&f_alpha0_aux_F), &params };

    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real interval(M_PI / (a - sigma));
    Real low(i * interval + std::numeric_limits<Real>::epsilon());
    Real high((i + 1) * interval);

    //assert(GSL_FN_EVAL(&F, low) * GSL_FN_EVAL(&F, high) < 0.0);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const uint maxIter(100);
    for (uint j(0);;++j)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 0.0, 1e-15));

        if (status != GSL_CONTINUE) break;
        if (j >= maxIter) throw std::runtime_error("GreensFunction3DRadAbs: alpha0_i: failed to converge");
    }

    const Real alpha(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);
    return alpha;
}

void GreensFunction3DRadAbs::updateAlphaTable0(const Real t) const
{
    RealVector& alphaTable_0(getAlphaTable(0));
    alphaTable_0.clear();
    alphaTable_0.reserve(MAX_ALPHA_SEQ);

    const Real alpha0_0(alpha0_i(0));
    alphaTable_0.push_back(alpha0_0);

    const Real Dt(D * t);

    //    const Real alpha_cutoff(sqrt((- log(TOLERANCE * 1e-2) / Dt)
    //                                 + alpha0_0 * alpha0_0));
    const Real alpha_cutoff(sqrt((-log(TOLERANCE * 1e-3) / Dt)));
    for (uint i(1);;)
    {
        const Real alpha0_i(this->alpha0_i(i));
        alphaTable_0.push_back(alpha0_i);

        if (alpha0_i > alpha_cutoff && i >= 10) // make at least 10 terms
            break;

        ++i;

        if (i >= MAX_ALPHA_SEQ)
            break;
    }
}

Real GreensFunction3DRadAbs::f_alpha(Real alpha, int n) const
{
    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real hSigma(h * sigma);
    const Real hSigma_m_n(hSigma - n);

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n, sigmaAlpha));
    const Real ys1(s.y(n, sigmaAlpha));
    const Real js2(s.j(n + 1, sigmaAlpha));
    const Real ys2(s.y(n + 1, sigmaAlpha));
    const Real ja(s.j(n, aAlpha));
    const Real ya(s.y(n, aAlpha));

    const Real term1((hSigma_m_n * js1 + sigmaAlpha * js2) * ya);
    const Real term2((hSigma_m_n * ys1 + sigmaAlpha * ys2) * ja);

    const Real factor(2.0 * alpha * sqrt(a * sigma) * M_1_PI);

    return (term1 - term2) * factor;
}

static inline Real G(const uint n, const uint k)
{
    return factorial(n + k) * (factorial_r(k) * factorial_r(n - k));
}

static Real P(int n, Real x)
{
    Real result(0.0);

    Real sx2(1.0);
    int term1(1);

    const Real x2sq_r(1.0 / gsl_pow_2(x + x));
    const uint maxm(n / 2);
    for (uint m(0); m <= maxm; ++m)
    {
        const Real value(term1 * sx2 * G(n, 2 * m));
        result += value;
        term1 = -term1;
        sx2 *= x2sq_r;
    }

    return result;
}

static std::pair<Real, Real> P2(int n, Real x)
{
    Real result(0.0);
    Real resultp(0.0);

    Real sx2(1.0);
    int term1(1);

    const Real x2sq_r(1.0 / gsl_pow_2(x + x));
    const uint np1(n + 1);
    const uint maxm(n / 2);
    for (uint m(0); m <= maxm; ++m)
    {
        const Real sx2p(term1 * sx2);
        const uint m2(2 * m);
        const Real value(sx2p * G(n, m2));
        result += value;

        const Real valuep(sx2p * G(np1, m2));
        resultp += valuep;

        term1 = -term1;
        sx2 *= x2sq_r;
    }

    if (n % 2)
    {
        resultp += term1 * sx2 * G(np1, np1);
    }

    return std::make_pair(result, resultp);
}

static Real Q(int n, Real x)
{
    Real result(0.0);

    Real sx2(1.0 / (x + x));
    int term1(1);

    const Real x2sq(sx2 * sx2);
    const uint maxm((n + 1) / 2); // sum_(0)^((n-1)/2)
    for (uint m(0); m < maxm; ++m)
    {
        const Real value(term1 * sx2 * G(n, 2 * m + 1));
        result += value;

        term1 = -term1;  // (-1)^m
        sx2 *= x2sq;
    }

    return result;
}

static std::pair<Real, Real> Q2(int n, Real x)
{
    Real result(0.0);
    Real resultp(0.0);

    Real sx2(1.0 / (x + x));
    int term1(1);  // (-1)^m

    const Real x2sq(sx2 * sx2);
    const uint np1(n + 1);
    const uint maxm((n + 1) / 2); // sum_(0)^((n-1)/2)
    for (uint m(0); m < maxm; ++m)
    {
        const Real sx2p(term1 * sx2);
        const uint m2p1(2 * m + 1);
        const Real value(sx2p * G(n, m2p1));
        result += value;

        const Real valuep(sx2p * G(np1, m2p1));
        resultp += valuep;

        term1 = -term1; // (-1)^m
        sx2 *= x2sq;
    }

    if (!(n % 2))
    {
        resultp += term1 * sx2 * G(np1, np1);
    }

    return std::make_pair(result, resultp);
}

Real GreensFunction3DRadAbs::f_alpha_aux(Real alpha, int n) const
{
    if (alpha == 0.0) return -1.0;

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real n_m_hSigma(n - h * sigma);

    /*(a - s) u -
      ArcTan[(P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
      Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))/
      (Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) +
      P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))]
      */

    const Real Pa(P(n, aAlpha));
    const Real Qa(Q(n, aAlpha));

    Real Ps;
    Real Psp;
    // boost::tie(Ps, Psp) = P2(n, sigmaAlpha);
    {
        std::pair<Real, Real> res(P2(n, sigmaAlpha));
        Ps = res.first;
        Psp = res.second;
    }

    Real Qs;
    Real Qsp;
    // boost::tie(Qs, Qsp) = Q2(n, sigmaAlpha);
    {
        std::pair<Real, Real> res(Q2(n, sigmaAlpha));
        Qs = res.first;
        Qsp = res.second;
    }

    const Real n_m_hSigmaPs(n_m_hSigma * Ps);
    const Real n_m_hSigmaQs(n_m_hSigma * Qs);
    const Real sigmaAlphaPsp(sigmaAlpha * Psp);
    const Real sigmaAlphaQsp(sigmaAlpha * Qsp);

    const Real Qa_Pa(Qa / Pa);
    const Real A(sigmaAlphaQsp - n_m_hSigmaPs);
    const Real B(sigmaAlphaPsp + n_m_hSigmaQs);

    // this form, dividing all terms by Pa, prevents overflow.
    const Real angle((A - Qa_Pa * B) / (Qa_Pa * A + B));
    const Real term1((a - sigma) * alpha);
    const Real term2(std::atan(angle));
    return term1 - term2;
}

struct f_alpha_aux_params
{
    GreensFunction3DRadAbs const* const gf;
    const int n;
    const Real value;
};

static Real f_alpha_aux_F(Real alpha, f_alpha_aux_params const* params)
{
    return params->gf->f_alpha_aux(alpha, params->n) - params->value;
}

Real GreensFunction3DRadAbs::alpha_i(int i, int n, gsl_root_fsolver* solver) const
{
    const Real target(M_PI * i + M_PI_2);
    const Real factor(1.0 / (a - sigma));
    Real low((target - M_PI_2) * factor);
    Real high((target + M_PI_2) * factor);

    f_alpha_aux_params params = { this, n, target };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&f_alpha_aux_F), &params };
    gsl_root_fsolver_set(solver, &F, low, high);

    const uint maxIter(100);
    for (uint k(0);;++k)
    {
        gsl_root_fsolver_iterate(solver);

        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 1e-6, 1e-15));

        if (status != GSL_CONTINUE) break;
        if (k >= maxIter) throw std::runtime_error("GreensFunction3DRadAbs: alpha_i: failed to converge");
    }
    return gsl_root_fsolver_root(solver);
}

uint GreensFunction3DRadAbs::alphaOffset(uint n) const
{
    assert(alphaOffsetTable.size() > n);

    if (alphaOffsetTable[n] >= 0)
        return alphaOffsetTable[n];

    uint offset(alphaOffsetTable[n - 1]);
    const Real factor(1.0 / (a - sigma));

    Real target(offset * M_PI + M_PI_2);
    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real alphaMid(target * factor);
    const Real alphaHalfRange(M_PI_2 * factor);
    Real low(alphaMid - alphaHalfRange * (1.0 - 1e-3)); // avoid zero.
    Real high(alphaMid + alphaHalfRange);

    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // (Pi * offset + Pi/2) +- Pi/2 / (a - sigma)
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / (a - sigma).

    Real lowvalue(f_alpha(low, n));
    Real highvalue(f_alpha(high, n));

    for (;;) // this can be much faster if better initial guess is given.
    {
        if (lowvalue * highvalue < 0) // low and high straddle?
            break;

        ++offset;
        target = M_PI * offset + M_PI_2;
        low = (target - M_PI_2) * factor;
        high = (target + M_PI_2) * factor;

        lowvalue = highvalue;
        highvalue = f_alpha(high, n);
    }

    alphaOffsetTable[n] = offset;
    return offset;
}

void GreensFunction3DRadAbs::updateAlphaTable(const uint n, const Real t) const
{
    if (n >= MAX_ORDER)
        throw std::range_error((boost::format("GreensFunction3DRadAbs: n >= 0 && n < MAX_ORDER : n=%d, MAX_ORDER=%d") % n % MAX_ORDER).str());

    if (n == 0)
    {
        updateAlphaTable0(t);
        return;
    }

    const uint offset(alphaOffset(n));

    RealVector& alphaTable_n(getAlphaTable(n));
    alphaTable_n.clear();
    alphaTable_n.reserve(MAX_ALPHA_SEQ);

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

    const Real alphan_0(alpha_i(offset, n, solver));
    const Real alphan_0_sq(alphan_0 * alphan_0);

    alphaTable_n.push_back(alphan_0);

    const Real Dt(D * t);

    const Real threshold(TOLERANCE * 1e-2 * alphan_0_sq * exp(-Dt * alphan_0_sq));

    const uint end(offset + MAX_ALPHA_SEQ);
    for (uint i(offset + 1);;)
    {
        const Real alpha_i(this->alpha_i(i, n, solver));
        alphaTable_n.push_back(alpha_i);

        // cutoff
        const Real alpha_i_sq(alpha_i * alpha_i);
        if (alpha_i_sq * exp(-Dt * alpha_i_sq) < threshold)
            break;

        ++i;

        if (i >= end)
        {
            log_.info("alphaTable (%d): didn't converge. t = %.16g, %s", n, t, dump().c_str());
            break;
        }
    }
    gsl_root_fsolver_free(solver);
}

Real GreensFunction3DRadAbs::p_0_i(Real alpha, Real r) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real angle_r(alpha * (r - sigma));
    Real num1 = alpha * sigma * cos(angle_r) + hsigma_p_1 * sin(angle_r);
    const Real num2(num_r0(alpha));
    const Real den(2 * M_PI * r * r0 * ((a - sigma) * sigmasq * alphasq + hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));
    const Real result(num1 * num2 / den);
    return result;
}

Real GreensFunction3DRadAbs::p_survival_i(Real alpha) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));
    const Real num1(h * sigmasq * hsigma_p_1 - a * (hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq) * cos_a);
    const Real num2(num_r0(alpha));
    const Real den(r0 * hsigma_p_1 * alpha * (-hsigma_p_1 * (a + a * h * sigma - h * sigmasq) + (sigma - a) * sigmasq * alphasq));
    const Real result(-2.0 * num1 * num2 / den);
    return result;
}

Real GreensFunction3DRadAbs::dp_survival_i(Real alpha) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));
    const Real num1(alpha * (h * sigmasq * hsigma_p_1 - (a * (hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq)) * cos_a));
    const Real num2(num_r0(alpha));
    const Real den(r0 * hsigma_p_1 * (-hsigma_p_1 * (a + a * h * sigma - h * sigmasq)) + (sigma - a) * sigmasq * alphasq);
    const Real result(2.0 * D * num1 * num2 / den);
    return result;
}

Real GreensFunction3DRadAbs::leavea_i(Real alpha) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));
    const Real num1(alpha * (hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq) * cos_a);
    const Real num2(num_r0(alpha));
    const Real den(2 * a * M_PI * r0 * hsigma_p_1 * (hsigma_p_1 * (a + a * h * sigma - h * sigmasq) + (a - sigma) * sigmasq * alphasq));
    const Real result(D * num1 * num2 / den);
    return result;
}

Real GreensFunction3DRadAbs::leaves_i(Real alpha) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real num(h * alpha * num_r0(alpha));
    const Real den(2 * M_PI * r0 * ((a - sigma) * sigmasq * alphasq + hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));
    const Real result(-D * num / den);
    return result;
}

Real GreensFunction3DRadAbs::p_leavea_i(Real alpha, Real pleave_factor) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real angle_a(alpha * (a - sigma));
    const Real cos_a(cos(angle_a));
    const Real num1((hsigma_p_1 * hsigma_p_1 + sigmasq * alphasq) * cos_a);
    const Real result(-2.0 * a * num1 * pleave_factor / hsigma_p_1);
    return result;
}

Real GreensFunction3DRadAbs::p_leaves_i(Real alpha, Real pleave_factor) const
{
    const Real num(h * sigma * sigma);
    const Real result(2.0 * num * pleave_factor);
    return result;
}

Real GreensFunction3DRadAbs::p_survival_den(Real alpha) const
{
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real den(r0 * alpha * ((a - sigma) * sigmasq * alphasq + hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));
    return den;
}

Real GreensFunction3DRadAbs::num_r0(Real alpha) const
{
    const Real angle_r0(alpha * (r0 - sigma));
    const Real result(alpha * sigma * cos(angle_r0) + hsigma_p_1 * sin(angle_r0));
    return result;
}

Real GreensFunction3DRadAbs::pleaveFactor(Real alpha) const
{
    return num_r0(alpha) / p_survival_den(alpha);
}

Real GreensFunction3DRadAbs::p_int_r_i(Real r, Real alpha, Real num_r0) const
{
    const Real angle_r(alpha * (r - sigma));
    const Real sigmasq(sigma * sigma);
    const Real alphasq(alpha * alpha);
    const Real hsigma(h * sigma);
    const Real num1(alpha * (hsigma * sigma - hsigma * r * cos(angle_r) - (r - sigma) * cos(angle_r)) + (hsigma_p_1 + r * sigma * alphasq) * sin(angle_r));
    const Real num2(num_r0);
    const Real den(r0 * alphasq * ((a - sigma) * sigmasq * alphasq + hsigma_p_1 * (a + a * h * sigma - h * sigmasq)));
    const Real result(2 * num1 * num2 / den);
    return result;
}

void GreensFunction3DRadAbs::createPsurvTable(RealVector& table) const
{
    const RealVector& alphaTable_0(getAlphaTable(0));
    table.clear();
    table.reserve(alphaTable_0.size());
    std::transform(alphaTable_0.begin(), alphaTable_0.end(), std::back_inserter(table), boost::bind(&GreensFunction3DRadAbs::p_survival_i, this, _1));
}

void GreensFunction3DRadAbs::createNum_r0Table(RealVector& table) const
{
    const RealVector& alphaTable_0(alphaTable[0]);
    table.clear();
    table.reserve(alphaTable_0.size());
    std::transform(alphaTable_0.begin(), alphaTable_0.end(), std::back_inserter(table), boost::bind(&GreensFunction3DRadAbs::num_r0, this, _1));
}

void GreensFunction3DRadAbs::createPleaveFactorTable(RealVector& table) const
{
    const RealVector& alphaTable_0(alphaTable[0]);
    table.clear();
    table.reserve(alphaTable_0.size());
    std::transform(alphaTable_0.begin(), alphaTable_0.end(), std::back_inserter(table), boost::bind(&GreensFunction3DRadAbs::pleaveFactor, this, _1));
}

void GreensFunction3DRadAbs::createPleavesTable(RealVector& table, RealVector const& pleaveFactorTable) const
{
    const RealVector& alphaTable_0(alphaTable[0]);
    assert(pleaveFactorTable.size() >= alphaTable_0.size());
    table.clear();
    table.reserve(alphaTable_0.size());

    for (uint i(0); i < alphaTable_0.size(); ++i)
    {
        const Real alpha(alphaTable_0[i]);
        table.push_back(p_leaves_i(alpha, pleaveFactorTable[i]));
    }
}

void GreensFunction3DRadAbs::createPleaveaTable(RealVector& table, RealVector const& pleaveFactorTable) const
{
    const RealVector& alphaTable_0(alphaTable[0]);
    assert(pleaveFactorTable.size() >= alphaTable_0.size());
    table.clear();
    table.reserve(alphaTable_0.size());

    for (uint i(0); i < alphaTable_0.size(); ++i)
    {
        const Real alpha(alphaTable_0[i]);
        table.push_back(p_leavea_i(alpha, pleaveFactorTable[i]));
    }
}

Real GreensFunction3DRadAbs::p_0_i_exp(uint i, Real t, Real r) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * p_0_i(alpha, r);
}

Real GreensFunction3DRadAbs::p_survival_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    return p_survival_i_alpha(alpha, t);
}

Real GreensFunction3DRadAbs::p_survival_i_alpha(Real alpha, Real t) const
{
    return std::exp(-D * t * alpha * alpha) * p_survival_i(alpha);
}

Real GreensFunction3DRadAbs::p_survival_2i_exp(uint i, Real t) const
{
    const Real Dt(D * t);
    const Real alpha0(getAlpha0(2 * i));
    const Real p0(std::exp(-Dt * alpha0 * alpha0) * p_survival_i(alpha0));
    const Real alpha1(getAlpha0(2 * i + 1));
    const Real p1(std::exp(-Dt * alpha1 * alpha1) * p_survival_i(alpha1));
    return p0 + p1;
}

Real GreensFunction3DRadAbs::p_survival_i_exp_table(uint i, Real t, RealVector const& table) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * table[i];
}

Real GreensFunction3DRadAbs::p_leave_i_exp_table(uint i, Real t, RealVector const& table) const
{
    const Real alpha(getAlpha0(i));
    return expm1(-D * t * alpha * alpha) * table[i];
}

Real GreensFunction3DRadAbs::dp_survival_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * dp_survival_i(alpha);
}

Real GreensFunction3DRadAbs::leavea_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * leavea_i(alpha);
}

Real GreensFunction3DRadAbs::leaves_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * leaves_i(alpha);
}

Real GreensFunction3DRadAbs::p_leavea_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    const Real num_r0(this->num_r0(alpha));
    const Real den(p_survival_den(alpha));
    return exp(-D * t * alpha * alpha) * p_leavea_i(alpha, num_r0 / den);
}

Real GreensFunction3DRadAbs::p_leaves_i_exp(uint i, Real t) const
{
    const Real alpha(getAlpha0(i));
    const Real num_r0(this->num_r0(alpha));
    const Real den(p_survival_den(alpha));
    return exp(-D * t * alpha * alpha) * p_leaves_i(alpha, num_r0 / den);
}

Real GreensFunction3DRadAbs::p_int_r_i_exp(uint i, Real t, Real r) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * p_int_r_i(r, alpha, num_r0(alpha));
}

Real GreensFunction3DRadAbs::p_int_r_i_exp_table(uint i, Real t, Real r, RealVector& num_r0Table) const
{
    const Real alpha(getAlpha0(i));
    return std::exp(-D * t * alpha * alpha) * p_int_r_i(r, alpha, num_r0(alpha));      //num_r0Table[i]);
}

Real GreensFunction3DRadAbs::p_0(Real t, Real r) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::p_0_i_exp, this, _1, t, r), MAX_ALPHA_SEQ);
}

uint GreensFunction3DRadAbs::guess_maxi(Real t) const
{
    const uint safety(2);
    if (!std::isfinite(t)) return safety;

    const Real alpha0(getAlpha0(0));
    const Real Dt(D * t);
    const Real thr(exp(-Dt * alpha0 * alpha0) * TOLERANCE * 1e-1);

    if (thr <= 0.0) return MAX_ALPHA_SEQ;

    const Real max_alpha(sqrt(alpha0 * alpha0 - log(thr) / Dt));
    const uint maxi(safety + static_cast<uint>(max_alpha * (a - sigma) / M_PI));
    return std::min(maxi, MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::p_survival(Real t) const
{
    RealVector psurvTable;
    return p_survival_table(t, psurvTable);
}

Real GreensFunction3DRadAbs::p_survival_table(Real t, RealVector& psurvTable) const
{
    Real p;
    const Real distToa(a - r0);
    const Real distTos(r0 - sigma);
    const Real H(6.0); // a fairly strict criterion for safety.
    const Real maxDist(H * sqrt(6.0 * D * t));

    if (distToa > maxDist)
    {
        if (distTos > maxDist) // far from anything; it'll survive.
        {
            p = 1.0;
        }
        else // close only to s, ignore a
        {
            p = p_survival_irr(t, r0, kf, D, sigma);
        }
    }
    else
    {
        if (distTos > maxDist)  // close only to a.
        {
            p = p_survival_nocollision(t, r0, D, a);
        }
        else  // close to both boundaries.  do the normal calculation.
        {
            const uint maxi(guess_maxi(t));

            if (psurvTable.size() < maxi + 1)
            {
                getAlpha0(maxi);  // this updates the table
                createPsurvTable(psurvTable);
            }

            p = funcSum_all(boost::bind(&GreensFunction3DRadAbs::p_survival_i_exp_table, this, _1, t, psurvTable), maxi);
        }
    }

    return p;
}

Real GreensFunction3DRadAbs::p_leave_table(Real t, RealVector const& table) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::p_leave_i_exp_table, this, _1, t, table), table.size());
}

Real GreensFunction3DRadAbs::dp_survival(Real t) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::dp_survival_i_exp, this, _1, t), MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::leaves(Real t) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::leaves_i_exp, this, _1, t), MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::leavea(Real t) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::leavea_i_exp, this, _1, t), MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::p_leaves(Real t) const
{
    return funcSum_all(boost::bind(&GreensFunction3DRadAbs::p_leaves_i_exp, this, _1, t), guess_maxi(t));
}

Real GreensFunction3DRadAbs::p_leavea(Real t) const
{
    return funcSum_all(boost::bind(&GreensFunction3DRadAbs::p_leavea_i_exp, this, _1, t), guess_maxi(t));
}

Real GreensFunction3DRadAbs::p_int_r(Real r, Real t) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::p_int_r_i_exp, this, _1, t, r), MAX_ALPHA_SEQ);
}

Real GreensFunction3DRadAbs::p_int_r_table(Real r, Real t, RealVector const& num_r0Table) const
{
    return funcSum(boost::bind(&GreensFunction3DRadAbs::p_int_r_i_exp_table, this, _1, t, r, num_r0Table), num_r0Table.size());
}

struct p_survival_table_params
{
    GreensFunction3DRadAbs const* const gf;
    RealVector& table;
    const Real rnd;
};

Real p_survival_table_F(Real t, p_survival_table_params const* params)
{
    return params->rnd - params->gf->p_survival_table(t, params->table);
}

struct p_survival_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real rnd;
};

//NOT USED static Real p_survival_F(Real t, p_survival_params const* params)
//{
//    return params->rnd - params->gf->p_survival(t);
//}

struct p_survival_2i_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real t;
};

//NOT USED static Real p_survival_2i_F(Real ri, p_survival_2i_params const* params)
//{
//    return params->gf->p_survival_2i_exp(static_cast<uint>(ri), params->t);
//}

struct p_survival_i_alpha_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real t;
};

//NOT USED static Real p_survival_i_alpha_F(Real alpha, p_survival_i_alpha_params const* params)
//{
//    return params->gf->p_survival_i_alpha(alpha, params->t);
//}

struct p_leave_params
{
    GreensFunction3DRadAbs const* const gf;
    RealVector const& table;
    const Real rnd;
};

Real p_leave_F(Real t, p_leave_params const* params)
{
    return -params->gf->p_leave_table(t, params->table) - params->rnd;
}

struct p_int_r_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real t;
    const Real rnd;
};

static Real p_int_r_F(Real r, p_int_r_params const* params)
{
    return params->gf->p_int_r(r, params->t) - params->rnd;
}

Real GreensFunction3DRadAbs::drawTime(Real rnd) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 <= a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 <= a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (r0 == a || a == sigma) return 0.0;

    Real t_guess;
    Real dist;

    if (kf != 0)
    {
        dist = std::min(a - r0, r0 - sigma);
    }
    else
    {
        dist = a - r0;
    }

    t_guess = dist * dist / (6.0 * D);
    t_guess *= .1;

    const Real minT(std::min(sigma * sigma / D * MIN_T_FACTOR, t_guess * 1e-6));
    RealVector psurvTable;
    p_survival_table_params params = { this, psurvTable, rnd };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_survival_table_F), &params };
    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
    const Real value(GSL_FN_EVAL(&F, t_guess));
    if (value < 0.0)
    {
        high *= 10;
        for (;;)
        {
            const Real high_value(GSL_FN_EVAL(&F, high));

            if (high_value >= 0.0)
                break;

            if (fabs(high) >= 1e10)
                throw std::runtime_error((boost::format("GreensFunction3DRadAbs: couldn't adjust high. F(%.16g) = %.16g; r0 = %.16g, %s") % high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());

            high *= 10;
        }
    }
    else
    {
        Real low_value_prev(value);
        low *= .1;
        for (;;)
        {
            const Real low_value(GSL_FN_EVAL(&F, low));

            if (low_value <= 0.0)
                break;

            // FIXME: 
            if (fabs(low) <= minT || fabs(low_value - low_value_prev) < TOLERANCE)
            {
                log_.info("GreensFunction3DRadAbs: couldn't adjust low. F(%.16g) = %.16g; r0 = %.16g, %s", low, GSL_FN_EVAL(&F, low), r0, dump().c_str());
                return low;
            }
            low_value_prev = low_value;
            low *= .1;
        }
    }

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    const Real t(findRoot(F, solver, low, high, 0.0, TOLERANCE, "drawTime"));
    gsl_root_fsolver_free(solver);
    return t;
}

GreensFunction::EventKind GreensFunction3DRadAbs::drawEventType(Real rnd, Real t) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t > 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t > 0.0 : t=%.16g") % t).str());
    }

    if (kf == 0) return EventKind::IV_ESCAPE;

    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always IV_ESCAPE or 
    // IV_REACTION, respectively. This avoids numerical instability in 
    // calculating leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    const uint H(6);
    const Real max_dist(H * sqrt(6.0 * D * t));
    const Real a_dist(a - r0);
    const Real s_dist(r0 - sigma);

    if (a_dist > max_dist)
    {
        if (s_dist < max_dist)
        {
            return EventKind::IV_REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if (s_dist > max_dist)
        {
            return EventKind::IV_ESCAPE;
        }
    }

    const Real reaction(leaves(t) * 4.0 * M_PI * sigma * sigma);
    const Real escape(leavea(t) * 4.0 * M_PI * a * a);
    const Real value(reaction / (reaction + escape));
    return rnd <= value ? EventKind::IV_REACTION : EventKind::IV_ESCAPE;
}

Real GreensFunction3DRadAbs::drawPleavea(gsl_function const& F, gsl_root_fsolver* solver, Real t_guess, RealVector& pleaveFactorTable, RealVector& pleaveaTable) const
{
    const Real minT(1e-12);

    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
    const Real value(GSL_FN_EVAL(&F, t_guess));
    if (value < 0.0)
    {
        high *= 10;
        for (;;)
        {
            const Real high_value(GSL_FN_EVAL(&F, high));

            if (high_value >= 0.0)
                break;

            if (fabs(high) >= 1e10)
                throw std::runtime_error((boost::format("GreensFunction3DRadAbs: couldn't adjust high. Fa(%.16g) = %.16g; r0 = %.16g, %s") % high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());

            log_.info("drawTime2: adjusting high: %.16g Fa = %.16g", high, high_value);
            high *= 10;
        }
    }
    else
    {
        Real low_value_prev(value);
        low *= .1;
        for (;;)
        {
            updateAlphaTable0(low);
            createPleaveFactorTable(pleaveFactorTable);
            createPleaveaTable(pleaveaTable, pleaveFactorTable);

            const Real low_value(GSL_FN_EVAL(&F, low));

            if (low_value <= 0.0)
                break;

            // FIXME: 
            if (fabs(low) <= minT || fabs(low_value - low_value_prev) < TOLERANCE)
            {
                log_.info("couldn't adjust low. Fa(%.16g) = %.16g; r0 = %.16g, %s", low, GSL_FN_EVAL(&F, low), r0, dump().c_str());
                return minT;
            }
            low_value_prev = low_value;

            log_.info("drawTime2: adjusting low: %.16g, Fa = %.16g", low, low_value);
            low *= .1;
        }
    }

    const Real t(findRoot(F, solver, low, high, 0., TOLERANCE, "drawTime2: a"));

    return t;
}

Real GreensFunction3DRadAbs::drawPleaves(gsl_function const& F, gsl_root_fsolver* solver, Real t_guess, RealVector& pleaveFactorTable, RealVector& pleavesTable) const
{
    const Real minT(1e-12);

    Real low(t_guess);
    Real high(t_guess);

    // adjust high and low to make sure that f(low) and f(high) straddle.
    const Real value(GSL_FN_EVAL(&F, t_guess));
    if (value < 0.0)
    {
        high *= 10;
        for (;;)
        {
            const Real high_value(GSL_FN_EVAL(&F, high));

            if (high_value >= 0.0) break;
            if (fabs(high) >= 1e10)
                throw std::runtime_error((boost::format("GreensFunction3DRadAbs: couldn't adjust high. Fs(%.16g) = %.16g; r0 = %.16g, %s") % high % GSL_FN_EVAL(&F, high) % r0 % dump()).str());

            log_.info("drawTime2: adjusting high: %.16g Fs = %.16g", high, high_value);
            high *= 10;
        }
    }
    else
    {
        //Real low_value_prev(value);
        low *= .1;
        for (;;)
        {
            updateAlphaTable0(low);
            createPleaveFactorTable(pleaveFactorTable);
            createPleavesTable(pleavesTable, pleaveFactorTable);
            const Real low_value(GSL_FN_EVAL(&F, low));

            if (low_value <= 0.0) break;
            // FIXME: 
            if (fabs(low) <= minT) //|| fabs(low_value - low_value_prev) < TOLERANCE) 
            {
                log_.info("couldn't adjust low.  returning minT (=%.16g);" "Fs(%.16g) = %.16g; r0 = %.16g, %s", minT, low, GSL_FN_EVAL(&F, low), r0, dump().c_str());
                return minT;
            }
            //low_value_prev = low_value;

            log_.info("drawTime2: adjusting low: %.16g, Fs = %.16g", low, low_value);
            low *= .1;
        }
    }

    return findRoot(F, solver, low, high, 0., TOLERANCE, "drawTime2: s");
}

Real GreensFunction3DRadAbs::drawR(Real rnd, Real t) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (t == 0.0)return r0;

    const Real psurv(p_survival(t));

    //    RealVector num_r0Table;
    //    createNum_r0Table(num_r0Table, r0);

    p_int_r_params params = { this, t, /*num_r0Table,*/ rnd * psurv };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&p_int_r_F), &params };

    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    Real low(r0);
    Real high(r0);

    const Real sqrt6Dt(sqrt(6.0 * D * t));
    if (GSL_FN_EVAL(&F, r0) < 0.0)
    {
        for (uint H(3);;)
        {
            high = r0 + H * sqrt6Dt;
            if (high > a)
            {
                if (GSL_FN_EVAL(&F, a) < 0.0)
                {
                    log_.info("drawR: p_int_r(a) < 0.0. returning a");
                    return a;
                }

                high = a;
                break;
            }

            const Real value(GSL_FN_EVAL(&F, high));
            if (value > 0.0) break;
            ++H;
        }

    }
    else
    {
        for (uint H(3);;)
        {
            low = r0 - H * sqrt6Dt;
            if (low < sigma)
            {
                if (GSL_FN_EVAL(&F, sigma) > 0.0)
                {
                    log_.info("drawR: p_int_r(sigma) > 0.0. returning sigma");
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value(GSL_FN_EVAL(&F, low));
            if (value < 0.0) break;
            ++H;
        }
    }

    // root finding by iteration.

    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, low, high);

    const uint maxIter(100);
    for (uint i(0);;++i)
    {
        gsl_root_fsolver_iterate(solver);
        low = gsl_root_fsolver_x_lower(solver);
        high = gsl_root_fsolver_x_upper(solver);
        const int status(gsl_root_test_interval(low, high, 1e-15, TOLERANCE));

        if (status != GSL_CONTINUE) break;
        if (i >= maxIter) throw std::runtime_error("GreensFunction3DRadAbs: drawR: failed to converge");
    }

    const Real r(gsl_root_fsolver_root(solver));
    gsl_root_fsolver_free(solver);
    return r;
}

Real GreensFunction3DRadAbs::p_n_alpha(uint i, uint n, Real r, Real t) const
{
    const Real mDt(-D * t);

    const Real alpha(getAlpha(n, i));
    const Real alphasq(alpha * alpha);

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real hSigma(h * sigma);
    const Real hSigma_m_n(hSigma - n);

    const Real term1(alphasq * alphasq * exp(mDt * alphasq));

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n, sigmaAlpha));
    const Real js2(s.j(n + 1, sigmaAlpha));
    const Real ja(s.j(n, aAlpha));
    const Real ya(s.y(n, aAlpha));
    const Real jr(s.j(n, r * alpha));
    const Real yr(s.y(n, r * alpha));
    const Real jr0(s.j(n, r0 * alpha));
    const Real yr0(s.y(n, r0 * alpha));

    const Real J(hSigma_m_n * js1 + sigmaAlpha * js2);
    const Real Jsq(J * J);

    const Real JY1(ja * yr - ya * jr);
    const Real JY2(ja * yr0 - ya * jr0);

    const Real num(Jsq * JY1 * JY2);
    const Real den1(a * (n + n * n - sigma * (h + h * h * sigma + sigma * alphasq)) * ja * ja);
    const Real den2(sigma * Jsq);
    const Real den(den1 + den2);
    const Real result(term1 * num / den);
    return result;
}

Real GreensFunction3DRadAbs::p_n(int n, Real r, Real t, Real max_alpha) const
{
    const uint min_i(2);
    Real p(0.0);

    for (uint i(0);;)
    {
        const Real alpha(getAlpha(n, i));
        const Real p_i(p_n_alpha(i, n, r, t));
        p += p_i;

        if (alpha >= max_alpha && i >= min_i) break;
        if (i == MAX_ALPHA_SEQ) break;

        ++i;
    }

    return p;
}

void GreensFunction3DRadAbs::makep_nTable(RealVector& p_nTable, Real r, Real t) const
{
    p_nTable.clear();

    const Real factor(a * sigma / (M_PI * 2));
    const Real Dt(D * t);
    const Real alpha00(getAlpha(0, 0));
    const Real max_alpha(sqrt(alpha00 * alpha00 - log(THETA_TOLERANCE * 1e-1) / Dt));

    const Real p_0(p_n(0, r, t, max_alpha) * factor);

    p_nTable.push_back(p_0);

    if (p_0 == 0) return;
    const Real threshold(fabs(THETA_TOLERANCE * p_0));
    Real p_n_prev_abs(fabs(p_0));
    for (uint n(1);;)
    {
        if (getAlpha(n, 0) >= max_alpha)
            break;

        Real p_n(this->p_n(n, r, t, max_alpha) * factor);

        p_nTable.push_back(p_n);
        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if (p_n_abs < threshold && p_n_prev_abs < threshold && p_n_abs <= p_n_prev_abs)
            break;

        n++;
        if (n >= MAX_ORDER) break;
        p_n_prev_abs = p_n_abs;
    }
}

Real GreensFunction3DRadAbs::dp_n_alpha_at_a(uint i, uint n, Real t) const
{
    const Real mDt(-D * t);
    const Real alpha(getAlpha(n, i));

    const Real alphasq(alpha * alpha);

    const Real aAlpha(a * alpha);
    const Real sigmaAlpha(sigma * alpha);
    const Real hSigma(h * sigma);
    const Real hSigma_m_n(hSigma - n);

    const Real term1(alphasq * alpha * exp(mDt * alphasq));

    const SphericalBesselGenerator& s(SphericalBesselGenerator::instance());

    const Real js1(s.j(n, sigmaAlpha));
    const Real js2(s.j(n + 1, sigmaAlpha));
    const Real ja(s.j(n, aAlpha));
    const Real ya(s.y(n, aAlpha));
    const Real jr0(s.j(n, r0 * alpha));
    const Real yr0(s.y(n, r0 * alpha));

    const Real J(hSigma_m_n * js1 + sigmaAlpha * js2);
    const Real Jsq(J * J);

    const Real JY(-jr0 * ya + ja * yr0);

    const Real num(Jsq * JY);

    const Real den1(a * (n + n * n - sigma * (h + h * h * sigma + sigma * alphasq)) * ja * ja);
    const Real den2(sigma * Jsq);
    const Real den(den1 + den2);
    const Real result(term1 * num / den);
    return result;
}

Real GreensFunction3DRadAbs::dp_n_at_a(int n, Real t, Real max_alpha) const
{
    const uint min_i(2);
    Real p(0.0);
    for (uint i(0);;)
    {
        const Real alpha(getAlpha(n, i));
        const Real p_i(dp_n_alpha_at_a(i, n, t));

        p += p_i;
        if (alpha >= max_alpha && i >= min_i) break;
        if (i == MAX_ALPHA_SEQ) break;
        ++i;
    }
    return p;
}

void GreensFunction3DRadAbs::makedp_n_at_aTable(RealVector& p_nTable, Real t) const
{
    p_nTable.clear();

    const Real factor(D * sigma / (a * M_PI * 2));
    const Real Dt(D * t);
    const Real alpha00(getAlpha(0, 0));
    const Real max_alpha(sqrt(Dt * alpha00 * alpha00 - log(THETA_TOLERANCE * 1e-1) / Dt));
    const Real p_0(dp_n_at_a(0, t, max_alpha) * factor);

    p_nTable.push_back(p_0);

    if (p_0 == 0) return;

    const Real threshold(fabs(THETA_TOLERANCE * p_0));

    Real p_n_prev_abs(fabs(p_0));
    for (uint n(1);;)
    {
        if (getAlpha(n, 0) >= max_alpha) break;

        Real p_n(dp_n_at_a(n, t, max_alpha) * factor);

        p_nTable.push_back(p_n);

        const Real p_n_abs(fabs(p_n));
        // truncate when converged enough.
        if (p_n_abs < threshold && p_n_prev_abs < threshold && p_n_abs <= p_n_prev_abs)
            break;

        n++;
        if (n >= MAX_ORDER) break;
        p_n_prev_abs = p_n_abs;
    }
}

Real GreensFunction3DRadAbs::p_theta(Real theta, Real r, Real t) const
{
    if (!(theta >= 0.0 && theta <= M_PI))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
    }

    // r \in (sigma, a);  not defined at r == sigma and r == a.
    if (!(r >= sigma && r < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r >= sigma && r < a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t >= 0.0 : t=%.16g") % t).str());
    }

    if (t == 0.0) return 0.0;

    RealVector p_nTable;
    makep_nTable(p_nTable, r, t);
    const Real p(p_theta_table(theta, r, t, p_nTable));
    return p;
}

Real GreensFunction3DRadAbs::dp_theta(Real theta, Real r, Real t) const
{
    if (!(theta >= 0.0 && theta <= M_PI))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
    }

    // r \in [ sigma, a ]  ;  unlike p_theta
    // defined at r == sigma and r == a.
    if (!(r >= sigma && r <= a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r >= sigma && r <= a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t >= 0.0 : t=%.16g") % t).str());
    }

    if (t == 0.0)return 0.0;

    RealVector p_nTable;
    makedp_n_at_aTable(p_nTable, t);
    const Real p(p_theta_table(theta, r, t, p_nTable));
    return p;
}

static Real p_theta_n(uint n, RealVector const& p_nTable, RealVector const& lgndTable)
{
    return p_nTable[n] * lgndTable[n] * (2 * n + 1);
}

Real GreensFunction3DRadAbs::p_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const
{
    const uint tableSize(p_nTable.size());
    RealVector lgndTable(tableSize);
    gsl_sf_legendre_Pl_array(tableSize - 1, cos(theta), &lgndTable[0]);
    return funcSum_all(boost::bind(&p_theta_n, _1, p_nTable, lgndTable), tableSize) * sin(theta);
}

void GreensFunction3DRadAbs::make_p_thetaTable(RealVector& pTable, Real r, Real t, uint n, RealVector const& p_nTable) const
{
    const Real thetaStep(M_PI / n);

    pTable.push_back(0.0);

    Real p_prev(0.0);
    for (uint i(1);;)
    {
        const Real theta(thetaStep * i);

        Real p(p_theta_table(theta, r, t, p_nTable));

        if (p < 0.0)
        {
            log_.info("drawTheta: p<0 %.16g", p);
            p = 0.0;
        }

        const Real value((p_prev + p) * 0.5);
        pTable.push_back(*(pTable.end() - 1) + value);

        if (/* value < pTable[i] * std::numeric_limits<Real>::epsilon() || */
            i >= n - 1)
        {
            break;   // pTable is valid in [0,i].
        }

        p_prev = p;
        ++i;
    }

}

Real GreensFunction3DRadAbs::ip_theta(Real theta, Real r, Real t) const
{
    if (!(theta >= 0.0 && theta <= M_PI))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
    }

    // r \in (sigma, a)
    if (!(r >= sigma && r < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r >= sigma && r < a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t >= 0.0 : t=%.16g") % t).str());
    }
    if (t == 0.0 || theta == 0.0) return 0.0;

    RealVector p_nTable;
    makep_nTable(p_nTable, r, t);
    const Real p(ip_theta_table(theta, r, t, p_nTable));
    return p;
}

Real GreensFunction3DRadAbs::idp_theta(Real theta, Real r, Real t) const
{
    if (!(theta >= 0.0 && theta <= M_PI))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: theta >= 0.0 && theta <= M_PI : theta=%.16g, M_PI=%.16g") % theta % M_PI).str());
    }

    // r \in [ sigma, a ]
    if (!(r >= sigma && r <= a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r >= sigma && r <= a : r=%.16g, sigma=%.16g, a=%.16g") % r % sigma % a).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t >= 0.0 : t=%.16g") % t).str());
    }

    if (t == 0.0 || theta == 0.0) return 0.0;

    RealVector p_nTable;
    makedp_n_at_aTable(p_nTable, t);
    const Real p(ip_theta_table(theta, r, t, p_nTable));
    return p;
}

static Real ip_theta_n(uint n, RealVector const& p_nTable, RealVector const& lgndTable1)
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1(lgndTable1[n]);   // n-1
    const Real lgnd_n_p1(lgndTable1[n + 2]); // n+1

    // the term (1 + 2 n) is canceled out.
    return p_nTable[n] * (lgnd_n_m1 - lgnd_n_p1);
}

Real GreensFunction3DRadAbs::ip_theta_table(Real theta, Real r, Real t, RealVector const& p_nTable) const
{
    const uint tableSize(p_nTable.size());
    const Real cos_theta(cos(theta));

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable(tableSize + 2);
    lgndTable[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array(tableSize, cos_theta, &lgndTable[1]);

    return funcSum_all(boost::bind(&ip_theta_n, _1, p_nTable, lgndTable), tableSize);
}

struct GreensFunction3DRadAbs::ip_theta_params
{
    GreensFunction3DRadAbs const* const gf;
    const Real r;
    const Real t;
    RealVector const& p_nTable;
    const Real value;
};

Real GreensFunction3DRadAbs::ip_theta_F(Real theta, ip_theta_params const* params)
{
    return params->gf->ip_theta_table(theta, params->r, params->t, params->p_nTable) - params->value;
}

Real GreensFunction3DRadAbs::drawTheta(Real rnd, Real r, Real t) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: rnd < 1.0 && rnd >= 0.0 : rnd=%.16g") % rnd).str());
    }

    if (!(r0 >= sigma && r0 < a))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r0 >= sigma && r0 < a : r0=%.16g, sigma=%.16g, a=%.16g") % r0 % sigma % a).str());
    }

    if (!(r >= sigma))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: r >= sigma : r=%.16g, sigma=%.16g") % r % sigma).str());
    }

    if (!(t >= 0.0))
    {
        throw std::invalid_argument((boost::format("GreensFunction3DRadAbs: t >= 0.0 : t=%.16g") % t).str());
    }

    if (t == 0.0) return 0.0;

    const Real high(M_PI);
    RealVector p_nTable;

    if (r >= a)
    {
        //puts("dp");
        makedp_n_at_aTable(p_nTable, t);
    }
    else
    {
        makep_nTable(p_nTable, r, t);
    }

    const Real ip_theta_pi(ip_theta_table(high, r, t, p_nTable));
    ip_theta_params params = { this, r, t, p_nTable, rnd * ip_theta_pi };
    gsl_function F = { reinterpret_cast<double(*)(double, void*)>(&ip_theta_F), &params };
    const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
    gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));
    gsl_root_fsolver_set(solver, &F, 0.0, high);

    const uint maxIter(100);
    for (uint i(0);;++i)
    {
        gsl_root_fsolver_iterate(solver);
        const Real low(gsl_root_fsolver_x_lower(solver));
        const Real high(gsl_root_fsolver_x_upper(solver));
        const int status(gsl_root_test_interval(low, high, 1e-11, THETA_TOLERANCE));

        if (status != GSL_CONTINUE) break;
        if (i >= maxIter) throw std::runtime_error("GreensFunction3DRadAbs: drawTheta: failed to converge");
    }

    Real theta = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);
    return theta;
}

std::string GreensFunction3DRadAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << D << ", r0 = " << r0 << ", sigma = " << sigma << ", a = " << a << ", kf = " << kf << ", h = " << h << std::endl;
    return ss.str();
}

