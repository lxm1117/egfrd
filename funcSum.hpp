#ifndef FUNCSUM_HPP
#define FUNCSUM_HPP

#include <functional>
#include "Defs.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

static const Real TOLERANCE(1e-8);

// --------------------------------------------------------------------------------------------------------------------------------

Real funcSum_all(std::function<Real(uint i)> f, std::size_t max_i);

Real funcSum(std::function<Real(uint i)> f, std::size_t max_i, Real tolerance = TOLERANCE);

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* FUNCSUM_HPP */
