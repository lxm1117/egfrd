#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <vector>
#include <cmath>
#include <boost/bind.hpp>
#include <gsl/gsl_sum.h>

#include "Logger.hpp"
#include "funcSum.hpp"

typedef std::vector<Real> RealVector;

static Logger& _log(Logger::get_logger("funcSum"));

Real 
funcSum_all(boost::function<Real(unsigned int i)> f, size_t max_i)
{
    Real sum(0.0);

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    sum = p_0;

    RealVector::size_type i(1); 
    while(i < max_i)
    {
        const Real p_i(f(i));
        sum += p_i;

        ++i;
    }

    return sum;
}


Real 
funcSum_all_accel(boost::function<Real(unsigned int i)> f,
                  size_t max_i, Real tolerance)
{
    RealVector pTable;
    pTable.reserve(max_i);

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    pTable.push_back(p_0);

    RealVector::size_type i(1);
    for(;  i < max_i; ++i)
    {
        const Real p_i(f(i));
        pTable.push_back(p_i);
    }

    Real sum;
    Real error;
    gsl_sum_levin_utrunc_workspace* 
        workspace(gsl_sum_levin_utrunc_alloc(i));
    gsl_sum_levin_utrunc_accel(&pTable[0], pTable.size(), workspace, 
                                &sum, &error);
    if (fabs(error) >= fabs(sum * tolerance))
    {
        _log.error("series acceleration error: %.16g"
                  " (rel error: %.16g), terms_used = %d (%d given)",
                  fabs(error), fabs(error / sum),
                  workspace->terms_used, pTable.size());
        // TODO look into this crashing behaviour
    }

    gsl_sum_levin_utrunc_free(workspace);

    return sum;
}


Real 
funcSum(boost::function<Real(unsigned int i)> f, size_t max_i, Real tolerance)
// funcSum
// ==
// Will simply calculate the sum over a certain function f, until it converges
// (i.e. the sum > tolerance*current_term for a CONVERGENCE_CHECK number of 
// terms), or a maximum number of terms is summed (usually 2000).
//
// Input:
// - f: A       function object
// - max_i:     maximum number of terms it will evaluate
// - tolerance: convergence condition, default value 1-e8 (see .hpp)
// 
// About Boost::function
// ===
// Boost::function doesn't do any type checking: It will take any object 
// and any signature you provide in its template parameter, and create an 
// object that's callable according to your signature and calls the object. 
// If that's impossible, it's a compile error.
// (From: http://stackoverflow.com/questions/527413/how-boostfunction-and-boostbind-work)
{
    // DEFAULT = 4
    const unsigned int CONVERGENCE_CHECK(4);    

    Real sum(0.0);
    RealVector pTable;

    const Real p_0(f(0));
    if (p_0 == 0.0)
    {
        return 0.0;
    }

    pTable.push_back(p_0);
    sum = p_0;

    bool extrapolationNeeded(true);

    unsigned int convergenceCounter(0);

    RealVector::size_type i(1); 
    while(i < max_i)
    {
        const Real p_i(f(i));
        pTable.push_back(p_i);
        sum += p_i;

        ++i;

        if (fabs(sum) * tolerance >= fabs(p_i)) // '=' is important
        {
            ++convergenceCounter;          
        }
        // this screws it up; why?
        else
        {
            convergenceCounter = 0;
        }


        if (convergenceCounter >= CONVERGENCE_CHECK)
        {
            extrapolationNeeded = false;
            break;
        }
        
    }

    if (extrapolationNeeded)
    {
        Real error;
        gsl_sum_levin_utrunc_workspace* 
            workspace(gsl_sum_levin_utrunc_alloc(i));
        gsl_sum_levin_utrunc_accel(&pTable[0], pTable.size(), workspace, 
            &sum, &error);
        if (fabs(error) >= fabs(sum * tolerance * 10))
        {
            _log.error("series acceleration error: %.16g"
                      " (rel error: %.16g), terms_used = %d (%d given)",
                      fabs(error), fabs(error / sum),
                      workspace->terms_used, pTable.size());
        }

        gsl_sum_levin_utrunc_free(workspace);
    }

    return sum;
}

#if defined(ENABLE_GF_TESTFUNCTIONS)

#include <gsl/gsl_rng.h>
#include "GreensFunction.hpp"
#include "GreensFunction3D.hpp"
#include "GreensFunction3DAbs.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "GreensFunction3DRadInf.hpp"

double GF_CLASS DrawGreensFunction(unsigned int gfn, unsigned int gfc, double rnd, double D, double r0, double a, double kf, double sigma, double t, double r)
{
   switch (gfn)
   {
   case 0:
   { auto gf = GreensFunction3DAbsSym(D, a);
   try
   {
      switch (gfc)
      {
      case 0: return gf.drawTime(rnd);
      case 1: return gf.drawR(rnd, t);
      case 2: return 0.0; // (double)(int)gf.drawEventType(rnd, t); 
      case 3: return 0.0; // gf.drawTheta(rnd, r, t); 
      }
   }
   catch (...)
   {
      return NAN;
   }
   }break;

   case 1: { auto gf = GreensFunction3D(D, r0);
      try
      {
         switch (gfc)
         {
         case 0: return gf.drawTime(rnd);
         case 1: return gf.drawR(rnd, t);
         case 2: return 0.0; //(double)(int)gf.drawEventType(rnd, t); 
         case 3: return gf.drawTheta(rnd, r, t);
         }
      }
      catch (...)
      {
         return NAN;
      }
   } break;
   case 2: { auto gf = GreensFunction3DAbs(D, r0, a);
      try
      {
         switch (gfc)
         {
         case 0: return gf.drawTime(rnd);
         case 1: return gf.drawR(rnd, t);
         case 2: return (double)(int)gf.drawEventType(rnd, t);
         case 3: return gf.drawTheta(rnd, r, t);
         }
      }
      catch (...)
      {
         return NAN;
      }
   } break;
   case 3: { auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
      try
      {
         switch (gfc)
         {
         case 0: return gf.drawTime(rnd);
         case 1: return gf.drawR(rnd, t);
         case 2: return (double)(int)gf.drawEventType(rnd, t);
         case 3: return gf.drawTheta(rnd, r, t);
         }
      }
      catch (...)
      {
         return -9999.0;
      }
   } break;
   case 4: { auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);
      try
      {
         switch (gfc)
         {
         case 0: return gf.drawTime(rnd);
         case 1: return gf.drawR(rnd, t);
         case 2: return 0.0; //(double)(int)gf.drawEventType(rnd, t); 
         case 3: return gf.drawTheta(rnd, r, t);
         }
      }
      catch (...)
      {
         return NAN;
      }
   }break;
   }
}


void GF_CLASS TestGreensFunction(unsigned int gfn, unsigned int gfc, unsigned int size, double *buffer, double D, double r0, double a, double kf, double sigma, double t, double r)
{
   gsl_set_error_handler([](const char*, const char*, int, int) { throw  std::runtime_error("GSL"); });

   switch (gfn)
   {
   case 0:
   { auto gf = GreensFunction3DAbsSym(D, a);
   for (unsigned int i = 0; i < size; ++i)
   {
      double rnd = (double)i / size;
      try
      {
         switch (gfc)
         {
         case 0: buffer[i] = gf.drawTime(rnd); break;
         case 1: buffer[i] = gf.drawR(rnd, t); break;
         case 2: buffer[i] = 0.0; // (double)(int)gf.drawEventType(rnd, t); break;
         case 3: buffer[i] = 0.0; // gf.drawTheta(rnd, r, t); break;
         }
      }
      catch (...)
      {
         buffer[i] = NAN;
      }
   }
   }break;

   case 1: { auto gf = GreensFunction3D(D, r0);
      for (unsigned int i = 0; i < size; ++i)
      {
         double rnd = (double)i / size;
         try
         {
            switch (gfc)
            {
            case 0: buffer[i] = gf.drawTime(rnd); break;
            case 1: buffer[i] = gf.drawR(rnd, t); break;
            case 2: buffer[i] = 0.0; //(double)(int)gf.drawEventType(rnd, t); break;
            case 3: buffer[i] = gf.drawTheta(rnd, r, t); break;
            }
         }
         catch (...)
         {
            buffer[i] = NAN;
         }
      }
   } break;
   case 2: { auto gf = GreensFunction3DAbs(D, r0, a);
      for (unsigned int i = 0; i < size; ++i)
      {
         double rnd = (double)i / size;
         try
         {
            switch (gfc)
            {
            case 0: buffer[i] = gf.drawTime(rnd); break;
            case 1: buffer[i] = gf.drawR(rnd, t); break;
            case 2: buffer[i] = (double)(int)gf.drawEventType(rnd, t); break;
            case 3: buffer[i] = gf.drawTheta(rnd, r, t); break;
            }
         }
         catch (...)
         {
            buffer[i] = NAN;
         }
      }
   } break;
   case 3: { auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
      for (unsigned int i = 0; i < size; ++i)
      {
         double rnd = (double)i / size;
         try
         {
            switch (gfc)
            {
            case 0: buffer[i] = gf.drawTime(rnd); break;
            case 1: buffer[i] = gf.drawR(rnd, t); break;
            case 2: buffer[i] = (double)(int)gf.drawEventType(rnd, t); break;
            case 3: buffer[i] = gf.drawTheta(rnd, r, t); break;
            }
         }
         catch (...)
         {
            buffer[i] = -9999.0;
         }
      }
   } break;
   case 4: { auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);
      for (unsigned int i = 0; i < size; ++i)
      {
         double rnd = (double)i / size;
         try
         {
            switch (gfc)
            {
            case 0: buffer[i] = gf.drawTime(rnd); break;
            case 1: buffer[i] = gf.drawR(rnd, t); break;
            case 2: buffer[i] = 0.0; //(double)(int)gf.drawEventType(rnd, t); break;
            case 3: buffer[i] = gf.drawTheta(rnd, r, t); break;
            }
         }
         catch (...)
         {
            buffer[i] = NAN;
         }
      }
   }break;
   }
}

#endif
