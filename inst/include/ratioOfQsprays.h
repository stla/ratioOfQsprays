#ifndef __ROQHEADER__
#define __ROQHEADER__

#include <Rcpp.h>
#include "qspray.h"
#include <boost/multiprecision/gmp.hpp>
typedef std::vector<signed int>                             powers;
typedef boost::multiprecision::mpq_rational                 gmpq;


// ---------------------------------------------------------------------------//
using namespace QSPRAY;

static Qspray<gmpq> gcdQsprays(const Qspray<gmpq>& Q1, const Qspray<gmpq>& Q2) {
  return scalarQspray<gmpq>(1);
}

static Qspray<gmpq> QuotientOfQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B) {
  return qsprayDivision(A, B).first;
}

namespace RATIOOFQSPRAYS {

  template<typename T>
  class RatioOfQsprays {

    Qspray<T> numerator;
    Qspray<T> denominator;
    int       dimension;

  public:
    // constructors ---------------
    RatioOfQsprays()
      : numerator(scalarQspray<T>(T(0))), 
        denominator(scalarQspray<T>(T(1))),
        dimension(0)
        {}

    RatioOfQsprays(Qspray<T> numerator_, Qspray<T> denominator_) 
      : numerator(numerator_), 
        denominator(denominator_),
        dimension(
          std::max<int>(
            numerator_.numberOfVariables(), denominator_.numberOfVariables()
          )
        )
        {}

    RatioOfQsprays(int k)
      : numerator(scalarQspray<T>(T(k))), 
        denominator(scalarQspray<T>(T(1))),
        dimension(0)
        {}
    
    // methods --------------------
    Qspray<T> getNumerator() {
      return numerator;
    }

    Qspray<T> getDenominator() {
      return denominator;
    }

    void simplify() {
    	Qspray<T> G = gcdQsprays(numerator, denominator);
    	numerator   = QuotientOfQsprays(numerator, G);
    	denominator = QuotientOfQsprays(denominator, G);
      if(denominator.isConstant()) {
        Qspray<T> d = scalarQspray<T>(T(1) / denominator.constantTerm());
        numerator   *= d;
        denominator *= d;
      }
    }

    RatioOfQsprays<T> operator+=(const RatioOfQsprays<T>& ROQ2) {
    	numerator   = numerator * ROQ2.denominator + denominator * ROQ2.numerator;
    	denominator = denominator * ROQ2.denominator;
    	RatioOfQsprays ROQ(numerator, denominator);
    	ROQ.simplify();
    	return ROQ;
    }

    RatioOfQsprays<T> operator+(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ += ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator-=(const RatioOfQsprays<T>& ROQ2) {
      numerator   = numerator * ROQ2.denominator - denominator * ROQ2.numerator;
      denominator = denominator * ROQ2.denominator;
      RatioOfQsprays ROQ(numerator, denominator);
      ROQ.simplify();
      return ROQ;
    }

    RatioOfQsprays<T> operator-(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ -= ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator*=(const RatioOfQsprays<T>& ROQ2) {
      numerator   = numerator * ROQ2.numerator;
      denominator = denominator * ROQ2.denominator;
      RatioOfQsprays ROQ(numerator, denominator);
      ROQ.simplify();
      return ROQ;
    }

    RatioOfQsprays<T> operator*(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ *= ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator/=(const RatioOfQsprays<T>& ROQ2) {
      numerator   = numerator * ROQ2.denominator;
      denominator = denominator * ROQ2.numerator;
      RatioOfQsprays ROQ(numerator, denominator);
      ROQ.simplify();
      return ROQ;
    }

    RatioOfQsprays<T> operator/(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ /= ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> power(int n) {
      RatioOfQsprays<T> Result(1);
      RatioOfQsprays ROQ(numerator, denominator);
      if(n >= 0) {
        while(n) {
          if(n & 1) {
            Result *= ROQ;
          }
          n >>= 1;
          ROQ *= ROQ;
        }
      } else {
        ROQ.power(-n);
      }
      return Result;
    }

    bool operator==(const RatioOfQsprays<T>& ROQ2) {
      Qspray<T> Q = numerator * ROQ2.denominator + denominator * ROQ2.numerator; 
      return Q.isConstant() && Q.constantTerm() == T(0);
    }

  };

  // -------------------------------------------------------------------------- //
  static Rcpp::List returnRatioOfQsprays(RatioOfQsprays<gmpq> ROQ) {
    return Rcpp::List::create(
      Rcpp::Named("numerator")   = returnQspray(ROQ.getNumerator()),
      Rcpp::Named("denominator") = returnQspray(ROQ.getDenominator())
    );
  }

  // -------------------------------------------------------------------------- //
  static RatioOfQsprays<gmpq> makeRatioOfQsprays(
    const Rcpp::List& Numerator, 
    const Rcpp::List& Denominator
  ) {
    Rcpp::List Powers1 = Numerator["powers"];
    Rcpp::List Powers2 = Denominator["powers"];
    Rcpp::StringVector coeffs1 = Numerator["coeffs"];
    Rcpp::StringVector coeffs2 = Denominator["coeffs"];
    qspray S1;
    for(int i = 0; i < Powers1.size(); i++) {
      Rcpp::IntegerVector Exponents = Powers1(i);
      gmpq coeff(Rcpp::as<std::string>(coeffs1(i)));
      powers pows(Exponents.begin(), Exponents.end());
      S1[pows] = coeff;
    }
    qspray S2;
    for(int i = 0; i < Powers2.size(); i++) {
      Rcpp::IntegerVector Exponents = Powers2(i);
      gmpq coeff(Rcpp::as<std::string>(coeffs2(i)));
      powers pows(Exponents.begin(), Exponents.end());
      S2[pows] = coeff;
    }
    Qspray<gmpq> Q1(S1);
    Qspray<gmpq> Q2(S2);
    RatioOfQsprays<gmpq> ROQ(Q1, Q2);
    return ROQ;
  }

} // end namespace RATIOOFQSPRAYS


// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
using namespace RATIOOFQSPRAYS;

typedef Qspray<RatioOfQsprays<gmpq>>                                   SymbolicQspray;
typedef std::unordered_map<powers, RatioOfQsprays<gmpq>, PowersHasher> symbolicQspray;


static Rcpp::List returnSymbolicQspray(SymbolicQspray SQ) { // used to return a list to R
  symbolicQspray S = SQ.get();
  if(S.size() == 0) {
    return Rcpp::List::create(Rcpp::Named("powers") = R_NilValue,
                              Rcpp::Named("coeffs") = R_NilValue);
  } else {
    Rcpp::List Powers(S.size());
    powers pows;
    unsigned int row = 0, col = 0;
    Rcpp::List Coeffs(S.size());
    unsigned int i = 0;
    for(auto it = S.begin(); it != S.end(); ++it) {
      pows = it->first;
      Rcpp::IntegerVector Exponents(pows.size());
      col = 0;
      for(auto ci = pows.begin(); ci != pows.end(); ++ci) {
        Exponents(col++) = *ci;
      }
      Powers(row++) = Exponents;
      Coeffs(i++) = returnRatioOfQsprays(it->second);
    }
    return Rcpp::List::create(Rcpp::Named("powers") = Powers,
                              Rcpp::Named("coeffs") = Coeffs);
  }
}

// -------------------------------------------------------------------------- //
static SymbolicQspray makeSymbolicQspray(
  const Rcpp::List& Powers, const Rcpp::List& Coeffs
) {
  symbolicQspray S;
  for(int i = 0; i < Powers.size(); i++) {
    Rcpp::IntegerVector Exponents = Powers(i);
    powers pows(Exponents.begin(), Exponents.end());
    Rcpp::List coeff = Coeffs(i);
    Rcpp::List numerator   = coeff["numerator"];
    Rcpp::List denominator = coeff["denominator"];
    S[pows] = makeRatioOfQsprays(numerator, denominator);
  }
  return SymbolicQspray(S);
}



#endif