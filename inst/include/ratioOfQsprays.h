#ifndef __ROQHEADER__
#define __ROQHEADER__

#define CGAL_EIGEN3_ENABLED 1
#include "gmp.h"
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/polynomial_utils.h>
#include <CGAL/Polynomial/Monomial_representation.h>

#include "qspray.h"

typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,1>::Type Poly1;
typedef CGAL::Polynomial_traits_d<Poly1>                     PT1;
typedef std::pair<CGAL::Exponent_vector, PT1::Innermost_coefficient_type> Monomial1;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 2>::Type Poly2;
typedef CGAL::Polynomial_traits_d<Poly2>                     PT2;
typedef std::pair<CGAL::Exponent_vector, PT2::Innermost_coefficient_type> Monomial2;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 3>::Type Poly3;
typedef CGAL::Polynomial_traits_d<Poly3>                     PT3;
typedef std::pair<CGAL::Exponent_vector, PT3::Innermost_coefficient_type> Monomial3;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 4>::Type Poly4;
typedef CGAL::Polynomial_traits_d<Poly4>                     PT4;
typedef std::pair<CGAL::Exponent_vector, PT4::Innermost_coefficient_type> Monomial4;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 5>::Type Poly5;
typedef CGAL::Polynomial_traits_d<Poly5>                     PT5;
typedef std::pair<CGAL::Exponent_vector, PT5::Innermost_coefficient_type> Monomial5;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 6>::Type Poly6;
typedef CGAL::Polynomial_traits_d<Poly6>                     PT6;
typedef std::pair<CGAL::Exponent_vector, PT6::Innermost_coefficient_type> Monomial6;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 7>::Type Poly7;
typedef CGAL::Polynomial_traits_d<Poly7>                     PT7;
typedef std::pair<CGAL::Exponent_vector, PT7::Innermost_coefficient_type> Monomial7;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 8>::Type Poly8;
typedef CGAL::Polynomial_traits_d<Poly8>                     PT8;
typedef std::pair<CGAL::Exponent_vector, PT8::Innermost_coefficient_type> Monomial8;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq, 9>::Type Poly9;
typedef CGAL::Polynomial_traits_d<Poly9>                     PT9;
typedef std::pair<CGAL::Exponent_vector, PT9::Innermost_coefficient_type> Monomial9;


// ---------------------------------------------------------------------------//
using namespace QSPRAY;

static Qspray<gmpq> gcdQsprays(const Qspray<gmpq>& Q1, const Qspray<gmpq>& Q2) {
  return scalarQspray<gmpq>(1);
}

static Qspray<gmpq> QuotientOfQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B) {
  return qsprayDivision(A, B).first;
}

namespace RATIOOFQSPRAYS {

  namespace utils {

    // -------------------------------------------------------------------------- //
    static std::string Gmpq2str(CGAL::Gmpq r) {
      CGAL::Gmpz numer = r.numerator();
      CGAL::Gmpz denom = r.denominator();
      size_t n = mpz_sizeinbase(numer.mpz(), 10) + 2;
      size_t d = mpz_sizeinbase(denom.mpz(), 10) + 2;
      char* cnumer = new char[n];
      char* cdenom = new char[d];
      cnumer = mpz_get_str(cnumer, 10, numer.mpz());
      cdenom = mpz_get_str(cdenom, 10, denom.mpz());
      std::string snumer = cnumer;
      std::string sdenom = cdenom;
      delete[] cnumer;
      delete[] cdenom;
      return snumer + "/" + sdenom;
    }

    template <typename PolyX, typename PTX, typename MonomialX>
    std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients(Qspray<gmpq> Q1, Qspray<gmpq> Q2) {
      typename PTX::Construct_polynomial constructPolynomial;
      typename std::list<MonomialX> terms1;
      qspray S1 = Q1.get();
      for(const auto& term : S1) {
        powers expnts = term.first;
        terms1.push_back(
          std::make_pair(
            CGAL::Exponent_vector(expnts.begin(), expnts.end()),
            CGAL::Gmpq(QSPRAY::utils::q2str(term.second))
          )
        );
      }
      typename std::list<MonomialX> terms2;
      qspray S2 = Q2.get();
      for(const auto& term : S2) {
        powers expnts = term.first;
        terms2.push_back(
          std::make_pair(
            CGAL::Exponent_vector(expnts.begin(), expnts.end()),
            CGAL::Gmpq(QSPRAY::utils::q2str(term.second))
          )
        );
      }
      typename PTX::Monomial_representation mrepr;
      // first polynomial
      PolyX P1 = constructPolynomial(terms1.begin(), terms1.end());  
      std::list<MonomialX> monomials1;
      mrepr(P1, std::back_inserter(monomials1));
      typename std::list<MonomialX>::iterator it_monoms;
      qspray SS1;
      for(it_monoms = monomials1.begin(); it_monoms != monomials1.end(); it_monoms++) {
        CGAL::Exponent_vector exponents = (*it_monoms).first;
        powers expnts(exponents.begin(), exponents.end());
        gmpq coeff(Gmpq2str((*it_monoms).second));
        SS1[expnts] = coeff;
      }
      // second polynomial
      PolyX P2 = constructPolynomial(terms2.begin(), terms2.end());  
      std::list<MonomialX> monomials2;
      mrepr(P1, std::back_inserter(monomials2));
      typename std::list<MonomialX>::iterator it_monoms;
      qspray SS2;
      for(it_monoms = monomials2.begin(); it_monoms != monomials2.end(); it_monoms++) {
        CGAL::Exponent_vector exponents = (*it_monoms).first;
        powers expnts(exponents.begin(), exponents.end());
        gmpq coeff(Gmpq2str((*it_monoms).second));
        SS2[expnts] = coeff;
      }
      //
      return std::pair<Qspray<gmpq>,Qspray<gmpq>>(Qspray<gmpq>(SS1), Qspray<gmpq>(SS2));
    }

    static std::pair<Qspray<gmpq>,Qspray<gmpq>> simplify(Qspray<gmpq> Q1, Qspray<gmpq> Q2) {
      int d1 = 0;
      qspray S1 = Q1.get();
      for(const auto& term : S1) {
        int n = term.first.size();
        if(n > d1) {
          d1 = n;
        }
      }
      int d2 = 0;
      qspray S2 = Q2.get();
      for(const auto& term : S2) {
        int n = term.first.size();
        if(n > d2) {
          d2 = n;
        }
      }
      const int X = std::max<int>(d1, d2);
      if(X == 4) {
        return getQuotients<Poly4, PT4, Monomial4>(Q1, Q2);
      } else {
        Rcpp::stop("");
      }
    }


  } // end of namespace RATIOOFQSPRAYS::utils

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