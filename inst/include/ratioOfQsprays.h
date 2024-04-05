#ifndef ___RATIO_OF_QSPRAYS___
#define ___RATIO_OF_QSPRAYS___

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
typedef CGAL::Polynomial_traits_d<Poly1>                    PT1;
typedef std::pair<CGAL::Exponent_vector, PT1::Innermost_coefficient_type> Monomial1;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,2>::Type Poly2;
typedef CGAL::Polynomial_traits_d<Poly2>                    PT2;
typedef std::pair<CGAL::Exponent_vector, PT2::Innermost_coefficient_type> Monomial2;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,3>::Type Poly3;
typedef CGAL::Polynomial_traits_d<Poly3>                    PT3;
typedef std::pair<CGAL::Exponent_vector, PT3::Innermost_coefficient_type> Monomial3;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,4>::Type Poly4;
typedef CGAL::Polynomial_traits_d<Poly4>                    PT4;
typedef std::pair<CGAL::Exponent_vector, PT4::Innermost_coefficient_type> Monomial4;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,5>::Type Poly5;
typedef CGAL::Polynomial_traits_d<Poly5>                    PT5;
typedef std::pair<CGAL::Exponent_vector, PT5::Innermost_coefficient_type> Monomial5;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,6>::Type Poly6;
typedef CGAL::Polynomial_traits_d<Poly6>                    PT6;
typedef std::pair<CGAL::Exponent_vector, PT6::Innermost_coefficient_type> Monomial6;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,7>::Type Poly7;
typedef CGAL::Polynomial_traits_d<Poly7>                    PT7;
typedef std::pair<CGAL::Exponent_vector, PT7::Innermost_coefficient_type> Monomial7;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,8>::Type Poly8;
typedef CGAL::Polynomial_traits_d<Poly8>                    PT8;
typedef std::pair<CGAL::Exponent_vector, PT8::Innermost_coefficient_type> Monomial8;
typedef CGAL::Polynomial_type_generator<CGAL::Gmpq,9>::Type Poly9;
typedef CGAL::Polynomial_traits_d<Poly9>                    PT9;
typedef std::pair<CGAL::Exponent_vector, PT9::Innermost_coefficient_type> Monomial9;


// ---------------------------------------------------------------------------//
using namespace QSPRAY;


namespace RATIOOFQSPRAYS {

  namespace utils {

    // -------------------------------------------------------------------------- //
    static inline std::string Gmpq2str(CGAL::Gmpq r) {
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

    // -------------------------------------------------------------------------- //
    template <typename PolyX, typename PTX, typename MonomialX>
    static Qspray<gmpq> getGCD(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      // CGAL polynomial constructor
      typename PTX::Construct_polynomial constructPolynomial;
      // converts the first Qspray to a CGAL polynomial
      typename std::list<MonomialX> terms1;
      qspray S1 = Q1.get();
      for(const auto& term : S1) {
        powers expnts = 
          QSPRAY::utils::growPowers(term.first, term.first.size(), 4);
        terms1.push_back(
          std::make_pair(
            CGAL::Exponent_vector(expnts.begin(), expnts.end()),
            CGAL::Gmpq(QSPRAY::utils::q2str(term.second))
          )
        );
      }
      PolyX P1 = constructPolynomial(terms1.begin(), terms1.end()); 
      // converts the second Qspray to a CGAL polynomial
      typename std::list<MonomialX> terms2;
      qspray S2 = Q2.get();
      for(const auto& term : S2) {
        powers expnts = 
          QSPRAY::utils::growPowers(term.first, term.first.size(), 4);
        terms2.push_back(
          std::make_pair(
            CGAL::Exponent_vector(expnts.begin(), expnts.end()),
            CGAL::Gmpq(QSPRAY::utils::q2str(term.second))
          )
        );
      }
      PolyX P2 = constructPolynomial(terms2.begin(), terms2.end());  
      // take the CGAL GCD up to a constant factor
      typename PTX::Gcd_up_to_constant_factor gcd_utcf;
      // (to get the 'ordinary' GCD we woule use typename PTX::Gcd gcd)
      PolyX D = gcd_utcf(P1, P2);
      // extract the monomials of the GCD
      std::list<MonomialX> monomials;
      typename PTX::Monomial_representation mrepr;
      mrepr(D, std::back_inserter(monomials));
      // now make the Qspray corresponding to the GCD
      qspray S;
      typename std::list<MonomialX>::iterator itmons;
      for(itmons = monomials.begin(); itmons != monomials.end(); itmons++) {
        CGAL::Exponent_vector exponents = (*itmons).first;
        powers expnts(exponents.begin(), exponents.end());
        gmpq coeff(Gmpq2str((*itmons).second));
        S[expnts] = coeff;
      }
      return Qspray<gmpq>(S);
    }

    static Qspray<gmpq> getGCD1(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly1, PT1, Monomial1>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD2(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly2, PT2, Monomial2>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD3(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly3, PT3, Monomial3>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD4(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly4, PT4, Monomial4>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD5(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly5, PT5, Monomial5>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD6(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly6, PT6, Monomial6>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD7(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly7, PT7, Monomial7>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD8(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly8, PT8, Monomial8>(Q1, Q2);
    }
    static Qspray<gmpq> getGCD9(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getGCD<Poly9, PT9, Monomial9>(Q1, Q2);
    }

    static Qspray<gmpq> callGCD(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      int d1 = Q1.numberOfVariables();
      int d2 = Q2.numberOfVariables();
      // qspray S1 = Q1.get();
      // for(const auto& term : S1) {
      //   int n = term.first.size();
      //   if(n > d1) {
      //     d1 = n;
      //   }
      // }
      // int d2 = 0;
      // qspray S2 = Q2.get();
      // for(const auto& term : S2) {
      //   int n = term.first.size();
      //   if(n > d2) {
      //     d2 = n;
      //   }
      // }
      const int X = std::max<int>(1, std::max<int>(d1, d2));
      if(X == 1) {
        return getGCD1(Q1, Q2);
      } else if (X == 2) {
        return getGCD2(Q1, Q2);
      } else if (X == 3) {
        return getGCD3(Q1, Q2);
      } else if (X == 4) {
        return getGCD4(Q1, Q2);
      } else if (X == 5) {
        return getGCD5(Q1, Q2);
      } else if (X == 6) {
        return getGCD6(Q1, Q2);
      } else if (X == 7) {
        return getGCD7(Q1, Q2);
      } else if (X == 8) {
        return getGCD8(Q1, Q2);
      } else if (X == 9) {
        return getGCD9(Q1, Q2);
      } else {
        Rcpp::stop("Cannot deal with more than nine variables.");
      }
    }

    static Qspray<gmpq> QuotientOfQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B) {
      return qsprayDivision(A, B).first;
    }

    template<typename T>
    static void simplifyFraction(Qspray<T> &A, Qspray<T> &B) {
      Qspray<T> G = callGCD(A, B);
      G.clean();
      A = QuotientOfQsprays(A, G);
      B = QuotientOfQsprays(B, G);
      if(B.isConstant()) {
        Qspray<T> d = scalarQspray<T>(T(1) / B.constantTerm());
        A *= d;
        B *= d;
      }
    }

  } // end of namespace RATIOOFQSPRAYS::utils


  // ------------------------------------------------------------------------ //
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

    RatioOfQsprays<T> operator+=(const RatioOfQsprays<T>& ROQ2) {
    	numerator    = numerator * ROQ2.denominator + denominator * ROQ2.numerator;
    	denominator *= ROQ2.denominator;
      utils::simplifyFraction(numerator, denominator);
    	RatioOfQsprays ROQ(numerator, denominator);
    	return ROQ;
    }

    RatioOfQsprays<T> operator+(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ += ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator-=(const RatioOfQsprays<T>& ROQ2) {
      numerator    = numerator * ROQ2.denominator - denominator * ROQ2.numerator;
      denominator *= ROQ2.denominator;
      utils::simplifyFraction(numerator, denominator);
      RatioOfQsprays ROQ(numerator, denominator);
      return ROQ;
    }

    RatioOfQsprays<T> operator-(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ -= ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator*=(const RatioOfQsprays<T>& ROQ2) {
      numerator   *= ROQ2.numerator;
      denominator *= ROQ2.denominator;
      utils::simplifyFraction(numerator, denominator);
      RatioOfQsprays ROQ(numerator, denominator);
      return ROQ;
    }

    RatioOfQsprays<T> operator*(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      ROQ *= ROQ2;
      return ROQ;
    }

    RatioOfQsprays<T> operator/=(const RatioOfQsprays<T>& ROQ2) {
      numerator   *= ROQ2.denominator;
      denominator *= ROQ2.numerator;
      utils::simplifyFraction(denominator, numerator);
      RatioOfQsprays ROQ(numerator, denominator);
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

  // ------------------------------------------------------------------------ //
  static Rcpp::List returnRatioOfQsprays(RatioOfQsprays<gmpq> ROQ) {
    return Rcpp::List::create(
      Rcpp::Named("numerator")   = returnQspray(ROQ.getNumerator()),
      Rcpp::Named("denominator") = returnQspray(ROQ.getDenominator())
    );
  }

  // ------------------------------------------------------------------------ //
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


static Rcpp::List returnSymbolicQspray(SymbolicQspray SQ) { // to return a list to R
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