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
#include <CGAL/number_utils.h>

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
    template <typename PolyX, typename PTX, typename MonomialX, int X>
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {

      CGAL::IO::set_pretty_mode(std::cout);
      
      // CGAL polynomial constructor
      typename PTX::Construct_polynomial constructPolynomial;

      // converts the first Qspray to a CGAL polynomial
      typename std::list<MonomialX> terms1;
      Polynomial<gmpq> S1 = Q1.get();
      for(const auto& term : S1) {
        powers expnts = 
          QSPRAY::utils::growPowers(term.first, term.first.size(), X);
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
      Polynomial<gmpq> S2 = Q2.get();
      for(const auto& term : S2) {
        powers expnts = 
          QSPRAY::utils::growPowers(term.first, term.first.size(), X);
        terms2.push_back(
          std::make_pair(
            CGAL::Exponent_vector(expnts.begin(), expnts.end()),
            CGAL::Gmpq(QSPRAY::utils::q2str(term.second))
          )
        );
      }
      PolyX P2 = constructPolynomial(terms2.begin(), terms2.end());  

      // take the CGAL GCD up to a constant factor

      Rcpp::Rcout << "--------------------\n";
      Rcpp::Rcout << "NUMERATOR: \n";
      std::cout << P1 << "\n";
      Rcpp::Rcout << "DENOMINATOR:\n";
      std::cout << P2 << "\n";
      Rcpp::Rcout << "--------------------\n";

      typename PTX::Gcd_up_to_constant_factor gcd_utcf;
      // (to get the 'ordinary' GCD we would use typename PTX::Gcd gcd)

      Rcpp::Rcout << "--------------------\n";
      Rcpp::Rcout << "GCD\n";

      PolyX D = gcd_utcf(P1, P2);

      std::cout << D << "\n";
      // Rcpp::Rcout << "--------------------\n";
      // Rcpp::Rcout << "DIVISION 1\n";

      // divisions
      PolyX QA = CGAL::integral_division(P1, D);
      PolyX QB = CGAL::integral_division(P2, D);

      // std::cout << QA << "\n";
      Rcpp::Rcout << "_______________________________\n";
      // Rcpp::Rcout << "--------------------\n";
      // Rcpp::Rcout << "DIVISION 2\n";
      // std::cout << QB << "\n";
      Rcpp::Rcout << "_______________________________\n";

      // now make the Qspray corresponding to the QA
      std::list<MonomialX> monomialsA;
      typename PTX::Monomial_representation mrepr;
      mrepr(QA, std::back_inserter(monomialsA));
      Polynomial<gmpq> SA;
      typename std::list<MonomialX>::iterator itmons;
      for(itmons = monomialsA.begin(); itmons != monomialsA.end(); itmons++) {
        CGAL::Exponent_vector exponents = (*itmons).first;
        powers expnts(exponents.begin(), exponents.end());
        QSPRAY::utils::simplifyPowers(expnts);
        gmpq coeff(Gmpq2str((*itmons).second));
        SA[expnts] = coeff;
      }
      // now make the Qspray corresponding to QB
      std::list<MonomialX> monomialsB;
      mrepr(QB, std::back_inserter(monomialsB));
      Polynomial<gmpq> SB;
      for(itmons = monomialsB.begin(); itmons != monomialsB.end(); itmons++) {
        CGAL::Exponent_vector exponents = (*itmons).first;
        powers expnts(exponents.begin(), exponents.end());
        QSPRAY::utils::simplifyPowers(expnts);
        gmpq coeff(Gmpq2str((*itmons).second));
        SB[expnts] = coeff;
      }
      //
      return std::pair<Qspray<gmpq>,Qspray<gmpq>>(Qspray<gmpq>(SA), Qspray<gmpq>(SB));
      // // extract the monomials of the GCD
      // std::list<MonomialX> monomials;
      // typename PTX::Monomial_representation mrepr;
      // mrepr(D, std::back_inserter(monomials));
      // // now make the Qspray corresponding to the GCD
      // qspray S;
      // typename std::list<MonomialX>::iterator itmons;
      // for(itmons = monomials.begin(); itmons != monomials.end(); itmons++) {
      //   CGAL::Exponent_vector exponents = (*itmons).first;
      //   powers expnts(exponents.begin(), exponents.end());
      //   gmpq coeff(Gmpq2str((*itmons).second));
      //   S[expnts] = coeff;
      // }
      // return Qspray<gmpq>(S);
    }

    // static Qspray<gmpq> getGCD1(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly1, PT1, Monomial1, 1>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD2(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly2, PT2, Monomial2, 2>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD3(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly3, PT3, Monomial3, 3>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD4(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly4, PT4, Monomial4, 4>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD5(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly5, PT5, Monomial5, 5>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD6(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly6, PT6, Monomial6, 6>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD7(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly7, PT7, Monomial7, 7>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD8(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly8, PT8, Monomial8, 8>(Q1, Q2);
    // }
    // static Qspray<gmpq> getGCD9(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
    //   return getGCD<Poly9, PT9, Monomial9, 9>(Q1, Q2);
    // }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients1(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly1, PT1, Monomial1, 1>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients2(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly2, PT2, Monomial2, 2>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients3(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly3, PT3, Monomial3, 3>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients4(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly4, PT4, Monomial4, 4>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients5(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly5, PT5, Monomial5, 5>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients6(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly6, PT6, Monomial6, 6>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients7(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly7, PT7, Monomial7, 7>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients8(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly8, PT8, Monomial8, 8>(Q1, Q2);
    }
    static std::pair<Qspray<gmpq>,Qspray<gmpq>> getQuotients9(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      return getQuotients<Poly9, PT9, Monomial9, 9>(Q1, Q2);
    }

    static std::pair<Qspray<gmpq>,Qspray<gmpq>> callGCD(Qspray<gmpq>& Q1, Qspray<gmpq>& Q2) {
      int d1 = Q1.numberOfVariables();
      int d2 = Q2.numberOfVariables();
      const int X = std::max<int>(1, std::max<int>(d1, d2));
      if(X == 1) {
        return getQuotients1(Q1, Q2);
      } else if (X == 2) {
        return getQuotients2(Q1, Q2);
      } else if (X == 3) {
        return getQuotients3(Q1, Q2);
      } else if (X == 4) {
        return getQuotients4(Q1, Q2);
      } else if (X == 5) {
        return getQuotients5(Q1, Q2);
      } else if (X == 6) {
        return getQuotients6(Q1, Q2);
      } else if (X == 7) {
        return getQuotients7(Q1, Q2);
      } else if (X == 8) {
        return getQuotients8(Q1, Q2);
      } else if (X == 9) {
        return getQuotients9(Q1, Q2);
      } else {
        Rcpp::stop("Cannot deal with more than nine variables.");
      }
    }

    // static inline Qspray<gmpq> QuotientOfQsprays(Qspray<gmpq>& A, Qspray<gmpq>& B) {
    //   return qsprayDivision(A, B).first;
    // }

    // template<typename T>
    static inline void simplifyFraction(Qspray<gmpq> &A, Qspray<gmpq> &B) {
      if(!A.isConstant() && !B.isConstant()) {
        std::pair<Qspray<gmpq>,Qspray<gmpq>> QAQB = callGCD(A, B);
        A = QAQB.first;  // is clean
        B = QAQB.second; // is clean
      }
      if(B.isConstant()) {
        gmpq b = B.constantTerm();
        if(b == gmpq(0)) {
          Rcpp::stop("division by zero");
        }
        gmpq lambda = gmpq(1) / b;
        A.scale(lambda);
        B.scale(lambda);
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
      : numerator(Qspray<T>(T(0))), 
        denominator(Qspray<T>(T(1))),
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

    RatioOfQsprays(T x)
      : numerator(Qspray<T>(x)), 
        denominator(Qspray<T>(T(1))),
        dimension(0)
        {}

    RatioOfQsprays(int k)
      : numerator(Qspray<T>(T(k))), 
        denominator(Qspray<T>(T(1))),
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

      // Rcpp::Rcout << "ADDITION\n";

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

    void unsafeMultiply(const RatioOfQsprays<T>& ROQ2) {
      numerator   *= ROQ2.numerator;
      denominator *= ROQ2.denominator;
    }

    RatioOfQsprays<T> operator*=(const RatioOfQsprays<T>& ROQ2) {
      numerator   *= ROQ2.numerator;
      denominator *= ROQ2.denominator;

      // Rcpp::Rcout << "MULTIPLICATION\n";

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
      utils::simplifyFraction(numerator, denominator);
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
        int n0 = n, b = 1, p = 0;
        while(n) {
          if(n & 1) {
            Result.unsafeMultiply(ROQ);
            p += b;
            if(p == n0) {
              break;
            }
          }
          n >>= 1;
          ROQ.unsafeMultiply(ROQ);
          b *= 2;
        }
      } else {
        Result = ROQ.power(-n);
      }
      return Result;
    }

    bool operator==(const RatioOfQsprays<T>& ROQ2) {
      Qspray<T> Q = numerator * ROQ2.denominator - denominator * ROQ2.numerator; 
      return Q.isConstant() && Q.constantTerm() == T(0);
    }

    bool operator==(RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      const RatioOfQsprays<T> ROQ3 = ROQ2;
      return ROQ == ROQ3;
    }

    bool operator!=(const RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      return !(ROQ == ROQ2);
    }

    bool operator!=(RatioOfQsprays<T>& ROQ2) {
      RatioOfQsprays<T> ROQ(numerator, denominator);
      const RatioOfQsprays<T> ROQ3 = ROQ2;
      return ROQ != ROQ3;
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
    Polynomial<gmpq> S1;
    for(int i = 0; i < Powers1.size(); i++) {
      Rcpp::IntegerVector Exponents = Powers1(i);
      gmpq coeff(Rcpp::as<std::string>(coeffs1(i)));
      powers pows(Exponents.begin(), Exponents.end());
      S1[pows] = coeff;
    }
    Polynomial<gmpq> S2;
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

#endif