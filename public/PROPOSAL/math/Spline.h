#pragma once

#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include "PROPOSAL/math/Function.h"

namespace PROPOSAL {

class Spline {
   public:
    Spline(std::vector<double>, std::vector<double>);
    Spline(std::vector<Polynom>, std::vector<double>);

    Spline(std::string spline_path);

    Spline(const Spline&);

    virtual Spline* clone() const = 0;

    virtual double evaluate(double x);

    virtual void Derivative();
    virtual void Antiderivative(double c);

    // virtual bool save(std::string, bool) = 0;

    std::pair<double, double> GetDomain();
    std::vector<Polynom> GetFunctions() { return splines_; }
    std::vector<double> GetSubintervalls() { return subintervall_; };
    unsigned int GetNumberIntervalls() { return n_subintervalls_; };

   protected:
    virtual void calculate_splines(std::vector<double> x,
                                   std::vector<double> y) = 0;

    std::vector<Polynom> splines_;
    std::vector<double> subintervall_;
    unsigned int n_subintervalls_;
    std::vector<double> x_;
    std::vector<double> y_;
};

}  // namespace PROPOSAL

struct spline_container {
   public:
    unsigned int degree;
    std::vector<double> coeff;
    double h_i;

    friend std::fstream& operator<<(std::fstream& stream, spline_container& s) {
        stream << s.degree << " " << s.h_i << " " << std::endl;
        for (auto a_i : s.coeff)
            stream << a_i << " ";
        return stream;
    }

    friend std::fstream& operator>>(std::fstream& stream, spline_container& s) {
        stream >> s.degree >> s.h_i;
        double a_i;
        for (unsigned int i = 0; i < s.degree; ++i) {
            stream >> a_i;
            s.coeff.push_back(a_i);
        }
        return stream;
    }
};

//---------------------------------------------------------------------------//
//------------------------------- Linear Spline -----------------------------//
//---------------------------------------------------------------------------//

namespace PROPOSAL {

class Linear_Spline : public Spline {
   public:
    Linear_Spline(std::vector<double>, std::vector<double>);
    Linear_Spline(std::vector<Polynom>, std::vector<double>);
    Linear_Spline(std::string spline_path);
    Linear_Spline(const Linear_Spline&);

    Spline* clone() const override { return new Linear_Spline(*this); };

    // bool save(std::string, bool) override;

   private:
    void calculate_splines(std::vector<double> x,
                           std::vector<double> y) override;
    // std::vector<linear_spline_container> GetSplineContainer();
};

}  // namespace PROPOSAL

//----------------------------------------------------------------------------//
//-------------------------------- Cubic Spline
//------------------------------//
//----------------------------------------------------------------------------//

// struct cubic_spline_container  // : linear_spline_container
// {
//     double h_i;
//     double a_0;
//     double a_1;
//     double a_2;
//     double a_3;

//     friend std::fstream& operator<<(std::fstream& stream,
//                                     cubic_spline_container& s) {
//         stream << s.h_i << " " << s.a_0 << " " << s.a_1 << " " << s.a_2
//         << "
//         "
//                << s.a_3;

//         return stream;
//     }

//     friend std::fstream& operator>>(std::fstream& stream,
//                                     cubic_spline_container& s) {
//         stream >> s.h_i >> s.a_0 >> s.a_1 >> s.a_2 >> s.a_3;
//         return stream;
//     }
// };

namespace PROPOSAL {
class Cubic_Spline : public Spline {
   public:
    Cubic_Spline(std::vector<double>, std::vector<double>);
    Cubic_Spline(std::vector<Polynom>, std::vector<double>);
    Cubic_Spline(std::string spline_path, bool binary);
    Cubic_Spline(const Cubic_Spline&);

    Spline* clone() const override { return new Cubic_Spline(*this); };

    // bool save(std::string, bool) override;

    // std::vector<spline_container> GetSplineContainer();

   private:
    void calculate_splines(std::vector<double> x,
                           std::vector<double> y) override;
};

}  // namespace PROPOSAL
