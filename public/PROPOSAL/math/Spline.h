#pragma once

#include <vector>
#include <string>
#include "PROPOSAL/math/Function.h"

namespace PROPOSAL {

class Spline 
{
    public:

        Spline(std::vector<double>, std::vector<double>);
        Spline(std::vector<Polynom>, std::vector<double>);

        Spline(std::string spline_path);

        Spline(const Spline&);

        virtual double evaluate(double x);

        virtual bool save(std::string, bool) = 0;

    protected:
        virtual void calculate_splines(std::vector<double> x, std::vector<double> y) = 0;

        std::vector<Polynom> splines_;
        std::vector<double> subintervall_;
        unsigned int n_subintervalls_;
        std::vector<double> x_;
        std::vector<double> y_;
};

struct spline_container {
    double a_0;
    double h_i;

    friend std::fstream& operator <<(std::fstream& stream, spline_container& s) {
        stream << s.h_i << " " << s.a_0 << std::endl;
        return stream;
    }

    friend std::fstream& operator >>(std::fstream& stream, spline_container& s) {
        stream >> s.h_i >> s.a_0;
        return stream;
    }
};

} // namespace PROPOSAL

////----------------------------------------------------------------------------//
////------------------------------- Linear Spline ------------------------------//
////----------------------------------------------------------------------------//

//namespace PROPOSAL {

//struct linear_spline_container: spline_container
//{
//    double a_1;
    
//    friend std::fstream& operator <<(std::fstream& stream, linear_spline_container& s) {
//        stream << s.h_i << " " << s.a_0 << " " << s.a_1 << std::endl;
//        return stream;
//    }

//    friend std::fstream& operator >>(std::fstream& stream, linear_spline_container& s) {
//        stream >> s.h_i >> s.a_0 >> s.a_1;
//        return stream;
//    }
//};

//class Linear_Spline : Spline
//{
//    public:

//        Linear_Spline(std::vector<double>, std::vector<double>);
//        Linear_Spline(std::vector<Polynom>, std::vector<double>);

//        Linear_Spline(std::string spline_path);

//        Linear_Spline(const Linear_Spline&);
   
//        bool save(std::string, bool);

//    private:
//        void calculate_splines(std::vector<double> x, std::vector<double> y) override;
//        std::vector<linear_spline_container> GetSplineContainer();

//};

//} // namespace PROPOSAL

//----------------------------------------------------------------------------//
//-------------------------------- Cubic Spline ------------------------------//
//----------------------------------------------------------------------------//

namespace PROPOSAL {

struct cubic_spline_container // : linear_spline_container
{
    double h_i;
    double a_0;
    double a_1;
    double a_2;
    double a_3;
    
    friend std::fstream& operator <<(std::fstream& stream, cubic_spline_container& s) {
        stream << s.h_i << " " << s.a_0 << " " << s.a_1 << " " << s.a_2 << " " << s.a_3;

        return stream;
    }

    friend std::fstream& operator >>(std::fstream& stream, cubic_spline_container& s) {
        stream >> s.h_i >> s.a_0 >> s.a_1 >> s.a_2 >> s.a_3;
        return stream;
    }
};

class Cubic_Spline : public Spline
{
    public:

        Cubic_Spline(std::vector<double>, std::vector<double>);
        Cubic_Spline(std::vector<Polynom>, std::vector<double>);

        Cubic_Spline(std::string spline_path, bool binary);

        Cubic_Spline(const Cubic_Spline&);
        
        bool save(std::string, bool);
    
    private:
        void calculate_splines(std::vector<double> x, std::vector<double> y) override;
        std::vector<cubic_spline_container> GetSplineContainer();

};

} // namespace PROPOSAL
