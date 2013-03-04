#include "generation/gen/CteqPDF.h"
#include <algorithm>
#include <cmath>
#include "PROPOSAL/Output.h"
#include <iostream>
#include <sstream>

using namespace std;


extern "C" {
   double ctq6pdf_(const int&, double&, double&);
   void setctq6_(int&);
   }


   std::string CteqPDF::tdir="";
   CteqPDF* CteqPDF::Ctq=NULL;

/**
 * Parton distribution functions necessary for neutrino interaction cross sections are calculated here.
 * cteq library is called if the parameterization table ".cteqPDF_raw.data" is not found.
 */


 //   static native void SetCtq6(int set);
 //   static native double Ctq6Pdf(int p, double x, double q);


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize cteq library. Must be compiled and located within the ld library path.
     */

    void CteqPDF::CteqLoad(){
        if(!loaded){
            loaded=true;


                setctq6_(pdf);


        }
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with the default PDF set.
     */

    CteqPDF::CteqPDF(){

        init();
        Ctq=this;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with PDF set specified with i.
     */

    CteqPDF::CteqPDF(int i){

        init();
        pdf=i;
        Ctq=this;
    }

    void CteqPDF::init(){

        num=3;
        xn=1.e-6;
        xm=1.;
        xq1=xn, xq2=xn*2;
        loaded=false;
        pdf=1;
        tdir="";

        Qn=0.2262;
        Qm=5.7471e4;
        jt=false;

    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates 3 linear combinations of PDFs. This is the function called by neutrino interaction classes.
     */

    double CteqPDF::PDF(int f, double x, double q){
        if(q<Qn) return 0;
        if(x<xq1){
            double pdf1=PDF(f, xq1, q), pdf2=PDF(f, xq2, q);
            if(pdf1<=0 || pdf2<=0) return 0;
            return pdf1*pow(x/xq1, log(pdf2/pdf1)/log(xq2/xq1));
        }
        if(jt) return max(J.at(f-1)->interpolate(x, q), 0.);
        this->f=f;
        return functionInt(x, q);
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Creates or reads the ".cteqPDF_raw.data" table.
     */

    void CteqPDF::interpolate(){
        int g=4;
        int n1=100;
        int n2=20;
        bool raw=true;

        bool rsv;
        rsv=Output::raw;
        Output::raw=raw;

        bool flag;
        stringstream ss;
        ss<<tdir<<".cteqPDF";

        if(Output::raw) ss<<"_raw"; else ss<<"_ascii";
        ss<<".data";

        string name=ss.str();
        ss.clear();
        ss.str("");
        do {
            if(Output::texi) return;
            flag=false;
            try{
                Output::open(name);

                jt=false;
                cerr<<"Parameterizing ctqpdF ... "<<endl;

                J.resize(num);
                for(int i=0; i<num; i++){
                    f=i+1;
                    J.at(i) = new Interpolate(n1, xn, xm, n2, Qn, Qm, this, g, false, false, true, g, false, false, true, g, false, false, false);
                }
                jt=true;
                cerr<<"done"<<endl;
                cerr<<"Finished parameterizations"<<endl;

                Output::close();
            }catch(int a){
                throw 0;
                flag=true;
                Output::Delete(name);
            }
        } while(flag);

        Output::raw=rsv;


        double step=log(xm/xn)/n1;
        xq1=xn*exp(step/2);
        xq2=xq2*exp(step);
    }

    //----------------------------------------------------------------------------------------------------//



    /**
     * 2d parametrization - interface to Interpolate
     */

    double CteqPDF::functionInt(double x, double q){
        if(q<Qn) return 0;
        int i;
        double sum;
        CteqLoad();

        switch(f){
        case 1:    // q-q~
            sum=0;
            for(i=1; i<6; i++) sum+=ctq6pdf_(i, x, q);
            for(i=1; i<6; i++) sum-=ctq6pdf_(-i, x, q);
            break;
        case 2:    // q+q~
            sum=0;
            for(i=1; i<6; i++) sum+=ctq6pdf_(i, x, q);
            for(i=1; i<6; i++) sum+=ctq6pdf_(-i, x, q);
            break;
        case 3:    // 2(s-c+b)
        default:
            sum=2*(ctq6pdf_(3, x, q)-ctq6pdf_(4, x, q)+ctq6pdf_(5, x, q));
        }
        return x*sum;
    }

    //----------------------------------------------------------------------------------------------------//

    double CteqPDF::CTQ2PDF(int i, double x, double q){

        return ctq6pdf_(i, x, q);

    }

    /**
     * Creates  the cteq table by calling the cteq library.
     */

//    static void main(string[] args){
//        CteqPDF F = new CteqPDF();
//        F.interpolate();
//    }

