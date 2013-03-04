#ifndef CTEQPDF_H
#define CTEQPDF_H

#include <string>
#include "PROPOSAL/PhysicsModel.h"
#include "PROPOSAL/Interpolate.h"
#include <vector>

/**
 * Parton distribution functions necessary for neutrino interaction cross sections are calculated here.
 * cteq library is called if the parameterization table ".cteqPDF_raw.data" is not found.
 */

class CteqPDF : public PhysicsModel{


private:


    int num;
    double xn;
    double xm;
    double xq1, xq2;



   // static native void SetCtq6(int set);
   // static native double Ctq6Pdf(int p, double x, double q);

    //static
    bool loaded;
    //static
    int pdf;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize cteq library. Must be compiled and located within the ld library path.
     */

   // static void CteqLoad();

public:
    static std::string tdir;
    static CteqPDF *Ctq;
    double Qn;
    double Qm;

    int f;

    std::vector<Interpolate*> J;
    bool jt;


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with the default PDF set.
     */

    CteqPDF();

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with PDF set specified with i.
     */

    CteqPDF(int i);

    void CteqLoad();

    double CTQ2PDF(int i, double x, double q);

    void init();

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates 3 linear combinations of PDFs. This is the function called by neutrino interaction classes.
     */

    double PDF(int f, double x, double q);


    //----------------------------------------------------------------------------------------------------//

    /**
     * Creates or reads the ".cteqPDF_raw.data" table.
     */

    void interpolate();

    //----------------------------------------------------------------------------------------------------//


    /**
     * 2d parametrization - interface to Interpolate
     */

    double functionInt(double x, double q);

    //----------------------------------------------------------------------------------------------------//

    /**
     * Creates the cteq table by calling the cteq library.
     */

    //static void main(string[] args);
};


#endif // CTEQPDF_H
