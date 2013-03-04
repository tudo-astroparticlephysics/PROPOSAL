#include "generation/gen/NeutrinoInt.h"
#include "generation/gen/AtmFlux.h"
#include <iostream>
#include "generation/gen/CteqPDF.h"
#include <sstream>
#include "generation/gen/AtmFlux.h"

using namespace std;


int main(){

    string mmcOpts;
    stringstream ss;


    ss << "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont ";
    ss << "-tdir="<<"/home/jkoehne/nobackup/IceCubeSim/V03-02-00PROPOSAL_final2/build/PROPOSAL ";
    ss << "-mediadef="<<"/home/jkoehne/nobackup/IceCubeSim/V03-02-00PROPOSAL_final2/build/PROPOSAL-icetray/resources/mediadef  ";
    ss << "-radius=600 -length=1200 ";
    ss << "-seed=1";
    ss << "-prop";

    mmcOpts=ss.str();
    AtmFlux *test = new AtmFlux();
    test->setup(mmcOpts);

    return 0;

}

//    stringstream name;
//    NeutrinoInt *in = new NeutrinoInt();

//    ofstream out;
//    double E;
//    double q;
//    double x[]={
//        0.65,
//        0.55,
//        0.35,
//        0.275,
//        0.175,
//        0.125,
//        0.045,
//        0.015
//    };

//    double y;
//    double wq;


//     E=45.;
//     for(int j =0; j<8; j++){
//         name<<"/home/koehne/Daten/"<<"diff_nu"<<j<<"_"<<E<<".txt";
//         out.open(name.str().c_str());

//        for(int i =1 ;i<11;i++){

//            q=sqrt(x[j]*2*1.e-3*in->m->get_MM()*y*E);
//            y=0.1*i;
//            wq=in->dS2dxdy(q, y ,E,true,true);
//            if(i!=1)
//            out<<wq/E<<"\t"<<y<<endl;



//        }
//        out.close();
//        name.str("");
//        name.clear();

//        name<<"/home/koehne/Daten/"<<"diff_antinu"<<j<<"_"<<E<<".txt";
//        out.open(name.str().c_str());

//       for(int i =1 ;i<11;i++){

//           q=sqrt(x[j]*2*1.e-3*in->m->get_MM()*y*E);
//           y=0.1*i;
//           wq=in->dS2dxdy(q, y ,E,true,false);
//           if(i!=1)
//           out<<wq/E<<"\t"<<y<<endl;



//       }
//       out.close();
//       name.str("");
//       name.clear();

//     }

//     E=150.;

//     for(int j =0; j<8; j++){
//         name<<"/home/koehne/Daten/"<<"diff_nu"<<j<<"_"<<E<<".txt";
//         out.open(name.str().c_str());

//        for(int i =1 ;i<11;i++){

//            q=sqrt(x[j]*2*1.e-3*in->m->get_MM()*y*E);
//            y=0.1*i;
//            wq=in->dS2dxdy(q, y ,E,true,true);
//            if(i!=1)
//            out<<wq/E<<"\t"<<y<<endl;



//        }
//        out.close();
//        name.str("");
//        name.clear();

//        name<<"/home/koehne/Daten/"<<"diff_antinu"<<j<<"_"<<E<<".txt";
//        out.open(name.str().c_str());

//       for(int i =1 ;i<11;i++){

//           q=sqrt(x[j]*2*1.e-3*in->m->get_MM()*y*E);
//           y=0.1*i;
//           wq=in->dS2dxdy(q, y ,E,true,false);
//           if(i!=1)
//           out<<wq/E<<"\t"<<y<<endl;



//       }
//       out.close();
//       name.str("");
//       name.clear();

//     }
//    return 0;
//}


//int main(){


//    CteqPDF *pdf = new CteqPDF();
//    pdf->CteqLoad();
//    //double q;
//    ofstream out;
//    out.open("/home/koehne/Daten/pdf_test.txt");
//    ofstream out2;
//    out2.open("/home/koehne/Daten/pdf_test2.txt");

//    out.precision(16);
//    out2.precision(16);
//    double E;

//    double x[]={
//        0.85,
//        0.75,
//        0.65,
//        0.4,
//        0.25,
//        0.18,
//        0.13,
//        0.08,
//        0.05,
//        0.032,
//        0.021,
//        0.013,
//        0.008,
//        0.005,
//        0.0032,
//        0.00253,
//        0.0021,
//        0.00161,
//        0.0013,
//        0.00102,
//        0.0008,
//        0.000632,
//        0.0005,
//        0.0004,
//        0.000253,
//        0.000162,
//        0.000102,
//        0.000063

//    };
//    double F2;

//    ofstream out3;
//    stringstream outname;
//    double q;


//    for(int i=0; i<28;i++){
//        outname<<"/home/koehne/Daten/"<<"pdf"<<i<<".txt";
//        out3.open(outname.str().c_str());

//        for(int j =0;j<1e4;j++){

//            q=1000*1e-4*(j+1);
//            F2=pdf->PDF(2,x[i],q);
//            out3<<F2<<"\t"<<q<<endl;

//        }
//        outname.str("");
//        outname.clear();
//        out3.close();

//    }

   // for(int i = )

//    NeutrinoTot * tot = new NeutrinoTot();


//    tot->interpolate("");
//    for(int j =0;j<1e4;j++){

//            E= 1. + 349.*1e-4*j;
//            if(j!=0){
//            out<<E<<"\t"<<tot->dSdy(E,true,true)/E<<endl;
//            out2<<E<<"\t"<<tot->dSdy(E,true,false)/E<<endl;
//            }

//    }
//    return 0;
//}
