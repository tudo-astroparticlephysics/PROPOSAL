#include "MMCexception.h"
#include "MathModel.h"
#include "typeinfo"
#include "Output.h"
#include <cmath>
#include "Mpairproduction.h"
#include "MpairContinuous.h"
#include "PhysicsModel.h"
#include "Amanda.h"
#include "Photonuclear.h"
#include "BremsContinuous.h"
#include "EpairContinuous.h"
#include "PhotoContinuous.h"
#include "IonizContinuous.h"

#include "Photonuclear.h"
#include "Ionizationloss.h"

#include "EpairStochastic.h"
#include "IonizStochastic.h"
#include "BremsStochastic.h"
#include "PhotoStochastic.h"
#include "Mpairproduction.h"
#include "MpairContinuous.h"
#include "PhysicsModel.h"
#include "Amanda.h"
#include <sstream>
#include <fstream>

using namespace std;




std::ifstream randomIn;
//double PhysicsModel::ebig_ =BIGENERGY;

#include <time.h>

int main(){
    cout.precision(16);
    randomIn.open("randomNumbers.txt");

    string opts;
    //    opts="-mediadef=/net/nas.e5.physik.tu-dortmund.de/home/tfuchs/MasterArbeit/IcecubeMMC++/build/src/mediadef -scat -cont -lpm -raw -romb=5 -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2";
    opts="-mediadef=/net/nas.e5.physik.tu-dortmund.de/home/tfuchs/MasterArbeit/IcecubeMMC++/build/src/mediadef         -cont -lpm -raw -romb=5 -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -seed=0";
    double energy;




    int StutzstellenProDekade = 30;
    double MultFak = pow(10,1.0/StutzstellenProDekade);
    int TeilchenProEnergie = 50;
    int Dekaden = 13;


    Propagate* Prop;



    energy = 1e14;
    int Romb=5;


    Propagate::g = Romb;
    Propagate* Prop_pol = new Propagate("ice",500000.,-1,"mu");
    //Prop_pol->contiCorr = true;
    Prop_pol->interpolate("all","");
    Prop = Prop_pol;


    Prop->set_seed(123);
    int init = clock();
    int final;

    double CheckValue=0;
//    for(int k=0;k<StutzstellenProDekade*Dekaden;k++){
//        energy /= MultFak;
//        for(int j = 0;j<TeilchenProEnergie;j++){
                for(int i=0;i<10;i++){CheckValue+= Prop->propagateTo(1e20,energy);
                cout<<i<<endl;}
//                CheckValue+= Prop->propagateTo(1e20,10*energy);
//        }
//    }




//    final = clock()-init;
//    cout << endl << "Zeit:" << ( (double)final / ( (double)CLOCKS_PER_SEC) )/60;
//    cout << endl << "Alles fertig!";
//    cout << endl << "Check: " << CheckValue;

};
