#include "gtest/gtest.h"
#include "PROPOSAL/StandardNormal.h"
#include <iostream>
#include <string>

using namespace std;

class RndFromFile{
private:
    double rnd_;
    string Path_;
    ifstream in_;

public:
    RndFromFile(string Path){
        Path_ = Path;
        in_.open(Path_.c_str());
        in_>>rnd_;
        if(!in_.good())cout << "less than one rnd_number!" << endl;
    }

    double rnd(){
        in_>>rnd_;
        if(!in_.good()){
            in_.close();
            in_.clear();
            in_.open(Path_.c_str());
            in_>>rnd_;
        }
        return rnd_;
    }
};


TEST(StandardNormal , StandardNormalRandomNumber ) {

    ifstream in;
    in.open("bin/TestFiles/StandardNormal.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double rnd;
    double average;
    double sigma;
    double xmin;
    double xmax;
    bool   cutoff;
    double StandardNormalRnd;
    double StandardNormalRnd_new;
    double abs_StandardNormalRnd;

    cout.precision(16);

    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    StandardNormal *normal = new StandardNormal(IROMB, IMAXS, IPREC);

    while(in.good())
    {

        rnd = Rand->rnd();

        in>>cutoff>>average>>sigma>>xmin>>xmax>>StandardNormalRnd;
        if(!in.good())break;
        StandardNormalRnd_new = normal->StandardNormalRandomNumber(rnd,average,sigma,xmin,xmax,cutoff);

        abs_StandardNormalRnd=fabs(StandardNormalRnd);

        ASSERT_NEAR(StandardNormalRnd_new, StandardNormalRnd,1e-14*abs_StandardNormalRnd);

    }
}

TEST(StandardNormal , StandardNormalRandomNumber_interpol ) {

    ifstream in;
    in.open("bin/TestFiles/StandardNormal_interpol.txt");

    char firstLine[256];
    in.getline(firstLine,256);

    double rnd;
    double average;
    double sigma;
    double xmin;
    double xmax;
    bool   cutoff;
    double StandardNormalRnd;
    double StandardNormalRnd_new;
    double abs_StandardNormalRnd;


    cout.precision(16);


    RndFromFile* Rand = new RndFromFile("bin/TestFiles/rnd.txt");
    StandardNormal *normal = new StandardNormal(IROMB, IMAXS, IPREC);

    while(in.good())
    {
        normal->EnableInterpolation();

        rnd = Rand->rnd();

        in>>cutoff>>average>>sigma>>xmin>>xmax>>StandardNormalRnd;
        if(!in.good())break;
        StandardNormalRnd_new = normal->StandardNormalRandomNumber(rnd,average,sigma,xmin,xmax,cutoff);

        abs_StandardNormalRnd=fabs(StandardNormalRnd);

        ASSERT_NEAR(StandardNormalRnd_new, StandardNormalRnd,1e-7*abs_StandardNormalRnd);
    }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
