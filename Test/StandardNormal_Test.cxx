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
    in.open("bin/StandardNormal.txt");

    double rnd;
    double average;
    double sigma;
    double xmin;
    double xmax;
    bool   cutoff;
    double StandardNormalRnd;
    double StandardNormalRnd_new;

    cout.precision(16);

    StandardNormal *normal = new StandardNormal(IROMB, IMAXS, IPREC);
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    while(in.good())
    {
        rnd = Rand->rnd();

        in>>cutoff>>average>>sigma>>xmin>>xmax>>StandardNormalRnd;

        StandardNormalRnd_new = normal->StandardNormalRandomNumber(rnd,average,sigma,xmin,xmax,cutoff);

        ASSERT_NEAR(StandardNormalRnd_new, StandardNormalRnd, 1e-7*StandardNormalRnd);

    }
}

TEST(StandardNormal , StandardNormalRandomNumber_interpol ) {

    ifstream in;
    in.open("bin/StandardNormal_interpol.txt");

    double rnd;
    double average;
    double sigma;
    double xmin;
    double xmax;
    bool   cutoff;
    double StandardNormalRnd;
    double StandardNormalRnd_new;

    cout.precision(16);

    StandardNormal *normal = new StandardNormal(IROMB, IMAXS, IPREC);
    normal->EnableInterpolation();
    RndFromFile* Rand = new RndFromFile("bin/rnd.txt");

    while(in.good())
    {
        rnd = Rand->rnd();

        in>>cutoff>>average>>sigma>>xmin>>xmax>>StandardNormalRnd;

        StandardNormalRnd_new = normal->StandardNormalRandomNumber(rnd,average,sigma,xmin,xmax,cutoff);

        ASSERT_NEAR(StandardNormalRnd_new, StandardNormalRnd, 1e-7*StandardNormalRnd);

    }
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
