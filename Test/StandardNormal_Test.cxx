#include "gtest/gtest.h"
#include "PROPOSAL/StandardNormal.h"
#include <iostream>
#include <string>
#include "PROPOSAL/Output.h"

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
        if(!in_.good())log_warn("less than one rnd_number!");
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


TEST(Comparison , Comparison_equal ) {
    StandardNormal A;
    StandardNormal B;
    EXPECT_TRUE(A==B);
    StandardNormal* C = new StandardNormal(4,19,1e-3);
    StandardNormal* D = new StandardNormal(4,19,1e-3);
    EXPECT_TRUE(*C==*D);
    StandardNormal* E = new StandardNormal(5,20,1e-6);
    EXPECT_TRUE(A==*E);

}

TEST(Comparison , Comparison_not_equal ) {
    StandardNormal A;
    StandardNormal B(4,19,1e-3);
    EXPECT_TRUE(A!=B);
    StandardNormal* C = new StandardNormal(4,19,1e-3);
    StandardNormal* D = new StandardNormal(4,12,1e-3);
    EXPECT_TRUE(*C!=*D);
    StandardNormal* E = new StandardNormal(5,20,1e-6);
    EXPECT_TRUE(A==*E);
    E->SetVal1(0.7);
    EXPECT_TRUE(A!=*E);


}

TEST(Assignment , Copyconstructor ) {
    StandardNormal A;
    StandardNormal B =A;

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Copyconstructor2 ) {
    StandardNormal A(10,30,1e-3);
    StandardNormal B(A);

    EXPECT_TRUE(A==B);

}

TEST(Assignment , Operator ) {
    StandardNormal A;
    StandardNormal B(5,10,0.4);

    EXPECT_TRUE(A!=B);

    B=A;

    EXPECT_TRUE(A==B);
}

TEST(Assignment , Swap ) {
    StandardNormal A;
    StandardNormal B;
    EXPECT_TRUE(A==B);
    StandardNormal* C = new StandardNormal(5,20,0.3);
    StandardNormal* D = new StandardNormal(5,20,0.3);
    EXPECT_TRUE(*C==*D);

    A.swap(*C);
    EXPECT_TRUE(A==*D);
    EXPECT_TRUE(B==*C);


}


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

        ASSERT_NEAR(StandardNormalRnd_new, StandardNormalRnd,1e-13*abs_StandardNormalRnd);

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
