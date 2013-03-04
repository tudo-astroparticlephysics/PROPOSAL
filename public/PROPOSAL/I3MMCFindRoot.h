// I3MMCFindRoot header file

#ifndef I3MMC_I3MMCFINDROOT_H_INCLUDED
#define I3MMC_I3MMCFINDROOT_H_INCLUDED

class I3MMCFindRoot:public I3MMCMathModel{
public:
	I3MMCFindRoot(void);
	I3MMCFindRoot(int maxSteps, double precision);
	~I3MMCFindRoot(void);

	double FindRoot(double function2use(double x), double dfunction2use(double x),double min, double max, double startX, double rightSide);

private:
	double _precision;
    int _maxSteps;
    };

#endif // I3MMC_I3MMCFINDROOT_H_INCLUDED
