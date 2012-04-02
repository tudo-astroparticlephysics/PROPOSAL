/*
 * Output.h
 *
 *  Created on: 29.09.2010
 *      Author: schmitz
 */

#include "CrossSections.h"



#ifndef OUTPUT_H_
#define OUTPUT_H_



#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "MMCexception.h"
#include "deque"



/**
 * Class containing the Output Methods. This Class will provide just the fundamentel Methods so far. You will need to implement more methods later. Look in the Java source file to do so.
 */

class PROPOSALParticle;

 class Output{

    

        private:
		// Members
		

			PROPOSALParticle *p; // PROPOSALParticle which needs to be printed somehow

                         // File to write interpolationtables

                        static int OUTNUM;

                        static std::vector<char> c; // needed to format the double

                        const static char z0='0', zm='-', zd='.', ze='e'; // needed to format the double. Set to 0, -, ., e in Constructor.

                        static long POWOUT; // need to format the double;


			   //?!?!


	public:
		// Members
                        static std::ofstream output_P;
                        std::stringstream history; // file to be written in - Maybe fstream is not good here..
                        static std::fstream ReadWrite;
                        int igen, gens; // ?!

                        static std::string version;

                        const static int HISTSIZE=1000;
                        int HIST; // output "styler", preset is -1
                        static bool AMASIM; // needed in Decay

                        std::streambuf *stdout_backup;

			std::fstream err; // Error file - Maybe fstream is not good here..
                        static bool outf;
                        static bool inf;
                        static bool I3flag;
                        static bool RecDec;

			bool DEBUG; // Debug flag

                        static bool raw; // set to false in Constructor

                        static bool texi; // set to false in Constructor

                        std::vector<PROPOSALParticle*> I3hist;

                        std::stringstream out;
		// Constructors
	

			Output(PROPOSALParticle *part);


			Output();
	
		// Functions

            //----------------------------------------------------------------------------------------------------//

              /**
               * enables redirection of the stderr output into a file
               */

                        bool setStderr(std::string name) throw (MMCexception);
            //----------------------------------------------------------------------------------------------------//


               /**
                * enables redirection of the stdout output into a string
                */

                        void setStdout() ;//throw (MMCexception);

            //----------------------------------------------------------------------------------------------------//

                /**
                 * gets the stdout output string, resets the stdout output
                 */

                        static void particleToFileAscii(PROPOSALParticle *p,std::string filename, bool isPrimary);

                        static void particleToFileRaw(PROPOSALParticle *p, std::string filename, bool isPrimary);

                        std::string getStdout();// throw (MMCexception);

                        static double read();

                        static void write(double g);

			// format the double

                        static std::string f(double d);

			// needed in Decay

			void output(int wint, std::string comp, double de, double ef);

                        // ?!

                        void init(std::string name);

                        // ?!

                        static void open(std::string name);

                        // ?!

                        static void close();

                        // ?!

                        static bool FileExist(std::string path);

                        static void Delete(std::string name){};

                        static std::deque<std::string> splitString(std::string args, std::string Delimiters);

                        // initialize for regular history output

                        void initDefault(std::string type);

                        void initDefault(int igen, int gens, std::string type, double time, double x, double y, double z, double theta, double phi);

                        //initialize for F2000 history output

                        void initF2000(int igen, int gens, std::string type, double time, double x, double y, double z, double theta, double phi);

		// Getter
	
		bool get_DEBUG(){ return DEBUG; }
		bool get_AMASIM(){return AMASIM;}
                static void set_I3flag(bool new_flag){I3flag = new_flag;}



		// Setter

		void set_DEBUG(bool newDebug){ DEBUG = newDebug; }

		
		// Maybe later to implement (exsistent in Java, not needed yet in c++)		
		/*	


			void open(std::string name); // opens the table // ** Not yet implemeted **
 
			void close(); // closes the table // ** Not yet implemeted **

			void Delete(std::string name);  // ** Not yet implemeted **
			
			 	
		










		


			// boolean setStderr(std::string name); // enables redirection of the stdout output into a string

			// void setStdout(); // gets the stdout output string, resets the stdout output

			// void init(String type); // initialize for no history output - no need to call this function

			// void initDefault(String type); // initialize for regular history output




			// boolean exists(String name); // See if the file "name" exists
			

			
			

			// Reader reader(String name) throws Exception // Opens reader Stream

			// public static String[] splitString(String args) // Splits the string into the array
		*/
			
			




};

#endif /* OUTPUT_H_ */
