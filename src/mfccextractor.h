/*
This file is part of SpkVer_iVector - speaker verification software base
on Total Variability and Projection matrixes

SpkVer_iVector is free software: you can redistribute it and/or modify
it under the terms of the MIT License.

SpkVer_iVector is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
MIT License for more details.

The SpkVer_iVector project shows the
limits of voice authentication in a forensic context.
The "Person Authentification by Voice: A Need of Caution" paper
proposes a good overview of this point (cf. "Person
Authentification by Voice: A Need of Caution", Bonastre J.F.,
Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-
chagnolleau I., Eurospeech 2003, Genova].
The conclusion of the paper of the paper is proposed bellow:
[Currently, it is not possible to completely determine whether the
similarity between two recordings is due to the speaker or to other
factors, especially when: (a) the speaker does not cooperate, (b) there
is no control over recording equipment, (c) recording conditions are not
known, (d) one does not know whether the voice was disguised and, to a
lesser extent, (e) the linguistic content of the message is not
controlled. Caution and judgment must be exercised when applying speaker
recognition techniques, whether human or automatic, to account for these
uncontrolled factors. Under more constrained or calibrated situations,
or as an aid for investigative purposes, judicious application of these
techniques may be suitable, provided they are not considered as infallible.
At the present time, there is no scientific process that enables one to
uniquely characterize a persones voice or to identify with absolute
certainty an individual from his or her voice.]

Copyright (C) 2004-2023
Waldek Maciejko
*/

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

namespace ava{

class MfccExtractor{
    private:
        // ustawiono parametry domyslne
    std::string SOURCEFORMAT = "WAV";
    int WINDOWSIZE=200000; //      = 200000.0 // dlugosc okna
    int TARGETRATE=100000; //      = 100000.0 // przesuniecie okna    
    std::string SOURCEKIND="WAVEFORM"; //	= WAVEFORM
    std::string TARGETFORMAT="HTK";	//= HTK
    //std::string TARGETKIND="MFCC_E"; //    	= MFCC_E // ustalenie parametrow ktore zostan¹ obliczone
    std::string TARGETKIND="MFCC_D_A";
    char SAVECOMPRESSED='F'; //  = F
    char SAVEWITHCRC='F'; //     = F
    char USEHAMMING='T'; //     	= T
    int NSTREAMS=1; //	= 1                 // liczba strumieni
    float PREEMCOEF=0.97; //      	= 0.0mfcc, 0.97plp # 1-az^-1   // preemfaza
    int LOFREQ=300; //		= 300 // dolan czestotliwosc odciecia
    int HIFREQ=4000; //		= 8000 // gorna czestotliwosc odciecia
    int NUMCHANS=24; //        = 24 dla mfcc, =20 dla plp // liczba kanalow w banku filtrow
    int CEPLIFTER=0; //	= 0 dla mfcc, 22 dla plp // wspolczynnik liftracji
    int NUMCEPS=20; //         = 20 // liczba wspolczynnikow
    char FORCEOUT='T'; //        = T
    char NATURALWRITEORDER='T'; // = T // zapis jako little endian - ta wartosc jest zmieniana T->F aby odycztac w HList
    char NATURALREADORDER='T'; // = T // odczyt jako little endian
//---------------------------------- interfejs
    public:
    MfccExtractor(void);
    //---------------------------------- konstruktor dla mfcc

    MfccExtractor(std::string,	// ="WAV",
                int,	// = 200000.0,
                int,	// = 100000.0,
                std::string,		// = "WAVEFORM",
                std::string,		// = "HTK",
                std::string,		// = "MFCC_E",
                char,		// = 'F',
                char,		// = 'F',
                char,		// = 'T',
                int,		 //= 1,
                float,		 //= 0.0,
                int,		// = 300,
                int,		// = 8000,
                int,		// = 24,
                int,		// = 0,
                int,		// = 16,
                char,		// = 'T',
                char,		// = 'T',
                char);	 	// = 'T'

//-------------------------------------
    std::string generujConfigMFCC(std::string SciezkaConfig, std::string nazwa);  // funkcja generuje plik config
    void wywolajHCopy(std::string path_f_hcopy_bin,
                      std::string pathToMfccConfig,
                      std::string path_d_workSpace,
                      std::string path_f_flt_fea_list_str);
/*
    int numcepsReturn();   				// zwraca liczbe NUMCEPS + 1

    string zwrocSOURCEFORMAT();
    //int ZwrocSOURCERATE(); //      = 625
    int zwrocWINDOWSIZE(); //      = 200000.0
    int zwrocTARGETRATE(); //      = 100000.0
    string zwrocSOURCEKIND(); //	= WAVEFORM
    string zwrocTARGETFORMAT();	//= HTK
    string zwrocTARGETKIND(); //    	= MFCC_E
    char zwrocSAVECOMPRESSED(); //  = F
    char zwrocSAVEWITHCRC(); //     = F
    char zwrocUSEHAMMING(); //     	= T
    int zwrocNSTREAMS(); //	= 1
    float zwrocPREEMCOEF(); //      	= 0.0           # 1-az^-1
    int zwrocLOFREQ(); //		= 300
    int zwrocHIFREQ(); //		= 8000
    int zwrocNUMCHANS(); //        = 24
    int zwrocCEPLIFTER(); //	= 0
    int zwrocNUMCEPS(); //         = 16
    char zwrocFORCEOUT(); //        = T
    int zwrocLPCORDER();
    char zwrocNATURALWRITEORDER(); // = T
    char zwrocNATURALREADORDER(); // = T
    */

};
}

