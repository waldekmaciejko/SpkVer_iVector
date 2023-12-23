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

#include "mfccextractor.h"

namespace ava{

MfccExtractor::MfccExtractor()
{

}

MfccExtractor::MfccExtractor(std::string _SOURCEFORMAT,
        //				int _SOURCERATE, // samodzielny odczyt f probkowania pliku
                        int _WINDOWSIZE,
                        int _TARGETRATE,
                        std::string _SOURCEKIND,
                        std::string _TARGETFORMAT,
                        std::string _TARGETKIND,
                        char _SAVECOMPRESSED,
                        char _SAVEWITHCRC,
                        char _USEHAMMING,
                        int _NSTREAMS,
                        float _PREEMCOEF,
                        int _LOFREQ,
                        int _HIFREQ,
                        int _NUMCHANS,
                        int _CEPLIFTER,
                        int _NUMCEPS,
                        char _FORCEOUT,
                        char _NATURALWRITEORDER,
                        char _NATURALREADORDER){

                        SOURCEFORMAT=_SOURCEFORMAT;
      //					SOURCERATE=_SOURCERATE; // samodzielny odczyt f probkowania pliku
                        WINDOWSIZE=_WINDOWSIZE;
                        TARGETRATE=_TARGETRATE;
                        SOURCEKIND=_SOURCEKIND;
                        TARGETFORMAT=_TARGETFORMAT;
                        TARGETKIND=_TARGETKIND;
                        SAVECOMPRESSED=_SAVECOMPRESSED;
                        SAVEWITHCRC=_SAVEWITHCRC;
                        USEHAMMING=_USEHAMMING;
                        NSTREAMS=_NSTREAMS;
                        PREEMCOEF=_PREEMCOEF;
                        LOFREQ=_LOFREQ;
                        HIFREQ=_HIFREQ;
                        NUMCHANS=_NUMCHANS;
                        CEPLIFTER=_CEPLIFTER;
                        NUMCEPS=_NUMCEPS;
                        FORCEOUT=_FORCEOUT;
                        NATURALWRITEORDER=_NATURALWRITEORDER;
                        NATURALREADORDER=_NATURALREADORDER;

}

//-------------------generuje plik konfig gla obiektu mfcc
std::string MfccExtractor::generujConfigMFCC(std::string SciezkaConfig, std::string nazwa){

    std::string __SOURCEFORMAT=this->SOURCEFORMAT;
    //int __SOURCERATE = this->SOURCERATE;             //samodzielny odczyt f probkowania pliku
    int __WINDOWSIZE = this->WINDOWSIZE;
    int __TARGETRATE = this->TARGETRATE;
    std::string __SOURCEKIND = this->SOURCEKIND;
    std::string __TARGETFORMAT = this->TARGETFORMAT;
    std::string __TARGETKIND = this->TARGETKIND;
    char __SAVECOMPRESSED = this->SAVECOMPRESSED;
    char __SAVEWITHCRC = this->SAVEWITHCRC;
    char __USEHAMMING	= this->USEHAMMING;
    int __NSTREAMS = this->NSTREAMS;
    float __PREEMCOEF	= this->PREEMCOEF;           // 1-az^-1
    int __LOFREQ = this->LOFREQ;
    int __HIFREQ = this->HIFREQ;
    int __NUMCHANS = this->NUMCHANS;
    int __CEPLIFTER = this->CEPLIFTER;
    int __NUMCEPS = this->NUMCEPS;
    char __FORCEOUT = this->FORCEOUT;
    char __NATURALWRITEORDER = this->NATURALWRITEORDER;
    char __NATURALREADORDER = this->NATURALREADORDER;

    std::string _SciezkaConfig=SciezkaConfig+"\\"+nazwa;
    //------------------------- tworz plik tekstowy
    std::ofstream DoConfig;
    DoConfig.open(_SciezkaConfig.c_str());

    DoConfig<<"SOURCEFORMAT = "<<__SOURCEFORMAT<<std::endl
        //<<"SOURCERATE = "<<__SOURCERATE<<endl	//samodzielny odczyt f probkowania pliku
        <<"WINDOWSIZE = "<<__WINDOWSIZE<<std::endl
        <<"TARGETRATE = "<<__TARGETRATE<<std::endl
        <<"SOURCEKIND = "<<__SOURCEKIND<<std::endl
        <<"TARGETFORMAT = "<<__TARGETFORMAT<<std::endl
        <<"TARGETKIND = "<<__TARGETKIND<<std::endl
        <<"SAVECOMPRESSED = "<<__SAVECOMPRESSED<<std::endl
        <<"SAVEWITHCRC = "<<__SAVEWITHCRC<<std::endl
        <<"USEHAMMING = "<<__USEHAMMING<<std::endl
        <<"NSTREAMS = "<<__NSTREAMS<<std::endl
        <<"PREEMCOEF = "<<__PREEMCOEF<<std::endl           // 1-az^-1
        <<"LOFREQ = "<<__LOFREQ<<std::endl
        <<"HIFREQ = "<<__HIFREQ<<std::endl
        <<"NUMCHANS = "<<__NUMCHANS<<std::endl
        <<"CEPLIFTER = "<<__CEPLIFTER<<std::endl
        <<"NUMCEPS = "<<__NUMCEPS<<std::endl
        <<"FORCEOUT = "<<__FORCEOUT<<std::endl
        <<"NATURALWRITEORDER = "<<__NATURALWRITEORDER<<std::endl
        <<"NATURALREADORDER = "<<__NATURALREADORDER<<std::endl;

        DoConfig.close();
        return _SciezkaConfig;
}

/*!
 * \brief MfccExtractor::wywolajHCopy - call the Hcopy binary
 * \param path_f_hcopy_bin
 * \param pathToMfccConfig
 * \param path_d_workSpace
 * \param path_f_flt_fea_list_str
 */
void MfccExtractor::wywolajHCopy(std::string path_f_hcopy_bin,
                                  std::string pathToMfccConfig,
                                  std::string path_d_workSpace,
                                  std::string path_f_flt_fea_list_str)
{
    std::string mode_a="-C";
    std::string mode_b="-S";
    std::string call=path_f_hcopy_bin+" "+mode_a+" "+pathToMfccConfig +" "+mode_b+" "+ path_d_workSpace+"\\"+path_f_flt_fea_list_str;
    //std::cout<<call<<std::endl;
    system(call.c_str());
}
//-------------------------------------------------------------------------

}
