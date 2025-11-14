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

#include "filelists.h"

namespace ava{
/*!
 * \brief FileLists::FileLists
 * \param _workSpace    -   path to WORKSPACE
 * \param _flt_fea_list_str -   path to HCopy list pcm mfc
 * \param _v_to_korpus  - table with paths to PCM
 * \param _flt_fea_mfcc_str - path to MFCC
 */
FileLists::FileLists(std::string _workSpace,
                     std::string _flt_fea_list_str,
                     std::vector<std::string> _v_to_korpus,
                     std::string _flt_fea_mfcc_str){

    this->workSpace=_workSpace;
    this->flt_fea_list_str=_flt_fea_list_str;
    this->v_to_korpus=_v_to_korpus;
    this->flt_fea_mfcc_str=_flt_fea_mfcc_str;
}

//----------------------------- funkcja generujaca listy dla biblioteki HCopy

//void FileLists::listaMfcc(string daneIn,
//                                string pathOut,
//                                string pathFeatureOut,
//                                string daneOutRozszerzenie,
//                                string nazwaListy,
//                                string featureId,
//                                int i){
//using namespace std;

// *
// * naneIn	-	sciezka do plikow wejsciowych (wav lub htk)
// * pathOut	-	sciezka do folderu gdzie zostan¹ umieszczone listy
// * string PathFeatureOut	-
// * daneOutRozszerzenie	-	rozszerzenie plikow wyjsciowych
// * count		-		mowi ile plikow znajduje sie w analizowanym folderze
// * nazwaListy	-	do ktorej zostana wpisane scieki do plikow
// * featureId	-		na podstawie Id oraz "i" zostana nadane nazwy plikow cech
// */

//    ofstream doPliku;
//    doPliku.open(pathOut.c_str(),ios_base::app);
//    string doWav=winPathToPosixPath(daneIn.c_str());
//    string doFeature=winPathToPosixPath(pathFeatureOut.c_str());
//    doPliku<<doWav<<" "<<doFeature<<"/"<<featureId<<i<<"."<<daneOutRozszerzenie<<endl;
//    doPliku.close();
//}

/*!
 * \brief Generate lists of files PCM ->MFCC
 *  Create files list for HTK HCopy binry in the form:
 *      - path to wav          path to mfcc
 *      - path to wav          path to mfcc
 * \brief Organization of the files of corpora have to be:
 * \brief \M01\M1_*_1.wav
 * \brief \M01\M1_*_2.wav
 * \brief \M01\M1_*_3.wav ...
 * \brief -------------------
 * \brief \F01\F01_*_1.wav
 * \brief \F01\F01_*_2.wav ...
 * \brief \F02\F02_*_3.wav ... etc
 * \param[out] multimapMowcyMFC - multimap in the form e.g. "F01" "D:\..\M01.mfcc"
 * \param[out] numFilesUBM      - number of PCM files in Corpora
 * \param[out] numSpkUBM        - number of speakers in Corpora, the differentiation
 *                              is make using M01, M02.. labels on files structure
 * \param[out] setSpk           - set of unique labels e.g. M01, M02...
 *
 * ATTENTION!!! variables numFilesUBM and numSpkUBM are critical
 * todo: develop tool to check this values-->PCM files->MFC  files->list
 *
 *
 * \return void
 *
 */
void FileLists::listWavToMfccUBM(std::multimap<std::string, std::string>& multimapMowcyMFC,
                                  unsigned int& numFilesUBM,
                                  unsigned int& numSpkUBM,
                                  std::set<std::string>& setSpk){

    std::ofstream flt_fea_mfcc_stream(winPathToPosixPath((this->workSpace+this->flt_fea_list_str).c_str()));
    numFilesUBM = 0;
    for(std::string s : this->v_to_korpus){       

       // std::ofstream flt_fea_mfcc_stream(winPathToPosixPath((this->workSpace+this->flt_fea_list_str).c_str()));
        std::string tmp=(this->workSpace+this->flt_fea_list_str).c_str();

        //QDirIterator it(QString::fromStdString(this->p_to_korpus), QStringList()<<"*.wav", QDir::Files, QDirIterator::Subdirectories);
        QDirIterator it(QString::fromStdString(s), QStringList()<<"*.wav", QDir::Files, QDirIterator::Subdirectories);
        std::string str_tmp={};
        std::string wavExt=".wav";
        QString qpathToWav={};
        std::string  spathToWav={};
        std::string flt_fea_fea=winPathToPosixPath((this->workSpace+this->flt_fea_mfcc_str).c_str());

        std::string s1={};
        std::string s2={};

        while(it.hasNext())
        {
            qpathToWav=it.next();
            QFileInfo fi(qpathToWav);
            spathToWav=qpathToWav.toStdString();
            str_tmp = fi.fileName().toStdString();
            size_t pos=str_tmp.find(wavExt);
            str_tmp.replace(pos, wavExt.length(), ".mfc");
            s2 = flt_fea_fea +"/" + str_tmp;
            //flt_fea_mfcc_stream<<spathToWav<<" "<<flt_fea_fea<<"/"<<str_tmp<<std::endl;
            flt_fea_mfcc_stream<<spathToWav<<" "<<s2<<std::endl;
            numFilesUBM=numFilesUBM+1;
            s1 = ava::takeSpkLabel(spathToWav);
            if(s1!="0"){
                multimapMowcyMFC.insert(std::make_pair(s1, s2));
                if(setSpk.count(s1)==0){
                    setSpk.insert(s1);
                    numSpkUBM=numSpkUBM+1;
                }
            }
        }        
    }
    flt_fea_mfcc_stream.close();
    //return multimapMowcyMFC;
}

void FileLists::listWavToMfccEnrollTest(std::multimap<std::string, std::string>& multimapEnrollMFC,
                                std::multimap<std::string, std::string>& multimapTestMFC,
                                std::set<std::string>& setSpkEnrollTest,
                                float test_size){

    std::multimap<std::string, std::string> multimapEnrollPCM{};
    std::multimap<std::string, std::string> multimapTestPCM{};

    //1.
    //create map anf fullfill: label of spk  -   path to pcm
    std::multimap<std::string, std::string> multimapMowcyPCM{};

    for(std::string s : this->v_to_korpus){

        std::string tmp=(this->workSpace+"/"+this->flt_fea_list_str).c_str();

        QDirIterator it(QString::fromStdString(s),
                        QStringList()<<"*.wav",
                        QDir::Files,
                        QDirIterator::Subdirectories);

        std::string wavExt=".wav";
        QString qpathToWav={};
        std::string  spathToWav={};
        std::string s1={};

        while(it.hasNext())
        {
            qpathToWav=it.next();
            QFileInfo fi(qpathToWav);
            spathToWav=qpathToWav.toStdString();

            // take label of spk from path
            s1 = ava::takeSpkLabel(spathToWav);
            if(s1!="0"){
                multimapMowcyPCM.insert(std::make_pair(s1, spathToWav));
              if(setSpkEnrollTest.count(s1)==0){
                      setSpkEnrollTest.insert(s1);
              }
            }
        }        
    }

    uint numSpkEnrollTestFiles{};
    std::vector<uint> ranTestNumVector{};
    std::vector<uint> ranEnrollNumVec{};

    uint counter = 0;

    //2. iteration PCMs of each spk, take randoom vectors and insetrt into
    //  suitable multimaps
    using mmapIterator=std::multimap<std::string, std::string>::iterator;
    for(auto s : setSpkEnrollTest){

        //take all file paths with spk label s
        std::pair<mmapIterator, mmapIterator> res = multimapMowcyPCM.equal_range(s);

        // take amount of files in each spk (label) using equal_range built-in function
        numSpkEnrollTestFiles=std::distance(res.first, res.second);

        // random index between test and enroll
        ava::randomTestSplit(ranTestNumVector,
                             ranEnrollNumVec,
                             numSpkEnrollTestFiles,
                             test_size);

        // find test index and insert into multimapTestPCM else insert into multimapEnrollPCM
        for(mmapIterator it = res.first; it != res.second; it++){
            if(std::find(ranTestNumVector.begin(),
                         ranTestNumVector.end(),
                         counter)!=ranTestNumVector.end()){
                multimapTestPCM.insert(std::make_pair(it->first, it->second));
            }else multimapEnrollPCM.insert(std::make_pair(it->first, it->second));
            counter=counter+1;
        }

        counter=0;
    ranTestNumVector.clear();
    ranEnrollNumVec.clear();
    }

    //3. create lists for Hcopy and multimaps mfc
    std::string tmp{};
    std::string str_tmp{};
    std::string wavExt=".wav";
    QString qpathToWav={};
    std::string s2={};
    std::string flt_fea_fea=winPathToPosixPath((this->workSpace+this->flt_fea_mfcc_str).c_str());

    std::ofstream flt_fea_mfcc_stream(winPathToPosixPath((this->workSpace+this->flt_fea_list_str).c_str()));

    //write test date to list for Hcopy and multimap mfc
    for(mmapIterator it = multimapTestPCM.begin(); it!=multimapTestPCM.end(); it++){
        qpathToWav=QString::fromStdString(it->second);
        QFileInfo fi(qpathToWav);
        str_tmp = fi.fileName().toStdString();
        size_t pos=str_tmp.find(wavExt);
        str_tmp.replace(pos, wavExt.length(), ".mfc");
        s2 = flt_fea_fea +"//" + str_tmp;
        flt_fea_mfcc_stream<<it->second<<" "<<s2<<std::endl;
        multimapTestMFC.insert(std::make_pair(it->first, s2));

    }

    //write enroll date to list for Hcopy and multimap mfc
    for(mmapIterator it = multimapEnrollPCM.begin(); it!=multimapEnrollPCM.end(); it++){
        qpathToWav=QString::fromStdString(it->second);
        QFileInfo fi(qpathToWav);
        str_tmp = fi.fileName().toStdString();
        size_t pos=str_tmp.find(wavExt);
        str_tmp.replace(pos, wavExt.length(), ".mfc");
        s2 = flt_fea_fea +"//" + str_tmp;
        flt_fea_mfcc_stream<<it->second<<" "<<s2<<std::endl;
        multimapEnrollMFC.insert(std::make_pair(it->first, s2));
    }
    flt_fea_mfcc_stream.close();
}

} //ava



