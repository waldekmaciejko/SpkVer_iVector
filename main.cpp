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

#include <cstdlib>
#include <iostream>
#include <string>
#include <string.h>

#include "config-fileparser.h"

#include <algorithm>
#include <cstdlib>
#include <string>

#include "mfccefeatures.h"
#include "mfccextractor.h"
#include "funhelpers.h"
#include "filelists.h"
#include "gpldamodel.h"
#include "trainTVPMiVector.h"
#include "Log4AVA.h"
#include "enrollFRRFAR.h"
#include "timeheleprs.h"
#include "createworkspace.h"
#include "unittests.h"

#include <typeinfo>
#include <sys/sysinfo.h>


int main(int argc, char* argv[]){

    //avatest::test1(); -OK
    //avatest::test2(); -OK
    //avatest::test3(); -OK
    //avatest::test4(); -OK
    //avatest::test5(); -OK
    //avatest::test6(); //-OK
    //avatest::test7(); //-OK
    //avatest::test8();
    //avatest::test9(); //-OK
    //avatest::test10();//-OK
    //avatest::test11();//-OK
    //avatest::test12();//-OK
    //avatest::test13();//-OK
    //avatest::test14();

    //FileParser fpconfig;
    //FileParser fpworkspace;
    //FileParser fpdata;

    FileParserCfg fpconfig;
    FileParserWrkSpace fpworkspace;
    FileParserData fpdata;

    std::string pthcfg{};
    std::string pthwrkspc{};
    std::string pthdata{};

    bool verbose = false;

    std::stringstream ss{};

    // Parse config file
    if(!(strcmp(argv[1], "-h"))){
        std::cout << " ===========================================\n"
                  << " This is ivector speaker comparing system. \n"
                  << "\n"
                  << " Compiled on : "<< ava::readFile("/proc/version")<<"\n"
                  << "\n"
                  << " Example: ./SpkrVer_iVector -c config.cfg -w workspace.cfg -d data.cfg -v\n"
                  << "        ./SpkrVer_iVector -h\n"
                  << "\n"
                  << " For details read README.md. \n"
                  << " -c   path to textfile, configuration file \n"
                  << " -w   path to textfile, list of paths of workspace \n"
                  << " -d   path to textfile, paths to data: train, enroll, model, etc \n"
                  << " -v   use verbose mode, otherwise use silent mode (all logs to raport)"
                  << " ===========================================\n";
        return -1;

    }else if(argc > 1){
        for(int ndx{}; ndx < argc; ++ndx){
            if (!(strcmp(argv[ndx], "-c"))){
                fpconfig.parseConfigFile(argv, ndx, "-c");
                pthcfg = argv[ndx+1];
            }
            if(!(strcmp(argv[ndx], "-w"))){
                fpworkspace.parseConfigFile(argv, ndx, "-w");
                pthwrkspc = argv[ndx+1];
            }
            if(!(strcmp(argv[ndx], "-d"))){
                fpdata.parseConfigFile(argv, ndx, "-d");
                pthdata = argv[ndx+1];
            }
            if(!(strcmp(argv[ndx], "-v"))){
                verbose = true;
            }
        }
    }

    //std::cout << fpconfig.numEigenVoices << '\n' ;
    //std::cout << "test120\n";


    //std::cout << "========" <<fpdata.pWork_path << '\n';

    // for(uint i=0; i<train.size(); ++i){
    //     std::cout <<train[i]<<'\n';
    // }
    // for(uint i=0; i<test.size(); ++i){
    //     std::cout <<test[i]<<'\n';
    // }

    // std::cout << fpdata.pToModel<< '\n';
    // std::cout << fpworkspace.pTrainMfcc_path<< '\n';
    // std::cout << fpworkspace.pTrainList_file<< '\n';
    // std::cout << fpdata.pToModel << "\n";
    // std::cout << fpworkspace.pEnrollList_file << "\n";
    // std::cout << fpworkspace.pTrainList_file << "\n";

    CreateWorkspace *cw = new CreateWorkspace(fpdata.pWork_path);
    cw->validateDirs();
    cw->validateSrtct();

    fpdata.pWork_path =cw->getWorkspace();

    std::cout << fpdata.pWork_path  << std::endl;

    ava::Log4AVA logger(fpdata.pWork_path,
                        fpworkspace.pToSaveWSpace+"/",
                        ava::currentTime(),
                        fpconfig.numOfComponents,
                        fpconfig.numTdim,
                        verbose);

    std::string s = logger.returnSavePath();

    auto test = fpdata.getList_test();
    auto train = fpdata.getList_train();
    std::string pToModelFull = fpdata.pWork_path +fpdata.pToModel;


    logger.save("\nStart time:              ", ava::currentTime());

    logger.save("\n*************************", "RAPORT START\n");
    logger.save("Used data.cfg:           ", pthcfg);
    logger.save("Used workspace.cfg:      ", pthwrkspc);
    logger.save("Used data:               ", pthdata);
    logger.save("Train data: ", train);
    logger.save("Test data: ", test);
    logger.save("\n*************************", "CONFIG PARAMETERS\n");
    logger.save("numOfComponents:         ", fpconfig.numOfComponents);
    logger.save("numIterationsEM:         ", fpconfig.numIterationsEM);
    logger.save("numIterationsKmean_init: ", fpconfig.numIterationsKmean_init);
    logger.save("numTdim:                 ", fpconfig.numTdim);
    logger.save("numIterationsTV:         ", fpconfig.numIterationsTV);
    logger.save("performLDA:              ", fpconfig.performLDA);
    logger.save("performWCCN:             ", fpconfig.performWCCN);
    logger.save("numEigenVectors:         ", fpconfig.numEigenVectors);
    logger.save("zmiennWhitening:         ", fpconfig.zmiennWhitening);
    logger.save("numEigenVoices:          ", fpconfig.numEigenVoices);
    logger.save("numIterGPLDA:            ", fpconfig.numIter);
    logger.save("minimumDivergence:       ", fpconfig.minimumDivergence);
    logger.save("SAVETRAINMODEL:          ", fpconfig.SAVETRAINMODEL);
    logger.save("train UBM, TV, Proj Matr:", fpconfig.train);
    logger.save("test:                    ", fpconfig.enroll);
    logger.save("\n*************************", "EXECUTION\n");

    //*************** END TRAIN PARAMETERS ***************

    if(fpconfig.train){

        auto start = std::chrono::high_resolution_clock::now();

        logger.save("Traininng Total Variability, Projection Matrcis, iVectors ", "STARTS");

        ava::trainTVPMiVector(train,
                              fpdata.pWork_path,
                              fpworkspace.pTrainList_file,
                              fpworkspace.pTrainMfcc_path,
                              fpworkspace.pCfg_path,
                              fpworkspace.pToSaveWSpace,
                              fpworkspace.pTmp_file,
                              fpconfig.pHCopy_bin,
                              fpconfig.numOfComponents,
                              fpconfig.numIterationsEM,
                              fpconfig.numIterationsKmean_init,
                              fpconfig.numTdim,
                              fpconfig.numIterationsTV,
                              fpconfig.performLDA,
                              fpconfig.performWCCN,
                              fpconfig.numEigenVectors,
                              fpconfig.zmiennWhitening,
                              fpconfig.numEigenVoices,
                              fpconfig.numIter,
                              fpconfig.minimumDivergence,
                              fpconfig.SAVETRAINMODEL,
                              logger);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss << " main "<< __LINE__ << " (duration sec): "<<duration.count()<<std::endl;
        logger.save("Execution time in point ", ss);
        ss.str("");
        logger.save("\nTraininng Total Variability, Projection Matrics, iVectors ", "FINISHED\n");
    }
    if(fpconfig.enroll){

        //*************** START ENROLL ***************
        auto start = std::chrono::high_resolution_clock::now();

        logger.save("Enroll ", "STARTS");

        ava::FileLists fileListEnroll(fpdata.pWork_path,
                                      fpworkspace.pEnrollList_file,
                                      test,
                                      fpworkspace.pEnrollMfcc_path);

        std::multimap<std::string, std::string> multimapEnrollMFC{};
        std::multimap<std::string, std::string> multimapTestMFC{};
        std::set<std::string> setSpkEnrollTest{};

        fileListEnroll.listWavToMfccEnrollTest(multimapEnrollMFC,
                                               multimapTestMFC,
                                               setSpkEnrollTest,
                                               fpconfig.test_size);

        //----------------------------------- utworz lancuch do przekazania do
        //----------------------------------- konsoli wywylujacy HCopy dla ENROLL
        ava::MfccExtractor *ME;                    // parametry MFCC zdefiniowano w pliku mfccExtractor
        ME=new ava::MfccExtractor();
        std::string pCfgHcopy_file=ME->generujConfigMFCC(fpdata.pWork_path+
                                                               fpworkspace.pCfg_path, "conf.mfcc");    //-- wygeneruj cfg dla HCopy
        ME->wywolajHCopy(fpconfig.pHCopy_bin,
                         pCfgHcopy_file,
                         fpdata.pWork_path,
                         fpworkspace.pEnrollList_file);


        //     //*************** END ENROLL ***************

        //     //---------------- READ MODELS FOR FAR AND FRR

        arma::mat normMean{};
        arma::mat normStd{};
        arma::gmm_diag ubmModel;
        arma::mat T={};
        arma::mat TS={};
        arma::mat TSi={};
        arma::mat projectionMatrix={};
        std::multimap<std::string, arma::mat> enrolledSpeakersByIdx{};

        normMean.load(pToModelFull + "/normMean.csv");
        normStd.load(pToModelFull+"/normStd.csv");
        ubmModel.load(pToModelFull+"/ubmModel.gmm");
        T.load(pToModelFull+"/T.csv");
        TS.load(pToModelFull+"/TS.csv");
        TSi.load(pToModelFull+"/TSi.csv");
        projectionMatrix.load(pToModelFull+"/projectionMatrix.csv");

        //     //---------------- END READ MODELS

        ava::enrollII(setSpkEnrollTest,
                      multimapEnrollMFC,
                      fpdata.pWork_path,
                      fpworkspace.pToLogs,
                      normMean,
                      normStd,
                      ubmModel,
                      projectionMatrix,
                      fpconfig.numOfComponents,
                      fpconfig.numTdim,
                      T,
                      TS,
                      TSi,
                      enrolledSpeakersByIdx,
                      logger);

        ava::FRR(multimapTestMFC,
                 normMean,
                 normStd,
                 fpconfig.numOfComponents,
                 ubmModel,
                 fpconfig.numTdim,
                 projectionMatrix,
                 enrolledSpeakersByIdx,
                 fpdata.pWork_path,
                 fpworkspace.pToLogs,
                 fpconfig.scoringMethod,
                 T,
                 TS,
                 TSi,
                 logger);

        ava::FAR(multimapTestMFC,
                 normMean,
                 normStd,
                 fpconfig.numOfComponents,
                 ubmModel,
                 fpconfig.numTdim,
                 projectionMatrix,
                 enrolledSpeakersByIdx,
                 fpdata.pWork_path,
                 fpworkspace.pToLogs,
                 fpconfig.scoringMethod,
                 T,
                 TS,
                 TSi,
                 logger);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss << " (duration sec): "<<duration.count()<<std::endl;
        logger.save("Enroll FINISHED", ss);
    }
    if(fpconfig.enroll == false && fpconfig.train == false){
        std::cout << "We have nothing to do.";
        logger.save("We have nothing to do. ", "Please switch on train or enfroll\n");
    }
    logger.save("Finish time:              ", ava::currentTime());

    return 0;
}


