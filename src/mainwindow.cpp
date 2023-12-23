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

#include <QApplication>
#include <QCoreApplication>
#include <QStringList>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QDebug>

#include <algorithm>
#include <cstdlib>

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "mfccefeatures.h"
#include "mfccextractor.h"
#include "funhelpers.h"
#include "filelists.h"
#include "gpldamodel.h"
#include "trainTVPMiVector.h"
#include "Log4AVA.h"
#include "enrollFRRFAR.h"

#include "unittests.h"

#include <typeinfo>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    MainWindow::ava300();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::ava300()
{    
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
    //---------------------------------
    //--------------------- SET PARAMETERS
    //---------------------------------
    std::string corpora_dir = "j:\\__Korpusy__";
    std::string workspace_dir = "j:\\__Badania__\\ava300_QtGUI_GPLDA\\ava300\\WORKSPACE";
    //*************** TRAIN KORPUS ***************

    std::vector<std::string> pTrain_path{};
    pTrain_path.push_back(corpora_dir + "\\PTDB-TUG\\SPEECH_DATA_silence_removed\\FEMALE\\MIC");
    pTrain_path.push_back(corpora_dir + "\\PTDB-TUG\\SPEECH_DATA_silence_removed\\MALE\\MIC");
    //pTrain_path.push_back("j:\\__Korpusy__\\korpus_roboczy\\train");

    //*************** ENROLL TEST KORPUS *********

    std::vector<std::string> pEnrollAndVerify_path{};
    pEnrollAndVerify_path.push_back(corpora_dir + "\\PTDB-TUG\\SPEECH_DATA_silence_removed\\FEMALE_ENROLLVERIFY\\MIC");
    pEnrollAndVerify_path.push_back(corpora_dir + "\\PTDB-TUG\\SPEECH_DATA_silence_removed\\MALE_ENROLLVERIFY\\MIC");
    //pEnrollAndVerify_path.push_back("j:\\__Korpusy__\\korpus_roboczy\\enrollverify");

    //*************** READ MODELS TO ENROLL AND TEST

    std::string pToModel = workspace_dir + "\\calc\\21-12-2023_13h29m30s_gmm128_TV100";

    //*************** WORKING DIRECTORYS

    std::string pWork_path      = "j:\\__Badania__\\ava300_QtGUI_GPLDA\\ava300\\WORKSPACE\\";
    std::string pTrainList_file = "\\list\\train_mfcc.lst";
    std::string pTrainMfcc_path = "\\feat\\train";
    std::string pCfg_path       = "\\cfg";
    std::string pToSaveWSpace   = "\\calc";
    std::string pTmp_file       = "\\tmp";
    std::string pToLogs         = "\\logs";

    std::string pHCopy_bin      = "j:\\__Badania__\\ava300_QtGUI_GPLDA\\ava300\\libs\\hcopy\\HCopy.exe";

    std::string pEnrollList_file   = "list\\enroll_mfcc.lst";
    std::string pEnrollMfcc_path   = "\\feat\\enroll\\";

    //--GMM

    const unsigned int numOfComponents = 128;     // 1024
    unsigned int numIterationsEM = 3;            // test: 3
    unsigned int numIterationsKmean_init = 10;

    //--Total Variability

    unsigned int numTdim = 100; // test:100 powinna wynosic około 1000
    unsigned int numIterationsTV = 10;   //optymalna wartosc wynosi 20

    //--Projection Matrix

    bool performLDA = true;
    bool performWCCN = true;
    uint numEigenVectors = 16;

    //--G-PLDA

    std::string zmiennWhitening = "ZCA"; //PCA, 000 - brak whitening
    uint numEigenVoices = 16;
    uint numIter = 10;
    bool minimumDivergence = true;
    bool SAVETRAINMODEL = true;  //-- SAVE TRAINED MODEL?
    bool train = true;           //-- TRAIN TOTAL VARIABILITY, PROJECTION MATRIX, iVECTOR?

    float test_size = 0.9; //-- BALANCE ENROLL/TEST
    uint scoringMethod = 1; // 1 - CSS, 0 - GPLDA

    //---------------------------------
    //--------------------- PARAMETERS END
    //---------------------------------

    //-------------- logging    

    ava::Log4AVA logger(pWork_path,
                        pToSaveWSpace+"\\",
                        ava::currentTime(),
                        numOfComponents,
                        numTdim);

    std::string s = logger.returnSavePath();

    logger.save("*************************", "RAPORT START");
    logger.save("start time:              ", ava::currentTime());
    logger.save("pTrain_path:             ", pTrain_path);
    logger.save("numOfComponents:         ", numOfComponents);
    logger.save("numIterationsEM:         ", numIterationsEM);
    logger.save("numIterationsKmean_init: ", numIterationsKmean_init);
    logger.save("numTdim:                 ", numTdim);
    logger.save("numIterationsTV:         ", numIterationsTV);
    logger.save("performLDA:              ", performLDA);
    logger.save("performWCCN:             ", performWCCN);
    logger.save("numEigenVectors:         ", numEigenVectors);
    logger.save("zmiennWhitening:         ", zmiennWhitening);
    logger.save("numEigenVoices:          ", numEigenVoices);
    logger.save("numIterGPLDA:            ", numIter);
    logger.save("minimumDivergence:       ", minimumDivergence);
    logger.save("SAVETRAINMODEL:          ", SAVETRAINMODEL);
    logger.save("train UBM, TV, Proj Matr:", train);

    //*************** END TRAIN PARAMETERS ***************

    if(train){
        ava::trainTVPMiVector(pTrain_path,
                                pWork_path,
                                pTrainList_file,
                                pTrainMfcc_path,
                                pCfg_path,
                                pToSaveWSpace,
                                pTmp_file,
                                pHCopy_bin,
                                numOfComponents,
                                numIterationsEM,
                                numIterationsKmean_init,
                                numTdim,
                                numIterationsTV,
                                performLDA,
                                performWCCN,
                                numEigenVectors,
                                zmiennWhitening,
                                numEigenVoices,
                                numIter,
                                minimumDivergence,
                                SAVETRAINMODEL,
                                logger);
    }else{

        //*************** START ENROLL ***************

        auto start = std::chrono::high_resolution_clock::now();

        ava::FileLists fileListEnroll(pWork_path,
                                      pEnrollList_file,
                                      pEnrollAndVerify_path,
                                      pEnrollMfcc_path);

        std::multimap<std::string, std::string> multimapEnrollMFC{};
        std::multimap<std::string, std::string> multimapTestMFC{};
        std::set<std::string> setSpkEnrollTest{};

        fileListEnroll.listWavToMfccEnrollTest(multimapEnrollMFC,
                                               multimapTestMFC,
                                               setSpkEnrollTest,
                                               test_size);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        std::cout<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration sec: "<<duration.count()<<std::endl;


        //----------------------------------- utworz lancuch do przekazania do
        //----------------------------------- konsoli wywylujacy HCopy dla ENROLL

        start = std::chrono::high_resolution_clock::now();

        ava::MfccExtractor *ME;                    // parametry MFCC zdefiniowano w pliku mfccExtractor
        ME=new ava::MfccExtractor();
        std::string pCfgHcopy_file=ME->generujConfigMFCC(pWork_path+
                                                            pCfg_path, "conf.mfcc");    //-- wygeneruj cfg dla HCopy
        ME->wywolajHCopy(pHCopy_bin,
                         pCfgHcopy_file,
                         pWork_path,
                         pEnrollList_file);

        stop = std::chrono::high_resolution_clock::now();
        duration=stop-start;
        std::cout<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration sec: "<<duration.count()<<std::endl;

        //*************** END ENROLL ***************

        //---------------- READ MODELS FOR FAR AND FRR

        arma::mat normMean{};
        arma::mat normStd{};
        arma::gmm_diag ubmModel;
        arma::mat T={};
        arma::mat TS={};
        arma::mat TSi={};
        arma::mat projectionMatrix={};
        std::multimap<std::string, arma::mat> enrolledSpeakersByIdx{};

        normMean.load(pToModel + "\\normMean.csv");
        normStd.load(pToModel+"\\normStd.csv");
        ubmModel.load(pToModel+"\\ubmModel.gmm");
        T.load(pToModel+"\\T.csv");
        TS.load(pToModel+"\\TS.csv");
        TSi.load(pToModel+"\\TSi.csv");
        projectionMatrix.load(pToModel+"\\projectionMatrix.csv");

        //---------------- END READ MODELS

        ava::enrollII(setSpkEnrollTest,
                    multimapEnrollMFC,
                    pWork_path,
                    pToLogs,
                    normMean,
                    normStd,
                    ubmModel,
                    projectionMatrix,
                    numOfComponents,
                    numTdim,
                    T,
                    TS,
                    TSi,
                    enrolledSpeakersByIdx,
                    logger);

        ava::FRR(multimapTestMFC,
                 normMean,
                 normStd,
                 numOfComponents,
                 ubmModel,
                 numTdim,
                 projectionMatrix,
                 enrolledSpeakersByIdx,
                 pWork_path,
                 pToLogs,
                 scoringMethod,
                 T,
                 TS,
                 TSi,
                 logger);

        ava::FAR(multimapTestMFC,
                 normMean,
                 normStd,
                 numOfComponents,
                 ubmModel,
                 numTdim,
                 projectionMatrix,
                 enrolledSpeakersByIdx,
                 pWork_path,
                 pToLogs,
                 scoringMethod,
                 T,
                 TS,
                 TSi,
                 logger);


        logger.save("finish time:              ", ava::currentTime());
    }
}


