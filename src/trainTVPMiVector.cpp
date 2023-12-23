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

#include "trainTVPMiVector.h"

namespace ava{


    void trainTVPMiVector(std::vector<std::string> pTrain_path,
                            std::string pWork_path,
                            std::string pTrainList_file,
                            std::string pTrainMfcc_path,
                            std::string pCfg_path,
                            std::string pToSaveWSpace,
                            std::string pTmp_file,
                            std::string pHCopy_bin,
                            const unsigned int numOfComponents,
                            unsigned int numIterationsEM,
                            unsigned int numIterationsKmean_init,
                            unsigned int numTdim,
                            unsigned int numIterationsTV,
                            bool performLDA,
                            bool performWCCN,
                            uint numEigenVectors,
                            std::string zmiennWhitening,
                            uint numEigenVoices,
                            uint numIter,
                            bool minimumDivergence,
                            bool SAVETRAINMODEL,
                            ava::Log4AVA& logger){

        logger.save("*************************", "FUNCTION trainTVPMiVector() START");

        auto start = std::chrono::high_resolution_clock::now();

        ava::FileLists fileListUBM(pWork_path,
                                       pTrainList_file,
                                       pTrain_path,
                                       pTrainMfcc_path);

        std::multimap<std::string, std::string> multimapUBMSpkMFCC{};
        std::set<std::string> setSpk{};
        unsigned int numFilesUBM{};
        unsigned int numSpkUBM{};

        fileListUBM.listWavToMfccUBM(multimapUBMSpkMFCC,
                                      numFilesUBM,
                                      numSpkUBM,
                                      setSpk);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        std::stringstream ss{};
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------- generate config for HTK HCopy.exe

        start = std::chrono::high_resolution_clock::now();

        ava::MfccExtractor *ME;                    // parametry MFCC zdefiniowano w pliku mfccExtractor
        ME=new ava::MfccExtractor();
        std::string pCfgHcopy_file=ME->generujConfigMFCC(pWork_path+
                                                            pCfg_path, "conf.mfcc");    //-- conf.mfcc file config

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------------- create string to send to console and call HCopy
        //----------------------------------- result: mfcc files
        start = std::chrono::high_resolution_clock::now();

        ME->wywolajHCopy(pHCopy_bin,
                         pCfgHcopy_file,
                         pWork_path,
                         pTrainList_file);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------------- building UBM
        //  1. read generated MFCC and concatenate (dozen files - one speaker)

        arma::mat normMean{};
        arma::mat normStd{};
        arma::mat normYtrans{};

        start = std::chrono::high_resolution_clock::now();

        ava::readRecursivlyMFCCHcopyFile(pWork_path,
                                              pTrainMfcc_path,
                                              normMean,
                                              normStd,
                                              normYtrans);

        unsigned int liczbaRamekY = normYtrans.n_rows;
        unsigned int numFeatures = normYtrans.n_cols;

        //normYtrans.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\workspace_maw200\\tmp\\normYtrans.csv", arma::csv_ascii);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        start = std::chrono::high_resolution_clock::now();
        arma::gmm_diag ubmModel;

        /*
         *
        (armadill docs) learn the model parameters via multi-threaded k-means and/or EM algorithms;
        return a bool value, with true indicating success, and false indicating failure;
        the parameters have the following meanings:
        */
        bool status = ubmModel.learn(normYtrans.t(),
                                      numOfComponents,
                                      arma::eucl_dist,
                                      arma::random_subset,
                                      numIterationsKmean_init,
                                      numIterationsEM,
                                      1e-10,
                                      false);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        normYtrans.clear();

        unsigned int prodNumOfComponentsNumFeatures=numOfComponents*numFeatures;

        std::vector<arma::mat> N{};
        N.reserve(numFilesUBM);
        std::vector<arma::mat> F{};
        F.reserve(numFilesUBM);
        std::vector<arma::mat> Nc{};
        Nc.reserve(numFilesUBM);

        ava::statBaumWelchVector(N,
                                 F,
                                 Nc,
                                 pWork_path,
                                 pTrainMfcc_path,
                                 normMean,
                                 normStd,
                                 ubmModel,
                                 numOfComponents,
                                 numFilesUBM,
                                 numFeatures,
                                 prodNumOfComponentsNumFeatures);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------------- Total Variability
        /*
         *
         * Dehak. N, Kenny P, Dehak R, Dumouchel P, Ouellet P, "Front-end factor analysis for speaker verification"
         * Kenny P. Ouellet P. Dehak N. Gupta V. Dumouchel P. "A study of inter-speaker variability in speaker verification"
         *
         */
        start = std::chrono::high_resolution_clock::now();

        arma::mat sigma = arma::reshape(ubmModel.dcovs, prodNumOfComponentsNumFeatures, 1);
        arma::mat repSigma = arma::repmat(sigma, 1, numTdim);

        arma::mat T = ava::funkTotalVariabilityVector(numTdim,
                                                       sigma,
                                                       numSpkUBM,
                                                       numFilesUBM,
                                                       repSigma,
                                                       N,
                                                       prodNumOfComponentsNumFeatures,
                                                       numFeatures,
                                                       numOfComponents,
                                                       F,
                                                       Nc,
                                                       numIterationsTV);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------------- i-Vector
        /*
         *
         * Dehak. N, Kenny P, Dehak R, Dumouchel P, Ouellet P, "Front-end factor analysis for speaker verification"
         * Kenny P. Ouellet P. Dehak N. Gupta V. Dumouchel P. "A study of inter-speaker variability in speaker verification"
         *
         */

        arma::mat ubmMu = ubmModel.means;
        arma::mat I=arma::eye(numTdim, numTdim);
        std::map<std::string, arma::mat> ivectorPerSpk{};
        start = std::chrono::high_resolution_clock::now();
        arma::mat ivectorsTrain{}; // matrix put on funkiVectorFromUBM as ref and fullfiled as [numTdim, numFiles]
        arma::mat TS =T/repSigma;
        arma::mat TSi = T.t();

        ivectorPerSpk = ava::funkiVectorFromUBM(repSigma,
                                                    T,
                                                    ubmMu,
                                                    numTdim,
                                                    setSpk,
                                                    multimapUBMSpkMFCC,
                                                    numOfComponents,
                                                    normMean,
                                                    normStd,
                                                    numFeatures,
                                                    ubmModel,
                                                    ivectorsTrain,
                                                    TS,
                                                    TSi);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //----------------------------------- Projection Matrix
        /*
         *
         * Dehak. N, Kenny P, Dehak R, Dumouchel P, Ouellet P, "Front-end factor analysis for speaker verification"
         *
         * LDA-WCCN combination
         *
         *
         */

        arma::mat projectionMatrix=arma::eye(numTdim, numTdim);
        std::vector<arma::mat> w{};
        std::vector<uint> utterancePerSpeaker{};

        if(performLDA){
            start = std::chrono::high_resolution_clock::now();

            ava::performaLDAfunk(numEigenVectors,
                                 performLDA,
                                 numTdim,
                                 ivectorPerSpk,
                                 ivectorsTrain,
                                 projectionMatrix,
                                 w,
                                 utterancePerSpeaker);

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
            ss.str("");
            ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
            logger.save("log:", ss);
        }

        /* WCCN
         *
         * Dehak Najim, Reda Dehak, James R. Glass, Douglas A. Reynolds, Patric Kenny
         * "Cosine Similarity Scoring without Score Normalization Techniques" Odyssey 2010.
         *
         */

        if(performWCCN ){
            start = std::chrono::high_resolution_clock::now();

            ava::performaWCCNfunk(numEigenVectors,
                                  projectionMatrix,
                                  w);

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
            ss.str("");
            ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
            logger.save("log:", ss);
        }

        /* G-PLDA
         *
         * Sizov A, Kong A.L., Kinnunen T, Unifying probabilistic linear discriminant
         * analysis variants in biometric authentication
         *
         */
        arma::mat ivectorMatrix{};
        arma::mat mu{};
        arma::mat W{};

        start = std::chrono::high_resolution_clock::now();

        ava::performWhitening(ivectorPerSpk,
                                ivectorMatrix,
                                projectionMatrix,
                                zmiennWhitening,
                                mu,
                                W);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);

        //Test12
        arma::mat ivectorMatrixAfterNorm(ivectorMatrix.n_rows, ivectorMatrix.n_cols, arma::fill::randn);

        ava::divdeMatrixRBRVectorNorm(ivectorMatrix,
                                      ivectorMatrixAfterNorm);

        arma::mat S = ivectorMatrixAfterNorm*ivectorMatrixAfterNorm.t();
        std::vector<arma::mat> ivector{};

        ava::divideVectorIntoSpeakerCells(utterancePerSpeaker,
                                          ivectorMatrixAfterNorm,
                                          ivector,
                                          numEigenVectors);

        uint D = numEigenVectors;
        uint K = numSpkUBM;

        // moment pierwszego rzedu defniujacy cechy danego mowcy
        arma::mat ff(D, K, arma::fill::zeros);
        std::vector<std::vector<arma::mat>> ivectorSorted{};

        ava::sortIVectors(ivectorSorted,
                          ff,
                          utterancePerSpeaker,
                          ivector);
    //-------------------------
        arma::mat V(D, numEigenVoices, arma::fill::randn);
        arma::mat lambda =arma::pinv(S/numFilesUBM);

        //lambda.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\lambda2.csv", arma::csv_ascii);
        //S.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\S2.csv", arma::csv_ascii);
        //ivectorSorted.at(0).at(0).save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ivectorSorted00.csv", arma::csv_ascii);
        //ff.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ff.csv", arma::csv_ascii);

        //ivectorSorted.at(1).at(1).save("j:\\__Badania__\\QtGUI_ava300\\ava300\\workspace_maw200\\tmp\\ivectorSorted_1_1.csv", arma::csv_ascii);

        start = std::chrono::high_resolution_clock::now();

        ava::trainGPLDA(numEigenVoices,
                            V,
                            lambda,
                            numIter,
                            minimumDivergence,
                            utterancePerSpeaker,
                            K,
                            ivectorSorted,
                            ff,
                            numFilesUBM,
                            S);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        ss.str("");
        ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
        logger.save("log:", ss);
        //std::cout<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration sec: "<<duration.count()<<std::endl;

        arma::mat lam=arma::pinv(lambda);
        ava::gpldaModel gpldaMod(mu, W, V, lam);

        if(SAVETRAINMODEL){

            projectionMatrix.save(logger.returnSavePath()+"\\projectionMatrix.csv", arma::csv_ascii);
            mu.save(logger.returnSavePath()+"\\mu.csv", arma::csv_ascii);
            W.save(logger.returnSavePath()+"\\W.csv", arma::csv_ascii);
            V.save(logger.returnSavePath()+"\\V.csv", arma::csv_ascii);
            lam.save(logger.returnSavePath()+"\\lam.csv", arma::csv_ascii);
            ubmMu.save(logger.returnSavePath()+"\\ubmMu.csv", arma::csv_ascii);
            ubmModel.save(logger.returnSavePath()+"\\ubmModel.gmm");
            T.save(logger.returnSavePath()+"\\T.csv", arma::csv_ascii);
            TS.save(logger.returnSavePath()+"\\TS.csv", arma::csv_ascii);
            TSi.save(logger.returnSavePath()+"\\TSi.csv", arma::csv_ascii);
            normMean.save(logger.returnSavePath()+"\\normMean.csv", arma::csv_ascii);
            normStd.save(logger.returnSavePath()+"\\normStd.csv", arma::csv_ascii);

            logger.save("MODEL SAVED IN: ", logger.returnSavePath());

        }
        logger.save("*************************", "FUNCTION trainTVPMiVector() STOP");
    }
}
