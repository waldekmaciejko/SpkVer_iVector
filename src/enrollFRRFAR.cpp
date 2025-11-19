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

#include "enrollFRRFAR.h"

namespace ava {

/*!
 * \brief enroll - build the enroll models
 * \param setSpkEnrollTest
 * \param multimapEnrollMFC
 * \param pWork_path
 * \param pToLogs
 * \param normMean
 * \param normStd
 * \param ubmModel
 * \param projectionMatrix
 * \param numOfComponents
 * \param numTdim
 * \param T
 * \param TS
 * \param TSi
 * \param enrolledSpeakersByIdx
 * \param logger
 */
void enroll(std::set<std::string>& setSpkEnrollTest,
            std::multimap<std::string, std::string>& multimapEnrollMFC,
            std::string pWork_path,
            std::string pToLogs,
            arma::mat& normMean,
            arma::mat& normStd,
            arma::gmm_diag& ubmModel,
            arma::mat& projectionMatrix,
            uint numOfComponents,
            uint numTdim,
            arma::mat& T,
            arma::mat& TS,
            arma::mat& TSi,
            std::multimap<std::string, arma::mat>& enrolledSpeakersByIdx,
            ava::Log4AVA& logger){



    auto start = std::chrono::high_resolution_clock::now();

    arma::mat tmp{};
    arma::mat Yenroll{};
    arma::mat logLikelihood = {};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat amax={};
    arma::mat ivector{};
    arma::mat n={};
    arma::mat f={};
    arma::mat ivectorMat{};
    arma::mat I=arma::eye(numTdim, numTdim);

    for(auto speakerIdx : setSpkEnrollTest){
        for(auto para : multimapEnrollMFC){
            if(para.first==speakerIdx){

                ava::saveEnrolFRRFARlList(speakerIdx,
                                            para.second,
                                            pWork_path+pToLogs,
                                            "enroll");

                tmp = ava::readMFCCHCopyFile(para.second);
                uint numFeatures=tmp.n_cols;
                Yenroll = ava::normZ(tmp.t(),
                                          normMean,
                                          normStd,
                                          tmp.n_rows,
                                          numFeatures);

                for(unsigned int comp=0; comp<numOfComponents; comp++)
                {
                    logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Yenroll.t(), comp));
                }

                // Compute a posteriori normalized probability
                amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
                logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
                gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

                /*
                  *  statystyki Bauma-Welcha
                 */

                 n = arma::sum(gamma, 0);
                 //f = arma::operator*(Yenroll.t(), gamma)-arma::repelem(n, numFeatures, 1);
                 f = arma::operator*(Yenroll.t(), gamma)- arma::repmat(n, numFeatures, 1)%ubmModel.means;

                 //i-vector Extraction
                 ivector = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);

                 // Intersession Compensation
                 ivector = projectionMatrix * ivector;

                 ivectorMat = arma::join_rows(ivectorMat, ivector);

                 logLikelihoodSum.clear();
                 logLikelihood.clear();
            }
        }
        enrolledSpeakersByIdx.insert({speakerIdx, arma::mean(ivectorMat, 1)});
    }


    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    // std::stringstream ss{};
    // ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
    // logger.save("log:", ss);

    std::string verbose_str = "enroll module - line: "+ std::to_string(__LINE__) + ", (duration sec): " + std::to_string(duration.count());
    logger.save("Execution in point ", verbose_str);

}

void enrollII(std::set<std::string>& setSpkEnrollTest,
            std::multimap<std::string, std::string>& multimapEnrollMFC,
            std::string pWork_path,
            std::string pToLogs,
            arma::mat& normMean,
            arma::mat& normStd,
            arma::gmm_diag& ubmModel,
            arma::mat& projectionMatrix,
            uint numOfComponents,
            uint numTdim,
            arma::mat& T,
            arma::mat& TS,
            arma::mat& TSi,
            std::multimap<std::string, arma::mat>& enrolledSpeakersByIdx,
            ava::Log4AVA& logger){



    auto start = std::chrono::high_resolution_clock::now();

    arma::mat tmp{};
    arma::mat Yenroll{};
    arma::mat logLikelihood = {};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat amax={};
    arma::mat ivector{};
    arma::mat n={};
    arma::mat f={};
    arma::mat ivectorMat{};
    arma::mat I=arma::eye(numTdim, numTdim);
    uint numFeatures={};

    for(auto speakerIdx : setSpkEnrollTest){
        for(auto para : multimapEnrollMFC){
            if(para.first==speakerIdx){

                ava::saveEnrolFRRFARlList(speakerIdx,
                                            para.second,
                                            pWork_path+pToLogs,
                                            "enroll");

                tmp = ava::readMFCCHCopyFile(para.second);
                //tmp = arma::join_cols(tmp, ava::readMFCCHCopyFile(para.second));
                //uint numFeatures=tmp.n_cols;
            //}
            uint numFeatures=tmp.n_cols;
            Yenroll = ava::normZ(tmp.t(),
                                      normMean,
                                      normStd,
                                      tmp.n_rows,
                                      numFeatures);

            for(unsigned int comp=0; comp<numOfComponents; comp++)
            {
                logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Yenroll.t(), comp));
            }

            // Compute a posteriori normalized probability
            amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
            logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
            gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

            /*
              *  statystyki Bauma-Welcha
             */

             n = arma::sum(gamma, 0);
             //f = arma::operator*(Yenroll.t(), gamma)-arma::repelem(n, numFeatures, 1);
             f = arma::operator*(Yenroll.t(), gamma)- arma::repmat(n, numFeatures, 1)%ubmModel.means;

             //i-vector Extraction
             ivector = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);

             // Intersession Compensation
             ivector = projectionMatrix * ivector;

             ivectorMat = arma::join_rows(ivectorMat, ivector);

             logLikelihoodSum.clear();
             logLikelihood.clear();
            }
        }
        enrolledSpeakersByIdx.insert({speakerIdx, arma::mean(ivectorMat, 1)});
        ivectorMat.clear();
    }


    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    // std::stringstream ss{};
    // ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
    // logger.save("log:", ss);

    std::string verbose_str = "enroll module - line: "+ std::to_string(__LINE__) + ", (duration sec): " + std::to_string(duration.count());
    logger.save("Execution in point ", verbose_str);
}


void FRR(std::multimap<std::string, std::string>& multimapTestMFC,
         arma::mat& normMean,
         arma::mat& normStd,
         uint numOfComponents,
         arma::gmm_diag& ubmModel,
         uint numTdim,
         arma::mat& projectionMatrix,
         std::multimap<std::string, arma::mat>& enrolledSpeakersByIdx,
         std::string pWork_path,
         std::string pToLogs,
         uint scoringMethod,
         arma::mat& T,
         arma::mat& TS,
         arma::mat& TSi,
         ava::Log4AVA& logger){

    // the assumption: there is only ona key for each speaker
    // label in enrolledSpeakersByIdx

    arma::mat ivectorToTest{};
    arma::mat scoreFRR{};
    arma::mat score{};

    arma::mat tmp{};
    arma::mat Yenroll{};
    arma::mat logLikelihood = {};
    arma::mat amax={};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat n={};
    arma::mat f={};
    arma::mat I=arma::eye(numTdim, numTdim);
    arma::mat ivectorMat{};
    arma::mat ivector{};

    //*************** FALSE REJECTION RATE

    auto start = std::chrono::high_resolution_clock::now();

         for(auto para : multimapTestMFC){

                tmp = ava::readMFCCHCopyFile(para.second);
                uint numFeatures=tmp.n_cols;
                Yenroll = ava::normZ(tmp.t(),
                                          normMean,
                                          normStd,
                                          tmp.n_rows,
                                          numFeatures);

                for(unsigned int comp=0; comp<numOfComponents; comp++){
                    logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Yenroll.t(), comp));
                }

                // Compute a posteriori normalized probability
                amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
                logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
                gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

                /*
                  *  statystyki Bauma-Welcha
                */

                 n = arma::sum(gamma, 0);
                 //f = arma::operator*(Yenroll.t(), gamma)-arma::repelem(n, numFeatures, 1);
                 f = arma::operator*(Yenroll.t(), gamma)- arma::repmat(n, numFeatures, 1)%ubmModel.means;

                 //i-vector Extraction
                 ivector = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);

                 // Intersession Compensation
                 ivector = projectionMatrix * ivector;

                 ava::saveEnrolFRRFARlList(para.first,
                                             enrolledSpeakersByIdx.find(para.first)->first,
                                             pWork_path+pToLogs,
                                             "FRR");

                 //std::cout<<"test  : "<<para.first<<std::endl;
                 //std::cout<<"enroll:"<<enrolledSpeakersByIdx.find(para.first)->first<<std::endl;

                 // 1. in multimapTestMFC find label of the speaker (para.first)
                 // 2. than in (multimap) enrolledSpeakersByIdx find label from 1.
                 // 3. and put the matrix into the arma::mat  ivectorToTest
                 ivectorToTest = enrolledSpeakersByIdx.find(para.first)->second;

                 if(scoringMethod == 1){
                    //CSS
                     score = arma::dot(ivectorToTest, ivector)/(arma::norm(ivector)*arma::norm(ivectorToTest));
                     scoreFRR = arma::join_cols(scoreFRR, score);

                  }
                  else if (scoringMethod == 0) {
                    //GPLDA
                    std::cout<<"TODO"<<std::endl;
                  }
                 logLikelihoodSum.clear();
                 logLikelihood.clear();
    }

         logger.saveScore(scoreFRR, 0);
         scoreFRR.clear();
         ivectorToTest.clear();
         score.clear();

         auto stop = std::chrono::high_resolution_clock::now();
         std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
         // std::stringstream ss{};
         // ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
         // logger.save("log:", ss);

         std::string verbose_str = "FRR was finished - line: "+ std::to_string(__LINE__) + ", (duration sec): " + std::to_string(duration.count());
         logger.save("Execution in point ", verbose_str);
}

void FAR(std::multimap<std::string, std::string>& multimapTestMFC,
         arma::mat& normMean,
         arma::mat& normStd,
         uint numOfComponents,
         arma::gmm_diag& ubmModel,
         uint numTdim,
         arma::mat& projectionMatrix,
         std::multimap<std::string, arma::mat>& enrolledSpeakersByIdx,
         std::string pWork_path,
         std::string pToLogs,
         uint scoringMethod,
         arma::mat& T,
         arma::mat& TS,
         arma::mat& TSi,
         ava::Log4AVA& logger){

    arma::mat ivectorToTest{};
    arma::mat scoreFAR{};
    arma::mat score{};

    arma::mat tmp{};
    arma::mat Yenroll{};
    arma::mat logLikelihood = {};
    arma::mat amax={};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat n={};
    arma::mat f={};
    arma::mat I=arma::eye(numTdim, numTdim);
    arma::mat ivectorMat{};
    arma::mat ivector{};

    //*************** FALSE ACCEPTANCE RATE

         auto start = std::chrono::high_resolution_clock::now();

         for(auto para : multimapTestMFC){

                tmp = ava::readMFCCHCopyFile(para.second);
                uint numFeatures=tmp.n_cols;
                Yenroll = ava::normZ(tmp.t(),
                                          normMean,
                                          normStd,
                                          tmp.n_rows,
                                          numFeatures);

                for(unsigned int comp=0; comp<numOfComponents; comp++){
                    logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Yenroll.t(), comp));
                }

                // Compute a posteriori normalized probability
                amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
                logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
                gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

                /*
                  *  statystyki Bauma-Welcha
                */

                 n = arma::sum(gamma, 0);
                 //f = arma::operator*(Yenroll.t(), gamma)-arma::repelem(n, numFeatures, 1);
                 f = arma::operator*(Yenroll.t(), gamma)- arma::repmat(n, numFeatures, 1)%ubmModel.means;

                 //i-vector Extraction
                 ivector = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);

                 // Intersession Compensation
                 ivector = projectionMatrix * ivector;

                 for(auto paraEnroll : enrolledSpeakersByIdx){
                     if(paraEnroll.first != para.first){

                     //std::cout<<"test  : "<<para.first<<std::endl;
                     //std::cout<<"enroll:"<<enrolledSpeakersByIdx.find(para.first)->first<<std::endl;

                     ava::saveEnrolFRRFARlList(para.first,
                                                 paraEnroll.first,
                                                 pWork_path+pToLogs,
                                                 "FAR");

                     // 1. in multimapTestMFC find label of the speaker (para.first)
                     // 2. than in (multimap) enrolledSpeakersByIdx find label from 1.
                     // 3. and put the matrix into the arma::mat  ivectorToTest
                     ivectorToTest = paraEnroll.second;

                     if(scoringMethod == 1){
                        //CSS
                         score = arma::dot(ivectorToTest, ivector)/(arma::norm(ivector)*arma::norm(ivectorToTest));
                         scoreFAR = arma::join_cols(scoreFAR, score);

                      }
                      else if (scoringMethod == 0) {
                        //GPLDA
                        std::cout<<"TODO"<<std::endl;
                      }
                     }
                 }
                 logLikelihoodSum.clear();
                 logLikelihood.clear();
         }

     auto stop = std::chrono::high_resolution_clock::now();
     std::chrono::duration<double>  duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
     // std::stringstream ss{};
     // ss<<"file: "<<__FILE__<<" | line: "<<__LINE__<<" | duration since previous point [sec]: "<<duration.count();
     // logger.save("log:", ss);

     std::string verbose_str = "FAR was finished - line: "+ std::to_string(__LINE__) + ", (duration sec): " + std::to_string(duration.count());
     logger.save("Execution in point ", verbose_str);

    //*************** END FALSE ACCEPTANCE RATE

    logger.saveScore(scoreFAR, 1);
    scoreFAR.clear();
    ivectorToTest.clear();
    score.clear();
}
}
