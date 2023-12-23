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

#include <ftw.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <Eigen/Dense>
#include <armadillo>
#include <chrono>
#include <regex>
#include <ctime>
#include <set>

#include <QDirIterator>

#include "mfccefeatures.h"

namespace ava{

//---------------- zestaw funkcji tworzacych liste plikow z podanej lokalizacji workSpace i zapisujacych te liste w pliku flt_fea_mfcc.lst
//---------------- brak dzialajacej funkccji bibliotecznej filesystem - dlatego to toporne rozwiazanie

int analizeDirectoryElement (const char *fpath,
                            const struct stat *sb,
                            int tflag,
                            struct FTW *ftwbuf);

void WalkDirectoryTree (const char * pchFileName);

//----------------------------- funkcja zmieniajaca \ na /

std::string winPathToPosixPathAnsi(std::string winSciezka);

//-----------------------------  zamienia sciezke POSIX na WinNT
std::string winPathToPosixPath(const char* winSciezka);

//----------------------------- funkcja generujaca listy dla biblioteki HCopy
//naneIn	-	sciezka do plikow wejsciowych (wav lub htk)
//pathOut	-	sciezka do folderu gdzie zostanĂ‚Ä… umieszczone listy
//string PathFeatureOut	-
//daneOutRozszerzenie	-	rozszerzenie plikow wyjsciowych
//count		-		mowi ile plikow znajduje sie w analizowanym folderze
//nazwaListy	-	do ktorej zostana wpisane scieki do plikow
//featureId	-		na podstawie Id oraz "i" zostana nadane nazwy plikow cech

void ListaMfcc(std::string daneIn,
                    std::string pathOut,
                    std::string pathFeatureOut,
                    std::string daneOutRozszerzenie,
                    std::string nazwaListy,
                    std::string featureId,
                    int i);

//----------------------------- funkcja badajaca rozmiar pliku
int get_file_size(std::string filename);


//----------------------------- funkcja zapisujaca macierz Egien do pliku

void save_egien_matrix_to_txt_file(std::string workSpace,
                                                       std::string flt_tmp,
                                                       std::string nazwa,
                                                       Eigen::MatrixXf& matrix);

void save_egien_matrix_to_txt_file(std::string workSpace,
                                                std::string flt_tmp,
                                                std::string nazwa,
                                                Eigen::VectorXf& matrix);
unsigned int listWavToMfcc(std::string workSpace,
                              std::string flt_fea_list_str,
                              std::string p_to_korpus,
                              std::string flt_fea_mfcc_str);

arma::mat funkTotalVariability(unsigned int numTdim,
                               arma::mat& sigma,
                               unsigned int numSpkUBM,
                               arma::mat& repSigma,
                               arma::cube& N,
                               unsigned int prodNumOfComponentsNumFeatures,
                               unsigned int numFeatures,
                               unsigned int numOfComponents,
                               arma::cube& F,
                               arma::cube& Nc,
                               unsigned int numIterationsTV);

arma::mat funkTotalVariabilityVector(unsigned int numTdim,
                                       arma::mat& sigma,
                                       unsigned int numSpkUBM,
                                       unsigned int numFilesUBM,
                                       arma::mat& repSigma,
                                        std::vector<arma::mat>& NN,
                                       unsigned int prodNumOfComponentsNumFeatures,
                                       unsigned int numFeatures,
                                       unsigned int numOfComponents,
                                       std::vector<arma::mat>& FF,
                                       std::vector<arma::mat>& NNcc,
                                       unsigned int numIterationsTV);

std::map<std::string, arma::mat> funkiVectorFromUBM(arma::mat& repSigma,
                                                    arma::mat& T,
                                                    arma::mat& ubmMu,
                                                    uint16_t numTdim,
                                                    std::set<std::string>& setSpk,
                                                    std::multimap<std::string, std::string>& multimapMowcyMFC,
                                                    unsigned int numOfComponents,
                                                    arma::mat& normMean,
                                                    arma::mat& normStd,
                                                    unsigned int numFeatures,
                                                    arma::gmm_diag& ubmModel,
                                                    arma::mat& ivectorsTrain,
                                                    arma::mat& TS,
                                                    arma::mat& TSi);

/*!
 * \brief takeSpkUBM_MLabel -   is a key function in parsing UBM PCM WAV files paths
 * \param pathToWav
 * \return
 */
std::string takeSpkLabel(std::string pathToFile);

arma::mat normZ(arma::mat Y,
                  arma::mat normMean,
                  arma::mat normStd,
                  unsigned int liczbaRamekY,
                  unsigned int numFeatures);

arma::mat multipicColByCol(arma::mat& A, arma::mat& B);

arma::mat multipicColByCol(arma::mat& A, arma::mat&& B);

arma::mat takeMostSignificantEigenVecLDA(std::map<std::string, arma::mat>& ivectorPerSpk,
                                         unsigned int numTdim,
                                         arma::mat& ivectorsTrain,
                                         uint numEigenvectors,
                                         std::vector<uint>& utterancePerSpeaker);

arma::mat takeMostSignificantEigenVecLDA2(std::map<std::string, arma::mat>& ivectorPerSpk,
                                              unsigned int numTdim,
                                              arma::mat& ivectorsTrain,
                                              uint numEigenvectors,
                                              std::vector<uint>& utterancePerSpeaker);

void performaLDAfunk(uint numEigenvectors,
                     bool performLDA,
                     uint numTdim,
                     std::map<std::string, arma::mat>& ivectorPerSpk,
                     arma::mat& ivectorsTrain,
                     arma::mat& projectionMatrix,
                     std::vector<arma::mat>& w,
                     std::vector<uint>& utterancePerSpeaker);

void performaWCCNfunk(uint numEigenvectors,
                      arma::mat& projectionMatrix,
                      std::vector<arma::mat>& w);

void performWhitening(std::map<std::string, arma::mat>& ivectorPerSpk,
                      arma::mat& ivectorMatrix,
                      arma::mat& projectionMatrix,
                      std::string zmiennWhitening,
                      arma::mat& mu,
                      arma::mat& W);

void normArmaMatrix(arma::mat& matToNorm,
                    arma::mat& vecAfterNorm);

void divdeMatrixRBRVectorNorm(arma::mat& matrixToDivide,
                              arma::mat& matrixAfterDivide);

void divideVectorIntoSpeakerCells(std::vector<uint>& utterancePerSpeaker,
                                  arma::mat& ivectorMatrixAfterNorm,
                                  std::vector<arma::mat>& ivector,
                                  uint numEigenVectors);

void sortIVectors(std::vector<std::vector<arma::mat>>& ivectorsSorted,
                  arma::mat& ff,
                  std::vector<uint>& utterancePerSpeaker,
                  std::vector<arma::mat>& ivector);

void trainGPLDA(uint numEigenVoices,
                arma::mat& V,
                arma::mat& lambda,
                uint numIter,
                bool minimumDivergence,
                std::vector<uint>& utterancePerSpeaker,
                uint K,
                std::vector<std::vector<arma::mat>>& ivectorSorted,
                arma::mat& ff,
                uint numFilesUBM,
                arma::mat& S);

arma::mat readMFCCHCopyFile(std::string path_to_mfc);


void readRecursivlyMFCCHcopyFile(std::string pWork_path,
                                      std::string pFeaMfcc_path,
                                      arma::mat& normMean,
                                      arma::mat& normStd,
                                      arma::mat& normYtrans);

void statBaumWelch(arma::cube& N,
                   arma::cube& F,
                   arma::cube& Nc,
                   std::string pWork_path,
                   std::string pFeaMfcc_path,
                   arma::mat& normMean,
                   arma::mat& normStd,
                   arma::gmm_diag& ubmModel,
                   const unsigned int numOfComponents,
                   unsigned int numFilesUBM,
                   unsigned int numFeatures,
                   unsigned int prodNumOfComponentsNumFeatures);

void statBaumWelchVector(std::vector<arma::mat>& NN,
                           std::vector<arma::mat>& FF,
                           std::vector<arma::mat>& NNcc,
                           std::string pWork_path,
                           std::string pFeaMfcc_path,
                           arma::mat& normMean,
                           arma::mat& normStd,
                           arma::gmm_diag& ubmModel,
                           const unsigned int numOfComponents,
                           unsigned int numFilesUBM,
                           unsigned int numFeatures,
                           unsigned int prodNumOfComponentsNumFeatures);

void randomTestSplit(std::vector<uint>& ranTestNumVector,
                     std::vector<uint>& ranEnrollNumVec,
                     uint numSpkEnrollTestFiles,
                     float test_size);

arma::mat normZprim(arma::mat Y,
                  arma::mat normMean,
                  arma::mat normStd,
                  unsigned int liczbaRamekY,
                  unsigned int numFeatures);

std::string currentTime();

std::string extractFNameFrPath(std::string pathToMFCFile);

void saveEnrolFRRFARlList(std::string spkLabel,
                             std::string ss,
                             std::string pathToSave,
                             std::string nameOfTheList);

}


