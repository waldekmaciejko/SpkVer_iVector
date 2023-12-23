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

#ifndef ENROLLFRRFAR_H
#define ENROLLFRRFAR_H

#include "mfccefeatures.h"
#include "mfccextractor.h"
#include "funhelpers.h"
#include "filelists.h"
#include "gpldamodel.h"
#include "Log4AVA.h"

namespace ava{

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
            ava::Log4AVA& logger);

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
            ava::Log4AVA& logger);


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
         ava::Log4AVA& logger);

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
         ava::Log4AVA& logger);
}

#endif // ENROLLFRRFAR_H
