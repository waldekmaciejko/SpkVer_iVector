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

#ifndef LOG4AVA_H
#define LOG4AVA_H

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>

#include "funhelpers.h"

namespace ava {

class Log4AVA{

private:
    std::string fileToLog;      // path to file with logs
    std::string pathDirToLog;   // path to dir with logs

public:
    Log4AVA(std::string dir,
                 std::string subPath,
                 std::string currentTime,
                 unsigned int gmmNumComponents,
                 unsigned int TotalVnumTdim);

    void saveScore(arma::mat score, bool frr_far);

    void save(std::string variableName, std::string variableValue);

    void save(const char* s1, const char* s2);

    void save(std::string variableName, unsigned int variableValue);

    void save(std::string variableName, float variableValue);

    void save(std::string variableName, bool variableValue);

    void save(std::string variableName,  std::vector<std::string> pTrain_path);

    void save(std::string variableName,
               std::stringstream& ss);

    std::string returnSavePath();

};

}

#endif // LOG4AVA_H