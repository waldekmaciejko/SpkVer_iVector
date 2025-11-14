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

#include "Log4AVA.h"


namespace ava {

Log4AVA::Log4AVA(std::string dir,
                 std::string subPath,
                 std::string currentTime,
                 unsigned int gmmNumComponents,
                 unsigned int TotalVnumTdim,
                 bool verbose){

    // create directory to raport and calc
    // std::string pDirectory = dir +
    //                          subPath +
    //                          ava::currentTime()+
    //                          "_gmm"+std::to_string(gmmNumComponents)+
    //                          "_TV"+std::to_string(TotalVnumTdim);

    //std::string pMkdir = "mkdir " + pDirectory;
    //std::system(pMkdir.c_str());

    //this->fileToLog = pDirectory+"/"+currentTime+".log";
    //this->pathDirToLog = pDirectory;

    std::filesystem::path pDirectory = dir +
                                       subPath +
                                       ava::currentTime() +
                                       "_gmm"+std::to_string(gmmNumComponents)+
                                       "_TV"+std::to_string(TotalVnumTdim);

    std::string pMkdir = "mkdir " + std::string(pDirectory);
    std::system(pMkdir.c_str());

    this->fileToLog =std::string(pDirectory) +
                                    "/"+
                                    ava::currentTime() +
                                    ".log";

    this->pathDirToLog = std::string(pDirectory);
    this->verbose = verbose;
}

void Log4AVA::saveScore(arma::mat score, bool frr_far){
    if(frr_far==0)  score.save(this->pathDirToLog+"/scoreFRR", arma::csv_ascii);
    if(frr_far==1)  score.save(this->pathDirToLog+"/scoreFAR", arma::csv_ascii);
}

std::string Log4AVA::returnSavePath(){
    return this->pathDirToLog;
}

void Log4AVA::save(std::string variableName, std::string variableValue){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<variableValue<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<variableValue<<std::endl;
    }
}

void Log4AVA::save(const char* s1, const char* s2){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<std::string(s1)<<" "<<std::string(s2)<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << std::string(s1)<<" "<<std::string(s2)<<std::endl;
    }

}

void Log4AVA::save(std::string variableName, unsigned int variableValue){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<variableValue<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<variableValue<<std::endl;
    }
}

void Log4AVA::save(std::string variableName, int variableValue){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<variableValue<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<variableValue<<std::endl;
    }
}

void Log4AVA::save(std::string variableName, float variableValue){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<variableValue<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<variableValue<<std::endl;
    }
}

void Log4AVA::save(std::string variableName, bool variableValue){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<std::boolalpha<<variableValue<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<std::boolalpha<<variableValue<<std::endl;
    }
}

void Log4AVA::save(std::string variableName,  std::vector<std::string> pTrain_path){
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<std::endl;
    for(auto i : pTrain_path){
        logToFile<<i<<std::endl;
    }
    logToFile.close();

    std::cout << variableName << std::endl;

    if (this->verbose == true){

        for(auto i : pTrain_path){
            std::cout << i << std::endl;
        }
    }
}

void Log4AVA::save(std::string variableName,
                   std::stringstream& ss){
    std::string raport =ss.str();
    std::cout<<raport<<std::endl;
    std::ofstream logToFile;
    logToFile.open(this->fileToLog, std::ios::app);
    logToFile<<variableName<<" "<<raport<<std::endl;
    logToFile.close();

    if (this->verbose == true){
        std::cout << variableName<<" "<<raport<<std::endl;
    }
}
}
