#ifndef CONFIG-FILEPARSER_H
#define CONFIG-FILEPARSER_H

#include <iostream>
#include <filesystem>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

/*
 * This class contains tools for parsing configuration files
 * from command line, like config.cfg and paths.cfg
 */

class FileParser
{
public:
    FileParser();

    void parseConfigFile(char **argv, int ndx, std::string flag);
    bool strToBool(std::string inp);
    std::string extractValueFromLine(std::string line);


};

// ==============================================================

class FileParserCfg : private FileParser{

public:
    FileParserCfg();
    void parseConfigFile(char **argv, int ndx, std::string flag);

    //---------------------------------
    //--------------------- SET PARAMETERS
    //---------------------------------

    std::string pHCopy_bin{};

    std::string pEnrollList_file{};
    std::string pEnrollMfcc_path{};

    //--GMM

    int numOfComponents{};     // 1024
    int numIterationsEM{};            // test: 3
    int numIterationsKmean_init{};

    //--Total Variability

    int numTdim{}; // test:100 powinna wynosic około 1000
    int numIterationsTV{};   //optymalna wartosc wynosi 20

    //--Projection Matrix

    bool performLDA{};
    bool performWCCN{};
    int numEigenVectors{};

    //--G-PLDA

    std::string zmiennWhitening{}; //PCA, 000 - brak whitening
    int numEigenVoices{};
    int numIter{};
    bool minimumDivergence{};
    bool SAVETRAINMODEL{};  //-- SAVE TRAINED MODEL?
    bool train{};           //-- TRAIN TOTAL VARIABILITY, PROJECTION MATRIX, iVECTOR?
    bool enroll{};

    float test_size{}; //-- BALANCE ENROLL/TEST
    int scoringMethod{}; // 1 - CSS, 0 - GPLDA

private:
    std::vector<std::string> cfg_list{};
    std::map<std::string, std::string> parserCfg{};

};

class FileParserWrkSpace : private FileParser{

public:
    FileParserWrkSpace();
    void parseConfigFile(char **argv, int ndx, std::string flag);

    //*************** WORKING DIRECTORYS

    std::string pTrainList_file{};
    std::string pTrainMfcc_path{};
    std::string pCfg_path{};
    std::string pToSaveWSpace{};
    std::string pTmp_file{};
    std::string pToLogs{};

    std::string pEnrollList_file{};
    std::string pEnrollMfcc_path{};

private:
    std::vector<std::string> cfg_list{};
    std::map<std::string, std::string> parserCfg{};

};

class FileParserData : private FileParser{

public:
    FileParserData();
    void parseConfigFile(char **argv, int ndx, std::string flag);

    std::string pWork_path{};
    std::string pToModel{};

    std::vector<std::string> getList_train();
    std::vector<std::string> getList_test();

private:
    std::vector<std::string> cfg_list{};
    std::map<std::string, std::string> parserCfg{};

    std::vector<std::string> list_train{};
    std::vector<std::string> list_test{};

};

#endif // CONFIG-FILEPARSER_H
