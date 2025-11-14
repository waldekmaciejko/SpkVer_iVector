#include "config-fileparser.h"


FileParser::FileParser()
{
    int i {0};
}

void FileParser::parseConfigFile(char **argv, int ndx, std::string flag)
{
    int i {0};
}

bool FileParser::strToBool(std::string inp) {
        bool op;
        std::istringstream(inp) >> std::boolalpha >> op;
        return op;
}

std::string FileParser::extractValueFromLine(std::string line)
{
    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
    // determine start of value
    uint pos_start = line.find("=");
    // determine end of value
    uint pos_end = line.find(";");
    std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);
    return p1;
}


FileParserCfg::FileParserCfg()
{
    this->cfg_list.push_back("pHCopy_bin");

    //--GMM

    this->cfg_list.push_back("numOfComponents");           // 1024
    this->cfg_list.push_back("numIterationsEM");            // test: 3
    this->cfg_list.push_back("numIterationsKmean_init");

    //--Total Variability

    this->cfg_list.push_back("numTdim"); // test:100 powinna wynosic około 1000
    this->cfg_list.push_back("numIterationsTV"); //optymalna wartosc wynosi 20

    //--Projection Matrix

    this->cfg_list.push_back("performLDA");
    this->cfg_list.push_back("performWCCN");
    this->cfg_list.push_back("numEigenVectors");

    //--G-PLDA

    this->cfg_list.push_back("zmiennWhitening"); //PCA, 000 - brak whitening
    this->cfg_list.push_back("numEigenVoices");
    this->cfg_list.push_back("numIter");
    this->cfg_list.push_back("minimumDivergence");
    this->cfg_list.push_back("SAVETRAINMODEL");      //-- SAVE TRAINED MODEL?
    this->cfg_list.push_back("train");               //-- TRAIN TOTAL VARIABILITY, PROJECTION MATRIX, iVECTOR?
    this->cfg_list.push_back("enroll");

    this->cfg_list.push_back("test_size");          //-- BALANCE ENROLL/TEST
    this->cfg_list.push_back("scoringMethod");
}

void FileParserCfg::parseConfigFile(char **argv, int ndx, std::string flag){

    if(!std::filesystem::exists(argv[ndx+1])){
        std::cerr << "File: "<<argv[ndx+1]<<" not found \n";
        throw;
    }else{
        std::fstream s{argv[ndx+1], s.in};
        if(!s.is_open()){
            std::cerr << " Problem with open the file: "<<argv[ndx+1]<<'\n';
            throw;
        }if(flag == "-c"){
            for(std::string line; std::getline(s, line);){
                for(std::vector<std::string>::iterator itr = this->cfg_list.begin();
                     itr != this->cfg_list.end();
                     ++itr){
                    if(!line.find('#'))
                        continue;
                    if(!line.find(*itr)){
                        // remove spaces in string line has read from file
                        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
                        // determine start of value
                        uint pos_start = line.find("=");
                        // determine end of value
                        uint pos_end = line.find(";");
                        std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);
                        this->parserCfg[*itr] = p1;
                    }
                }
            }

            this->pHCopy_bin = this->parserCfg.at("pHCopy_bin");
            //--GMM

            this->numOfComponents = std::stoi(this->parserCfg.at("numOfComponents"));          // 1024
            this->numIterationsEM = std::stoi(this->parserCfg.at("numIterationsEM"));            // test: 3
            this->numIterationsKmean_init = std::stoi(this->parserCfg.at("numIterationsKmean_init"));
            //--Total Variability

            this->numTdim = std::stoi(this->parserCfg.at("numTdim")); // test:100 powinna wynosic około 1000
            this->numIterationsTV = std::stoi(this->parserCfg.at("numIterationsTV")); //optymalna wartosc wynosi 20

            //--Projection Matrix

            this->performLDA = strToBool(this->parserCfg.at("performLDA"));
            this->performWCCN = strToBool(this->parserCfg.at("performWCCN"));
            this->numEigenVectors = std::stoi(this->parserCfg.at("numEigenVectors"));

            //--G-PLDA

            this->zmiennWhitening = this->parserCfg.at("zmiennWhitening"); //PCA, 000 - brak whitening
            this->numEigenVoices = std::stoi(this->parserCfg.at("numEigenVoices"));
            this->numIter = std::stoi(this->parserCfg.at("numIter"));
            this->minimumDivergence = strToBool(this->parserCfg.at("minimumDivergence"));
            this->SAVETRAINMODEL = strToBool(this->parserCfg.at("SAVETRAINMODEL"));      //-- SAVE TRAINED MODEL?
            this->train = strToBool(this->parserCfg.at("train"));               //-- TRAIN TOTAL VARIABILITY, PROJECTION MATRIX, iVECTOR?
            this->enroll = strToBool(this->parserCfg.at("enroll"));

            this->test_size = std::stof(this->parserCfg.at("test_size"));          //-- BALANCE ENROLL/TEST
            this->scoringMethod = std::stoi(this->parserCfg.at("scoringMethod"));       // 1 - CSS, 0 - GPLDA
        }
    }
}


FileParserWrkSpace::FileParserWrkSpace()
{
    this->cfg_list.push_back("pTrainMfcc_path");
    this->cfg_list.push_back("pTrainList_file");
    this->cfg_list.push_back("pCfg_path");
    this->cfg_list.push_back("pToSaveWSpace");
    this->cfg_list.push_back("pTmp_file");
    this->cfg_list.push_back("pToLogs");

    this->cfg_list.push_back("pEnrollList_file");
    this->cfg_list.push_back("pEnrollMfcc_path");

}

void FileParserWrkSpace::parseConfigFile(char **argv, int ndx, std::string flag){

    if(!std::filesystem::exists(argv[ndx+1])){
        std::cerr << "File: "<<argv[ndx+1]<<" not found \n";
        throw;
    }else{
        std::fstream s{argv[ndx+1], s.in};
        if(!s.is_open()){
            std::cerr << " Problem with open the file: "<<argv[ndx+1]<<'\n';
            throw;
        }if(flag == "-w"){
            for(std::string line; std::getline(s, line);){
                for(std::vector<std::string>::iterator itr = this->cfg_list.begin();
                     itr != this->cfg_list.end();
                     ++itr){
                    if(!line.find('#'))
                        continue;
                    if(!line.find(*itr)){
                        // remove spaces in string line has read from file
                        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
                        // determine start of value
                        uint pos_start = line.find("=");
                        // determine end of value
                        uint pos_end = line.find(";");
                        std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);

                        this->parserCfg[*itr] = p1;
                    }
                }
            }
            // for(std::map<std::string, std::string>::iterator itr=parserCfg.begin();
            //       itr!=parserCfg.end();
            //       ++itr)
            //        std::cout<< itr->first <<" "<<itr->second <<'\n';

            //std::cout<<parserCfg.at("numOfComponents") <<'\n';

        this->pTrainList_file = this->parserCfg.at("pTrainList_file");
        this->pTrainMfcc_path = this->parserCfg.at("pTrainMfcc_path");
        this->pCfg_path = this->parserCfg.at("pCfg_path");
        this->pToSaveWSpace = this->parserCfg.at("pToSaveWSpace");
        this->pTmp_file = this->parserCfg.at("pTmp_file");
        this->pToLogs = this->parserCfg.at("pToLogs");
        this->pEnrollMfcc_path = this->parserCfg.at("pEnrollMfcc_path");
        this->pEnrollList_file = this->parserCfg.at("pEnrollList_file");
        }
    }
}

FileParserData::FileParserData()
{

    this->cfg_list.push_back("pToModel");
    this->cfg_list.push_back("pWork_path");

}

void FileParserData::parseConfigFile(char **argv, int ndx, std::string flag){

    if(!std::filesystem::exists(argv[ndx+1])){
        std::cerr << "File: "<<argv[ndx+1]<<" not found \n";
        throw;
    }else{
        std::fstream s{argv[ndx+1], s.in};
        if(!s.is_open()){
            std::cerr << " Problem with open the file: "<<argv[ndx+1]<<'\n';
            throw;
        }if(flag == "-d"){
            for(std::string line; std::getline(s, line);){
                if(!line.find('#'))
                    continue;
                if(!line.find("train")){
                    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                    line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
                    //line.erase(std::remove_if(line.begin(), line.end(), '"'),
                    // determine start of value
                    uint pos_start = line.find("=");
                    // determine end of value
                    uint pos_end = line.find(";");
                    std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);
                    //std::string p1 = this->extractValueFromLine(line);
                    this->list_train.push_back(p1);
                    continue;
                }
                if(!line.find("test")){
                    line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                    line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
                    // determine start of value
                    uint pos_start = line.find("=");
                    // determine end of value
                    uint pos_end = line.find(";");
                    std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);
                    //std::string p1 = this->extractValueFromLine(line);
                    this->list_test.push_back(p1);
                    continue;
                }
                for(std::vector<std::string>::iterator itr = this->cfg_list.begin();
                     itr != this->cfg_list.end();
                     ++itr){
                    if(!line.find(*itr)){
                        // remove spaces in string line has read from file
                        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
                        line.erase(std::remove(line.begin(), line.end(), '"'), line.end());
                        // determine start of value
                        uint pos_start = line.find("=");
                        // determine end of value
                        uint pos_end = line.find(";");
                        std::string p1 = line.substr(pos_start+1, pos_end-pos_start-1);
                        //std::cout << p1 << "\n";
                        this->parserCfg[*itr] = p1;
                    }
                }
            }

        this->pToModel = this->parserCfg.at("pToModel");
        this->pWork_path = this->parserCfg.at("pWork_path");

        // for(uint i=0; i<this->list_train.size(); ++i){
        //     std::cout <<this->list_train[i]<<'\n';
        // }
        // for(uint i=0; i<this->list_test.size(); ++i){
        //     std::cout <<this->list_test[i]<<'\n';
        // }
    }
}
}

std::vector<std::string> FileParserData::getList_train()
{
    return this->list_train;
}

std::vector<std::string> FileParserData::getList_test()
{
    return this->list_test;
}


