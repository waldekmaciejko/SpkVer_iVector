#include "createworkspace.h"

CreateWorkspace::CreateWorkspace(std::string path){
    this->pathToWorkspace = path;
}

void CreateWorkspace::validateDirs(){


    if(this->pathToWorkspace.find("WORKSPACE")){
        this->pathToWorkspace += "/WORKSPACE";
        std::cout << "!!!!!!!!!!!!!" << this->pathToWorkspace  <<std::endl;
    }

    std::filesystem::path fsp = this->pathToWorkspace;


    // if path to workspace exists
    if (std::filesystem::exists(fsp)){

        std::cout << "Folder WORKSPACE exists\n";
        validateSrtct();

    }
    // if not create
    else{
        std::cerr << "WORKSPACE does not exists" << std::endl;

        // fsp += "/";
        // fsp += this->head;

        std::filesystem::create_directory(fsp);
        std::cerr << "WORKSPACE was created \n" << std::endl;
        validateSrtct();
    }
}

std::string CreateWorkspace::getWorkspace(){
    return this->pathToWorkspace;
}

void CreateWorkspace::validateSrtct(){

    for(std::size_t i=0;i<strct.size();i++){
        if(!(std::filesystem::exists(this->pathToWorkspace + "/" + this->strct[i]))) {
            std::filesystem::create_directory(this->pathToWorkspace + "/" + this->strct[i]);
        }
    }

    for(std::size_t i=0;i<strct_feat.size();i++){
        if(!(std::filesystem::exists(this->pathToWorkspace + "/feat/" + this->strct_feat[i]))) {
            std::filesystem::create_directory(this->pathToWorkspace + "/feat/" + this->strct_feat[i]);
        }
    }
}

