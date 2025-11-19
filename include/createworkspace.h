//
// Created by ai on 18.11.2025.
//

#ifndef CLION_TEST1_CREATEWORKSPACE_H
#define CLION_TEST1_CREATEWORKSPACE_H

#include <vector>
#include <string>
#include <filesystem>
#include <iostream>

class CreateWorkspace
{
public:
    CreateWorkspace(std::string path);
    void validateDirs();
    void validateSrtct();
    std::string getWorkspace();

private:
    std::string pathToWorkspace;
    std::string head = "WORKSPACE";
    std::vector<std::string> strct = {"calc", "cfg", "feat", "list", "logs"};
    std::vector<std::string> strct_feat = {"enroll", "train"};
};

#endif //CLION_TEST1_CREATEWORKSPACE_H
