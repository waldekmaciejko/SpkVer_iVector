#include "timeheleprs.h"

/*!
 * \brief currentTime - format string containes date time
 * \param
 * \return
 */
namespace ava {

std::string currentTime(){

    std::time_t now = std::time(0);
    std::tm *ltm = std::localtime(&now);

    std::string t1 = std::to_string(ltm->tm_sec);
    std::string t2 = std::to_string(ltm->tm_min);
    std::string t3 = std::to_string(ltm->tm_hour);
    std::string t4 = std::to_string(ltm->tm_mday);
    std::string t5 = std::to_string(1+ltm->tm_mon);
    std::string t6 = std::to_string(1900+ltm->tm_year);

    //std::string t = t3+"h"+t2+"m"+t1+"s"+"_"+t4+"-"+t5+"-"+t6;
    std::string t = t4+"-"+t5+"-"+t6+"_"+t3+"h"+t2+"m"+t1+"s";
    return t;
}
}
