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

#include "funhelpers.h"

namespace ava{
//[deprecated]
//---------------- zestaw funkcji tworzacych liste plikow z podanej lokalizacji workSpace i zapisujacych te liste w pliku flt_fea_mfcc.lst
//---------------- brak dzialajacej funkccji bibliotecznej filesystem - dlatego to toporne rozwiazanie
int analizeDirectoryElement (const char *fpath,
                            const struct stat *sb,
                            int tflag,
                            struct FTW *ftwbuf) {

  if (tflag == FTW_F) {
    std::string strFileName(fpath);
    //cout<<strFileName<<endl;

    std::string workSpace="j:\\__Badania__\\speakerRecognition_ivector_cpp\\maw200\\maw200\\workspace_maw200\\features";
    //std::string workSpace= sciezka_workSpace(path)

    std::ofstream saveList;
    std::string ws= workSpace + "/list/flt_fea_mfcc.lst";
    saveList.open(ws, std::ios::app);
    saveList<<strFileName<<std::endl;
    saveList.close();
  }
  return 0;
  }


void WalkDirectoryTree (const char * pchFileName) {

  int nFlags = 0;

  if (nftw(pchFileName, analizeDirectoryElement, 20, nFlags) == -1) {
    perror("nftw");
  }
  }

//----------------------------- funkcja zmieniajaca \ na /

std::string winPathToPosixPathAnsi(std::string winSciezka){

    std::stringstream s;
    std::string winSciezkaString;
    s<<winSciezka;
    s>>winSciezkaString;

    for(unsigned int i=0; i<winSciezkaString.length(); i++){
                if(winSciezkaString[i] == '\\')
                winSciezkaString[i]='/';
    }
        return winSciezkaString.c_str();

    }

//-----------------------------  zamienia sciezke POSIX na WinNT
std::string winPathToPosixPath(const char* winSciezka){

    std::stringstream s;
    std::string winSciezkaString;
    s<<winSciezka;
    s>>winSciezkaString;

    for(unsigned int i=0; i<winSciezkaString.length(); i++){
                if(winSciezkaString[i] == '\\')
                winSciezkaString[i]='/';
    }
        return winSciezkaString;
    }


//----------------------------- funkcja badajaca rozmiar pliku

int get_file_size(std::string filename){

    // argumenty
    // filename - sciezka do pliku

    // return
    // size - rozmiar pliku w bajtach


    FILE *p_file = NULL;
    p_file = fopen(filename.c_str(), "rb");
    fseek(p_file,0,SEEK_END);
    unsigned int  size = ftell(p_file);
    fclose(p_file);
    return size;
    }

// funkcja zapisujaca do pliku macierz Eigen
void save_egien_matrix_to_txt_file(std::string workSpace,
                                                       std::string flt_tmp,
                                                       std::string nazwa,
                                                       Eigen::MatrixXf& matrix){

    std::string p =workSpace+flt_tmp+nazwa;
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
    std::ofstream file(p.c_str());
    file<<matrix.format(CSVFormat);
    file.close();
    }

// funkcja zapisujaca do pliku wektor Eigen
void save_egien_matrix_to_txt_file(std::string workSpace,
                                                       std::string flt_tmp,
                                                       std::string nazwa,
                                                       Eigen::VectorXf& matrix){

    std::string p =workSpace+flt_tmp+nazwa;
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
    std::ofstream file(p.c_str());
    file<<matrix.format(CSVFormat);
    file.close();
    }

/*! [deprecated]
 * \brief funkTotalVariability
 * \param numTdim
 * \param sigma
 * \param numSpkUBM
 * \param repSigma
 * \param N
 * \param prodNumOfComponentsNumFeatures
 * \param numFeatures
 * \param numOfComponents
 * \param F
 * \param Nc
 * \param numIterationsTV
 * \return
 */
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
                               unsigned int numIterationsTV){

    arma::mat T=arma::randn(sigma.n_elem, numTdim);
    T=T/arma::norm(T);

    arma::mat I = arma::eye(numTdim, numTdim);

    arma::cube Ey(numTdim, 1, numSpkUBM);
    arma::cube Eyy(numTdim, numTdim, numSpkUBM);
    arma::cube Linv(numTdim, numTdim, numSpkUBM);
    arma::mat L(numTdim, numTdim);

    //unsigned int numIteratio/repnsTV = 20;   //optymalna wartosc wynosi 20

    arma::mat TtimeInvSSDiag = {};


    for(unsigned int i=0; i<numIterationsTV; i++)
    {

        //1.---- oblicz dystrybucje posteriori zmiennych ukrytych
            TtimeInvSSDiag =T/repSigma;
            for(unsigned int s=0; s<numSpkUBM; s++)
            {
                L = I + TtimeInvSSDiag.t()%arma::repmat((N.slice(s)).t(), 1, numTdim).t()*T;
                Linv.slice(s) =arma::pinv(L);
                Ey.slice(s) = Linv.slice(s) * TtimeInvSSDiag.t()*F.slice(s).t();
                //mat LL= Linv.slice(s) * TtimeInvSSDiag.t()*F.slice(s).t();
                Eyy.slice(s) = Linv(s) + Ey.slice(s)*Ey.slice(s).t();
            }

        //2.---- akumulacja statystyk w dla kazdego mowcy
        arma::mat Eymat(numTdim, numSpkUBM);
        arma::mat FFmat(prodNumOfComponentsNumFeatures, numSpkUBM);
        arma::mat Kt(prodNumOfComponentsNumFeatures, numTdim);
        arma::cube K(numFeatures, numTdim, numOfComponents);
        arma::mat AcLocal(numTdim, numTdim);
        arma::cube newT(numFeatures, numTdim, numOfComponents);

        // splaszczenie obiektow cube z wymiarem 1 do macierzy
        for(unsigned int s=0; s<numSpkUBM; s++)
        {
            Eymat.col(s)=(Ey.slice(s)).rows(0,int(numTdim-1));
            FFmat.col(s)=(F.slice(s)).cols(0,int(prodNumOfComponentsNumFeatures-1)).t();
        }

        Kt = FFmat*Eymat.t();

        for(unsigned int ii=0; ii<numOfComponents; ii++)
        {
            //Kt[prodNumOfComponentsNumFeatures, numTdim]->K[numFeatures, numTdim, numOfComponents]
            K.slice(ii)=Kt.rows(int(ii*numFeatures), int(ii*numFeatures+numFeatures-1));
        }

        for(unsigned int c=0; c<numOfComponents; c++)
        {
            AcLocal.zeros(numTdim, numTdim);
            for(unsigned int s=0; s<numSpkUBM; s++)
            {
                AcLocal = AcLocal + (float(Nc.slice(s)(0,c))*Eyy.slice(s));
            }

       //3.---- aktualizacjaprzestrzeni TV
           newT.slice(c) =  (arma::pinv(AcLocal)*(K.slice(c).t())).t();
        }

        arma::mat tmp={};

        for(unsigned int cc=0; cc<numOfComponents; cc++)
        {
            tmp = arma::join_cols(tmp, newT.slice(cc));
        }
        //T.save("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\T0.csv", arma::csv_ascii);
        T=tmp;
        //T.save("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\T1.csv", arma::csv_ascii);

      }
  return T;
  }

/*!
 * \brief funkTotalVariabilityVector - estimate total variability matrix (each audio as seperate spk)
 * \param numTdim
 * \param sigma
 * \param numSpkUBM
 * \param repSigma
 * \param NN
 * \param prodNumOfComponentsNumFeatures
 * \param numFeatures
 * \param numOfComponents
 * \param FF
 * \param NNcc
 * \param numIterationsTV
 * \return -[numTdim, number of files x gmm order]
 */

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
                                       unsigned int numIterationsTV,
                                       ava::Log4AVA& logger){

    arma::mat T=arma::randn(sigma.n_elem, numTdim);
    T=T/arma::norm(T);    

    arma::mat I = arma::eye(numTdim, numTdim);

    std::vector<arma::mat> Ey{};
    Ey.reserve(numFilesUBM);

    std::vector<arma::mat> Eyy{};
    Eyy.reserve(numFilesUBM);

    std::vector<arma::mat> Linv{};
    Linv.reserve(numFilesUBM);

    arma::mat L(numTdim, numTdim);

    //unsigned int numIteratio/repnsTV = 20;   //optymalna wartosc wynosi 20

    arma::mat TtimeInvSSDiag = {};

    std::vector<arma::mat> K{};
    K.reserve(numOfComponents);

    std::vector<arma::mat> newT{};
    newT.reserve(numOfComponents);

    arma::mat tmp={};

    // Test 4 - OK
    logger.save("Total Variability is training: ", ava::currentTime());

    for(unsigned int i=0; i<numIterationsTV; i++)
    {
        logger.save("  iteration " + std::to_string(i + 1),
                    " of " + std::to_string(numIterationsTV) + " at time: " + ava::currentTime());

        //1.---- estimate posteriori density of hidden variables
        //       (oblicz dystrybucje posteriori zmiennych ukrytych)
            TtimeInvSSDiag =T/repSigma;
            Linv.clear();
            Ey.clear();
            Eyy.clear();
            for(unsigned int s=0; s<numFilesUBM; s++){
                L = I + TtimeInvSSDiag.t()%arma::repmat((NN.at(s)).t(), 1, numTdim).t()*T;
                Linv.push_back(arma::pinv(L));
                Ey.push_back(Linv.at(s) * TtimeInvSSDiag.t()*FF.at(s).t());
                Eyy.push_back(Linv.at(s) + Ey.at(s)*Ey.at(s).t());
            }
        // Test 4 - END

        //2.---- accumulate statistics for each spk
        arma::mat Eymat(numTdim, numFilesUBM);
        arma::mat FFmat(prodNumOfComponentsNumFeatures, numFilesUBM);
        arma::mat Kt(prodNumOfComponentsNumFeatures, numTdim);
        arma::mat AcLocal(numTdim, numTdim);

        //arma::cube K(numFeatures, numTdim, numOfComponents);
        K.clear();

        //arma::cube newT(numFeatures, numTdim, numOfComponents);
        newT.clear();

        // splaszczenie obiektow cube z wymiarem 1 do macierzy
        for(unsigned int s=0; s<numFilesUBM; s++){
            Eymat.col(s)=(Ey.at(s)).rows(0,int(numTdim-1));
            FFmat.col(s)=(FF.at(s)).cols(0,int(prodNumOfComponentsNumFeatures-1)).t();
        }

        Kt = FFmat*Eymat.t();

        for(unsigned int ii=0; ii<numOfComponents; ii++) {
            K.push_back(Kt.rows(int(ii*numFeatures), int(ii*numFeatures+numFeatures-1)));
        }

        for(unsigned int c=0; c<numOfComponents; c++){
            AcLocal.zeros(numTdim, numTdim);
            for(unsigned int s=0; s<numFilesUBM; s++){
                AcLocal = AcLocal + (float(NNcc.at(s)(0,c))*Eyy.at(s));
            }
        //3.---- update total variability matrix (aktualizacjaprzestrzeni TV)
            newT.push_back((arma::pinv(AcLocal)*(K.at(c).t())).t());            
        }

        tmp.clear();

        for(unsigned int cc=0; cc<numOfComponents; cc++){
            //tmp = arma::join_cols(tmp, newT.slice(cc));
            tmp = arma::join_cols(tmp, newT.at(cc));
        }
        T=tmp;
       }

    //std::cout<<"T"<<std::endl;
    //T.print();
    return T;
  }


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
                                                    arma::mat& TSi){

    //arma::cube ivectorPerSpeaker
    //arma::mat TS =T/repSigma;
    //arma::mat TSi = T.t();
    //arma::mat ubmMu = ubmModel.means;
    arma::mat Yivector={};
    unsigned int liczbaRamekYivector {}; //= Y.n_rows;
    arma::mat logLikelihoodSum = {};
    arma::mat amax={};
    arma::mat gamma={};

    arma::mat I=arma::eye(numTdim, numTdim);
    arma::mat tmp6{};
    std::string m={};
    arma::mat ivectorPerSpk{};
    std::map<std::string, arma::mat> ivectorPerFile{};
    arma::mat logLikelihoodivector {};
    arma::mat n={};
    arma::mat f={};


    for(auto st : setSpk)
    {
        //std::cout<<"setSpk"<<st<<std::endl;
        for(auto para : multimapMowcyMFC)
        {
            //arma::mat YY{};
            if(para.first==st){

                ava::MfcceFeatures mfcefivector(para.second);
                liczbaRamekYivector = mfcefivector.returnArmaMat().n_rows;
                Yivector =(ava::normZ(mfcefivector.returnArmaMat().t(), normMean, normStd, liczbaRamekYivector, numFeatures)).t();
                //logLikelihoodivector = ubmModel.log_p(Yivector.t());

                for(unsigned int comp=0; comp<numOfComponents; comp++) {
                    logLikelihoodivector=arma::join_cols(logLikelihoodivector, ubmModel.log_p(Yivector, comp));
                }

                amax = max(logLikelihoodivector, 0);
                logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihoodivector-arma::repmat(amax, numOfComponents, 1))));
                gamma = arma::exp(logLikelihoodivector-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

                n = arma::sum(gamma, 0);
                //f = arma::operator*(Yivector, gamma) - arma::repmat((n*ubmMu.t()).t(), 1, numOfComponents);
                 f = arma::operator*(Yivector, gamma)- arma::repmat(n, numFeatures, 1)%ubmMu;

                tmp6 = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);

                ivectorPerSpk=arma::join_rows(ivectorPerSpk, tmp6);
                ivectorsTrain = arma::join_rows(ivectorsTrain, tmp6);
            }
            logLikelihoodivector.clear();
        }
        ivectorPerFile.insert(std::make_pair(st, ivectorPerSpk));        
        ivectorPerSpk.clear();
    }
    return ivectorPerFile;
    }

/*!
 * \brief takeSpkLabel
 *                      -   function to find label of spk in file name
 * \param pathToWav     -   path to file
 * \return              -   label of spk
 */
std::string takeSpkLabel(std::string pathToFile){

    std::regex reg1("[F|M]{1}[0-9]+");
    std::smatch m;
    bool found = regex_search(pathToFile,
                              m,
                              reg1);
    return m.str();
    }

/*
 *
 * normalizacja mean-std macierzy Y
 * Y[lRamekY, numFea] - [liczba ramek, liczba cech w ramce]
 *
 */

/*
arma::mat normZ(arma::mat Y,
                  arma::mat normMean,
                  arma::mat normStd,
                  unsigned int liczbaRamekY,
                  unsigned int numFeatures){

    arma::mat tmpY=(Y-arma::repelem(normMean, liczbaRamekY, 1));
    arma::mat Vec=arma::pow(normStd, -1);
    arma::mat YbyVec(liczbaRamekY, numFeatures);

    if(Y.n_cols!=Vec.n_cols)
        std::cerr<<"Error fail dimmensions"<<std::endl;
    else{
        for(unsigned int i = 0; i<tmpY.n_rows; i++)
        {
            YbyVec.row(i)=tmpY.row(i)%Vec;
        }
    }
    return YbyVec;
}
*/

/*!
 * \brief normZ
 * \param Y [nFeatures, nFrames]
 * \param normMean
 * \param normStd
 * \param liczbaRamekY
 * \param numFeatures
 * \return
 *
 * TODO - do cepstral normalization
 */
arma::mat normZ(arma::mat Y,
                  arma::mat normMean,
                  arma::mat normStd,
                  unsigned int liczbaRamekY,
                  unsigned int numFeatures){

    //arma::mat tmpY=arma::repelem(normMean, liczbaRamekY, 1);
    arma::mat tmpY=(Y-arma::repelem(normMean.t(), 1, liczbaRamekY));

    tmpY=tmpY.t();

    arma::mat Vec=arma::pow(normStd, -1);
    arma::mat YbyVec(liczbaRamekY, numFeatures);

    for(unsigned int i = 0; i<tmpY.n_rows; i++)
    {
        YbyVec.row(i)=tmpY.row(i)%Vec;
    }

    return YbyVec;
}





arma::mat multipicColByCol(arma::mat& A, arma::mat& B){
    uint rA=A.n_rows;
    uint cA=A.n_cols;

    uint rB=B.n_rows;
    uint cB=B.n_cols;

    uint flaga=1;

    arma::mat res(rA,cA);

    if (rA!=rB && rA!=cB){
        //std::cout<<"N_rows of A needs to be egual n_rows B"<<std::endl;
        flaga=0;
    }
    if (cB!=1){
        //std::cout<<"N_col of A needs to be egual n_rows B"<<std::endl;
        flaga=0;
    }
    if (rB==1){
        B=B.t();
        flaga=1;
    }

    if(flaga==1){
        for(uint i = 0; i<cA; i++){
            //std::cout<<A.col(i)%B.t()<<std::endl;
            res.col(i)=A.col(i)%B;
        }
    }else throw "false dimension of matrix";
    return res;

    }

arma::mat multipicColByCol(arma::mat& A, arma::mat&& B){
    uint rA=A.n_rows;
    uint cA=A.n_cols;

    uint rB=B.n_rows;
    uint cB=B.n_cols;

    uint flaga=1;

    arma::mat res(rA,cA);

    if (rA!=rB && rA!=cB){
        //std::cout<<"N_rows of A needs to be egual n_rows B"<<std::endl;
        flaga=0;
    }
    if (cB!=1){
        //std::cout<<"N_col of A needs to be egual n_rows B"<<std::endl;
        flaga=0;
    }
    if (rB==1){
        B=B.t();
        flaga=1;
    }

    if(flaga==1){
        for(uint i = 0; i<cA; i++){
            //std::cout<<A.col(i)%B.t()<<std::endl;
            res.col(i)=A.col(i)%B;
        }
    }else throw "false dimension of matrix";
    return res;
    }

void divdeMatrixRBRVectorNorm(arma::mat& matrixToDivide,
                              arma::mat& matrixAfterDivide){

    if(arma::size(matrixToDivide) != arma::size(matrixAfterDivide)){
        std::cout<<"Error: fail matrix dimmensions"<<std::endl;}
    else{
        arma::mat vecAfterNorm(1,matrixToDivide.n_cols, arma::fill::randu);
        normArmaMatrix(matrixToDivide, vecAfterNorm);
        for(uint i=0; i<matrixToDivide.n_rows; i++){
            matrixAfterDivide.row(i)=matrixToDivide.row(i)/vecAfterNorm;
        }
    }

    //}

}
/*!
 * \brief takeMostSignificantEigenVecLDA
 * \param ivectorPerFile - std::map<std::string, arma::mat> consists label of spk as a key
 *                          and arma::mat [numTdim, number of files of spk]
 * \param numTdim
 * \param ivectorsTrain
 * \param numEigenvectors
 * \param utterancePerSpeaker
 * \return
 */
arma::mat takeMostSignificantEigenVecLDA(std::map<std::string, arma::mat>& ivectorPerSpk,
                                              unsigned int numTdim,
                                              arma::mat& ivectorsTrain,
                                              uint numEigenvectors,
                                              std::vector<uint>& utterancePerSpeaker){
    //std::vector<uint> utterancePerSpeaker{};
    arma::mat Sw = arma::zeros(numTdim, numTdim); // between variability
    arma::mat Sb = arma::zeros(numTdim, numTdim); // within variability
    arma::mat ws={};
    arma::mat wsbar={};

    //mean vaule of each numTdim via all files
    arma::mat wbar = arma::mean(ivectorsTrain, 1);

    for(const auto& x : ivectorPerSpk){
        utterancePerSpeaker.push_back(x.second.n_cols);
        ws = x.second;
        wsbar = arma::mean(ws, 1);
        Sb = Sb + (wsbar - wbar)*((wsbar - wbar).t());
        Sw = Sw + arma::cov(ws.t());
    }

    //Sb.save("j:\\__Badania__\\QtGUI_ava300\\Sb.csv", arma::csv_ascii);
    //Sw.save("j:\\__Badania__\\QtGUI_ava300\\Sw.csv", arma::csv_ascii);

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_pair(eigval, eigvec, Sb, Sw);
    arma::mat A = {};

    arma::mat eigval_real = arma::real(eigval);
    arma::mat eigvec_real = arma::real(eigvec);

    eigval_real.replace(arma::datum::inf, 0);
    eigvec_real.replace(arma::datum::inf, 0);

    //arma::uvec ind = arma::sort_index(eigval_real, "descend");
    //numEigenvectors=numEigenvectors-1;

    //for(uint  i : ind(arma::span(0,numEigenvectors))){
    //    A=arma::join_cols(A, eigvec_real.row(i));
    //}
    // arma::inplace_trans(A);

    numEigenvectors = numEigenvectors - 1;

    A = eigvec_real.cols(0, numEigenvectors);
    //arma::inplace_trans(A);

    arma::mat Avecnorm={};
    arma::mat tmp={};
    numEigenvectors = numEigenvectors + 1;

    for(uint i=0; i<numEigenvectors; i++){
        tmp=arma::norm(A.col(i));
        Avecnorm= arma::join_cols(Avecnorm, tmp);
    }

    arma::mat Anorm(numTdim, numEigenvectors, arma::fill::zeros);

    for(uint ii=0; ii<numTdim; ii++){
        for(uint iii=0; iii<numEigenvectors; iii++){
            Anorm(ii, iii) = (A(ii, iii)/Avecnorm(iii,0));
        }
    }

    return Anorm;
    }

/*!
 * \brief takeMostSignificantEigenVecLDA2 - new version of takeMostSignificantEigenVecLDA
 *                  instead of armadillo c++ eigen_pair use eig_gen after the
 *                  prouct of iverted Sw matrix and Sb matrix.
 *
 * \param ivectorPerSpk
 * \param numTdim
 * \param ivectorsTrain
 * \param numEigenvectors
 * \param utterancePerSpeaker
 * \return
 */


arma::mat takeMostSignificantEigenVecLDA2(std::map<std::string, arma::mat>& ivectorPerSpk,
                                              unsigned int numTdim,
                                              arma::mat& ivectorsTrain,
                                              uint numEigenvectors,
                                              std::vector<uint>& utterancePerSpeaker){
    //std::vector<uint> utterancePerSpeaker{};
    arma::mat Sw = arma::zeros(numTdim, numTdim); // between variability
    arma::mat Sb = arma::zeros(numTdim, numTdim); // within variability
    arma::mat ws={};
    arma::mat wsbar={};

    //mean vaule of each numTdim via all files
    arma::mat wbar = arma::mean(ivectorsTrain, 1);

    for(const auto& x : ivectorPerSpk){
        utterancePerSpeaker.push_back(x.second.n_cols);
        ws = x.second;
        wsbar = arma::mean(ws, 1);
        Sb = Sb + (wsbar - wbar)*((wsbar - wbar).t());
        Sw = Sw + arma::cov(ws.t());
    }

    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::mat X1 = arma::inv(Sw);
    arma::mat X2 = X1*Sb;
    arma::eig_gen(eigval, eigvec, X2);
    arma::mat A = {};

    arma::mat eigval_real = arma::real(eigval);
    arma::mat eigvec_real = arma::real(eigvec);

    eigval_real.replace(arma::datum::inf, 0);
    eigvec_real.replace(arma::datum::inf, 0);

    numEigenvectors=numEigenvectors-1;

    A = eigvec_real.cols(0, numEigenvectors);
/*
    for(uint  i : ind(arma::span(0,numEigenvectors))){
        A=arma::join_cols(A, eigvec_real.row(i));
    }*/
    arma::inplace_trans(A);
    arma::mat Avecnorm={};
    arma::mat tmp={};

    numEigenvectors=numEigenvectors+1;

    for(uint i=0; i<numEigenvectors; i++){
        tmp=arma::norm(A.col(i));
        Avecnorm= arma::join_cols(Avecnorm, tmp);
    }

    arma::mat Anorm(numTdim, numEigenvectors, arma::fill::zeros);

    for(uint ii=0; ii<numTdim; ii++){
        for(uint iii=0; iii<numEigenvectors; iii++){
            Anorm(ii, iii) = (A(ii, iii)/Avecnorm(iii,0));
        }
    }
    return Anorm;
    }

void performaLDAfunk(uint numEigenvectors,
                     bool performLDA,
                     uint numTdim,
                     std::map<std::string, arma::mat>& ivectorPerSpk,
                     arma::mat& ivectorsTrain,
                     arma::mat& projectionMatrix,
                     std::vector<arma::mat>& w,
                     std::vector<uint>& utterancePerSpeaker){

    arma::mat A = {};

    if(performLDA){

        A = ava::takeMostSignificantEigenVecLDA(ivectorPerSpk,
                                                   numTdim,
                                                   ivectorsTrain,
                                                   numEigenvectors,
                                                   utterancePerSpeaker);

        ivectorsTrain =A.t() * ivectorsTrain;
        projectionMatrix = A.t() * projectionMatrix;
        uint tmp=0;
        arma::mat tmp3={};

        for(auto i : utterancePerSpeaker){
            i=tmp+i-1;
            tmp3=ivectorsTrain.cols(tmp,i);
            w.push_back(tmp3);
            tmp=i+1;
        }
    }
    }


/*!
 * \brief performaWCCNfunk - checked in test9
 * \param numEigenvectors
 * \param projectionMatrix
 * \param w
 */
void performaWCCNfunk(uint numEigenvectors,
                      arma::mat& projectionMatrix,
                      std::vector<arma::mat>& w){

    float alpha = 0.9;
    arma::mat B={};
    arma::mat W(numEigenvectors, numEigenvectors, arma::fill::zeros);

    for(auto ww : w){
        W = W + arma::cov(ww.t());
    }
    W=W/int(w.size());
    W=(1-alpha) * W + alpha * arma::eye(numEigenvectors, numEigenvectors);
    B=arma::chol(arma::pinv(W), "lower"); // Cholesky decomposition
    projectionMatrix = B * projectionMatrix;
}

void performWhitening(std::map<std::string, arma::mat>& ivectorPerSpk,
                      arma::mat& ivectorMatrix,
                      arma::mat& projectionMatrix,
                      std::string zmiennWhitening,
                      arma::mat& mu,
                      arma::mat& W){

    std::vector<arma::mat> ivector{};
    //uint tmp= 0;

    for(auto& i : ivectorPerSpk){
        ivector.push_back(projectionMatrix*i.second);
    }

    //uint numEigenVoices = 16;
    //uint K = numSpkUBM;
    //uint D = numEigenvectors;

    //arma::mat ivectorMatrix{};

    for(auto& ii : ivector){
        ivectorMatrix = arma::join_rows(ivectorMatrix, ii);
    }

    uint NN = ivectorMatrix.n_rows;
    //arma::mat mu =arma::mean(ivectorMatrix, 1);
    mu =arma::mean(ivectorMatrix, 1);

    for(uint ii=0; ii<ivectorMatrix.n_cols; ii++){
        ivectorMatrix.col(ii)=ivectorMatrix.col(ii)-mu;
    }

    //bool ZCA = false;
    //bool PCA = false;

    //arma::mat W{};

    if(zmiennWhitening=="ZCA"){
        arma::mat S{};
        arma::mat U_tmp{};
        arma::vec s_tmp{};
        arma::mat sD{};
        arma::mat sV{};

        S = arma::cov(ivectorMatrix.t());

        /*
         * SVD - Singular Value Decomposition
         * svd( mat U, vec s, mat V, mat X )
         * X = U*diagmat(s)*V.t()
         *
         */

        arma::svd(U_tmp, s_tmp, sV, S);
        //sD = arma::diagmat(s_tmp);
        W = arma::diagmat(1/(arma::sqrt(s_tmp) + arma::datum::eps))*sV.t();
        ivectorMatrix = W * ivectorMatrix;

    } else if(zmiennWhitening=="PCA"){
          arma::mat S{};
          arma::vec eigval;
          arma::mat eigvec;

          S = arma::cov(ivectorMatrix.t());
          arma::eig_sym(eigval, eigvec, S);

          // eigval.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\workspace_maw200\\tmp\\eigval.csv", arma::csv_ascii);

          W = arma::diagmat(1/(arma::sqrt(eigval) + arma::datum::eps))*eigvec.t();
          ivectorMatrix = W * ivectorMatrix;

      } else{
        W = arma::eye(ivectorMatrix.n_rows, ivectorMatrix.n_rows);
    }

}

void normArmaMatrix(arma::mat& matToNorm,
                    arma::mat& vecAfterNorm){

    for(uint i=0; i<matToNorm.n_cols; i++)
    {
        vecAfterNorm.col(i)=arma::norm(matToNorm.col(i));
    }

}

/*!
 * \brief divideVectorIntoSpeakerCells
 * \param [in] utterancePerSpeaker - vector with numbers of PCM files for each spk
 * \param [in] ivectorMatrixAfterNorm
 * \param [out] ivector
 */

void divideVectorIntoSpeakerCells(std::vector<uint>& utterancePerSpeaker,
                                  arma::mat& ivectorMatrixAfterNorm,
                                  std::vector<arma::mat>& ivector,
                                  uint numEigenVectors){
    uint r=0;
    uint k=0;
    arma::mat tmp2{};

    for(uint i=0; i<utterancePerSpeaker.size(); i++){
        k=r+utterancePerSpeaker[i]-1;
        tmp2 = ivectorMatrixAfterNorm(arma::span(0,numEigenVectors-1), arma::span(r,k));
        r=utterancePerSpeaker[i]+r;
        ivector.push_back(tmp2);
    }
}
/*!
 * \brief sortIVectors
 * \param ivectorsSorted
 * \param ff
 * \param utterancePerSpeaker       -   vector contains
 *                                  number o PCM files per spk
 *                                  in coropora for UBM and T
 * \param ivector
 */

void sortIVectors(std::vector<std::vector<arma::mat>>& ivectorsSorted,
                  arma::mat& ff,
                  std::vector<uint>& utterancePerSpeaker,
                  std::vector<arma::mat>& ivector){

    // convert std::vector to arma::mat and check unique values
    // arma::umat uniquePerSpk  - matrix contains
    //                            number o PCM files per spk
    //                            in coropora for UBM and T
    //                            (= utterancePerSpeaker)
    arma::umat uniquePerSpk = arma::conv_to<arma::umat>::from(utterancePerSpeaker);

    // matrix contains unique numbers of files
    // per spk for UBM and T
    arma::umat uniqueLengths = arma::unique(uniquePerSpk);

    // number of unique numbers of files
    uint numUniqueLengths = uniqueLengths.n_rows;

    arma::uvec idx{};
    std::vector<arma::mat> temp{};
    uint speakerIdx = 0;
    uint tmp100{};
    arma::mat rho{};

    for(uint i=0; i<numUniqueLengths; i++){

        idx=arma::find(uniquePerSpk == uniqueLengths(i,0));
        arma::mat ttt{};

        temp.clear();
        for(uint speakerIdxWithinUniqueLength=0;
            speakerIdxWithinUniqueLength<idx.n_rows;
            speakerIdxWithinUniqueLength++){

                tmp100=idx(speakerIdxWithinUniqueLength);
                rho=ivector[tmp100];
                temp.push_back(rho);             
                ff.col(speakerIdx)=arma::sum(rho, 1);
                speakerIdx=speakerIdx+1;
        }
        ivectorsSorted.push_back(temp);
    }
}

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
                arma::mat& S){

    // Train G-PLDA model using the EM algorithm described in [13]

    //uint numEigenVoices = 16;
    //arma::mat V(D, numEigenVoices, arma::fill::randn);
    //arma::mat lambda =arma::pinv(S/numFilesUBM);
    //uint numIter = 5;
    //bool minimumDivergence = true;

    arma::umat uniquePerSpk = arma::conv_to<arma::umat>::from(utterancePerSpeaker);
    arma::umat uniqueLengths = arma::unique(uniquePerSpk);
    uint numUniqueLengths = uniqueLengths.n_rows;

    arma::mat gamma(numEigenVoices, numEigenVoices, arma::fill::zeros);
    arma::mat EyTotal(numEigenVoices, K, arma::fill::zeros);
    arma::mat R(numEigenVoices, numEigenVoices, arma::fill::zeros);

    uint ivectroLength{};

    std::vector<arma::mat> iv{};

    arma::mat M{};

    for(uint iter=0; iter<numIter; iter++){
        //---arma::mat gamma(numEigenVoices, numEigenVoices, arma::fill::zeros);
        //---arma::mat EyTotal(numEigenVoices, K, arma::fill::zeros);
        //---arma::mat R(numEigenVoices, numEigenVoices, arma::fill::zeros);
        gamma.zeros();
        EyTotal.zeros();
        R.zeros();

        uint idx=0;
        for(uint lengthIndex=0; lengthIndex<numUniqueLengths; lengthIndex++){
            ivectroLength = uniqueLengths(lengthIndex);

            // izoluj wektory w tym samym rozmiarze
            //---std::vector<arma::mat> iv =ivectorSorted[lengthIndex];
            iv =ivectorSorted[lengthIndex];

            // oblicz wartosc M
            // Equation (A.7) in Unifying probabilistic linear discriminant
            // analysis variants in biometric authentication
            // M matrix should be symmetric
            M = arma::pinv(ivectroLength * (V.t()*(lambda*V)) + (arma::eye(numEigenVoices, numEigenVoices)));
            for(uint speakerIndex=0; speakerIndex<iv.size(); speakerIndex++){
                // First moment of latent variable for V
                // Equation (A.8) in [13]
                arma::mat Ey = M*V.t()*lambda*ff.col(idx);

                // Calculate second moment
                arma::mat Eyy = Ey*Ey.t();

                // Update Ryy
                // Equation (A.13) in [13]
                R = R + ivectroLength*(M + Eyy);

                // Append EyTotal
                EyTotal.col(idx)=Ey;
                idx = idx + 1;

                // If using minimum divergence, update gamma
                // Equation (A.18) in [13]
                if(minimumDivergence){                    
                    gamma = gamma + (M+Eyy);                   
                }
            }
        }
        // Calculate T
        // Equation (A.12) in [13]
        arma::mat TT = EyTotal *ff.t();

        // MAXIMIZATION
        // Equation (A.16) in [13]
        V = TT.t()*arma::pinv(R);
        lambda =arma::pinv((S-V*TT)/numFilesUBM);

        arma::mat tmp{};

        if(minimumDivergence){
            // Equation (A.18) in [13]
             gamma = gamma/K;
             // Equation (A.22) in [13]            
             V = V*arma::chol(gamma, "lower");
        }

    }

}

/*!
 * \brief readMFCCHCopyFile -    dunction that read HTK MFCC file
 * \param path_to_mfc       -    std::string path to HTK MFCC file
 * \return arma::mat [numOfFreames, numOfFeatures]
 */
arma::mat readMFCCHCopyFile(std::string path_to_mfc){

    int smap_period={}; // okres próbek (przesuniecie3) - 32 bity
    short int s_size={};    //rozmiar probki - 16 bitow
    short int htk_code={};   // rodzaj probek - 16 bitow
    float p1={};        // wartosc probki

    // dla macierzy typu Eigen
    //Eigen::MatrixXf m = Eigen::MatrixXf(2,2); // macierz danych MFCC, ktora pozniej dynamicznie zmieni rozmiar
    //Eigen::VectorXf v_e = Eigen::VectorXf(2); // wektor energii pozniej dynamicznie zmieni rozmiar

    arma::mat mm={};
    arma::mat vv_e={};
    unsigned int l_col={}; // liczba kolumn w macierzy cech (liczba cech)
    unsigned int n_o_f={};       // liczba ramek (wierszy w macierzy cech) - 32 bity

    path_to_mfc=winPathToPosixPath(path_to_mfc.c_str());
    std::ifstream fin(path_to_mfc, std::ios::in | std::ios::binary);

    //-------------- odczyt naglowka

    fin.read((char*) &(n_o_f), sizeof(int));
    fin.read((char*) &(smap_period), sizeof(int));
    fin.read((char*) &(s_size), sizeof(short int));
    fin.read((char*) &(htk_code), sizeof(short int));

    //-------------- okreslenie rozmiaru pliku

    unsigned int l_baj_in_file=get_file_size(path_to_mfc);

    unsigned int l_cech_w_pliku = (l_baj_in_file-12)/4; // liczba cech w pliku to iloczyn liczby kolumn i liczby wierszy
                                                        // (l_baj_in_file - 12 (naglowek)) podzilieć przez 4 bo 4 bajt yto jedna cecha

    l_col=l_cech_w_pliku/n_o_f;                   // liczba cech (kolumn w wektorze) w pliku

    //Eigen::MatrixXf mm(n_o_f, l_col);                 // dla macierzy typu Eigen
    //this->m.resize(n_o_f, l_col);

    mm.resize(n_o_f, l_col);                      // alokacja dla armadillo
    vv_e.resize(n_o_f, 1);

    for(unsigned int k=0; k<n_o_f; k++)                 // zapiesz do macierzy
    {
        for(unsigned  int w=0; w<l_col; w++)
        {
            fin.read((char*) &p1, sizeof(float));
        //    this->m(k,w) = p1; // dla macierzy typu Eigen
            mm(k,w)= p1;
        }
    }

    unsigned int id_e = l_col-1;                            // ustalenie indexu energii l_col-1
                                                            // poniewaz ostatni kolumna to E

    //this->v_e.resize(n_o_f);                              // dla macierzy typu Eigen
    //this->v_e=this->m.col(id_e);                          // dla macierzy typu Eigen

    vv_e.resize(n_o_f);                               // dla macierzy typu Armadillo - arma
    vv_e=mm.col(id_e);

    // uwaga sprawdzic wartosci v_e porownojac z danymi z HList
    return mm;
}

/*!
 * \brief readRecursivlyMFCCHcopyFile - read in loop htk mfcc files
 * \param pWork_path
 * \param pFeaMfcc_path
 * \param [in] normMean - mean matrix (vector) [1, nFeatures]
 * \param [out] normStd - st. deviation (vector) [1, nFeatures]
 * \param [in/out] normYtrans - retun concatanated matrix [nFrames, nFeatures]
 */
void readRecursivlyMFCCHcopyFile(std::string pWork_path,
                                      std::string pFeaMfcc_path,
                                      arma::mat& normMean,
                                      arma::mat& normStd,
                                      arma::mat& normYtrans){

    arma::mat Y{};
    std::string p{};

    QDirIterator it_mfcc(QString::fromStdString(pWork_path + "/" + pFeaMfcc_path),
                         QStringList()<<"*.mfc", QDir::Files,
                         QDirIterator::Subdirectories);

    // concatenation
    while(it_mfcc.hasNext())
    {
        //qInfo()<<it_mfcc.next()<<endl;
        //ava::MfcceFeatures mfcef(it_mfcc.next().toStdString());
        //Y=arma::join_cols(Y, mfcef.returnArmaMat());

        //[nFrames, nFeatures]
        //nFrames - number of hops
        //nFeatures - number of cepstral coefficients
        Y=arma::join_cols(Y, ava::readMFCCHCopyFile(it_mfcc.next().toStdString()));
    }

    unsigned int liczbaRamekY = Y.n_rows;
    unsigned int numFeatures = Y.n_cols;

    //[1, numFeatures]
    normMean = arma::mean(Y, 0);
    //[1, numFeatures]
    normStd = arma::stddev(Y, 0, 0);
   // normYtrans =(ava::normZ(Y, normMean, normStd, liczbaRamekY, numFeatures)).t();
    normYtrans =(ava::normZ(Y.t(), normMean, normStd, liczbaRamekY, numFeatures));   

   // it_mfcc.~QDirIterator();
}

/*! [deprecated]
 * \brief statBaumWelch
 * \param N
 * \param F
 * \param Nc
 * \param pWork_path
 * \param pFeaMfcc_path
 * \param normMean
 * \param normStd
 * \param ubmModel
 * \param numOfComponents
 * \param numFilesUBM
 * \param numFeatures
 * \param prodNumOfComponentsNumFeatures
 */
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
                   unsigned int prodNumOfComponentsNumFeatures){

    /*
     * Bauma-Welcha Statistics
     * Kenny P. Boulianne G, Dumouchel P "Eigenvoice modeling with sparse training data"
     *
     */

//----start = std::chrono::high_resolution_clock::now();

    //unsigned int n_o_f  = {};
    arma::mat logLikelihood = {};
    arma::mat amax={};
    arma::mat amax2={};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat n={};
    arma::mat f={};
    arma::mat Y_bw={};
//----    arma::cube Nc(1, numOfComponents, numFilesUBM);
    arma::cube Fc(numFeatures, numOfComponents, numFilesUBM);
    unsigned int idx_tw=0;
    arma::mat tmp{};

    QDirIterator it_mfcc2(QString::fromStdString(pWork_path + "/" + pFeaMfcc_path), QStringList()<<"*.mfc", QDir::Files, QDirIterator::Subdirectories);

    while(it_mfcc2.hasNext())
    {
/*
        ava::MfcceFeatures mfcef(it_mfcc2.next().toStdString());
        Y_bw = ava::normZ(mfcef.returnArmaMat(),
                          normMean,
                          normStd,
                          (mfcef.returnArmaMat()).n_rows,
                          numFeatures);
                          */

        tmp=ava::readMFCCHCopyFile(it_mfcc2.next().toStdString());

        Y_bw = ava::normZ(tmp,
                          normMean,
                          normStd,
                          tmp.n_rows,
                          tmp.n_cols);

                /*
         *  prawdopodobienstwo a posteriori
         */

        for(unsigned int comp=0; comp<numOfComponents; comp++)
        {
            logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Y_bw.t(), comp));
        }

        /*
         *  prawdopodobienstwo a posteriori znormalizowane
         */

        amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
        logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
        gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

       /*
         *  statystyki Bauma-Welcha
         */

        n = arma::sum(gamma, 0);
        f = arma::operator*(Y_bw.t(), gamma)-arma::repelem(n, numFeatures, 1);
        Nc.slice(idx_tw)=n;
        Fc.slice(idx_tw)=f;

        idx_tw=idx_tw+1;
        logLikelihood.clear();
    }

//----    stop = std::chrono::high_resolution_clock::now();
//----    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//----    std::cout<<"Line: "<<__LINE__<<" duration sec: "<<duration.count()<<std::endl;

    //----------------------------------- Przeksztalcenie statystyk Bauma-Welcha

    /*
     * Kenny P. Ouellet P. Dehak N. Gupta V. Dumouchel P. "A study of inter-speaker variability in speaker verification"
     * N(s) -> C F x C F macierz diagonalna, ktorej bloki to Nc(s)I (c=1, ...C)
     * F(s) -> C F x 1 supervector otrzymany przez kontatancje Fc(s)(c=1, ...C)
     * C to numer komponentu UBM-GMM
     */
//----    start = std::chrono::high_resolution_clock::now();

//----    unsigned int prodNumOfComponentsNumFeatures=numOfComponents*numFeatures;

//----    arma::cube N(1, prodNumOfComponentsNumFeatures, numFilesUBM);
//----    arma::cube F(1, prodNumOfComponentsNumFeatures, numFilesUBM);
    arma::mat muc = ubmModel.means;

    for(unsigned int i=0; i<numFilesUBM; i++)
    {
        N.slice(i)=repelem(Nc.slice(i), 1, numFeatures);
        F.slice(i)=arma::reshape((Fc.slice(i)-((arma::repmat(Nc.slice(i),numFeatures, 1))%muc)), 1, prodNumOfComponentsNumFeatures);
    }

}

/*!
 * \brief statBaumWelchVector
 *
 * \param [out] NN - std::vector<arma::mat>& NN, vect. size equals
 * number of UBM train files (general sum of all files of all spk)
 * matrixes in vector equal
 * [1, number of files x gmm order]
 *
 * \param FF - std::vector<arma::mat>& NN, vect. size equals
 * number of UBM train files (general sum of all files of all spk)
 * matrixes in vector equal
 * [1, number of files x gmm order]
 *
 * \param NNcc - - std::vector<arma::mat>& NN, vect. size equals
 * number of UBM train files (general sum of all files of all spk)
 * matrixes in vector equal
 * [1, gmm order]
 *
 * \param pWork_path
 * \param pFeaMfcc_path
 * \param normMean
 * \param normStd
 * \param ubmModel
 * \param numOfComponents
 * \param numFilesUBM
 * \param numFeatures
 * \param prodNumOfComponentsNumFeatures
 */

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
                           unsigned int prodNumOfComponentsNumFeatures){

    /*
     * The Baum-Welch statistics are the N (zeroth order) and F (first order)
     * statistics used in the EM algorithm, calculated using the final UBM.
     * Calculate the zeroth and first order Baum-Welch statistics over the training set.
    */

    /*
     * Statystyki Bauma-Welcha
     * Kenny P. Boulianne G, Dumouchel P "Eigenvoice modeling with sparse training data"
     *
     */

//----start = std::chrono::high_resolution_clock::now();

    //unsigned int n_o_f  = {};
    arma::mat logLikelihood = {};
    arma::mat amax={};
    arma::mat amax2={};
    arma::mat logLikelihoodSum = {};
    arma::mat gamma={};
    arma::mat n={};
    arma::mat f={};
    arma::mat Y_bw={};
//----    arma::cube Nc(1, numOfComponents, numFilesUBM);
    //arma::cube Fc(numFeatures, numOfComponents, numFilesUBM);
    std::vector<arma::mat> FFcc{};
    unsigned int idx_tw=0;
    arma::mat tmp{};

    QDirIterator it_mfcc2(QString::fromStdString(pWork_path + "/" + pFeaMfcc_path), QStringList()<<"*.mfc", QDir::Files, QDirIterator::Subdirectories);

    while(it_mfcc2.hasNext())
    {
/*
        ava::MfcceFeatures mfcef(it_mfcc2.next().toStdString());
        Y_bw = ava::normZ(mfcef.returnArmaMat(),
                          normMean,
                          normStd,
                          (mfcef.returnArmaMat()).n_rows,
                          numFeatures);
                          */

        tmp=ava::readMFCCHCopyFile(it_mfcc2.next().toStdString());

        Y_bw = ava::normZ(tmp.t(),
                          normMean,
                          normStd,
                          tmp.n_rows,
                          tmp.n_cols).t();

                /*
         *  prawdopodobienstwo a posteriori
         */

        for(unsigned int comp=0; comp<numOfComponents; comp++){
            //logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Y_bw.t(), comp));
            logLikelihood=arma::join_cols(logLikelihood, ubmModel.log_p(Y_bw, comp));
        }

        /*
         *  prawdopodobienstwo a posteriori znormalizowane
         */

        //test1 - OK/

        amax = max(logLikelihood, 0);                                                                                    // [1, liczba ramek]
        logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1)))); // [1, l.ramek]+[1,l.ramek]
        gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

        //test1 - END/

        /*
         *  statystyki Bauma-Welcha
         */

        //test2 - OK/

        n = arma::sum(gamma, 0);
        f = arma::operator*(Y_bw, gamma); //-arma::repelem(n, numFeatures, 1);

        //test2 - END/

        //Nc.slice(idx_tw)=n;
        NNcc.push_back(n);

        //Fc.slice(idx_tw)=f;
        FFcc.push_back(f);

        idx_tw=idx_tw+1;
        logLikelihood.clear();
    }

    //----    stop = std::chrono::high_resolution_clock::now();
    //----    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    //----    std::cout<<"Line: "<<__LINE__<<" duration sec: "<<duration.count()<<std::endl;

    //----------------------------------- Przeksztalcenie statystyk Bauma-Welcha

    /*
     * Kenny P. Ouellet P. Dehak N. Gupta V. Dumouchel P. "A study of inter-speaker
     * variability in speaker verification"
     * N(s) -> C F x C F diagonal matrix, blocks: Nc(s)I (c=1, ...C)
     * F(s) -> C F x 1 supervector after concatenation Fc(s)(c=1, ...C)
     * C to numer komponentu UBM-GMM
     */
    //----    start = std::chrono::high_resolution_clock::now();

    //----    unsigned int prodNumOfComponentsNumFeatures=numOfComponents*numFeatures;

    //----    arma::cube N(1, prodNumOfComponentsNumFeatures, numFilesUBM);
    //----    arma::cube F(1, prodNumOfComponentsNumFeatures, numFilesUBM);
    arma::mat muc = ubmModel.means;

    // test3 - OK

    for(unsigned int i=0; i<numFilesUBM; i++)
    {
        //N.slice(i)=repelem(Nc.slice(i), 1, numFeatures);
        NN.push_back(arma::repelem(NNcc[i], 1, numFeatures));

        //F.slice(i)=arma::reshape((Fc.slice(i)-((arma::repmat(Nc.slice(i),numFeatures, 1))%muc)), 1, prodNumOfComponentsNumFeatures);
        FF.push_back(arma::reshape((FFcc[i]-((arma::repmat(NNcc[i],numFeatures, 1))%muc)), 1, prodNumOfComponentsNumFeatures));

    }

    //test3 - END
}

/*!
 * \brief randomTestSplit - generate two vectors<int> of uint values in the
 *                          range <0,numSpkEnrollTestFiles>
 *                          random split bewteen test and enroll
 * \param [out] ranTestNumVector    - indexes of test wav files
 * \param [out] ranEnrollNumVec     - indexes of enroll files
 * \param [int] numSpkEnrollTestFiles   - amount of test and enroll files
 * \param [int] test_size   - proportion between test and enroll
 */
void randomTestSplit(std::vector<uint>& ranTestNumVector,
                     std::vector<uint>& ranEnrollNumVec,
                     uint numSpkEnrollTestFiles,
                     float test_size){

    // how many files going to test set
    uint numSpkTest = int(test_size*numSpkEnrollTestFiles);

    // how many files going to enroll set
    uint numSpkEnroll= numSpkEnrollTestFiles - numSpkTest;

    // algorithm of randomizing the indexes goind to test set
    std::vector<uint>::iterator it;
    std::mt19937 eng{std::random_device{}()};
    for(uint i =0; i<numSpkTest; i++){
        uint q = std::uniform_int_distribution<uint>{0, numSpkTest}(eng);
        do{
            q = std::uniform_int_distribution<uint>{0, numSpkEnrollTestFiles}(eng);
            it = std::find(ranTestNumVector.begin(), ranTestNumVector.end(), q);
        }while (it != ranTestNumVector.end());
        ranTestNumVector.push_back(q);
    }

    std::vector<uint>::iterator itt;

    // take remaing values from set <0,numSpkEnrollTestFiles> to enroll set
    for(uint i=0; i<numSpkEnrollTestFiles; i++){
        itt=std::find(ranTestNumVector.begin(), ranTestNumVector.end(), i);
        if(itt == ranTestNumVector.end()){
            ranEnrollNumVec.push_back(i);
            numSpkEnroll=numSpkEnroll+1;
        }
    }
}

// /*!
//  * \brief currentTime - format string containes date time
//  * \param
//  * \return
//  */
// std::string currentTime(){

//     std::time_t now = std::time(0);
//     std::tm *ltm = std::localtime(&now);

//     std::string t1 = std::to_string(ltm->tm_sec);
//     std::string t2 = std::to_string(ltm->tm_min);
//     std::string t3 = std::to_string(ltm->tm_hour);
//     std::string t4 = std::to_string(ltm->tm_mday);
//     std::string t5 = std::to_string(1+ltm->tm_mon);
//     std::string t6 = std::to_string(1900+ltm->tm_year);

//     //std::string t = t3+"h"+t2+"m"+t1+"s"+"_"+t4+"-"+t5+"-"+t6;
//     std::string t = t4+"-"+t5+"-"+t6+"_"+t3+"h"+t2+"m"+t1+"s";
//     return t;
// }

/*!
 * \brief extractFNameFrPath - function extract file name from path
 * \param pathToMFCFile
 * \return
 */
std::string extractFNameFrPath(std::string pathToMFCFile){

    // get filenam
    std::string base_filename = pathToMFCFile.substr(pathToMFCFile.find_last_of("/")+1);

    // remove extension
    std::string::size_type const p(base_filename.find_last_of('.'));
    std::string file_without_extension = base_filename.substr(0,p);

    return file_without_extension;
}

/*!
 * \brief saveEnrolFRRFARlList - function store the FRR and FAR
 * \param spkLabel  -   speaker label
 * \param ss        -   path to mfcc file
 * \param pathToSave    -   path to save file
 * \param nameOfTheList -   name of file list
 * \return
 */
void saveEnrolFRRFARlList(std::string spkLabel,
                             std::string ss,
                             std::string pathToSave,
                             std::string nameOfTheList){

    std::string base_file = ava::extractFNameFrPath(ss);

    std::ofstream logToFile;
    logToFile.open(pathToSave + "/" + nameOfTheList, std::ios::app);
    logToFile<<spkLabel<<" "<<ava::extractFNameFrPath(ss)<<std::endl;
}

std::string readFile(const std::string &file)
{
    std::ifstream is(file);
    if( !is.good() ){
        throw std::runtime_error("Error: stream has errors.");
    }
    std::stringstream ss;
    ss << is.rdbuf();
    std::string m;
    // Remove ending line character '\n' or '\r\n'.
    std::getline(ss, m);
    return m;

}


}
