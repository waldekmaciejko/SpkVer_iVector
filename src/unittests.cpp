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

#include"unittests.h"

namespace avatest {


void test1(){

    arma::mat logLikelihood = {{1,2,3},{1,2,3}};
    arma::mat amax;
    arma::mat logLikelihoodSum;
    arma::mat gamma;

    uint numOfComponents = 2;

    amax = max(logLikelihood, 0);
    logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1))));
    gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

    gamma.print();
}


void test2(){

    arma::mat logLikelihood = {{1,2,3},{1,2,3}};
    arma::mat Y_bw={{11,12,13},{14,15,16}};
    //arma::mat Y_bw={{11,12},{13,14},{15,16}};
    arma::mat amax;
    arma::mat logLikelihoodSum;
    arma::mat gamma;
    arma::mat n;
    arma::mat f;

    uint numOfComponents = 2;
    uint numFeatures = 2;

    amax = max(logLikelihood, 0);
    logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1))));
    gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();

    std::cout<<"gamma"<<std::endl;
    gamma.print();

    n = arma::sum(gamma, 0);
    f = arma::operator*(Y_bw, gamma); //-arma::repelem(n, numFeatures, 1);

    std::cout<<"Y"<<std::endl;
    Y_bw.print();

    //Nc.slice(idx_tw)=n;
    std::cout<<"n"<<std::endl;
    n.print();

    //Fc.slice(idx_tw)=f;
    std::cout<<"f"<<std::endl;
    f.print();
}

void test3(){

    uint numFilesUBM = 3;
    uint numFeatures = 4;
    uint prodNumOfComponentsNumFeatures = 8;
    uint numOfComponents = 2;

    arma::mat muc(numFeatures, 2);
    muc={{1,2},{1,2},{1,2},{1,2}};
    std::vector<arma::mat> NNcc{};
    std::vector<arma::mat> FFcc{};

    arma::mat n1={1,2};
    arma::mat n2={3,4};
    arma::mat n3={5,6};

    NNcc.push_back(n1);
    NNcc.push_back(n2);
    NNcc.push_back(n3);

    arma::mat f1={{1,2},{1,2},{1,2},{1,2}};
    arma::mat f2={{1,2},{1,2},{1,2},{1,2}};
    arma::mat f3={{1,2},{1,2},{1,2},{1,2}};

    FFcc.push_back(f1);
    FFcc.push_back(f2);
    FFcc.push_back(f3);

    for(unsigned int i=0; i<numFilesUBM; i++)
    {
        (arma::repelem(NNcc[i], 1, numFeatures)).print();
        (arma::reshape((FFcc[i]-((arma::repmat(NNcc[i],numFeatures, 1))%muc)), 1, prodNumOfComponentsNumFeatures)).print();

    }
}

void test4(){

    uint numFilesUBM = 3;
    uint numTdim = 2;

    arma::mat T={{1, 2}, {3, 4}, {5, 6}, {7, 8}, {9, 10}, {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}, {21, 22}, {23, 24}, {23, 24}, {23, 24}, {23, 24}, {23, 24}};
    arma::mat I = arma::eye(2, 2);

    std::vector<arma::mat> Ey{};
    Ey.reserve(numFilesUBM);

    std::vector<arma::mat> Eyy{};
    Eyy.reserve(numFilesUBM);

    std::vector<arma::mat> Linv{};
    Linv.reserve(numFilesUBM);

    arma::mat sigma = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    sigma=sigma.t();

    arma::mat repSigma = arma::repmat(sigma, 1, numTdim);

    arma::mat L{};

    arma::mat TtimeInvSSDiag = {};

    std::vector<arma::mat> NN{};
    std::vector<arma::mat> FF{};
    std::vector<arma::mat> NNcc{};

    arma::mat n1={1,2,3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    arma::mat n2={1,2,3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    arma::mat n3={1,2,3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    NN.push_back(n1);
    NN.push_back(n2);
    NN.push_back(n3);

    arma::mat f1={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    arma::mat f2={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    arma::mat f3={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    FF.push_back(f1);
    FF.push_back(f2);
    FF.push_back(f3);

    //1.---- oblicz dystrybucje posteriori zmiennych ukrytych
        TtimeInvSSDiag =T/repSigma;
        std::cout<<"TtimeInvSSDiag"<<std::endl;
        TtimeInvSSDiag.print();

        for(unsigned int s=0; s<numFilesUBM; s++){
            L = I + TtimeInvSSDiag.t()%arma::repmat((NN.at(s)).t(), 1, numTdim).t()*T;
            std::cout<<"L"<<std::endl;
            L.print();

            std::cout<<"Linv"<<std::endl;
            Linv.push_back(arma::pinv(L));
            (arma::pinv(L)).print();

            std::cout<<"Ey"<<std::endl;
            (Linv.at(s) * TtimeInvSSDiag.t()*FF.at(s).t()).print();
            Ey.push_back(Linv.at(s) * TtimeInvSSDiag.t()*FF.at(s).t());

            std::cout<<"Eyy"<<std::endl;
            (Linv.at(s) + Ey.at(s)*Ey.at(s).t()).print();

            }
      }


void test5(){

    uint numFilesUBM = 11;
    uint numTdim = 5;
    uint numOfComponents=2;
    uint numFeatures = 4;
    uint prodNumOfComponentsNumFeatures = 8;
    uint counter=0;

    arma::mat T = {{1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5},
                   {1,2,3,4,5}};

    arma::mat I = arma::eye(5, 5);

    std::vector<arma::mat> Ey{};
    Ey.reserve(numFilesUBM);

    std::vector<arma::mat> Eyy{};
    Eyy.reserve(numFilesUBM);

    std::vector<arma::mat> Linv{};
    Linv.reserve(numFilesUBM);

    arma::mat sigma = {1, 2, 3, 4, 5, 6, 7, 8};
    sigma=sigma.t();

    arma::mat repSigma = arma::repmat(sigma, 1, numTdim);

    arma::mat L{};

    arma::mat TtimeInvSSDiag = {};

    std::vector<arma::mat> NN{};
    std::vector<arma::mat> FF{};
    std::vector<arma::mat> NNcc{};

    arma::mat n1={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n2={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n3={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n4={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n5={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n6={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n7={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n8={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n9={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n10={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat n11={1, 2, 3, 4, 5, 6, 7, 8};

    NN.push_back(n1);
    NN.push_back(n2);
    NN.push_back(n3);
    NN.push_back(n4);
    NN.push_back(n5);
    NN.push_back(n6);
    NN.push_back(n7);
    NN.push_back(n8);
    NN.push_back(n9);
    NN.push_back(n10);
    NN.push_back(n11);

    arma::mat f1={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f2={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f3={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f4={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f5={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f6={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f7={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f8={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f9={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f10={1, 2, 3, 4, 5, 6, 7, 8};
    arma::mat f11={1, 2, 3, 4, 5, 6, 7, 8};

    FF.push_back(f1);
    FF.push_back(f2);
    FF.push_back(f3);
    FF.push_back(f4);
    FF.push_back(f5);
    FF.push_back(f6);
    FF.push_back(f7);
    FF.push_back(f8);
    FF.push_back(f9);
    FF.push_back(f10);
    FF.push_back(f11);

    NNcc.push_back(f1);
    NNcc.push_back(f2);
    NNcc.push_back(f3);
    NNcc.push_back(f4);
    NNcc.push_back(f5);
    NNcc.push_back(f6);
    NNcc.push_back(f7);
    NNcc.push_back(f8);
    NNcc.push_back(f9);
    NNcc.push_back(f10);
    NNcc.push_back(f11);

    uint numIterationsTV=2;

    std::vector<arma::mat> K{};
    K.reserve(numOfComponents);

    std::vector<arma::mat> newT{};
    newT.reserve(numOfComponents);
    arma::mat tmp={};

    //---------------------------------
    for(unsigned int i=0; i<numIterationsTV; i++)
    {
        //1.---- oblicz dystrybucje posteriori zmiennych ukrytych
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

        //2.---- akumulacja statystyk w dla kazdego mowcy
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
        //3.---- aktualizacjaprzestrzeni TV
            newT.push_back((arma::pinv(AcLocal)*(K.at(c).t())).t());
        }

        tmp.clear();

        for(unsigned int cc=0; cc<numOfComponents; cc++){
            //tmp = arma::join_cols(tmp, newT.slice(cc));
            tmp = arma::join_cols(tmp, newT.at(cc));
        }
        T=tmp;
       }
     std::cout<<"T"<<std::endl;
     T.print();
}

void test6(){

    arma::mat amax={};
    arma::mat gamma={};

    std::set<std::string> setSpk{};
    setSpk.insert("M1");
    setSpk.insert("M2");

    std::multimap<std::string, std::string> multimapMowcyMFC{};
    multimapMowcyMFC.insert(std::make_pair("M1", "M1"));
    multimapMowcyMFC.insert(std::make_pair("M2", "M2"));

    unsigned int liczbaRamekYivector = 6;

    // [numFeatures, numeFrames]
    arma::mat Yivector={ {1,2,3,4,5,6},
                         {1,2,3,4,5,6},
                         {1,2,3,4,5,6}};

    // [numComponents, numeFrames]
    arma::mat logLikelihood ={{1,2,3,4,5,6},
                             {1,2,3,4,5,6},
                             {1,2,3,4,5,6},
                             {1,2,3,4,5,6}};

    // [numFeatues, numComponents]
    arma::mat ubmMu = {{1, 2, 3, 4},
                       {1, 2, 3, 4},
                       {1, 2, 3, 4}};

    // [numComponents x numFeatures, numTdim]
    // [9, 4]
    arma::mat TS={{1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4},
                    {1,2,3,4}};

    // [numTdim, numComponents x numFeatures]
    // [4, 9]
    arma::mat TSi= {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                    {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
                    {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
                    {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};

    // [numComponents x numFeatures, numTdim]
    arma::mat T= {{1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4},
                     {1,2,3,4}};

    uint numTdim = 4;

    arma::mat I=arma::eye(numTdim, numTdim);

    arma::mat logLikelihoodSum = {};
    arma::mat n{};
    arma::mat f{};
    arma::mat tmp6;

    uint numOfComponents = 4;
    uint numFeatures = 3;

    for(auto st : setSpk)
    {
        //std::cout<<"setSpk"<<st<<std::endl;
        for(auto para : multimapMowcyMFC)
        {
            //arma::mat YY{};
            //if(para.first==st){

                //ava::MfcceFeatures mfcefivector(para.second);
                //liczbaRamekYivector = mfcefivector.returnArmaMat().n_rows;
                //Yivector =(ava::normZ(mfcefivector.returnArmaMat().t(), normMean, normStd, liczbaRamekYivector, numFeatures)).t();
                //logLikelihoodivector = ubmModel.log_p(Yivector.t());

               /* for(unsigned int comp=0; comp<numOfComponents; comp++) {
                    logLikelihoodivector=arma::join_cols(logLikelihoodivector, ubmModel.log_p(Yivector, comp));
                }
                */

                amax = max(logLikelihood, 0);
                //amax.print();
                logLikelihoodSum = amax + arma::log(arma::sum(arma::exp(logLikelihood-arma::repmat(amax, numOfComponents, 1))));
                //logLikelihoodSum.print();
                gamma = arma::exp(logLikelihood-(arma::repmat(logLikelihoodSum, numOfComponents, 1))).t();
                //gamma.print();

                n = arma::sum(gamma, 0);
                //n.print();
                f = arma::operator*(Yivector, gamma)- arma::repmat(n, numFeatures, 1)%ubmMu;
                f.print();
                std::cout<<"---"<<std::endl;

                tmp6 = (arma::pinv(I+(ava::multipicColByCol(TS, (arma::repmat(n, 1, numFeatures))).t()*T))*TSi)*arma::vectorise(f, 0);
                tmp6.print();
                std::cout<<"----"<<std::endl;

                //ivectorPerSpk=arma::join_rows(ivectorPerSpk, tmp6);
                //ivectorsTrain = arma::join_rows(ivectorsTrain, tmp6);
            //}
            //logLikelihoodivector.clear();
        }
        //ivectorPerFile.insert(std::make_pair(st, ivectorPerSpk));
        //ivectorPerSpk.clear();
    }
    //return ivectorPerFile;
    }

void test7(){
    /*
     * zweryfikowano wartosci: wbar, Sb oraz Sw
     */

    uint numTdim = 20;
    uint numEigenvectors = 4;

    arma::mat ivectorsTrain{};

   ivectorsTrain = {{500, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
           {400, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21},
           {300, 4, 5, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 22},
           {200, 5, 6, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 23},
           {100, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24},
           {500, 7, 8, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 25},
           {400, 8, 9, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 26},
           {300, 9, 10, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 27},
           {200, 10, 11, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 28},
           {100, 11, 12, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 29},
           {50, 12, 13, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 30},
           {500, 13, 14, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 31},
           {400, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 32},
           {300, 15, 16, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 33},
           {200, 16, 17, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 34},
           {100, 17, 18, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 35},
           {500, 18, 19, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 36},
           {400, 19, 20, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 37},
           {300, 20, 21, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 38}};

    arma::inplace_trans(ivectorsTrain);
    //ivectorsTrain.print();

    std::map<std::string, arma::mat> ivectorPerSpk{};
    ivectorPerSpk.insert(std::make_pair("M1", ivectorsTrain.submat(arma::span(0, 19),arma::span(0,4))));
    ivectorPerSpk.insert(std::make_pair("M2", ivectorsTrain.submat(arma::span(0, 19),arma::span(5,10))));
    ivectorPerSpk.insert(std::make_pair("M3", ivectorsTrain.submat(arma::span(0, 19),arma::span(11,15))));
    ivectorPerSpk.insert(std::make_pair("M4", ivectorsTrain.submat(arma::span(0, 19),arma::span(16,18))));

   // std::vector<uint> utterancePerSpeaker{};

    arma::mat Sw = arma::zeros(numTdim, numTdim); // between variability
    arma::mat Sb = arma::zeros(numTdim, numTdim); // within variability
    arma::mat ws={};
    arma::mat wsbar={};
    arma::mat wbar = arma::mean(ivectorsTrain, 1);

    std::cout<<"wbar"<<std::endl;
    wbar.print();

    for(const auto& x : ivectorPerSpk){
        ws = x.second;

        wsbar = arma::mean(ws, 1);
        //std::cout<<"wsbar"<<std::endl;
        //wsbar.t().print();

        Sb = Sb + (wsbar - wbar)*((wsbar - wbar).t());
        //std::cout<<"Sb"<<std::endl;
        //Sb.print();

        Sw = Sw + arma::cov(ws.t(),1);
        //std::cout<<"Sw"<<std::endl;
        //Sw.print();
    }

    std::cout<<"Sb"<<std::endl;
    Sb.print();
    std::cout<<"Sw"<<std::endl;
    Sw.print();

}

void test8(){
    /*
     * werify algorithm eigen vectors extractor to solve general equation
     * and normalization matrix A
     *
     */

    uint  numEigenvectors = 4;
    uint numTdim = 10;


    //arma::mat Sb = arma::randn(numEigenvectors, numEigenvectors);
    //arma::mat Sw = arma::randn(numEigenvectors, numEigenvectors);

    //Sw.save("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\Sw.csv", arma::csv_ascii);
    //Sb.save("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\Sb.csv", arma::csv_ascii);

    arma::mat Sw{};
    arma::mat Sb{};

    Sb.load("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\Sb.csv", arma::csv_ascii);
    Sw.load("j:\\__Badania__\\Qt_maw200\\workspace_maw200\\tmp\\Sw.csv", arma::csv_ascii);

    //arma::mat Sb{{10,20,40},{50,60,70},{80,90,100}};
    //arma::mat Sw{{10,20,40},{50,60,70},{80,90,100}};

    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_pair(eigval, eigvec, Sb, Sw);

    //eigvec.print();
    //eigval.print();

    arma::mat eigval_real = arma::real(eigval);
    arma::mat eigvec_real = arma::real(eigvec);

    eigval_real.replace(arma::datum::inf, 0);
    eigvec_real.replace(arma::datum::inf, 0);

    //eigval_real.print();

    numEigenvectors = numEigenvectors - 1;

    arma::mat A={};
    arma::mat Avecnorm{};
    arma::mat tmp{};
    A = eigvec_real.cols(0, numEigenvectors);

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

    Anorm.print();

    int stop=0;
}

void test9(){
    /*
     *
     * verify WCCN
     *
     */

    uint numEigenvectors=4;
    arma::mat projectionMatrix={{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
                                {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
                                {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
                                {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}};

    std::vector<arma::mat> w{};

    arma::mat w1={{500, 2, 3, 4},
                    {400, 3, 4, 5},
                    {300, 4, 5, 4},
                    {200, 5, 6, 4},
                    {100, 6, 7, 4}};

    arma::mat w2={{500, 7, 8, 4},
                    {400, 8, 9, 4},
                    {300, 9, 10, 4},
                    {200, 10, 11, 4},
                    {100, 11, 12, 4},
                    {50, 12, 13, 4}};

    arma::mat w3={{500, 13, 14, 4},
                    {400, 14, 15, 4},
                    {300, 15, 16, 4},
                    {200, 16, 17, 4},
                    {100, 17, 18, 4}};

    arma::mat w4={{500, 18, 19, 4},
                    {400, 19, 20, 4},
                    {300, 20, 21, 4}};

    w.push_back(w1.t());
    w.push_back(w2.t());
    w.push_back(w3.t());
    w.push_back(w4.t());

    float alpha = 0.9;
    arma::mat B={};
    arma::mat W(numEigenvectors, numEigenvectors, arma::fill::zeros);

    for(auto ww : w){
        W = W + arma::cov(ww.t());
    }

    W=W/int(w.size());
    W=(1-alpha) * W + alpha * arma::eye(numEigenvectors, numEigenvectors);
    B=arma::chol(arma::pinv(W), "lower");
    projectionMatrix = B * projectionMatrix;
    projectionMatrix.print();
}

void test10(){

    arma::mat S{};
    arma::mat U_tmp{};
    arma::vec s_tmp{};
    arma::mat sD{};
    arma::mat sV{};
    arma::mat W{};

    uint r = 4;
    uint c =10;

    //arma::mat ivectorMatrix = arma::randn(r, c);
    //ivectorMatrix.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ivectorMatrix.csv", arma::csv_ascii);

    arma::mat ivectorMatrix{};
    ivectorMatrix.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ivectorMatrix.csv", arma::csv_ascii);
    S = arma::cov(ivectorMatrix.t());

    /*
     * SVD - Singular Value Decomposition
     * svd( mat U, vec s, mat V, mat X )
     * X = U*diagmat(s)*V.t()
     *
     */

    arma::svd(U_tmp, s_tmp, sV, S);

    //sD = arma::diagmat(s_tmp);
    //sV.print();
    W = arma::diagmat(1/(arma::sqrt(s_tmp) + arma::datum::eps))*sV.t();
    //W.print();
    ivectorMatrix = W * ivectorMatrix;

    ivectorMatrix.print();

}

void test11(){

    arma::mat S{};
    //arma::cx_vec eigval;
    //arma::cx_mat eigvec;

    arma::vec eigval;
    arma::mat eigvec;

    arma::mat sD{};
    arma::mat sV{};
    arma::vec s_tmp{};
    arma::mat W{};

    arma::mat ivectorMatrix{};
    ivectorMatrix.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ivectorMatrix.csv", arma::csv_ascii);


    S = arma::cov(ivectorMatrix.t());

    arma::eig_sym(eigval, eigvec, S);

    eigval.print();
    eigvec.print();


    W = arma::diagmat(1/(arma::sqrt(eigval) + arma::datum::eps))*eigvec.t();
    ivectorMatrix = W * ivectorMatrix;
    ivectorMatrix.print();
}

// ivectorMatrixAfterNorm()

void test12(){

    arma::mat ivectorMatrix{};
    ivectorMatrix.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\ivectorMatrix.csv", arma::csv_ascii);

    arma::mat ivectorMatrixAfterNorm(ivectorMatrix.n_rows, ivectorMatrix.n_cols, arma::fill::randn);

    ava::divdeMatrixRBRVectorNorm(ivectorMatrix,
                                  ivectorMatrixAfterNorm);

    ivectorMatrixAfterNorm.print();
    arma::mat S = ivectorMatrixAfterNorm*ivectorMatrixAfterNorm.t();
    S.print();
}

void test13(){
//    ava::trainGPLDA(numEigenVoices,
//                        V,
//                        lambda,
//                        numIter,
//                        minimumDivergence,
//                        utterancePerSpeaker,
//                        K,
//                        ivectorSorted,
//                        ff,
//                        numFilesUBM,
//                        S);

    //- beginning parameters
    bool minimumDivergence = true;
    uint numEigenVoices = 10;
    uint numEigenvectors = 8;
    uint numIter=4;
    uint numFilesUBM = 4;
    uint K = numFilesUBM;
    uint numTdim = 5;

    //- test data initialization

    /*
     * unittest data: four speakers
     * 1. 3 utterances
     * 2. 4 utterances
     * 3. 4 utterances
     * 4. 3 utterances
     *
     */

    std::vector<uint> utterancePerSpeaker = {3, 4, 4, 3};

    /*
     * unittest ivectors for four speakers
     * 1. [numEigenvectors, number of utterances] - [numEigenvectors, 3]
     * 2. [numEigenvectors, number of utterances] - [numEigenvectors, 4]
     * 3. [numEigenvectors, number of utterances] - [numEigenvectors, 4]
     * 4. [numEigenvectors, number of utterances] - [numEigenvectors, 3]
     *
     */

    std::vector<std::vector<arma::mat>> ivectorSorted;

    arma::mat spk1 = arma::randn(numEigenvectors, 3);
    //spk1.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk1.csv", arma::csv_ascii);
    spk1.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk1.csv", arma::csv_ascii);

    arma::mat spk2 = arma::randn(numEigenvectors, 4);
    //spk2.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk2.csv", arma::csv_ascii);
    spk2.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk2.csv", arma::csv_ascii);

    arma::mat spk3 = arma::randn(numEigenvectors, 4);
    //spk3.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk3.csv", arma::csv_ascii);
    spk3.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk3.csv", arma::csv_ascii);

    arma::mat spk4 = arma::randn(numEigenvectors, 3);
    //spk4.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk4.csv", arma::csv_ascii);
    spk4.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\spk4.csv", arma::csv_ascii);

    std::vector<arma::mat> interU3{};
    interU3.push_back(spk1); //divided into number of utterances
    interU3.push_back(spk4);

    std::vector<arma::mat> interU4{};
    interU4.push_back(spk2);
    interU4.push_back(spk3);

    ivectorSorted.push_back(interU3);
    ivectorSorted.push_back(interU4);

    //--

    arma::mat lambda = arma::randn(numEigenvectors, numEigenvectors);
    //lambda.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\lambda.csv", arma::csv_ascii);
    lambda.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\lambda.csv", arma::csv_ascii);

    arma::mat V = arma::randn(numEigenvectors, numEigenVoices);
    //V.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\V.csv", arma::csv_ascii);
    V.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\V.csv", arma::csv_ascii);

    arma::mat ff = arma::randn(numEigenvectors, 4);
    //ff.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\ff.csv", arma::csv_ascii);
    ff.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\ff.csv", arma::csv_ascii);

    arma::mat S = arma::randn(numEigenvectors, numEigenvectors);
    //S.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\S.csv", arma::csv_ascii);
    S.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\gpldaUnitTest\\S.csv", arma::csv_ascii);

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
            iv = ivectorSorted[lengthIndex];

            // oblicz wartosc M
            // Equation (A.7) in Unifying probabilistic linear discriminant
            // analysis variants in biometric authentication
            // M matrix should be symmetric
            M = arma::pinv(ivectroLength * (V.t()*(lambda*V)) + (arma::eye(numEigenVoices, numEigenVoices)));
            //std::cout<<"M"<<std::endl;
            //M.print();
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
        lambda.print();

        arma::mat tmp{};

        if(minimumDivergence){
            // Equation (A.18) in [13]
             gamma = gamma/K;
             // Equation (A.22) in [13]
             // to pozniej V = V*arma::chol(gamma, "lower");
        }
    }
}

void test14(){

    arma::gmm_diag ubmModel;
    unsigned int numOfComponents=4;
    unsigned int numIterationsEM = 4;
    unsigned int numIterationsKmean_init = 10;

    arma::mat normYtrans = arma::randn(300,9);
    //normYtrans.save("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\Y.csv", arma::csv_ascii);
    normYtrans.load("j:\\__Badania__\\QtGUI_ava300\\ava300\\WORKSPACE\\tmp\\Y.csv", arma::csv_ascii);
/*
    bool status = ubmModel.learn(normYtrans.t(),
                                  numOfComponents,
                                  arma::eucl_dist,
                                  arma::random_subset,
                                  numIterationsKmean_init,
                                  numIterationsEM,
                                  1e-10,
                                  false);
                              */

    arma::rowvec weight={{0.2500, 0.2500, 0.2500, 0.2500}};

    arma::mat mu={{2.4558, 1.9005, -0.0701, 0.6122},
                   {-0.7783, 0.6381, 1.2347, -0.1949},
                    {1.8015, 2.6661, 0.1998, -0.6663},
                   {-0.5566, 0.7312, 0.0943, -1.5527},
                   {-0.4653, 0.1133, 0.5957, -1.5338},
                   {-1.1844, 1.9846, -0.7560, 0.3537},
                   {-1.0832, -2.4076, 0.0966, 1.7619},
                   {-0.0316, -1.5672, 0.2685, 0.9029},
                    {-1.2272, 0.4977, 0.0623, 0.5677}};

    arma::mat sigma={{0.4698, 0.6243, 0.1381, 0.4568},
                    {0.6441, 0.7638, 0.4073, 0.0442},
                    {0.4959, 0.3749, 0.2665, 0.9273},
                    {0.6643, 0.4567, 0.7901, 0.6557},
                    {0.2878, 0.9497, 0.6723, 0.4801},
                    {0.3855, 0.9073, 0.5667, 0.9538},
                    {0.7616, 0.1714, 0.7525, 0.5646},
                    {0.7702, 0.2603, 0.4011, 0.1949},
                    {0.6212, 0.4502, 0.8481, 0.3739}};

    ubmModel.set_params(mu, sigma, weight);

    arma::mat logLikelihoodivector{};

    std::cout<<"----"<<std::endl;
    ubmModel.means.print();
    std::cout<<"----"<<std::endl;
    ubmModel.dcovs.print();
    std::cout<<"----"<<std::endl;
    ubmModel.hefts.print();
    std::cout<<"----"<<std::endl;

    for(unsigned int comp=0; comp<numOfComponents; comp++) {
        logLikelihoodivector=arma::join_cols(logLikelihoodivector, ubmModel.log_p(normYtrans.t(),comp));
        logLikelihoodivector.print();
        std::cout<<"---"<<std::endl;
    }



}

void test15(){
// check ubmModel.log_p(Y, comp)



}

}


