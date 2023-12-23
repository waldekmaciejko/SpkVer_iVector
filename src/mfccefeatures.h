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

#pragma once

#include <Eigen/Dense>
#include <armadillo>

namespace  ava{

// do konstruktora klasy przekazywana jest sciezka dostepu do pliku HTK-Hcopy z cechami mfcc
// nastepnie nastpue odczyt tych parametrow do tablicy Eigen zgodnie z parametrami klasy
// umieszczonymi w sekcji publicznej

//using namespace arma;

class MfcceFeatures
{
    public:
        MfcceFeatures();
        // konstruktor otrzymuje sciezke do pliku HTK
        // plik ten jest dekodowany i doczytane wartosci
        // wspolczynnikow MFCC
        MfcceFeatures(std::string path_to_mfc);

        // dla macierzy typu Eigen
        /*
        Eigen::MatrixXf returnMatrixXf()
        {
            return this->m;
        }
        */

        arma::mat returnArmaMat()
        {
            return this->mm;

        }

        arma::mat returnArmaVec()
        {
            return this->vv_e;

        }

        unsigned int returnN_O_F()
        {
            return this->n_o_f;
        }

        unsigned int returnL_COL()
        {
            return this->l_col;
        }

        virtual ~MfcceFeatures();


    protected:

    public:


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
};
}

