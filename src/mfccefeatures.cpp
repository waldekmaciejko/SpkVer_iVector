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

#include "mfccefeatures.h"
#include "funhelpers.h"

namespace ava{

MfcceFeatures::MfcceFeatures()
{
    //ctor
}

MfcceFeatures::MfcceFeatures(std::string path_to_mfc)
{
    path_to_mfc=winPathToPosixPath(path_to_mfc.c_str());
    std::ifstream fin(path_to_mfc, std::ios::in | std::ios::binary);

    //-------------- odczyt naglowka

    fin.read((char*) &(this->n_o_f), sizeof(int));
    fin.read((char*) &(this->smap_period), sizeof(int));
    fin.read((char*) &(this->s_size), sizeof(short int));
    fin.read((char*) &(this->htk_code), sizeof(short int));

    //-------------- okreslenie rozmiaru pliku

    unsigned int l_baj_in_file=get_file_size(path_to_mfc);

    unsigned int l_cech_w_pliku = (l_baj_in_file-12)/4; // liczba cech w pliku to iloczyn liczby kolumn i liczby wierszy
                                                        // (l_baj_in_file - 12 (naglowek)) podzilieć przez 4 bo 4 bajt yto jedna cecha

    this->l_col=l_cech_w_pliku/n_o_f;                   // liczba cech (kolumn w wektorze) w pliku

    //Eigen::MatrixXf mm(n_o_f, l_col);                 // dla macierzy typu Eigen
    //this->m.resize(n_o_f, l_col);

    this->mm.resize(n_o_f, l_col);                      // alokacja dla armadillo
    this->vv_e.resize(n_o_f, 1);

    for(unsigned int k=0; k<n_o_f; k++)                 // zapiesz do macierzy
    {
        for(unsigned  int w=0; w<l_col; w++)
        {
            fin.read((char*) &p1, sizeof(float));
        //    this->m(k,w) = p1; // dla macierzy typu Eigen
            this->mm(k,w)= p1;
        }
    }

    unsigned int id_e = l_col-1;                            // ustalenie indexu energii l_col-1
                                                            // poniewaz ostatni kolumna to E

    //this->v_e.resize(n_o_f);                              // dla macierzy typu Eigen
    //this->v_e=this->m.col(id_e);                          // dla macierzy typu Eigen

    this->vv_e.resize(n_o_f);                               // dla macierzy typu Armadillo - arma
    this->vv_e=this->mm.col(id_e);

    // uwaga sprawdzic wartosci v_e porownojac z danymi z HList
}

MfcceFeatures::~MfcceFeatures()
{
    //dtor
}
}
